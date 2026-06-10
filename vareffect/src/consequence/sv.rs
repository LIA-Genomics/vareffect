//! Structural-variant (CNV) interval annotation.
//!
//! Where [`super::annotate`] consumes a single ref/alt allele pair, copy-number
//! variants are described as a genomic *interval* plus a direction (deletion or
//! duplication). This module annotates such intervals against every overlapping
//! transcript and assigns the SV-shaped SO consequence terms VEP uses for CNVs:
//! `transcript_ablation` / `feature_truncation` (deletions) and
//! `transcript_amplification` / `feature_elongation` (duplications).
//!
//! The emitted [`CnvConsequenceResult`] also carries per-transcript exon-overlap
//! geometry (5'/3'/last-exon coverage, overlapped exon count, CDS bases
//! affected) so downstream callers can drive region-centric ACMG/ClinGen CNV
//! scoring (Riggs 2020, PMID 31690835) without re-deriving exon math.

use super::complex::compute_cds_projected_length;
use super::{Consequence, Impact};
use crate::error::VarEffectError;
use crate::locate::LocateIndex;
use crate::transcript::TranscriptStore;
use crate::types::{Strand, TranscriptModel};

/// Direction of a structural copy-number event.
///
/// Riggs 2020 scores losses and gains with separate rubrics; the SV
/// consequence term also depends on direction (deletion -> ablation/truncation,
/// duplication -> amplification/elongation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SvKind {
    /// Copy-number loss (`<DEL>`, `g.A_Bdel`). A deletion of the interval.
    Deletion,
    /// Copy-number gain (`<DUP>`, `g.A_Bdup`). A duplication of the interval.
    Duplication,
}

/// Per-transcript annotation of a structural-variant interval.
///
/// One result is produced for each transcript whose `[tx_start, tx_end)`
/// overlaps the queried interval. Coordinates and exon numbering follow the
/// crate-wide 0-based / transcript-5'->3' conventions.
#[derive(Debug, Clone, PartialEq)]
pub struct CnvConsequenceResult {
    /// RefSeq transcript accession with version (e.g. `"NM_006772.2"`).
    pub transcript: String,
    /// HGNC gene symbol.
    pub gene_symbol: String,
    /// HGNC identifier including the `HGNC:` prefix, if present on the model.
    pub hgnc_id: Option<String>,
    /// Transcript strand.
    pub strand: Strand,
    /// SO consequence term assigned to this transcript for the SV.
    pub consequence: Consequence,
    /// VEP IMPACT for [`Self::consequence`].
    pub impact: Impact,
    /// Whether the SV interval fully spans `[tx_start, tx_end)`. `true` ->
    /// ablation/amplification; `false` -> truncation/elongation.
    pub spans_whole_transcript: bool,
    /// Fraction of the transcript span `[tx_start, tx_end)` covered by the SV
    /// interval, in `[0.0, 1.0]`.
    pub transcript_overlap_fraction: f64,
    /// Number of CDS bases of this transcript that fall inside the SV interval.
    /// `0` for non-coding transcripts or intronic/UTR-only overlaps.
    pub cds_bases_overlapped: u32,
    /// Number of exons (of any kind) that the SV interval overlaps.
    pub exons_overlapped: u16,
    /// `true` if the SV interval overlaps the transcript's most 5' exon
    /// (exon 1 in transcript order).
    pub overlaps_first_exon: bool,
    /// `true` if the SV interval overlaps the transcript's most 3' exon
    /// (the last exon in transcript order).
    pub overlaps_last_exon: bool,
    /// `true` if a breakpoint falls strictly inside the transcript span
    /// (i.e. the SV starts after `tx_start` and/or ends before `tx_end`).
    /// Always `false` when [`Self::spans_whole_transcript`] is `true`.
    pub breakpoint_within_transcript: bool,
    /// `true` if **both** SV breakpoints fall strictly inside the transcript
    /// span — i.e. `sv_start > tx_start && sv_end < tx_end`. This is the exact
    /// intragenic predicate (the SV is fully contained in `[tx_start, tx_end)`),
    /// correct for single-exon transcripts and for intragenic SVs that touch a
    /// terminal exon, where the first/last-exon flags would mislead. Drives the
    /// Riggs 2020 "complete containment within an HI gene" tests (2E loss / 2I
    /// gain, PMID 31690835).
    pub both_breakpoints_within_tx: bool,
}

/// Exon-overlap geometry of an SV interval against one transcript.
struct ExonOverlap {
    cds_bases_overlapped: u32,
    exons_overlapped: u16,
    overlaps_first_exon: bool,
    overlaps_last_exon: bool,
}

/// Compute exon-overlap geometry for `[sv_start, sv_end)` against `transcript`.
///
/// Exons in `transcript.exons` are ordered 5'->3' on the transcript, so the
/// most 5' exon is index `0` and the most 3' exon is index `len-1`, regardless
/// of strand. Genomic-coordinate overlap is therefore strand-agnostic.
fn compute_exon_overlap(sv_start: u64, sv_end: u64, transcript: &TranscriptModel) -> ExonOverlap {
    let last_idx = transcript.exons.len().saturating_sub(1);
    let mut exons_overlapped = 0u16;
    let mut overlaps_first_exon = false;
    let mut overlaps_last_exon = false;

    for (i, exon) in transcript.exons.iter().enumerate() {
        let lo = sv_start.max(exon.genomic_start);
        let hi = sv_end.min(exon.genomic_end);
        if lo < hi {
            exons_overlapped += 1;
            if i == 0 {
                overlaps_first_exon = true;
            }
            if i == last_idx {
                overlaps_last_exon = true;
            }
        }
    }

    ExonOverlap {
        cds_bases_overlapped: compute_cds_projected_length(sv_start, sv_end, transcript),
        exons_overlapped,
        overlaps_first_exon,
        overlaps_last_exon,
    }
}

/// Annotate one transcript against an SV interval.
///
/// `[sv_start, sv_end)` is the 0-based half-open SV footprint. The transcript
/// is assumed to overlap it (the caller filters via `query_overlap`).
fn annotate_transcript(
    sv_start: u64,
    sv_end: u64,
    sv_kind: SvKind,
    transcript: &TranscriptModel,
) -> CnvConsequenceResult {
    let spans_whole_transcript = sv_start <= transcript.tx_start && sv_end >= transcript.tx_end;
    let both_breakpoints_within_tx = sv_start > transcript.tx_start && sv_end < transcript.tx_end;

    let consequence = match (sv_kind, spans_whole_transcript) {
        (SvKind::Deletion, true) => Consequence::TranscriptAblation,
        (SvKind::Deletion, false) => Consequence::FeatureTruncation,
        (SvKind::Duplication, true) => Consequence::TranscriptAmplification,
        (SvKind::Duplication, false) => Consequence::FeatureElongation,
    };

    let tx_span = transcript.tx_end.saturating_sub(transcript.tx_start);
    let covered = sv_end
        .min(transcript.tx_end)
        .saturating_sub(sv_start.max(transcript.tx_start));
    let transcript_overlap_fraction = if tx_span == 0 {
        0.0
    } else {
        covered as f64 / tx_span as f64
    };

    let overlap = compute_exon_overlap(sv_start, sv_end, transcript);

    CnvConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        hgnc_id: transcript.hgnc_id.clone(),
        strand: transcript.strand,
        consequence,
        impact: consequence.impact(),
        spans_whole_transcript,
        transcript_overlap_fraction,
        cds_bases_overlapped: overlap.cds_bases_overlapped,
        exons_overlapped: overlap.exons_overlapped,
        overlaps_first_exon: overlap.overlaps_first_exon,
        overlaps_last_exon: overlap.overlaps_last_exon,
        breakpoint_within_transcript: !spans_whole_transcript,
        both_breakpoints_within_tx,
    }
}

/// Annotate a structural-variant interval against every overlapping transcript.
///
/// Returns one [`CnvConsequenceResult`] per transcript whose `[tx_start,
/// tx_end)` intersects `[start, end)` on `chrom`. An **empty** vec means the
/// interval overlapped no transcript — this is a valid "no overlap" answer, not
/// an error. Errors are reserved for malformed input (`start >= end`) so a
/// store/coordinate fault can never silently masquerade as "no overlap".
///
/// # Arguments
///
/// * `chrom` — UCSC-style chromosome name (`"chr17"`).
/// * `start` — 0-based inclusive start of the SV footprint.
/// * `end` — 0-based exclusive end of the SV footprint.
/// * `sv_kind` — [`SvKind::Deletion`] or [`SvKind::Duplication`].
/// * `store` — transcript store for the overlap query.
///
/// # Errors
///
/// Returns [`VarEffectError::Malformed`] if `start >= end` (a zero-length or
/// inverted interval is a caller bug, not bad reference data).
pub(crate) fn annotate_interval(
    chrom: &str,
    start: u64,
    end: u64,
    sv_kind: SvKind,
    store: &TranscriptStore,
) -> Result<Vec<CnvConsequenceResult>, VarEffectError> {
    if start >= end {
        return Err(VarEffectError::Malformed(format!(
            "SV interval start ({start}) must be < end ({end})"
        )));
    }

    let overlaps: Vec<(&TranscriptModel, &LocateIndex)> = store.query_overlap(chrom, start, end);
    let mut results = Vec::with_capacity(overlaps.len());
    for (tx, _idx) in overlaps {
        results.push(annotate_transcript(start, end, sv_kind, tx));
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{
        minus_strand_coding, noncoding_2_exon, plus_strand_coding, single_exon_coding,
    };
    use crate::transcript::TranscriptStore;

    /// Build a store over the four shared fixtures (chr1 plus, chr17 minus,
    /// chr2 non-coding, chr3 single-exon).
    fn fixture_store() -> TranscriptStore {
        TranscriptStore::from_transcripts(vec![
            plus_strand_coding(),
            minus_strand_coding(),
            noncoding_2_exon(),
            single_exon_coding(),
        ])
    }

    /// A deletion fully spanning a transcript -> transcript_ablation (HIGH).
    #[test]
    fn deletion_spanning_whole_transcript_is_ablation() {
        let store = fixture_store();
        // plus_strand_coding spans chr1 [1000, 5000); cover it entirely.
        let results = annotate_interval("chr1", 500, 6000, SvKind::Deletion, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.consequence, Consequence::TranscriptAblation);
        assert_eq!(r.impact, Impact::High);
        assert!(r.spans_whole_transcript);
        assert!(!r.breakpoint_within_transcript);
        assert_eq!(r.transcript_overlap_fraction, 1.0);
        assert!(r.overlaps_first_exon && r.overlaps_last_exon);
        assert_eq!(r.exons_overlapped, 3);
        assert_eq!(r.cds_bases_overlapped, 1500);
    }

    /// A duplication fully spanning a transcript -> transcript_amplification.
    #[test]
    fn duplication_spanning_whole_transcript_is_amplification() {
        let store = fixture_store();
        let results = annotate_interval("chr1", 500, 6000, SvKind::Duplication, &store).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].consequence, Consequence::TranscriptAmplification);
        assert_eq!(results[0].impact, Impact::High);
        assert!(results[0].spans_whole_transcript);
    }

    /// A deletion clipping only the 5' end -> feature_truncation; first exon
    /// covered, last exon not.
    #[test]
    fn partial_deletion_is_truncation_with_geometry() {
        let store = fixture_store();
        // Cover exon 0 [1000,2000) and part of intron; stop before exon 2.
        let results = annotate_interval("chr1", 900, 2500, SvKind::Deletion, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.consequence, Consequence::FeatureTruncation);
        assert_eq!(r.impact, Impact::High);
        assert!(!r.spans_whole_transcript);
        assert!(r.breakpoint_within_transcript);
        assert!(r.overlaps_first_exon);
        assert!(!r.overlaps_last_exon);
        assert_eq!(r.exons_overlapped, 1);
        // CDS in exon 0 is [1500,2000) = 500 bases.
        assert_eq!(r.cds_bases_overlapped, 500);
        assert!(r.transcript_overlap_fraction > 0.0 && r.transcript_overlap_fraction < 1.0);
    }

    /// A partial duplication -> feature_elongation.
    #[test]
    fn partial_duplication_is_elongation() {
        let store = fixture_store();
        let results = annotate_interval("chr1", 3200, 4200, SvKind::Duplication, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.consequence, Consequence::FeatureElongation);
        assert!(!r.spans_whole_transcript);
        // Exons 1 [3000,3500) and 2 [4000,5000) both touched; exon 0 not.
        assert!(!r.overlaps_first_exon);
        assert!(r.overlaps_last_exon);
        assert_eq!(r.exons_overlapped, 2);
    }

    /// Minus-strand transcript: exon 0 is the most 5' exon at the *highest*
    /// genomic coordinates. A deletion clipping the high-coordinate end must
    /// flag `overlaps_first_exon`, not last.
    #[test]
    fn minus_strand_first_exon_is_high_coordinate() {
        let store = fixture_store();
        // minus_strand_coding: exon 0 (5') = [18000,20000), exon 2 (3') =
        // [10000,12000). Clip just the 5' (high-coord) exon.
        let results = annotate_interval("chr17", 18500, 21000, SvKind::Deletion, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.strand, Strand::Minus);
        assert_eq!(r.consequence, Consequence::FeatureTruncation);
        assert!(r.overlaps_first_exon);
        assert!(!r.overlaps_last_exon);
    }

    /// A non-coding transcript still annotates, with `cds_bases_overlapped == 0`.
    #[test]
    fn noncoding_transcript_has_zero_cds_overlap() {
        let store = fixture_store();
        // noncoding_2_exon spans chr2 [100,600); cover it fully.
        let results = annotate_interval("chr2", 50, 700, SvKind::Deletion, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.consequence, Consequence::TranscriptAblation);
        assert!(r.spans_whole_transcript);
        assert_eq!(r.cds_bases_overlapped, 0);
        assert_eq!(r.exons_overlapped, 2);
    }

    /// No overlapping transcript -> empty vec, NOT an error.
    #[test]
    fn no_overlap_returns_empty_not_error() {
        let store = fixture_store();
        let results = annotate_interval("chr1", 100, 200, SvKind::Deletion, &store).unwrap();
        assert!(results.is_empty());
    }

    /// Unknown chromosome -> empty vec (mirrors `query_overlap`), not an error.
    #[test]
    fn unknown_chrom_returns_empty_not_error() {
        let store = fixture_store();
        let results = annotate_interval("chrZZ", 100, 5000, SvKind::Deletion, &store).unwrap();
        assert!(results.is_empty());
    }

    /// A fully-spanning deletion is not intragenic: `both_breakpoints_within_tx`
    /// is `false` (the SV reaches/passes both transcript bounds).
    #[test]
    fn full_span_is_not_both_breakpoints_within() {
        let store = fixture_store();
        let results = annotate_interval("chr1", 500, 6000, SvKind::Deletion, &store).unwrap();
        assert!(!results[0].both_breakpoints_within_tx);
    }

    /// An intragenic deletion confined to exon 1 + intron 1 on the plus-strand
    /// fixture (both breakpoints strictly inside `[1000, 5000)`) is genuinely
    /// intragenic even though it touches no terminal exon. Regression for the
    /// exact-bound predicate vs. the exon-flag proxy.
    #[test]
    fn intragenic_deletion_touching_no_terminal_exon() {
        let store = fixture_store();
        // [3200, 3800): inside exon 1 + intron 1, no first/last exon.
        let r = &annotate_interval("chr1", 3200, 3800, SvKind::Deletion, &store).unwrap()[0];
        assert!(!r.spans_whole_transcript);
        assert!(!r.overlaps_first_exon && !r.overlaps_last_exon);
        assert!(
            r.both_breakpoints_within_tx,
            "both breakpoints strictly inside the transcript"
        );
    }

    /// An intragenic deletion that touches the first exon but stays strictly
    /// inside the transcript span. The exon-flag proxy would mis-read this as
    /// 5'-partial (`overlaps_first_exon == true`); the exact bound test keeps it
    /// intragenic.
    #[test]
    fn intragenic_deletion_touching_first_exon_is_still_within() {
        let store = fixture_store();
        // [1500, 3200): touches exon 0 but starts after tx_start (1000) and ends
        // before tx_end (5000) -> intragenic.
        let r = &annotate_interval("chr1", 1500, 3200, SvKind::Deletion, &store).unwrap()[0];
        assert!(r.overlaps_first_exon);
        assert!(!r.overlaps_last_exon);
        assert!(
            r.both_breakpoints_within_tx,
            "first-exon-touching but bounded on both sides -> intragenic"
        );
    }

    /// Single-exon gene: the lone exon is both first and last, so the exon flags
    /// can never distinguish intragenic from terminal overlap. A partial,
    /// strictly-internal deletion must still be intragenic via the bound test.
    #[test]
    fn single_exon_intragenic_deletion() {
        let store = fixture_store();
        // [3000, 4000): strictly inside the single exon [2000, 5000).
        let r = &annotate_interval("chr3", 3000, 4000, SvKind::Deletion, &store).unwrap()[0];
        assert!(!r.spans_whole_transcript);
        assert!(
            r.overlaps_first_exon && r.overlaps_last_exon,
            "single exon is both first and last"
        );
        assert_eq!(r.exons_overlapped, 1);
        assert!(
            r.both_breakpoints_within_tx,
            "single-exon intragenic deletion is contained in the transcript"
        );
    }

    /// Single-exon gene fully spanned: not intragenic.
    #[test]
    fn single_exon_full_span_is_not_within() {
        let store = fixture_store();
        let r = &annotate_interval("chr3", 1000, 6000, SvKind::Deletion, &store).unwrap()[0];
        assert!(r.spans_whole_transcript);
        assert!(!r.both_breakpoints_within_tx);
    }

    /// An inverted/zero-length interval is a caller bug -> Malformed error,
    /// never a silent "no overlap".
    #[test]
    fn inverted_interval_errors() {
        let store = fixture_store();
        assert!(matches!(
            annotate_interval("chr1", 5000, 5000, SvKind::Deletion, &store),
            Err(VarEffectError::Malformed(_))
        ));
        assert!(matches!(
            annotate_interval("chr1", 5000, 1000, SvKind::Deletion, &store),
            Err(VarEffectError::Malformed(_))
        ));
    }
}
