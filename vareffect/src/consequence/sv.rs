//! Structural-variant consequence annotation (DEL / DUP / INV / INS / BND).
//!
//! Where [`super::annotate`] consumes a single ref/alt allele pair, structural
//! variants are described as a genomic *interval* plus a direction
//! ([`SvKind`]), or as a single breakpoint (breakends and symbolic insertions).
//! This module annotates such events against every overlapping transcript and
//! assigns the per-transcript Sequence Ontology term **set** Ensembl VEP emits
//! for structural variants (verified against the VEP REST API, GRCh38, RefSeq
//! transcripts):
//!
//! - **Sub-region terms** — which regions of the transcript the footprint
//!   overlaps: `coding_sequence_variant`, `5_prime_UTR_variant`,
//!   `3_prime_UTR_variant`, `intron_variant`, `non_coding_transcript_exon_variant`.
//! - **A copy-number headline term**, added only for unbalanced events:
//!   deletions → `transcript_ablation` (whole-transcript) or `feature_truncation`
//!   (partial); duplications → `transcript_amplification` / `feature_elongation`;
//!   a breakend whose breakpoint falls inside a transcript → `feature_truncation`.
//!   Inversions and symbolic insertions are balanced/positional and carry **no**
//!   headline term — VEP reports only the overlapped sub-regions (this is the
//!   resolution of Ensembl/ensembl-vep#79: an inversion is not a deletion).
//!
//! The emitted [`SvConsequenceResult`] also carries per-transcript exon-overlap
//! geometry (5'/3'/last-exon coverage, overlapped exon count, CDS bases
//! affected) so downstream callers can drive region-centric ACMG/ClinGen CNV
//! scoring (Riggs 2020, PMID 31690835) without re-deriving exon math. That
//! geometry is interval-oriented and most meaningful for DEL/DUP; for the
//! single-breakpoint events (BND/INS) it describes the 1 bp footprint.

use super::complex::compute_cds_projected_length;
use super::{Consequence, Impact};
use crate::error::VarEffectError;
use crate::locate::LocateIndex;
use crate::transcript::TranscriptStore;
use crate::types::{Strand, TranscriptModel};

/// Direction / shape of a structural variant.
///
/// Riggs 2020 scores losses and gains with separate rubrics, and the SO
/// consequence term also depends on direction. [`SvKind`] covers the three
/// *interval* SV shapes handled by [`crate::VarEffect::annotate_interval`];
/// breakends and symbolic insertions are single-breakpoint events with their
/// own entry points ([`crate::VarEffect::annotate_breakend`],
/// [`crate::VarEffect::annotate_sv_insertion`]).
///
/// `#[non_exhaustive]` so future interval shapes can be added without a SemVer
/// break for downstream `match` sites.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SvKind {
    /// Copy-number loss (`<DEL>`, `g.A_Bdel`). A deletion of the interval.
    Deletion,
    /// Copy-number gain (`<DUP>`, `g.A_Bdup`). A duplication of the interval.
    Duplication,
    /// Balanced inversion (`<INV>`). Copy-number-neutral: no ablation /
    /// truncation headline, only the overlapped sub-region terms.
    Inversion,
}

/// Internal SV shape, including the single-breakpoint events that are not part
/// of the public interval [`SvKind`]. Drives the headline-term rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SvShape {
    Deletion,
    Duplication,
    Inversion,
    /// Breakend breakpoint (`<BND>`): truncates any transcript it falls inside.
    Breakend,
    /// Symbolic insertion (`<INS>`): positional, no headline term.
    Insertion,
}

impl From<SvKind> for SvShape {
    fn from(kind: SvKind) -> Self {
        match kind {
            SvKind::Deletion => SvShape::Deletion,
            SvKind::Duplication => SvShape::Duplication,
            SvKind::Inversion => SvShape::Inversion,
        }
    }
}

/// Per-transcript annotation of a structural variant.
///
/// One result is produced for each transcript whose `[tx_start, tx_end)`
/// overlaps the queried footprint. Coordinates and exon numbering follow the
/// crate-wide 0-based / transcript-5'->3' conventions.
#[derive(Debug, Clone, PartialEq)]
pub struct SvConsequenceResult {
    /// RefSeq transcript accession with version (e.g. `"NM_006772.2"`).
    pub transcript: String,
    /// HGNC gene symbol.
    pub gene_symbol: String,
    /// HGNC identifier including the `HGNC:` prefix, if present on the model.
    pub hgnc_id: Option<String>,
    /// Transcript strand.
    pub strand: Strand,
    /// Transcript footprint start, `[tx_start, tx_end)` 0-based half-open.
    pub tx_start: u64,
    /// Transcript footprint end (exclusive).
    pub tx_end: u64,
    /// SO consequence terms VEP assigns to this transcript for the SV, ordered
    /// most-severe first (by IMPACT, then VEP's SO severity rank).
    pub consequences: Vec<Consequence>,
    /// Highest IMPACT among [`Self::consequences`].
    pub impact: Impact,
    /// Whether the SV footprint fully spans `[tx_start, tx_end)`.
    pub spans_whole_transcript: bool,
    /// Fraction of the transcript span `[tx_start, tx_end)` covered by the SV
    /// footprint, in `[0.0, 1.0]`.
    pub transcript_overlap_fraction: f64,
    /// Number of CDS bases of this transcript that fall inside the SV footprint.
    /// `0` for non-coding transcripts or intronic/UTR-only overlaps.
    pub cds_bases_overlapped: u32,
    /// Number of exons (of any kind) that the SV footprint overlaps.
    pub exons_overlapped: u16,
    /// `true` if the footprint overlaps the transcript's most 5' exon
    /// (exon 1 in transcript order).
    pub overlaps_first_exon: bool,
    /// `true` if the footprint overlaps the transcript's most 3' exon
    /// (the last exon in transcript order).
    pub overlaps_last_exon: bool,
    /// `true` if a breakpoint falls strictly inside the transcript span.
    /// Always `false` when [`Self::spans_whole_transcript`] is `true`.
    pub breakpoint_within_transcript: bool,
    /// `true` if **both** SV breakpoints fall strictly inside the transcript
    /// span — i.e. `sv_start > tx_start && sv_end < tx_end` (the SV is fully
    /// contained in `[tx_start, tx_end)`). Correct for single-exon transcripts
    /// and intragenic SVs touching a terminal exon, where the first/last-exon
    /// flags would mislead. Drives the Riggs 2020 "complete containment within
    /// an HI gene" tests (2E loss / 2I gain, PMID 31690835).
    pub both_breakpoints_within_tx: bool,
}

impl SvConsequenceResult {
    /// The most severe consequence term (first in the severity-ordered set),
    /// or `None` if no term was assigned.
    pub fn most_severe_consequence(&self) -> Option<Consequence> {
        self.consequences.first().copied()
    }
}

/// Exon-overlap geometry of an SV footprint against one transcript.
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

/// Append `c` to `out` if not already present (term sets are tiny, ≤ 5).
fn push_unique(out: &mut Vec<Consequence>, c: Consequence) {
    if !out.contains(&c) {
        out.push(c);
    }
}

/// Sub-region SO terms a footprint `[fs, fe)` overlaps within a transcript.
///
/// Coding transcripts split exonic overlap into CDS / 5'UTR / 3'UTR using the
/// genomic CDS bounds and strand (matching [`crate::locate::locate_variant`]);
/// non-coding transcripts emit `non_coding_transcript_exon_variant`. Any
/// intronic overlap adds `intron_variant`. The returned set is unordered.
fn sv_region_consequences(fs: u64, fe: u64, transcript: &TranscriptModel) -> Vec<Consequence> {
    let mut out: Vec<Consequence> = Vec::new();
    // A transcript is treated as coding only when it has CDS segments *and*
    // genomic CDS bounds; a malformed model missing the bounds degrades to
    // non-coding rather than panicking.
    let cds_bounds = match (transcript.cds_genomic_start, transcript.cds_genomic_end) {
        (Some(start), Some(end)) if !transcript.cds_segments.is_empty() => Some((start, end)),
        _ => None,
    };

    for exon in &transcript.exons {
        let lo = fs.max(exon.genomic_start);
        let hi = fe.min(exon.genomic_end);
        if lo >= hi {
            continue;
        }
        // Non-coding exon overlap. Split the coding case into CDS / UTR using
        // genomic CDS bounds; the UTR side follows strand (locate.rs: Plus
        // pos<cds_start = 5'UTR; Minus pos>=cds_end = 5'UTR).
        let Some((cds_start, cds_end)) = cds_bounds else {
            push_unique(&mut out, Consequence::NonCodingTranscriptExonVariant);
            continue;
        };

        if lo.max(cds_start) < hi.min(cds_end) {
            push_unique(&mut out, Consequence::CodingSequenceVariant);
        }
        if lo < hi.min(cds_start) {
            push_unique(
                &mut out,
                match transcript.strand {
                    Strand::Plus => Consequence::FivePrimeUtrVariant,
                    Strand::Minus => Consequence::ThreePrimeUtrVariant,
                },
            );
        }
        if lo.max(cds_end) < hi {
            push_unique(
                &mut out,
                match transcript.strand {
                    Strand::Plus => Consequence::ThreePrimeUtrVariant,
                    Strand::Minus => Consequence::FivePrimeUtrVariant,
                },
            );
        }
    }

    // Intronic overlap: any footprint base in a gap between genomically-sorted
    // exons. One `intron_variant` regardless of how many introns are touched.
    let mut spans: Vec<(u64, u64)> = transcript
        .exons
        .iter()
        .map(|e| (e.genomic_start, e.genomic_end))
        .collect();
    spans.sort_unstable_by_key(|&(start, _)| start);
    for pair in spans.windows(2) {
        let intron_start = pair[0].1;
        let intron_end = pair[1].0;
        if intron_start < intron_end && fs.max(intron_start) < fe.min(intron_end) {
            push_unique(&mut out, Consequence::IntronVariant);
            break;
        }
    }

    out
}

/// Order a consequence set most-severe first: IMPACT descending, then VEP's SO
/// severity rank ascending.
///
/// IMPACT dominates so `feature_truncation`/`feature_elongation` (HIGH since
/// Ensembl release 111) sort ahead of the MODIFIER sub-region terms, matching
/// VEP's `most_severe_consequence` selection.
fn order_by_severity(consequences: &mut [Consequence]) {
    consequences.sort_by(|a, b| {
        b.impact()
            .cmp(&a.impact())
            .then_with(|| a.severity_rank().cmp(&b.severity_rank()))
    });
}

/// Build the per-transcript result for one SV footprint against one transcript.
///
/// `[fs, fe)` is the 0-based half-open SV footprint (a single-base `[p, p+1)`
/// span for breakends and insertions). The transcript is assumed to overlap it.
fn build_sv_result(
    fs: u64,
    fe: u64,
    shape: SvShape,
    transcript: &TranscriptModel,
) -> SvConsequenceResult {
    let overlap = compute_exon_overlap(fs, fe, transcript);
    let spans_whole_transcript = fs <= transcript.tx_start && fe >= transcript.tx_end;
    let both_breakpoints_within_tx = fs > transcript.tx_start && fe < transcript.tx_end;
    // Contained = the footprint lies wholly inside `[tx_start, tx_end)`.
    let contained = fs >= transcript.tx_start && fe <= transcript.tx_end;
    let exonic = overlap.exons_overlapped > 0;

    // Headline-term rule, verified against the VEP REST API (RefSeq, GRCh38):
    // - whole-transcript DEL/DUP collapse to a single ablation/amplification term;
    // - `feature_truncation`: a non-whole deletion that overlaps an exon, or any
    //   breakend breakpoint (a break severs the transcript even from an intron);
    // - `feature_elongation`: a duplication/insertion contained within the
    //   transcript that overlaps an exon (an intronic-only or transcript-
    //   spanning event does not elongate the feature);
    // - inversions are balanced -> sub-region terms only (ensembl-vep#79).
    let mut consequences = match shape {
        SvShape::Deletion if spans_whole_transcript => vec![Consequence::TranscriptAblation],
        SvShape::Duplication if spans_whole_transcript => {
            vec![Consequence::TranscriptAmplification]
        }
        _ => {
            let mut c = sv_region_consequences(fs, fe, transcript);
            match shape {
                SvShape::Breakend => push_unique(&mut c, Consequence::FeatureTruncation),
                SvShape::Deletion if exonic => push_unique(&mut c, Consequence::FeatureTruncation),
                SvShape::Duplication | SvShape::Insertion if contained && exonic => {
                    push_unique(&mut c, Consequence::FeatureElongation)
                }
                _ => {}
            }
            c
        }
    };
    order_by_severity(&mut consequences);

    let impact = consequences
        .first()
        .map_or(Impact::Modifier, |c| c.impact());

    let tx_span = transcript.tx_end.saturating_sub(transcript.tx_start);
    let covered = fe
        .min(transcript.tx_end)
        .saturating_sub(fs.max(transcript.tx_start));
    let transcript_overlap_fraction = if tx_span == 0 {
        0.0
    } else {
        covered as f64 / tx_span as f64
    };

    SvConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        hgnc_id: transcript.hgnc_id.clone(),
        strand: transcript.strand,
        tx_start: transcript.tx_start,
        tx_end: transcript.tx_end,
        consequences,
        impact,
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

/// Annotate an interval structural variant (DEL / DUP / INV) against every
/// overlapping transcript.
///
/// Returns one [`SvConsequenceResult`] per transcript whose `[tx_start,
/// tx_end)` intersects `[start, end)` on `chrom`. An **empty** vec means the
/// interval overlapped no transcript — a valid "no overlap" answer, not an
/// error. Errors are reserved for malformed input (`start >= end`) so a
/// store/coordinate fault can never silently masquerade as "no overlap".
///
/// # Arguments
///
/// * `chrom` — UCSC-style chromosome name (`"chr17"`).
/// * `start` — 0-based inclusive start of the SV footprint.
/// * `end` — 0-based exclusive end of the SV footprint.
/// * `sv_kind` — interval SV shape ([`SvKind`]).
/// * `store` — transcript store for the overlap query.
///
/// # Errors
///
/// Returns [`VarEffectError::Malformed`] if `start >= end`.
pub(crate) fn annotate_interval(
    chrom: &str,
    start: u64,
    end: u64,
    sv_kind: SvKind,
    store: &TranscriptStore,
) -> Result<Vec<SvConsequenceResult>, VarEffectError> {
    if start >= end {
        return Err(VarEffectError::Malformed(format!(
            "SV interval start ({start}) must be < end ({end})"
        )));
    }

    let shape = SvShape::from(sv_kind);
    let overlaps: Vec<(&TranscriptModel, &LocateIndex)> = store.query_overlap(chrom, start, end);
    let mut results = Vec::with_capacity(overlaps.len());
    for (tx, _idx) in overlaps {
        results.push(build_sv_result(start, end, shape, tx));
    }
    Ok(results)
}

/// Annotate a single breakend breakpoint against every transcript it falls in.
///
/// A breakend at 0-based genomic `pos` is annotated as the 1 bp footprint
/// `[pos, pos + 1)`. A transcript the breakpoint falls inside is truncated by
/// the rearrangement → `feature_truncation`, alongside the sub-region term at
/// that base. A translocation is annotated by calling this once per mate
/// breakpoint (see [`crate::VarEffect::annotate_breakend_pair`]); mate pairing
/// (MATEID/EVENT) is intentionally not performed, matching VEP.
///
/// Returns one [`SvConsequenceResult`] per overlapping transcript; an empty vec
/// means the breakpoint fell in no transcript.
pub(crate) fn annotate_breakend(
    chrom: &str,
    pos: u64,
    store: &TranscriptStore,
) -> Vec<SvConsequenceResult> {
    let overlaps = store.query_overlap(chrom, pos, pos + 1);
    let mut results = Vec::with_capacity(overlaps.len());
    for (tx, _idx) in overlaps {
        results.push(build_sv_result(pos, pos + 1, SvShape::Breakend, tx));
    }
    results
}

/// Annotate a symbolic insertion (`<INS>`) at 0-based genomic `pos` against
/// every transcript it falls in.
///
/// A symbolic insertion carries no copy-number headline — VEP reports only the
/// sub-region the insertion point occupies (e.g. `intron_variant`,
/// `coding_sequence_variant`). The inserted sequence is unknown to the
/// annotator, so this is purely positional (the 1 bp footprint `[pos, pos + 1)`).
///
/// Returns one [`SvConsequenceResult`] per overlapping transcript.
pub(crate) fn annotate_sv_insertion(
    chrom: &str,
    pos: u64,
    store: &TranscriptStore,
) -> Vec<SvConsequenceResult> {
    let overlaps = store.query_overlap(chrom, pos, pos + 1);
    let mut results = Vec::with_capacity(overlaps.len());
    for (tx, _idx) in overlaps {
        results.push(build_sv_result(pos, pos + 1, SvShape::Insertion, tx));
    }
    results
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

    /// Sorted copy of a result's consequence terms for set comparison.
    fn terms(r: &SvConsequenceResult) -> Vec<&'static str> {
        let mut v: Vec<&'static str> = r.consequences.iter().map(|c| c.as_str()).collect();
        v.sort_unstable();
        v
    }

    // --- Deletion ----------------------------------------------------------

    /// A deletion fully spanning a transcript -> transcript_ablation only.
    #[test]
    fn deletion_spanning_whole_transcript_is_ablation_only() {
        let store = fixture_store();
        let results = annotate_interval("chr1", 500, 6000, SvKind::Deletion, &store).unwrap();
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.consequences, vec![Consequence::TranscriptAblation]);
        assert_eq!(r.impact, Impact::High);
        assert!(r.spans_whole_transcript);
        assert_eq!(r.transcript_overlap_fraction, 1.0);
        assert_eq!(r.cds_bases_overlapped, 1500);
    }

    /// A partial deletion clipping the plus-strand 5' end -> feature_truncation
    /// plus the overlapped sub-regions (5'UTR + CDS in exon 0, intron 0).
    #[test]
    fn partial_deletion_is_truncation_plus_subregions() {
        let store = fixture_store();
        // [900, 2500): exon 0 [1000,2000) (5'UTR [1000,1500)+CDS [1500,2000))
        // and intron 0 [2000,3000).
        let r = &annotate_interval("chr1", 900, 2500, SvKind::Deletion, &store).unwrap()[0];
        assert_eq!(
            terms(r),
            vec![
                "5_prime_UTR_variant",
                "coding_sequence_variant",
                "feature_truncation",
                "intron_variant",
            ]
        );
        assert_eq!(
            r.most_severe_consequence(),
            Some(Consequence::FeatureTruncation)
        );
        assert_eq!(r.impact, Impact::High);
        assert!(!r.spans_whole_transcript);
        assert!(r.overlaps_first_exon && !r.overlaps_last_exon);
    }

    // --- Duplication -------------------------------------------------------

    /// A duplication fully spanning a transcript -> transcript_amplification only.
    #[test]
    fn duplication_spanning_whole_transcript_is_amplification_only() {
        let store = fixture_store();
        let r = &annotate_interval("chr1", 500, 6000, SvKind::Duplication, &store).unwrap()[0];
        assert_eq!(r.consequences, vec![Consequence::TranscriptAmplification]);
        assert_eq!(r.impact, Impact::High);
    }

    /// A partial duplication -> feature_elongation plus sub-regions.
    #[test]
    fn partial_duplication_is_elongation_plus_subregions() {
        let store = fixture_store();
        // [3200, 4200): intron 1 [3500,4000), exon 1 [3000,3500) CDS,
        // exon 2 [4000,5000) CDS [4000,4500).
        let r = &annotate_interval("chr1", 3200, 4200, SvKind::Duplication, &store).unwrap()[0];
        assert!(r.consequences.contains(&Consequence::FeatureElongation));
        assert!(r.consequences.contains(&Consequence::CodingSequenceVariant));
        assert!(r.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(
            r.most_severe_consequence(),
            Some(Consequence::FeatureElongation)
        );
        assert!(!r.spans_whole_transcript);
    }

    // --- Inversion (balanced: NO ablation/truncation) ----------------------

    /// An inversion spanning the whole transcript lists every overlapped
    /// sub-region but carries NO ablation/truncation headline (ensembl-vep#79).
    #[test]
    fn inversion_whole_span_has_no_headline_term() {
        let store = fixture_store();
        let r = &annotate_interval("chr1", 500, 6000, SvKind::Inversion, &store).unwrap()[0];
        assert_eq!(
            terms(r),
            vec![
                "3_prime_UTR_variant",
                "5_prime_UTR_variant",
                "coding_sequence_variant",
                "intron_variant",
            ]
        );
        assert!(!r.consequences.contains(&Consequence::TranscriptAblation));
        assert!(!r.consequences.contains(&Consequence::FeatureTruncation));
        assert_eq!(r.impact, Impact::Modifier);
        assert_eq!(
            r.most_severe_consequence(),
            Some(Consequence::CodingSequenceVariant)
        );
    }

    /// A partial inversion clipping only the 5' end -> 5'UTR + CDS + intron,
    /// no headline.
    #[test]
    fn partial_inversion_is_subregions_only() {
        let store = fixture_store();
        let r = &annotate_interval("chr1", 900, 2500, SvKind::Inversion, &store).unwrap()[0];
        assert_eq!(
            terms(r),
            vec![
                "5_prime_UTR_variant",
                "coding_sequence_variant",
                "intron_variant"
            ]
        );
    }

    /// Minus-strand UTR assignment: a footprint clipping the high-coordinate 5'
    /// exon yields 5_prime_UTR_variant (not 3'), confirming strand-aware split.
    #[test]
    fn minus_strand_five_prime_utr_is_high_coordinate() {
        let store = fixture_store();
        // minus exon 0 [18000,20000): CDS [18000,19500) + 5'UTR [19500,20000).
        // [19600, 21000) hits only the 5'UTR portion (+ flank past tx_end).
        let r = &annotate_interval("chr17", 19600, 21000, SvKind::Inversion, &store).unwrap()[0];
        assert_eq!(r.strand, Strand::Minus);
        assert!(r.consequences.contains(&Consequence::FivePrimeUtrVariant));
        assert!(!r.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    // --- Non-coding --------------------------------------------------------

    /// A partial deletion of a non-coding transcript -> feature_truncation +
    /// non_coding_transcript_exon_variant (+ intron if touched).
    #[test]
    fn noncoding_partial_deletion() {
        let store = fixture_store();
        // noncoding_2_exon chr2: exon0 [100,300), intron [300,400), exon1 [400,600).
        let r = &annotate_interval("chr2", 250, 450, SvKind::Deletion, &store).unwrap()[0];
        assert_eq!(
            terms(r),
            vec![
                "feature_truncation",
                "intron_variant",
                "non_coding_transcript_exon_variant",
            ]
        );
        assert_eq!(r.cds_bases_overlapped, 0);
    }

    // --- Breakend ----------------------------------------------------------

    /// A breakend whose breakpoint falls in a CDS exon -> feature_truncation +
    /// coding_sequence_variant.
    #[test]
    fn breakend_in_cds_is_truncation_plus_cds() {
        let store = fixture_store();
        // pos 1700 is inside exon 0 CDS [1500,2000) of the plus fixture.
        let results = annotate_breakend("chr1", 1700, &store);
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(
            terms(r),
            vec!["coding_sequence_variant", "feature_truncation"]
        );
        assert_eq!(
            r.most_severe_consequence(),
            Some(Consequence::FeatureTruncation)
        );
        assert!(r.both_breakpoints_within_tx);
    }

    /// A breakend in an intron -> feature_truncation + intron_variant.
    #[test]
    fn breakend_in_intron_is_truncation_plus_intron() {
        let store = fixture_store();
        // pos 2500 is in intron 0 [2000,3000).
        let r = &annotate_breakend("chr1", 2500, &store)[0];
        assert_eq!(terms(r), vec!["feature_truncation", "intron_variant"]);
    }

    /// A breakend outside any transcript -> empty (no overlap).
    #[test]
    fn breakend_outside_transcript_is_empty() {
        let store = fixture_store();
        assert!(annotate_breakend("chr1", 100, &store).is_empty());
    }

    // --- Insertion ---------------------------------------------------------

    /// A symbolic insertion in an intron -> intron_variant only (no headline).
    #[test]
    fn insertion_in_intron_is_positional_only() {
        let store = fixture_store();
        let r = &annotate_sv_insertion("chr1", 2500, &store)[0];
        assert_eq!(terms(r), vec!["intron_variant"]);
        assert!(!r.consequences.contains(&Consequence::FeatureTruncation));
        assert_eq!(r.impact, Impact::Modifier);
    }

    /// A symbolic insertion in a CDS exon -> coding_sequence_variant +
    /// feature_elongation (a contained, exonic insertion lengthens the feature).
    #[test]
    fn insertion_in_cds_is_coding_plus_elongation() {
        let store = fixture_store();
        let r = &annotate_sv_insertion("chr1", 1700, &store)[0];
        assert_eq!(
            terms(r),
            vec!["coding_sequence_variant", "feature_elongation"]
        );
    }

    /// `feature_elongation` is exon-gated AND containment-gated: a duplication
    /// that extends beyond the transcript gets no headline, only sub-regions.
    #[test]
    fn partial_duplication_extending_out_has_no_elongation() {
        let store = fixture_store();
        // [3200, 6000): one breakpoint inside, one past tx_end (5000).
        let r = &annotate_interval("chr1", 3200, 6000, SvKind::Duplication, &store).unwrap()[0];
        assert!(!r.consequences.contains(&Consequence::FeatureElongation));
        assert!(!r.spans_whole_transcript);
        assert!(r.consequences.contains(&Consequence::CodingSequenceVariant));
    }

    /// `feature_truncation` is exon-gated: an intronic-only partial deletion
    /// removes no exonic sequence, so VEP emits only `intron_variant`.
    #[test]
    fn intronic_only_deletion_is_not_truncation() {
        let store = fixture_store();
        // [2200, 2800): wholly inside intron 0 [2000,3000).
        let r = &annotate_interval("chr1", 2200, 2800, SvKind::Deletion, &store).unwrap()[0];
        assert_eq!(terms(r), vec!["intron_variant"]);
        assert!(!r.consequences.contains(&Consequence::FeatureTruncation));
    }

    // --- Geometry / errors (carried over) ----------------------------------

    /// No overlapping transcript -> empty vec, NOT an error.
    #[test]
    fn no_overlap_returns_empty_not_error() {
        let store = fixture_store();
        assert!(
            annotate_interval("chr1", 100, 200, SvKind::Deletion, &store)
                .unwrap()
                .is_empty()
        );
    }

    /// Unknown chromosome -> empty vec, not an error.
    #[test]
    fn unknown_chrom_returns_empty_not_error() {
        let store = fixture_store();
        assert!(
            annotate_interval("chrZZ", 100, 5000, SvKind::Deletion, &store)
                .unwrap()
                .is_empty()
        );
    }

    /// A fully-spanning deletion is not intragenic.
    #[test]
    fn full_span_is_not_both_breakpoints_within() {
        let store = fixture_store();
        let r = &annotate_interval("chr1", 500, 6000, SvKind::Deletion, &store).unwrap()[0];
        assert!(!r.both_breakpoints_within_tx);
    }

    /// An intragenic deletion strictly inside the transcript, touching no
    /// terminal exon, is intragenic via the exact bound test.
    #[test]
    fn intragenic_deletion_touching_no_terminal_exon() {
        let store = fixture_store();
        let r = &annotate_interval("chr1", 3200, 3800, SvKind::Deletion, &store).unwrap()[0];
        assert!(!r.spans_whole_transcript);
        assert!(!r.overlaps_first_exon && !r.overlaps_last_exon);
        assert!(r.both_breakpoints_within_tx);
    }

    /// Single-exon gene fully spanned: not intragenic.
    #[test]
    fn single_exon_full_span_is_not_within() {
        let store = fixture_store();
        let r = &annotate_interval("chr3", 1000, 6000, SvKind::Deletion, &store).unwrap()[0];
        assert!(r.spans_whole_transcript);
        assert!(!r.both_breakpoints_within_tx);
        assert_eq!(r.consequences, vec![Consequence::TranscriptAblation]);
    }

    /// An inverted/zero-length interval is a caller bug -> Malformed error.
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
