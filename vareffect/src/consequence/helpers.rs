//! Internal helper functions shared across consequence submodules.

use super::{Consequence, ConsequenceResult, Impact, PtcStatus};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::locate::LocateIndex;
use crate::types::{Strand, TranscriptModel, TranscriptTier, UtrNumberingBreaks};

/// Fetch the 3-base reference codon from FASTA for a given CDS codon
/// start offset. Uses the batched [`fetch_cds_sequence`] internally so
/// a codon within a single CDS segment costs one FASTA seek instead of
/// three. For split codons spanning two CDS segments this costs two seeks.
pub(super) fn fetch_ref_codon(
    codon_start_offset: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<[u8; 3], VarEffectError> {
    let seq = fetch_cds_sequence(
        codon_start_offset,
        codon_start_offset + 3,
        chrom,
        transcript,
        index,
        fasta,
    )?;
    debug_assert_eq!(seq.len(), 3, "fetch_ref_codon: expected 3 bases");
    Ok([seq[0], seq[1], seq[2]])
}

/// Build the alternate codon by replacing the mutated position with the
/// coding-strand alt base.
pub(super) fn build_alt_codon(ref_codon: &[u8; 3], codon_position: u8, coding_alt: u8) -> [u8; 3] {
    let mut alt = *ref_codon;
    alt[codon_position as usize] = coding_alt;
    alt
}

/// Compute the cDNA position (1-based from transcript 5' end, includes UTR)
/// for a CDS variant at the given CDS offset.
///
/// Uses the precomputed [`LocateIndex::utr5_exonic_len`] for O(1) lookup.
/// cDNA counts transcript bases only, so the CDS pad of a 5'-incomplete CDS
/// is removed (VEP 116: ENST00000379099.3 c.4 is cDNA 2). An offset inside
/// the pad has no transcript base and maps to the first coding base.
pub(super) fn compute_cdna_position_for_cds(
    cds_offset: u32,
    index: &crate::locate::LocateIndex,
) -> u32 {
    index.utr5_exonic_len() + cds_offset.saturating_sub(index.cds_start_pad()) + 1
}

/// Compute the cDNA position (1-based from transcript 5' end) for any exonic
/// genomic position, including UTRs and non-coding exons.
///
/// Walks exons in transcript order (5' -> 3'), summing exonic bases until the
/// exon containing `pos` is reached, then adds the intra-exon offset.
///
/// Returns `None` if `pos` does not fall within any exon.
pub(crate) fn compute_cdna_position_exonic(pos: u64, transcript: &TranscriptModel) -> Option<u32> {
    let mut cdna = 0u32;
    for exon in &transcript.exons {
        if pos >= exon.genomic_start && pos < exon.genomic_end {
            let intra = match transcript.strand {
                Strand::Plus => (pos - exon.genomic_start) as u32,
                Strand::Minus => (exon.genomic_end - 1 - pos) as u32,
            };
            return Some(cdna + intra + 1);
        }
        cdna += (exon.genomic_end - exon.genomic_start) as u32;
    }
    None
}

/// Sort consequences by severity (most severe first) and compute the
/// maximum impact.
pub(super) fn finalize_consequences(consequences: &mut [Consequence]) -> Impact {
    consequences.sort_by_key(|c| c.severity_rank());
    consequences
        .iter()
        .map(|c| c.impact())
        .max()
        .unwrap_or(Impact::Modifier)
}

/// Build the common (non-coding) fields of a [`ConsequenceResult`].
pub(super) fn build_base_result(
    transcript: &TranscriptModel,
    mut consequences: Vec<Consequence>,
) -> ConsequenceResult {
    let impact = finalize_consequences(&mut consequences);
    ConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        protein_accession: transcript.protein_accession.clone(),
        consequences,
        impact,
        protein_start: None,
        protein_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        distance: None,
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
        is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
        is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
        hgvs_c: None,
        hgvs_p: None,
        predicts_nmd: false,
        ptc: PtcStatus::NotApplicable,
    }
}

/// Strip shared prefix (left-first) then shared suffix from REF/ALT alleles.
///
/// Returns `(trimmed_ref, trimmed_alt, pos_adjustment)` where
/// `pos_adjustment` is the number of prefix bases stripped. Add to the
/// original 0-based position to get the trimmed variant's start.
///
/// One (but not both) of the returned slices may be empty for pure
/// insertions / deletions. Both-empty means REF == ALT (no-op variant).
///
/// Left-first matches VEP's `--minimal` behavior.
pub(super) fn trim_alleles<'a>(
    ref_allele: &'a [u8],
    alt_allele: &'a [u8],
) -> (&'a [u8], &'a [u8], u64) {
    let mut r = ref_allele;
    let mut a = alt_allele;
    let mut pos_adj: u64 = 0;

    while !r.is_empty() && !a.is_empty() && r[0] == a[0] {
        r = &r[1..];
        a = &a[1..];
        pos_adj += 1;
    }

    while !r.is_empty() && !a.is_empty() && r[r.len() - 1] == a[a.len() - 1] {
        r = &r[..r.len() - 1];
        a = &a[..a.len() - 1];
    }

    (r, a, pos_adj)
}

/// Fetch coding-strand bases for CDS offsets `[start_offset, end_offset)`.
///
/// Batches consecutive offsets within the same CDS segment into a single
/// [`FastaReader::fetch_sequence`] call, reducing the number of mutex
/// acquires and file seeks from O(n) per base to O(m) per CDS segment
/// spanned. For minus-strand transcripts, each chunk is reverse-
/// complemented to yield coding-strand order.
pub(crate) fn fetch_cds_sequence(
    start_offset: u32,
    end_offset: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<Vec<u8>, VarEffectError> {
    let cumulative = index.cumulative_cds();
    if end_offset > index.total_cds_length() || start_offset > end_offset {
        return Err(VarEffectError::Malformed(format!(
            "{}: CDS offsets [{start_offset}, {end_offset}) outside the CDS of length {}",
            transcript.accession,
            index.total_cds_length(),
        )));
    }
    let len = (end_offset - start_offset) as usize;
    let mut seq = Vec::with_capacity(len);

    // Offsets inside the pad of a 5'-incomplete CDS have no genomic base;
    // VEP writes them as `N` (codons `nnT`), so the partial codon translates
    // to `X`.
    let pad = index.cds_start_pad();
    if start_offset < pad {
        let pad_end = end_offset.min(pad);
        seq.resize((pad_end - start_offset) as usize, b'N');
        if end_offset <= pad {
            return Ok(seq);
        }
    }
    let start_offset = start_offset.max(pad);

    // Find first and last CDS segments spanned by [start_offset, end_offset).
    let first_seg = cumulative.partition_point(|&c| c <= start_offset) - 1;
    let last_seg = cumulative.partition_point(|&c| c < end_offset) - 1;

    for seg_idx in first_seg..=last_seg {
        let seg = &transcript.cds_segments[seg_idx];
        let seg_cds_start = cumulative[seg_idx];
        let seg_cds_end = cumulative[seg_idx + 1];

        // Clamp to the requested range within this segment.
        let local_start = (start_offset.max(seg_cds_start) - seg_cds_start) as u64;
        let local_end = (end_offset.min(seg_cds_end) - seg_cds_start) as u64;

        let (gstart, gend) = match transcript.strand {
            Strand::Plus => (
                seg.genomic_start + local_start,
                seg.genomic_start + local_end,
            ),
            Strand::Minus => (seg.genomic_end - local_end, seg.genomic_end - local_start),
        };

        let chunk = fasta.fetch_sequence(chrom, gstart, gend)?;
        match transcript.strand {
            Strand::Plus => seq.extend_from_slice(&chunk),
            Strand::Minus => {
                // Plus-strand ascending genomic bytes → reverse complement
                // for coding strand in 5'→3' transcript order.
                let rc = crate::codon::reverse_complement(&chunk);
                seq.extend_from_slice(&rc);
            }
        }
    }
    Ok(seq)
}

/// Add the canonical splice terms that spill out of a one-base frameshift
/// intron (see [`crate::locate::helpers::exonic_splice_spill`]) and
/// recompute the impact.
pub(super) fn add_splice_spill(
    result: &mut ConsequenceResult,
    spill: Option<crate::locate::helpers::SpliceSpill>,
) {
    let Some(spill) = spill else {
        return;
    };
    for (hit, term) in [
        (spill.donor, Consequence::SpliceDonorVariant),
        (spill.acceptor, Consequence::SpliceAcceptorVariant),
    ] {
        if hit && !result.consequences.contains(&term) {
            result.consequences.push(term);
        }
    }
    // A base in a canonical window is not also a splice-region base (VEP 116
    // c.228C>A: `splice_acceptor_variant,synonymous_variant`).
    result
        .consequences
        .retain(|c| *c != Consequence::SpliceRegionVariant);
    result.impact = finalize_consequences(&mut result.consequences);
}

/// Clear everything numbered along the transcript on a transcript flagged
/// [`TranscriptModel::cds_alignment_gap`]: HGVS c. and p., cDNA/CDS/protein
/// positions, codons and amino acids; `ptc` becomes
/// [`PtcStatus::Indeterminate`] and `predicts_nmd` false. Terms, impact and
/// exon/intron numbers are kept (`VEP_DIVERGENCES.md` 3.3).
///
/// Past a [`TranscriptModel::utr_numbering_breaks`] bound only cDNA positions,
/// and HGVS c. or n. naming such a position, are cleared.
pub(super) fn without_gapped_numbering(
    transcript: &TranscriptModel,
    mut result: ConsequenceResult,
) -> ConsequenceResult {
    if transcript.cds_alignment_gap {
        result.hgvs_c = None;
        result.hgvs_p = None;
        result.cdna_position = None;
        result.cdna_position_end = None;
        result.cds_position = None;
        result.cds_position_end = None;
        result.protein_start = None;
        result.protein_end = None;
        result.codons = None;
        result.amino_acids = None;
        if result.ptc != PtcStatus::NotApplicable {
            result.ptc = PtcStatus::Indeterminate;
        }
        result.predicts_nmd = false;
    } else {
        let breaks = &transcript.utr_numbering_breaks;
        let past_cdna = |position: Option<u32>| {
            position
                .zip(breaks.cdna_from)
                .is_some_and(|(p, from)| p >= from)
        };
        if past_cdna(result.cdna_position) || past_cdna(result.cdna_position_end) {
            result.cdna_position = None;
            result.cdna_position_end = None;
        }
        if result
            .hgvs_c
            .as_deref()
            .is_some_and(|hgvs| names_broken_position(hgvs, breaks))
        {
            result.hgvs_c = None;
        }
    }
    result
}

/// `true` when HGVS c. or n. names a position at or past one of `breaks`:
/// `n.N` from `cdna_from`, `c.-N` (and `c.-N+M`) from `five_prime_from`,
/// `c.*N` from `three_prime_from`, at either end of a range. An unreadable
/// number counts as past.
///
/// # Arguments
///
/// * `hgvs` - HGVS c. or n., with or without the accession prefix
/// * `breaks` - The transcript's UTR numbering breaks
pub(super) fn names_broken_position(hgvs: &str, breaks: &UtrNumberingBreaks) -> bool {
    if !breaks.any() {
        return false;
    }
    let notation = hgvs.rsplit_once(':').map_or(hgvs, |(_, rest)| rest);
    let (coding, positions) = match (notation.strip_prefix("c."), notation.strip_prefix("n.")) {
        (Some(positions), _) => (true, positions),
        (None, Some(positions)) => (false, positions),
        (None, None) => return false,
    };
    positions
        .split(|c: char| !matches!(c, '0'..='9' | '+' | '-' | '*' | '_'))
        .next()
        .unwrap_or_default()
        .split('_')
        .any(|position| {
            let (bound, number) = match (coding, position.as_bytes().first()) {
                (true, Some(b'-')) => (breaks.five_prime_from, &position[1..]),
                (true, Some(b'*')) => (breaks.three_prime_from, &position[1..]),
                (false, Some(b'0'..=b'9')) => (breaks.cdna_from, position),
                _ => (None, position),
            };
            bound.is_some_and(|from| {
                number
                    .split(['+', '-'])
                    .next()
                    .and_then(|n| n.parse::<u32>().ok())
                    .is_none_or(|n| n >= from)
            })
        })
}

/// `true` when an alternate allele has a base other than A, C, G or T
/// (case-insensitive), typically an `N` from an unresolved stretch of an
/// inserted sequence.
pub(super) fn has_unknown_bases(alt: &[u8]) -> bool {
    alt.iter()
        .any(|base| !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// The site of an alternate allele, for [`without_unknown_alt_effects`].
pub(super) struct AltSite<'a> {
    /// The trimmed alternate allele.
    pub(super) alt: &'a [u8],
    /// `true` for an insertion or a length-changing delins, `false` for an
    /// SNV or MNV.
    pub(super) is_indel: bool,
    pub(super) chrom: &'a str,
    pub(super) transcript: &'a TranscriptModel,
    pub(super) index: &'a LocateIndex,
    pub(super) fasta: &'a FastaReader,
}

/// Apply VEP 116's handling of an alternate allele with unknown bases (an
/// `N`): it has no peptide, so no HGVS c., HGVS p. or amino acids, and every
/// peptide-derived term gives way to `coding_sequence_variant`
/// (`coding_unknown`). Codons and positions stay (NM_000518.5 `aCT/aAN`,
/// NM_003051.4 `Tct/Nct`). VEP's ambiguity-aware indel rules are kept
/// (NM_000546.6): an insertion or delins at the start codon is
/// `start_retained_variant` (plus `frameshift_variant` for a non-multiple of
/// 3), one overlapping the stop codon `stop_retained_variant` alone.
/// `coding_sequence_variant` is dropped next to `frameshift_variant` or
/// `start_retained_variant`. A frameshift's termination codon is
/// undetermined; otherwise there is none to report.
///
/// # Arguments
///
/// * `site` - The allele and the transcript it was annotated on
/// * `result` - The annotation computed from it
///
/// # Returns
///
/// `result` unchanged when the allele is all A, C, G or T.
///
/// # Errors
///
/// Returns an error if the stop codon cannot be read from the reference.
pub(super) fn without_unknown_alt_effects(
    site: &AltSite<'_>,
    mut result: ConsequenceResult,
) -> Result<ConsequenceResult, VarEffectError> {
    if !has_unknown_bases(site.alt) {
        return Ok(result);
    }
    result.hgvs_c = None;
    result.hgvs_p = None;
    result.amino_acids = None;
    let at_start = result.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::StartLost | Consequence::StartRetainedVariant
        )
    });
    let at_stop = site.is_indel && overlaps_stop_codon(site, &result)?;
    let mut consequences: Vec<Consequence> = Vec::with_capacity(result.consequences.len());
    for term in &result.consequences {
        let term = match term {
            Consequence::StartLost | Consequence::StartRetainedVariant if site.is_indel => {
                Consequence::StartRetainedVariant
            }
            Consequence::FrameshiftVariant if at_stop => Consequence::StopRetainedVariant,
            Consequence::MissenseVariant
            | Consequence::SynonymousVariant
            | Consequence::StopGained
            | Consequence::StopLost
            | Consequence::StopRetainedVariant
            | Consequence::StartLost
            | Consequence::StartRetainedVariant
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::ProteinAlteringVariant
            | Consequence::CodingSequenceVariant => {
                if at_stop {
                    Consequence::StopRetainedVariant
                } else {
                    Consequence::CodingSequenceVariant
                }
            }
            Consequence::TranscriptAblation
            | Consequence::TranscriptAmplification
            | Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::FrameshiftVariant
            | Consequence::FeatureTruncation
            | Consequence::FeatureElongation
            | Consequence::SpliceRegionVariant
            | Consequence::IncompleteTerminalCodonVariant
            | Consequence::FivePrimeUtrVariant
            | Consequence::ThreePrimeUtrVariant
            | Consequence::NonCodingTranscriptExonVariant
            | Consequence::IntronVariant
            | Consequence::NonCodingTranscriptVariant
            | Consequence::UpstreamGeneVariant
            | Consequence::DownstreamGeneVariant
            | Consequence::IntergenicVariant => *term,
        };
        if !consequences.contains(&term) {
            consequences.push(term);
        }
    }
    if site.is_indel && at_start && !consequences.contains(&Consequence::StartRetainedVariant) {
        consequences.push(Consequence::StartRetainedVariant);
    }
    if consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::FrameshiftVariant | Consequence::StartRetainedVariant
        )
    }) {
        consequences.retain(|c| *c != Consequence::CodingSequenceVariant);
    }
    result.impact = finalize_consequences(&mut consequences);
    let frameshift = consequences.contains(&Consequence::FrameshiftVariant);
    result.consequences = consequences;
    if frameshift {
        result.ptc = PtcStatus::Indeterminate;
    } else {
        result.ptc = PtcStatus::NotApplicable;
        result.predicts_nmd = false;
    }
    Ok(result)
}

/// `true` when the variant's CDS span (for an insertion, its two flanking
/// bases) touches the stop codon of a CDS that ends in one.
fn overlaps_stop_codon(
    site: &AltSite<'_>,
    result: &ConsequenceResult,
) -> Result<bool, VarEffectError> {
    let Some(first) = result.cds_position else {
        return Ok(false);
    };
    let last = result.cds_position_end.unwrap_or(first);
    let total = site.index.total_cds_length();
    if last + 2 < total || first > total {
        return Ok(false);
    }
    ends_in_stop_codon(site.chrom, site.transcript, site.index, site.fasta)
}

/// `true` when codon 1 is a complete start codon, the only case in which
/// VEP's start-codon predicates (`start_lost`, `start_retained_variant`)
/// apply. A 5'-incomplete CDS (`cds_start_NF`) has no start codon, even
/// when its first codon is whole and reads ATG:
/// `VariationEffect.pm::_overlaps_start_codon` returns 0 for it
/// (ensembl-variation release/115).
pub(super) fn has_complete_start_codon(transcript: &TranscriptModel, index: &LocateIndex) -> bool {
    transcript.cds_start_complete && index.cds_start_pad() == 0
}

/// `true` when the CDS ends in a stop codon on the genome: its length is a
/// whole number of codons and the last one translates to `*`.
///
/// False for a 3'-incomplete CDS (`cds_end_NF`) and for a stop that only
/// polyadenylation completes (NM_001424184.1, MT-ND4); VEP's
/// `_overlaps_stop_codon` returns 0 for `cds_end_NF` (`VEP_DIVERGENCES.md`
/// 3.1).
///
/// # Errors
///
/// Returns an error if the last codon cannot be read from the reference.
pub(super) fn ends_in_stop_codon(
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<bool, VarEffectError> {
    let total = index.total_cds_length();
    if !transcript.cds_end_complete
        || total < 3
        || !total.is_multiple_of(3)
        || total - 3 < index.cds_start_pad()
    {
        return Ok(false);
    }
    let codon = fetch_ref_codon(total - 3, chrom, transcript, index, fasta)?;
    let is_mito = transcript.chrom == "chrM";
    Ok(crate::codon::translate_codon_for_transcript(&codon, is_mito) == b'*')
}

/// End of the last whole codon of the CDS, in VEP's padded frame.
///
/// Equals [`LocateIndex::total_cds_length`] when the CDS is a whole number
/// of codons; otherwise it stops before the 1 or 2 bases of the incomplete
/// terminal codon. A window clamped to it always translates.
///
/// # Arguments
///
/// * `index` - Locate index of the transcript
///
/// # Returns
///
/// The padded CDS offset one past the last whole codon (0 for a
/// non-coding transcript).
pub(super) fn complete_codons_end(index: &LocateIndex) -> u32 {
    let total = index.total_cds_length();
    total - total % 3
}

/// Check whether a CDS offset falls in the incomplete terminal codon of a
/// transcript whose total CDS length is not divisible by 3.
///
/// Returns `true` if the CDS has a partial trailing codon (1 or 2 bases)
/// and the given offset falls within it. Very rare in MANE transcripts.
pub(super) fn is_incomplete_terminal_codon(cds_offset: u32, index: &LocateIndex) -> bool {
    let total = index.total_cds_length();
    let remainder = total % 3;
    remainder != 0 && cds_offset >= total - remainder
}
