//! Internal helper functions shared across consequence submodules.

use super::{Consequence, ConsequenceResult, Impact};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::locate::LocateIndex;
use crate::types::{Strand, TranscriptModel, TranscriptTier};

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
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
        is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
        is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
        hgvs_c: None,
        hgvs_p: None,
        predicts_nmd: false,
        ptc: super::PtcStatus::NotApplicable,
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

/// `true` when codon 1 is a complete start codon, the only case in which
/// VEP's start-codon predicates (`start_lost`, `start_retained_variant`)
/// apply. A 5'-incomplete CDS (`cds_start_NF`) has no start codon:
/// `VariationEffect.pm::_overlaps_start_codon` returns 0 for it
/// (ensembl-variation release/115).
pub(super) fn has_complete_start_codon(index: &LocateIndex) -> bool {
    index.cds_start_pad() == 0
}

/// `true` when the CDS ends in a stop codon on the genome: its length is a
/// whole number of codons and the last one translates to `*`.
///
/// False for a 3'-incomplete CDS and for a stop that only polyadenylation
/// completes (NM_001424184.1, MT-ND4). VEP 116 gives a deletion from such a
/// partial codon into the 3'UTR no stop-codon term (NM_001424184.1
/// c.658_*2del: `incomplete_terminal_codon_variant`,
/// `coding_sequence_variant`, `3_prime_UTR_variant`), and
/// `VariationEffect.pm::_overlaps_stop_codon` returns 0 for `cds_end_NF`.
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
    if total < 3 || !total.is_multiple_of(3) || total - 3 < index.cds_start_pad() {
        return Ok(false);
    }
    let codon = fetch_ref_codon(total - 3, chrom, transcript, index, fasta)?;
    let is_mito = transcript.chrom == "chrM";
    Ok(crate::codon::translate_codon_for_transcript(&codon, is_mito) == b'*')
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
