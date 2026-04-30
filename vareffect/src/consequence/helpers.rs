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
pub(super) fn compute_cdna_position_for_cds(
    cds_offset: u32,
    index: &crate::locate::LocateIndex,
) -> u32 {
    index.utr5_exonic_len() + cds_offset + 1
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
    let len = (end_offset - start_offset) as usize;
    let mut seq = Vec::with_capacity(len);

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

        let chunk = fasta.fetch_sequence_slice(chrom, gstart, gend)?;
        match transcript.strand {
            Strand::Plus => seq.extend_from_slice(chunk),
            Strand::Minus => {
                // Plus-strand ascending genomic bytes → reverse complement
                // for coding strand in 5'→3' transcript order.
                let rc = crate::codon::reverse_complement(chunk);
                seq.extend_from_slice(&rc);
            }
        }
    }
    Ok(seq)
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
