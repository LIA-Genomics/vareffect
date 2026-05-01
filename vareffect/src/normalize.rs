//! HGVS 3' normalization for indels.
//!
//! Implements the HGVS 3' rule: indels in repeated sequences are shifted to
//! the most 3' (coding-strand rightward) equivalent position before HGVS
//! formatting. The shift respects region boundaries: exonic shifts stay
//! within the exon, intronic shifts stay within the intron.
//!
//! The 3' direction is relative to the **coding strand** (transcript
//! direction):
//! - Plus-strand genes: 3' = increasing genomic coordinates (rightward)
//! - Minus-strand genes: 3' = decreasing genomic coordinates (leftward)
//!
//! Reference: HGVS Nomenclature v21.1 (2024), 3' rule.
//! Concordance target: VEP `--shift_hgvs` (ON by default since VEP 109).

use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::types::{Strand, TranscriptModel};

/// Compute the 3' shift offset (in genomic bases) for a deletion.
///
/// Starting from the deletion range `[del_start, del_end)` on the reference,
/// checks how many bases the deletion can slide toward the 3' end of the
/// coding strand while remaining equivalent. A deletion `[s, e)` can shift
/// right by 1 iff `ref[s] == ref[e]` (the base entering the window from the
/// 3' side equals the base leaving from the 5' side).
///
/// Stops at region boundaries: exonic shifts stop at exon/intron
/// junctions; intronic shifts stop at intron/exon junctions.
///
/// # Returns
///
/// Number of bases to shift (0 if no shift is possible).
pub(crate) fn compute_3prime_shift_deletion(
    chrom: &str,
    del_start: u64,
    del_end: u64,
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<u64, VarEffectError> {
    let del_len = del_end - del_start;
    if del_len == 0 {
        return Ok(0);
    }

    match transcript.strand {
        Strand::Plus => {
            // 3' = increasing genomic coordinate.
            let fetch_end = if let Some((_, exon_end)) =
                containing_exon_range(del_end - 1, transcript)
            {
                exon_end
            } else if let Some((_, intron_end)) = containing_intron_range(del_end - 1, transcript) {
                intron_end
            } else {
                return Ok(0);
            };
            if fetch_end <= del_end {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence_slice(chrom, del_start, fetch_end)?;
            let dl = del_len as usize;
            let mut shift = 0usize;
            // seq[shift] = ref[del_start + shift]; seq[dl + shift] = ref[del_end + shift].
            while dl + shift < seq.len() {
                if seq[shift] != seq[dl + shift] {
                    break;
                }
                shift += 1;
            }
            Ok(shift as u64)
        }
        Strand::Minus => {
            // 3' = decreasing genomic coordinate.
            let fetch_start = if let Some((exon_start, _)) =
                containing_exon_range(del_start, transcript)
            {
                exon_start
            } else if let Some((intron_start, _)) = containing_intron_range(del_start, transcript) {
                intron_start
            } else {
                return Ok(0);
            };
            if fetch_start >= del_start {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence_slice(chrom, fetch_start, del_end)?;
            let ds = (del_start - fetch_start) as usize;
            let de = (del_end - fetch_start) as usize;
            let mut shift = 0usize;
            // Compare ref[del_start - 1 - shift] with ref[del_end - 1 - shift]
            // in seq coords: seq[ds - 1 - shift] vs seq[de - 1 - shift].
            while shift < ds {
                if seq[ds - 1 - shift] != seq[de - 1 - shift] {
                    break;
                }
                shift += 1;
            }
            Ok(shift as u64)
        }
    }
}

/// Compute the 3' shift offset (in genomic bases) for an insertion.
///
/// Starting from the insertion point `ins_pos` (insertion is between
/// `pos - 1` and `pos` on the genome), checks how many bases the
/// insertion can slide toward the 3' end of the coding strand while
/// remaining equivalent. At each step, the reference base at the shift
/// position must match the corresponding inserted base (cycling through
/// the inserted sequence).
///
/// Both `fasta` bytes and `inserted_bases` are in **plus-strand**
/// orientation (VCF convention). The shift direction changes per strand
/// but the base comparison is always plus-strand vs plus-strand.
///
/// Stops at region boundaries: exonic shifts stop at exon/intron
/// junctions; intronic shifts stop at intron/exon junctions.
///
/// # Returns
///
/// Number of bases to shift (0 if no shift is possible).
pub(crate) fn compute_3prime_shift_insertion(
    chrom: &str,
    ins_pos: u64,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<u64, VarEffectError> {
    let ins_len = inserted_bases.len();
    if ins_len == 0 {
        return Ok(0);
    }

    match transcript.strand {
        Strand::Plus => {
            // 3' = increasing genomic coordinate.
            let fetch_end = if let Some(exon_range) = find_shift_exon_range(ins_pos, transcript) {
                exon_range.1
            } else if let Some((_, intron_end)) = containing_intron_range(ins_pos, transcript) {
                intron_end
            } else if ins_pos > 0 {
                if let Some((_, intron_end)) = containing_intron_range(ins_pos - 1, transcript) {
                    intron_end
                } else {
                    return Ok(0);
                }
            } else {
                return Ok(0);
            };
            if fetch_end <= ins_pos {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence_slice(chrom, ins_pos, fetch_end)?;
            let mut shift = 0usize;
            while shift < seq.len() {
                if seq[shift] != inserted_bases[shift % ins_len] {
                    break;
                }
                shift += 1;
            }
            Ok(shift as u64)
        }
        Strand::Minus => {
            // 3' = decreasing genomic coordinate.
            let fetch_start = if let Some(exon_range) =
                find_shift_exon_range_minus(ins_pos, transcript)
            {
                exon_range.0
            } else if ins_pos > 0 {
                if let Some((intron_start, _)) = containing_intron_range(ins_pos - 1, transcript) {
                    intron_start
                } else if let Some((intron_start, _)) = containing_intron_range(ins_pos, transcript)
                {
                    intron_start
                } else {
                    return Ok(0);
                }
            } else {
                return Ok(0);
            };
            if fetch_start >= ins_pos {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence_slice(chrom, fetch_start, ins_pos)?;
            let mut shift = 0usize;
            // Walk seq backward from ins_pos - 1, comparing against inserted_bases
            // cycled from its tail.
            while shift < seq.len() {
                let ref_base = seq[seq.len() - 1 - shift];
                let ins_idx = ins_len - 1 - (shift % ins_len);
                if ref_base != inserted_bases[ins_idx] {
                    break;
                }
                shift += 1;
            }
            Ok(shift as u64)
        }
    }
}

/// Find the exon whose `[genomic_start, genomic_end)` range contains `pos`.
///
/// Returns `None` if `pos` is intronic or outside all exons.
fn containing_exon_range(pos: u64, transcript: &TranscriptModel) -> Option<(u64, u64)> {
    transcript
        .exons
        .iter()
        .find(|e| pos >= e.genomic_start && pos < e.genomic_end)
        .map(|e| (e.genomic_start, e.genomic_end))
}

/// Find the intron whose `[genomic_start, genomic_end)` range contains `pos`.
///
/// Returns `None` if `pos` is exonic, outside the transcript, or the
/// transcript has fewer than two exons.
///
/// Intron boundaries are derived from adjacent exon pairs:
/// - Plus strand (exons ascending): `[exons[i].genomic_end, exons[i+1].genomic_start)`
/// - Minus strand (exons descending): `[exons[i+1].genomic_end, exons[i].genomic_start)`
fn containing_intron_range(pos: u64, transcript: &TranscriptModel) -> Option<(u64, u64)> {
    for pair in transcript.exons.windows(2) {
        let (intron_start, intron_end) = match transcript.strand {
            Strand::Plus => (pair[0].genomic_end, pair[1].genomic_start),
            Strand::Minus => (pair[1].genomic_end, pair[0].genomic_start),
        };
        if pos >= intron_start && pos < intron_end {
            return Some((intron_start, intron_end));
        }
    }
    None
}

/// Find the exon range for a plus-strand insertion shift.
///
/// The insertion point may sit exactly at an exon boundary
/// (`ins_pos == exon.genomic_end`). In that case, `ins_pos` itself is
/// not inside any exon, but `ins_pos - 1` (the last exonic base) is.
/// This fallback ensures boundary insertions get the correct exon
/// constraint while naturally producing shift = 0 (the shift window
/// `[ins_pos, exon_end)` is empty when `ins_pos == exon_end`).
fn find_shift_exon_range(pos: u64, transcript: &TranscriptModel) -> Option<(u64, u64)> {
    containing_exon_range(pos, transcript).or_else(|| {
        if pos > 0 {
            containing_exon_range(pos - 1, transcript)
        } else {
            None
        }
    })
}

/// Find the exon range for a minus-strand insertion shift.
///
/// For minus-strand, 3' = lower genomic coordinates, so we shift
/// leftward from `ins_pos`. The insertion point may sit at the start
/// of an exon (`ins_pos == exon.genomic_start`), in which case the
/// fetch window `[exon_start, ins_pos)` is empty and shift = 0.
/// Try `ins_pos - 1` first (leftward neighbor), then `ins_pos` itself.
fn find_shift_exon_range_minus(pos: u64, transcript: &TranscriptModel) -> Option<(u64, u64)> {
    if pos > 0
        && let Some(r) = containing_exon_range(pos - 1, transcript)
    {
        return Some(r);
    }
    containing_exon_range(pos, transcript)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::write_genome_binary;
    use crate::types::{Biotype, CdsSegment, Exon, TranscriptTier};
    use tempfile::TempDir;

    /// Build a tempdir containing a UCSC-style FASTA with known repeat
    /// sequences for shift testing.
    ///
    /// Layout (0-based):
    /// ```text
    /// >chr1                 (length 60)
    /// pos 0..4:   TCGA      (non-repeat prefix)
    /// pos 4..10:  AAAAAG    (mono-A repeat at 4..8, then AG)
    /// pos 10..22: ATATATATCCCG (di-AT repeat at 10..17, non-repeat suffix)
    /// pos 22..46: ATACGTGATGGCATACGTGATGGCNN (12bp tandem: ATACGTGATGGC x2 + NN)
    /// pos 46..60: GCGCGCGCGCGCGC (mono-GC tail)
    /// ```
    ///
    /// >chr17                (length 40)
    /// ```text
    /// pos 0..10:  TCGATCGATC
    /// pos 10..16: CCCCCC      (mono-C repeat)
    /// pos 16..30: ATGATGATGATGAT (tri-ATG repeat)
    /// pos 30..40: NNNNNNNNN\n
    /// ```
    fn write_test_fasta() -> (TempDir, FastaReader) {
        let tmp = TempDir::new().unwrap();

        let chr1_seq = b"TCGAAAAAAGATATATATCCCGATACGTGATGGCATACGTGATGGCNNGCGCGCGCGCGCGC";
        let chr17_seq = b"TCGATCGATCCCCCCCATGATGATGATGATNNNNNNNNNN";

        assert_eq!(chr1_seq.len(), 62);
        assert_eq!(chr17_seq.len(), 40);

        let bin_path = tmp.path().join("test.bin");
        let idx_path = tmp.path().join("test.bin.idx");
        let contigs: &[(&str, &[u8])] = &[
            ("chr1", chr1_seq.as_slice()),
            ("chr17", chr17_seq.as_slice()),
        ];
        write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();
        let reader = FastaReader::open_with_assembly(&bin_path, crate::Assembly::GRCh38).unwrap();
        (tmp, reader)
    }

    /// Plus-strand transcript on chr1 with two exons separated by an intron.
    ///
    /// Exon 0: [0, 25)    covers the repeats and start of tandem
    /// Exon 1: [30, 62)   covers end of tandem + GC tail
    /// Intron: [25, 30)
    fn chr1_plus_transcript() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_SHIFT_PLUS.1".into(),
            protein_accession: Some("NP_SHIFT_PLUS.1".into()),
            gene_symbol: "SHIFTPLUS".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr1".into(),
            strand: Strand::Plus,
            tx_start: 0,
            tx_end: 62,
            cds_genomic_start: Some(0),
            cds_genomic_end: Some(62),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: 0,
                    genomic_end: 25,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: 30,
                    genomic_end: 62,
                },
            ],
            cds_segments: vec![
                CdsSegment {
                    exon_index: 0,
                    genomic_start: 0,
                    genomic_end: 25,
                    phase: 0,
                },
                CdsSegment {
                    exon_index: 1,
                    genomic_start: 30,
                    genomic_end: 62,
                    phase: 0,
                },
            ],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 2,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    /// Minus-strand transcript on chr17 with one exon covering the repeats.
    ///
    /// Exon 0: [8, 30)   covers the C repeat and ATG repeat
    fn chr17_minus_transcript() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_SHIFT_MINUS.1".into(),
            protein_accession: Some("NP_SHIFT_MINUS.1".into()),
            gene_symbol: "SHIFTMINUS".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr17".into(),
            strand: Strand::Minus,
            tx_start: 8,
            tx_end: 30,
            cds_genomic_start: Some(8),
            cds_genomic_end: Some(30),
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: 8,
                genomic_end: 30,
            }],
            cds_segments: vec![CdsSegment {
                exon_index: 0,
                genomic_start: 8,
                genomic_end: 30,
                phase: 0,
            }],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 1,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    #[test]
    fn del_in_mono_repeat_plus() {
        // chr1 has a 5-base A-run at [4, 9). Deleting [4, 5) should shift
        // right until the 3' boundary of the run, then stop at ref[9]=G.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 4, 5, &tx, &fasta).unwrap();
        assert_eq!(shift, 4);
    }

    #[test]
    fn del_in_di_repeat_plus() {
        // chr1 [10, 18) is ATATATAT. Deleting AT at [10, 12) shifts through
        // 3 full periods before ref[18]=C breaks the dinucleotide cycle.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 10, 12, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_no_repeat_plus() {
        // chr1 [0, 4) = TCGA: no repeated context, no shift possible.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 0, 1, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_at_exon_boundary_plus() {
        // Deleting the last exonic base (exon [0, 25), del [24, 25)) cannot
        // shift because the 3' boundary equals del_end.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 24, 25, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_zero_length() {
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 5, 5, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_in_mono_repeat_minus() {
        // chr17 [9, 16) is a 7-base C-run. On the minus strand, 3' shifts
        // toward decreasing coordinate; the run extends back to pos 9 and
        // stops at pos 8 (T).
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_in_mono_repeat_plus() {
        // 5-base A-run at chr1 [4, 9): inserting A inside it shifts to the
        // 3' boundary of the run.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 5, b"A", &tx, &fasta).unwrap();
        assert_eq!(shift, 4);
    }

    #[test]
    fn ins_in_di_repeat_plus() {
        // chr1 [10, 18) ATATATAT: inserting AT inside cycles through 3 full
        // periods before ref[18]=C breaks the cycle.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 12, b"AT", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_no_repeat_plus() {
        // No repeated context at the insertion point: shift must be 0.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 1, b"G", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_at_exon_boundary_plus() {
        // ins_pos == exon_end: the boundary fallback resolves to the same
        // exon, leaving no room to shift.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 25, b"A", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_zero_length() {
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 5, b"", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_in_mono_repeat_minus() {
        // 6-base C-run at chr17 [10, 16); pos 9 is also C, extending the
        // minus-strand run back to pos 8 (T).
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 15, b"C", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_in_tri_repeat_minus() {
        // chr17 [16, 30) ATGATGATGATGAT — three-base period. Insert ATG on
        // the minus strand at ins_pos=25; the cycle terminates after 9 bases
        // when the upstream context (chr17[7]=C) breaks the period.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 25, b"ATG", &tx, &fasta).unwrap();
        assert_eq!(shift, 9);
    }

    #[test]
    fn ins_intronic_no_repeat() {
        // Insert in the intron [25, 30) at a position with no matching
        // reference base: shift is 0.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 27, b"A", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    /// Plus-strand transcript on chr1 with intron placed over GC repeat.
    ///
    /// Exon 0: [0, 46)    covers everything up to the NN pad
    /// Exon 1: [56, 62)   covers the last 6 bases of the GC tail
    /// Intron: [46, 56)   = NN GCGCGCGC (GC repeat at positions 48-55)
    fn chr1_plus_intron_repeat_transcript() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_INTRON_PLUS.1".into(),
            protein_accession: Some("NP_INTRON_PLUS.1".into()),
            gene_symbol: "INTRONPLUS".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr1".into(),
            strand: Strand::Plus,
            tx_start: 0,
            tx_end: 62,
            cds_genomic_start: Some(0),
            cds_genomic_end: Some(62),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: 0,
                    genomic_end: 46,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: 56,
                    genomic_end: 62,
                },
            ],
            cds_segments: vec![
                CdsSegment {
                    exon_index: 0,
                    genomic_start: 0,
                    genomic_end: 46,
                    phase: 0,
                },
                CdsSegment {
                    exon_index: 1,
                    genomic_start: 56,
                    genomic_end: 62,
                    phase: 0,
                },
            ],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 2,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    /// Minus-strand transcript on chr17 with two exons flanking a large intron.
    ///
    /// Exons ordered 5'→3' on the transcript (descending genomic):
    /// Exon 0: [30, 40)   (transcript 5' end, highest genomic coords)
    /// Exon 1: [0, 10)    (transcript 3' end, lowest genomic coords)
    /// Intron: [10, 30)   = CCCCCC ATGATGATGATGAT NNNNNNNNNN
    fn chr17_minus_intron_transcript() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_INTRON_MINUS.1".into(),
            protein_accession: Some("NP_INTRON_MINUS.1".into()),
            gene_symbol: "INTRONMINUS".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr17".into(),
            strand: Strand::Minus,
            tx_start: 0,
            tx_end: 40,
            cds_genomic_start: Some(0),
            cds_genomic_end: Some(40),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: 30,
                    genomic_end: 40,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: 0,
                    genomic_end: 10,
                },
            ],
            cds_segments: vec![
                CdsSegment {
                    exon_index: 0,
                    genomic_start: 30,
                    genomic_end: 40,
                    phase: 0,
                },
                CdsSegment {
                    exon_index: 1,
                    genomic_start: 0,
                    genomic_end: 10,
                    phase: 0,
                },
            ],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 2,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    #[test]
    fn del_intronic_gc_repeat_plus() {
        // Intron [46, 56) holds GCGCGCGC after a 2-base NN pad. Deleting GC
        // at [48, 50) shifts to the 3' intron boundary.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 48, 50, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_intronic_single_base_no_repeat_plus() {
        // chr1[48]=G, chr1[49]=C: no repeat, no shift.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 48, 49, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_intronic_stops_at_boundary_plus() {
        // The shift must stop at the intron/exon boundary instead of crossing
        // into the next exon, even when the repeat would otherwise continue.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 52, 54, &tx, &fasta).unwrap();
        assert_eq!(shift, 2);
    }

    #[test]
    fn ins_intronic_gc_repeat_plus() {
        // Inserting GC inside the intron's GC repeat shifts to the 3' intron
        // boundary.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 50, b"GC", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_intronic_c_repeat_minus() {
        // 6-base C-run at the 3' end of the minus-strand intron shifts back
        // 5 bases before hitting the intron start.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, &tx, &fasta).unwrap();
        assert_eq!(shift, 5);
    }

    #[test]
    fn del_intronic_stops_at_boundary_minus() {
        // Deletion already at the minus-strand 3' intron boundary cannot
        // shift further.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 10, 11, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_intronic_c_repeat_minus() {
        // Inserting C inside a minus-strand C-run shifts back to the intron
        // start (5 bases).
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 15, b"C", &tx, &fasta).unwrap();
        assert_eq!(shift, 5);
    }
}
