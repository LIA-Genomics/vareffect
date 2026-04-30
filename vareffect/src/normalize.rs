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
            // Find the region (exon or intron) containing the 3' end of
            // the deletion and use its boundary as the shift limit.
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
                return Ok(0); // no room to shift
            }
            let seq = fasta.fetch_sequence_slice(chrom, del_start, fetch_end)?;
            let dl = del_len as usize;
            let mut shift = 0usize;
            // seq[shift] = ref[del_start + shift]
            // seq[dl + shift] = ref[del_end + shift]
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
            // Find the region (exon or intron) containing the 3' coding
            // end (= lowest genomic coord) and use its boundary as limit.
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
                return Ok(0); // no room to shift
            }
            let seq = fasta.fetch_sequence_slice(chrom, fetch_start, del_end)?;
            let ds = (del_start - fetch_start) as usize;
            let de = (del_end - fetch_start) as usize;
            let mut shift = 0usize;
            // Compare ref[del_start - 1 - shift] with ref[del_end - 1 - shift].
            // In seq coords: seq[ds - 1 - shift] vs seq[de - 1 - shift].
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
            // Try exon first (with boundary fallback), then intron.
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
            // Try exon first (with boundary fallback), then intron.
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
            while shift < seq.len() {
                // Walk backward from ins_pos - 1: seq[seq.len() - 1 - shift]
                let ref_base = seq[seq.len() - 1 - shift];
                // Match from the end of inserted_bases, cycling.
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

    // -----------------------------------------------------------------------
    // Test FASTA helper
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Deletion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_in_mono_repeat_plus() {
        // chr1 pos 4..9: AAAAA, pos 9: G
        // Delete A at [4, 5). The run has 5 A's (pos 4..9).
        // ref[4]==ref[5]==A, ref[5]==ref[6]==A, ref[6]==ref[7]==A, ref[7]==ref[8]==A
        // ref[8] (A) != ref[9] (G) — but wait, pos 9 is 'G'.
        // Actually: TCGA AAAAA G ATATATATCCCG...
        //           0123 45678 9
        // Deleting [4,5): shift while ref[4+s]==ref[5+s].
        // s=0: ref[4]=A, ref[5]=A → match
        // s=1: ref[5]=A, ref[6]=A → match
        // s=2: ref[6]=A, ref[7]=A → match
        // s=3: ref[7]=A, ref[8]=A → match
        // s=4: ref[8]=A, ref[9]=G → stop
        // shift = 4
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 4, 5, &tx, &fasta).unwrap();
        assert_eq!(shift, 4);
    }

    #[test]
    fn del_in_di_repeat_plus() {
        // chr1 pos 10..18: ATATATAT
        // Delete AT at [10, 12).
        // s=0: ref[10]=A, ref[12]=A → match
        // s=1: ref[11]=T, ref[13]=T → match
        // s=2: ref[12]=A, ref[14]=A → match
        // s=3: ref[13]=T, ref[15]=T → match
        // s=4: ref[14]=A, ref[16]=A → match
        // s=5: ref[15]=T, ref[17]=T → match
        // s=6: ref[16]=A, ref[18]=C → stop
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 10, 12, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_no_repeat_plus() {
        // chr1 pos 0..4: TCGA — no repeats.
        // Delete T at [0, 1). ref[0]=T, ref[1]=C → no match.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 0, 1, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_at_exon_boundary_plus() {
        // chr1 exon 0 ends at 25. Delete A at [24, 25) — last base of exon.
        // ref[24] is the last exonic base. ref[25] is intronic (beyond exon).
        // Exon range is [0, 25), so fetch_end = 25. del_end = 25.
        // fetch_end <= del_end → return 0.
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

    // -----------------------------------------------------------------------
    // Deletion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_in_mono_repeat_minus() {
        // chr17 pos 9..16: CCCCCCC (7 C's — pos 9 is also C from "TCGATCGATC").
        // Minus strand: 3' = decreasing genomic.
        // Delete C at [15, 16). Exon [8, 30), so fetch [8, 16).
        // ds = 15 - 8 = 7, de = 16 - 8 = 8.
        // s=0: seq[6]=C, seq[7]=C → match
        // s=1: seq[5]=C, seq[6]=C → match
        // s=2: seq[4]=C, seq[5]=C → match
        // s=3: seq[3]=C, seq[4]=C → match
        // s=4: seq[2]=C, seq[3]=C → match
        // s=5: seq[1]=C, seq[2]=C → match (pos 9 = C)
        // s=6: seq[0]=T, seq[1]=C → stop (pos 8 = T)
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    // -----------------------------------------------------------------------
    // Insertion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_in_mono_repeat_plus() {
        // chr1 pos 4..9: AAAAA.
        // Insert A at ins_pos=5 (between pos 4 and 5).
        // Plus strand: check ref[5], ref[6], ref[7], ref[8] against inserted A.
        // ref[5]=A (match), ref[6]=A (match), ref[7]=A (match), ref[8]=A (match), ref[9]=G (stop)
        // shift = 4
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 5, b"A", &tx, &fasta).unwrap();
        assert_eq!(shift, 4);
    }

    #[test]
    fn ins_in_di_repeat_plus() {
        // chr1 pos 10..18: ATATATAT.
        // Insert AT at ins_pos=12 (between pos 11 and 12).
        // Check ref[12..] against [A,T] cycling:
        // ref[12]=A=ins[0] (match), ref[13]=T=ins[1] (match),
        // ref[14]=A=ins[0] (match), ref[15]=T=ins[1] (match),
        // ref[16]=A=ins[0] (match), ref[17]=T=ins[1] (match),
        // ref[18]=C!=ins[0]=A (stop)
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 12, b"AT", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_no_repeat_plus() {
        // chr1 pos 0: T. Insert G at pos 1 (between T and C).
        // ref[1] = C != G -> shift 0.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 1, b"G", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_at_exon_boundary_plus() {
        // Exon 0 ends at pos 25. Insert at ins_pos=25 (between pos 24 and 25).
        // find_shift_exon_range(25): containing_exon_range(25) finds exon 1 [30,62)?
        // No — pos 25 is in the intron [25,30). containing_exon_range returns None.
        // Fallback to pos-1=24 → exon 0 [0,25). fetch_end = 25.
        // fetch_end <= ins_pos (25 <= 25) → return 0.
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

    // -----------------------------------------------------------------------
    // Insertion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_in_mono_repeat_minus() {
        // chr17 pos 10..16: CCCCCC.
        // Minus strand: 3' = decreasing genomic.
        // Insert C at ins_pos=15 (between pos 14 and 15).
        // find_shift_exon_range_minus(15): pos-1=14 is in exon [8,30).
        // Fetch [8, 15) = "TCCCCCCC"... wait, let me re-examine.
        // chr17: TCGATCGATC CCCCCC ATGATGATGATGAT NNNNNNNNN
        //        0123456789 012345 0123456789...
        //                   10     16
        // Fetch [8, 15) = chr17[8..15] = "TCCCCCCC"
        // Wait: chr17[8]=T, chr17[9]=C, chr17[10..15]=CCCCC → "TCCCCCC" (7 bytes)
        // seq = [T, C, C, C, C, C, C]
        //
        // s=0: seq[6]=C, ins[0]=C → match (pos 14)
        // s=1: seq[5]=C, ins[0]=C → match (pos 13)
        // s=2: seq[4]=C, ins[0]=C → match (pos 12)
        // s=3: seq[3]=C, ins[0]=C → match (pos 11)
        // s=4: seq[2]=C, ins[0]=C → match (pos 10)
        // s=5: seq[1]=C, ins[0]=C → match (pos 9)
        // s=6: seq[0]=T, ins[0]=C → stop
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 15, b"C", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_in_tri_repeat_minus() {
        // chr17 pos 16..30: ATGATGATGATGAT
        // Minus strand: 3' = decreasing genomic.
        // Insert ATG (plus-strand) at ins_pos=25 (between pos 24 and 25).
        // Fetch [8, 25) = chr17[8..25].
        // chr17[8..25] = "TCCCCCCCATGATGATG" (17 bytes)
        //
        // ins = [A, T, G], ins_len = 3
        // s=0: seq[16]=G, ins[3-1-0]=ins[2]=G → match (pos 24)
        // s=1: seq[15]=T, ins[3-1-1]=ins[1]=T → match (pos 23)
        // s=2: seq[14]=A, ins[3-1-2]=ins[0]=A → match (pos 22)
        // s=3: seq[13]=G, ins[2]=G → match (pos 21)
        // s=4: seq[12]=T, ins[1]=T → match (pos 20)
        // s=5: seq[11]=A, ins[0]=A → match (pos 19)
        // s=6: seq[10]=G, ins[2]=G → match (pos 18)
        // s=7: seq[9]=T, ins[1]=T → match (pos 17)
        // s=8: seq[8]=A, ins[0]=A → match (pos 16)
        // s=9: seq[7]=C, ins[2]=G → stop
        // shift = 9
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 25, b"ATG", &tx, &fasta).unwrap();
        assert_eq!(shift, 9);
    }

    #[test]
    fn ins_intronic_no_repeat() {
        // chr1 pos 27 is in the intron [25, 30).
        // No repeat at this position: chr1[27] differs from inserted A.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 27, b"A", &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    // -----------------------------------------------------------------------
    // Intronic shift — test transcripts
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Intronic deletion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_intronic_gc_repeat_plus() {
        // Intron [46, 56). chr1[46..56] = NNGCGCGCGC.
        // Delete GC at [48, 50). 3' direction = rightward.
        // Shift window: [48, 56). seq = chr1[48..56] = GCGCGCGC (8 bytes).
        // dl = 2.
        // s=0: seq[0]=G == seq[2]=G, s=1: seq[1]=C == seq[3]=C,
        // s=2: G==G, s=3: C==C, s=4: G==G, s=5: C==C,
        // s=6: dl + 6 = 8 == seq.len() → stop.
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 48, 50, &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_intronic_single_base_no_repeat_plus() {
        // Delete G at [48, 49) in intron [46, 56).
        // chr1[48] = G, chr1[49] = C → no match. shift = 0.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 48, 49, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_intronic_stops_at_boundary_plus() {
        // Delete GC at [52, 54) in intron [46, 56).
        // Shift window: [52, 56). seq = GCGC (4 bytes), dl = 2.
        // s=0: G==G, s=1: C==C, s=2: dl+2 = 4 == seq.len() → stop.
        // shift = 2 (clamped at intron end, does NOT cross into exon 1).
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_deletion("chr1", 52, 54, &tx, &fasta).unwrap();
        assert_eq!(shift, 2);
    }

    // -----------------------------------------------------------------------
    // Intronic insertion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_intronic_gc_repeat_plus() {
        // Insert GC at ins_pos=50 in intron [46, 56).
        // chr1[50..56] = GCGCGC. Cycle through ins=[G,C].
        // s=0: G==G, s=1: C==C, s=2: G==G, s=3: C==C,
        // s=4: G==G, s=5: C==C, s=6: end of window → stop.
        // shift = 6
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr1_plus_intron_repeat_transcript();
        let shift = compute_3prime_shift_insertion("chr1", 50, b"GC", &tx, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    // -----------------------------------------------------------------------
    // Intronic deletion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_intronic_c_repeat_minus() {
        // chr17 intron [10, 30). chr17[10..16] = CCCCCC.
        // Delete C at [15, 16). Minus strand: 3' = decreasing genomic.
        // boundary_start = 10 (intron start). Fetch [10, 16) = CCCCCC.
        // ds = 15 - 10 = 5, de = 16 - 10 = 6.
        // s=0: seq[4]=C vs seq[5]=C → match
        // s=1: seq[3]=C vs seq[4]=C → match
        // s=2: seq[2]=C vs seq[3]=C → match
        // s=3: seq[1]=C vs seq[2]=C → match
        // s=4: seq[0]=C vs seq[1]=C → match
        // s=5: s == ds → stop
        // shift = 5
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, &tx, &fasta).unwrap();
        assert_eq!(shift, 5);
    }

    #[test]
    fn del_intronic_stops_at_boundary_minus() {
        // Delete C at [10, 11) in intron [10, 30). Minus strand: 3' = leftward.
        // boundary_start = 10. fetch_start (10) >= del_start (10) → return 0.
        // At the 3' intron boundary, no room to shift.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_deletion("chr17", 10, 11, &tx, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    // -----------------------------------------------------------------------
    // Intronic insertion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_intronic_c_repeat_minus() {
        // chr17 intron [10, 30). Insert C at ins_pos=15.
        // Minus strand: 3' = decreasing genomic. boundary_start = 10.
        // Fetch [10, 15) = CCCCC (5 bytes).
        // Walk backward: seq[4]=C==C, seq[3]=C==C, seq[2]=C==C,
        // seq[1]=C==C, seq[0]=C==C → shift = 5.
        let (_tmp, fasta) = write_test_fasta();
        let tx = chr17_minus_intron_transcript();
        let shift = compute_3prime_shift_insertion("chr17", 15, b"C", &tx, &fasta).unwrap();
        assert_eq!(shift, 5);
    }
}
