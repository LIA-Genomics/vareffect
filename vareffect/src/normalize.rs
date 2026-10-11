//! HGVS 3' normalization for indels.
//!
//! Implements the HGVS 3' rule: indels in repeated sequences are shifted to
//! the most 3' (coding-strand rightward) equivalent position before HGVS
//! formatting. Both follow VEP's genomic shift and cross exon/intron
//! junctions.
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
use crate::types::Strand;

/// Compute the 3' shift offset (in genomic bases) for a deletion.
///
/// Starting from the deletion range `[del_start, del_end)` on the reference,
/// counts how many bases the deletion can slide toward the 3' end of the
/// coding strand while remaining equivalent: a deletion `[s, e)` moves one
/// base 3' on the plus strand iff `ref[s] == ref[e]`, on the minus strand
/// iff `ref[s - 1] == ref[e - 1]`.
///
/// The shift is on the genomic sequence, as VEP 116's `_genomic_shift` and
/// `perform_shift` do for insertions and deletions alike: it crosses
/// exon/intron junctions (VEP 116 NM_000632.4 `c.1497+1del` for a deletion
/// of the last exon base before an identical intron base), reads at most
/// [`VEP_SHIFT_WINDOW`] bases past the deletion (clamped to the chromosome),
/// and moves at most as far as [`vep_max_shift`] allows.
///
/// # Returns
///
/// Number of bases to shift (0 if no shift is possible).
///
/// # Errors
///
/// [`VarEffectError::ChromNotFound`] if `chrom` is not in the genome index,
/// or any error from [`FastaReader::fetch_sequence`].
pub(crate) fn compute_3prime_shift_deletion(
    chrom: &str,
    del_start: u64,
    del_end: u64,
    strand: Strand,
    fasta: &FastaReader,
) -> Result<u64, VarEffectError> {
    let del_len = del_end.saturating_sub(del_start);
    if del_len == 0 {
        return Ok(0);
    }
    let chrom_len = fasta
        .chrom_length(chrom)
        .ok_or_else(|| VarEffectError::ChromNotFound {
            chrom: chrom.to_string(),
        })?;
    let max_shift = vep_max_shift(del_len, strand);
    let shift = match strand {
        // 3' = increasing genomic coordinate: compare the base leaving the
        // deletion on its 5' side with the one entering on its 3' side.
        Strand::Plus => {
            let end = del_end.saturating_add(VEP_SHIFT_WINDOW).min(chrom_len);
            if end <= del_end {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence(chrom, del_start, end)?;
            let deleted = del_len as usize;
            seq.iter()
                .zip(&seq[deleted..])
                .take(max_shift)
                .take_while(|(leaving, entering)| leaving == entering)
                .count()
        }
        // 3' = decreasing genomic coordinate.
        Strand::Minus => {
            let start = del_start.saturating_sub(VEP_SHIFT_WINDOW);
            if start >= del_start {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence(chrom, start, del_end)?;
            let before = (del_start - start) as usize;
            seq[..before]
                .iter()
                .rev()
                .zip(seq.iter().rev())
                .take(max_shift)
                .take_while(|(entering, leaving)| entering == leaving)
                .count()
        }
    };
    Ok(shift as u64)
}

/// Most bases VEP 116's `perform_shift` moves an indel of `len` bases: its
/// loop runs over the [`VEP_SHIFT_WINDOW`] flanking bases minus the indel,
/// so `VEP_SHIFT_WINDOW + 1 - len`; a longer indel gets the whole window,
/// except a minus-strand one exactly `VEP_SHIFT_WINDOW + 1` long, which the
/// reverse loop never shifts.
fn vep_max_shift(len: u64, strand: Strand) -> usize {
    (if len <= VEP_SHIFT_WINDOW {
        VEP_SHIFT_WINDOW + 1 - len
    } else if strand == Strand::Minus && len == VEP_SHIFT_WINDOW + 1 {
        0
    } else {
        VEP_SHIFT_WINDOW
    }) as usize
}

/// Bases VEP reads on each side of an indel for its HGVS 3' shift
/// (`$area_to_search` in `TranscriptVariationAllele::_genomic_shift`).
const VEP_SHIFT_WINDOW: u64 = 1000;

/// Compute the 3' shift offset (in genomic bases) for an insertion.
///
/// Starting from the insertion point `ins_pos` (insertion is between
/// `pos - 1` and `pos` on the genome), counts how many bases the insertion
/// can slide toward the 3' end of the coding strand while remaining
/// equivalent: at each step the reference base entering on the 3' side must
/// match the corresponding inserted base (cycling through the inserted
/// sequence).
///
/// Both `fasta` bytes and `inserted_bases` are in **plus-strand**
/// orientation (VCF convention). The shift direction changes per strand
/// but the base comparison is always plus-strand vs plus-strand.
///
/// The shift is on the genomic sequence and bounded as for deletions
/// ([`compute_3prime_shift_deletion`], [`vep_max_shift`]): it crosses
/// exon/intron junctions (VEP NM_000548.5 `c.336+1_336+2dup`). VEP 116 leaves
/// an 812-base ALMS1 insertion at NM_001378454.1 `c.2427_2428ins`, which an
/// uncapped shift would move to `c.2239_3050dup`. The bases at the shifted
/// point are [`rotate_inserted_bases`].
///
/// # Returns
///
/// Number of bases to shift (0 if no shift is possible).
///
/// # Errors
///
/// [`VarEffectError::ChromNotFound`] if `chrom` is not in the genome index,
/// or any error from [`FastaReader::fetch_sequence`].
pub(crate) fn compute_3prime_shift_insertion(
    chrom: &str,
    ins_pos: u64,
    inserted_bases: &[u8],
    strand: Strand,
    fasta: &FastaReader,
) -> Result<u64, VarEffectError> {
    let ins_len = inserted_bases.len() as u64;
    if ins_len == 0 {
        return Ok(0);
    }
    let chrom_len = fasta
        .chrom_length(chrom)
        .ok_or_else(|| VarEffectError::ChromNotFound {
            chrom: chrom.to_string(),
        })?;
    let max_shift = vep_max_shift(ins_len, strand);
    let shift = match strand {
        // 3' = increasing genomic coordinate.
        Strand::Plus => {
            let end = ins_pos.saturating_add(VEP_SHIFT_WINDOW).min(chrom_len);
            if end <= ins_pos {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence(chrom, ins_pos, end)?;
            seq.iter()
                .zip(inserted_bases.iter().cycle())
                .take(max_shift)
                .take_while(|(reference, inserted)| reference == inserted)
                .count()
        }
        // 3' = decreasing genomic coordinate: walk left from `ins_pos - 1`,
        // matching the inserted bases from their end.
        Strand::Minus => {
            let start = ins_pos.saturating_sub(VEP_SHIFT_WINDOW);
            if start >= ins_pos {
                return Ok(0);
            }
            let seq = fasta.fetch_sequence(chrom, start, ins_pos)?;
            seq.iter()
                .rev()
                .zip(inserted_bases.iter().rev().cycle())
                .take(max_shift)
                .take_while(|(reference, inserted)| reference == inserted)
                .count()
        }
    };
    Ok(shift as u64)
}

/// Rotate inserted bases to match an insertion shifted 3' by `shift` bases.
///
/// Shifting an insertion one base 3' moves its first coding-strand base to
/// its end (VEP `perform_shift`, `hgvs_allele_string`): `GT` inserted after
/// `...A|GTGG` and shifted 3 bases is `TG` after `...AGTG|G`. The bases are
/// plus-strand, so a minus-strand shift rotates them right.
///
/// # Returns
///
/// The plus-strand inserted bases at the shifted insertion point.
pub(crate) fn rotate_inserted_bases(inserted_bases: &[u8], shift: u64, strand: Strand) -> Vec<u8> {
    let mut rotated = inserted_bases.to_vec();
    if rotated.is_empty() {
        return rotated;
    }
    let steps = (shift % rotated.len() as u64) as usize;
    match strand {
        Strand::Plus => rotated.rotate_left(steps),
        Strand::Minus => rotated.rotate_right(steps),
    }
    rotated
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    #[test]
    fn vep_max_shift_follows_the_perform_shift_loop() {
        // VEP 116 `perform_shift`: 1001 - len steps; a longer indel gets the
        // window, except a minus-strand one exactly 1001 long.
        assert_eq!(vep_max_shift(1, Strand::Plus), 1000);
        assert_eq!(vep_max_shift(1000, Strand::Minus), 1);
        assert_eq!(vep_max_shift(1001, Strand::Plus), 1000);
        assert_eq!(vep_max_shift(1001, Strand::Minus), 0);
        assert_eq!(vep_max_shift(2000, Strand::Minus), 1000);
    }

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
        let reader = FastaReader::open(&bin_path).unwrap();
        (tmp, reader)
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
        let shift = compute_3prime_shift_deletion("chr1", 4, 5, Strand::Plus, &fasta).unwrap();
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
        let shift = compute_3prime_shift_deletion("chr1", 10, 12, Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_no_repeat_plus() {
        // chr1 pos 0..4: TCGA — no repeats.
        // Delete T at [0, 1). ref[0]=T, ref[1]=C → no match.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr1", 0, 1, Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_zero_length() {
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr1", 5, 5, Strand::Plus, &fasta).unwrap();
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
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, Strand::Minus, &fasta).unwrap();
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
        let shift = compute_3prime_shift_insertion("chr1", 5, b"A", Strand::Plus, &fasta).unwrap();
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
        let shift =
            compute_3prime_shift_insertion("chr1", 12, b"AT", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_no_repeat_plus() {
        // chr1 pos 0: T. Insert G at pos 1 (between T and C).
        // ref[1] = C != G -> shift 0.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_insertion("chr1", 1, b"G", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_crosses_donor_junction_plus() {
        // Exon 0 ends at 25; chr1[25..28] = CGT. CG inserted after the last
        // exon base slides two bases into the intron (VEP's genomic shift).
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr1", 25, b"CG", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 2);
    }

    #[test]
    fn ins_crosses_acceptor_junction_plus() {
        // Intron [25, 30), exon 1 from 30; chr1[28..33] = GATGG. GAT inserted
        // at 28 slides over GATG into the exon and stops at the second G.
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr1", 28, b"GAT", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 4);
    }

    #[test]
    fn ins_zero_length() {
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_insertion("chr1", 5, b"", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    // -----------------------------------------------------------------------
    // Insertion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_in_mono_repeat_minus() {
        // chr17[9..16] = CCCCCCC, chr17[8] = T. Minus strand: 3' = decreasing
        // genomic. C inserted at 15 slides over chr17[14] .. chr17[9]: shift 6.
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr17", 15, b"C", Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_in_tri_repeat_minus() {
        // chr17 pos 16..30: ATGATGATGATGAT
        // Minus strand: 3' = decreasing genomic.
        // Insert ATG (plus-strand) at ins_pos=25 (between pos 24 and 25).
        // seq = chr17[..25], walked from its end:
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
        let shift =
            compute_3prime_shift_insertion("chr17", 25, b"ATG", Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 9);
    }

    #[test]
    fn ins_intronic_no_repeat() {
        // chr1 pos 27 is in the intron [25, 30).
        // No repeat at this position: chr1[27] differs from inserted A.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_insertion("chr1", 27, b"A", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    // -----------------------------------------------------------------------
    // Intronic shift — test transcripts
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Intronic deletion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_intronic_gc_repeat_crosses_into_the_exon_plus() {
        // Delete GC at [48, 50) in intron [46, 56). chr1[48..62] is a GC
        // repeat to the chromosome end, so the genomic shift (VEP 116
        // `perform_shift`) crosses the junction at 56: 12 bases.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr1", 48, 50, Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 12);
    }

    #[test]
    fn del_intronic_single_base_no_repeat_plus() {
        // Delete G at [48, 49) in intron [46, 56).
        // chr1[48] = G, chr1[49] = C → no match. shift = 0.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr1", 48, 49, Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn del_intronic_crosses_the_junction_plus() {
        // Delete GC at [52, 54) in intron [46, 56): the repeat runs on past
        // the junction at 56 to the chromosome end (62), so shift = 8.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr1", 52, 54, Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 8);
    }

    // -----------------------------------------------------------------------
    // Intronic insertion shift — plus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_intronic_gc_repeat_plus() {
        // Insert GC at ins_pos=50 in intron [46, 56).
        // chr1[50..62] = GCGCGCGCGCGC. Cycle through ins=[G,C]; the shift
        // crosses the intron/exon junction at 56 and stops at the end of
        // chr1: shift = 12.
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr1", 50, b"GC", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 12);
    }

    // -----------------------------------------------------------------------
    // Intronic deletion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn del_intronic_c_repeat_minus() {
        // chr17 = TCGATCGAT CCCCCCC ATG...: a C run at [9, 16). Delete C at
        // [15, 16); minus strand, 3' = leftward. The genomic shift passes the
        // intron start at 10 and stops at chr17[8] = T: shift = 6.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr17", 15, 16, Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn del_intronic_crosses_the_intron_start_minus() {
        // Delete C at [10, 11), the 3' end of intron [10, 30). chr17[9] = C
        // continues the run into the exon, so the shift crosses it: 1.
        let (_tmp, fasta) = write_test_fasta();
        let shift = compute_3prime_shift_deletion("chr17", 10, 11, Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 1);
    }

    // -----------------------------------------------------------------------
    // Intronic insertion shift — minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn ins_crosses_acceptor_junction_minus() {
        // chr17 intron [10, 30), exon [0, 10); chr17[9..15] = CCCCCC.
        // Minus strand: 3' = decreasing genomic. C inserted at 15 slides over
        // the intronic Cs and the last exon base chr17[9]; chr17[8] = T stops it.
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr17", 15, b"C", Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 6);
    }

    #[test]
    fn ins_crosses_donor_junction_minus() {
        // Exon [30, 40) starts at 30 (its last base in coding order is 30);
        // chr17[27..30] = GAT. AT inserted at 30 slides over AT into the
        // intron and stops at G.
        let (_tmp, fasta) = write_test_fasta();
        let shift =
            compute_3prime_shift_insertion("chr17", 30, b"AT", Strand::Minus, &fasta).unwrap();
        assert_eq!(shift, 2);
    }

    // -----------------------------------------------------------------------
    // Shift window and rotation
    // -----------------------------------------------------------------------

    /// A chromosome `C` + 1100 `A` + `C` and a chromosome `G` + 20 `A`.
    fn write_poly_a_fasta() -> (TempDir, FastaReader) {
        let tmp = TempDir::new().unwrap();
        let mut long = vec![b'C'];
        long.extend(std::iter::repeat_n(b'A', 1100));
        long.push(b'C');
        let mut short = vec![b'G'];
        short.extend(std::iter::repeat_n(b'A', 20));
        let bin_path = tmp.path().join("poly.bin");
        let idx_path = tmp.path().join("poly.bin.idx");
        let contigs: &[(&str, &[u8])] = &[("chr2", long.as_slice()), ("chr3", short.as_slice())];
        write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();
        let reader = FastaReader::open(&bin_path).unwrap();
        (tmp, reader)
    }

    #[test]
    fn ins_shift_capped_at_vep_window() {
        // VEP moves an insertion at most `1001 - len` bases.
        let (_tmp, fasta) = write_poly_a_fasta();
        let shift = |pos, bases: &[u8], strand| {
            compute_3prime_shift_insertion("chr2", pos, bases, strand, &fasta).unwrap()
        };
        assert_eq!(shift(1, b"A", Strand::Plus), 1000);
        assert_eq!(shift(1, b"AA", Strand::Plus), 999);
        assert_eq!(shift(1101, b"A", Strand::Minus), 1000);
        assert_eq!(shift(1101, b"AAA", Strand::Minus), 998);
        let long = vec![b'A'; 1001];
        assert_eq!(shift(1, &long, Strand::Plus), 1000);
        assert_eq!(shift(1101, &long, Strand::Minus), 0);
        let longer = vec![b'A'; 1002];
        assert_eq!(shift(1101, &longer, Strand::Minus), 1000);
    }

    #[test]
    fn ins_shift_stops_at_chromosome_end() {
        let (_tmp, fasta) = write_poly_a_fasta();
        let shift = compute_3prime_shift_insertion("chr3", 1, b"A", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 20);
        let shift = compute_3prime_shift_insertion("chr3", 21, b"A", Strand::Plus, &fasta).unwrap();
        assert_eq!(shift, 0);
    }

    #[test]
    fn ins_shift_unknown_chromosome_is_an_error() {
        let (_tmp, fasta) = write_poly_a_fasta();
        let err = compute_3prime_shift_insertion("chr9", 1, b"A", Strand::Plus, &fasta);
        assert!(matches!(err, Err(VarEffectError::ChromNotFound { .. })));
    }

    #[test]
    fn rotate_inserted_bases_plus_moves_first_base_to_end() {
        assert_eq!(rotate_inserted_bases(b"GT", 3, Strand::Plus), b"TG");
        assert_eq!(rotate_inserted_bases(b"GAT", 4, Strand::Plus), b"ATG");
        assert_eq!(rotate_inserted_bases(b"GT", 4, Strand::Plus), b"GT");
    }

    #[test]
    fn rotate_inserted_bases_minus_moves_last_base_to_front() {
        assert_eq!(rotate_inserted_bases(b"AT", 1, Strand::Minus), b"TA");
        assert_eq!(rotate_inserted_bases(b"ATG", 2, Strand::Minus), b"TGA");
        assert_eq!(rotate_inserted_bases(b"ATG", 0, Strand::Minus), b"ATG");
    }

    #[test]
    fn rotate_inserted_bases_empty() {
        assert!(rotate_inserted_bases(b"", 5, Strand::Plus).is_empty());
    }
}
