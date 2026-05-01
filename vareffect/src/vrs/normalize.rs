//! VOCA "fully-justified" allele normalization.
//!
//! The Variant Overprecision Correction Algorithm (NCBI VOCA, also
//! called "expand-and-trim" or "fully-justified" normalization) is the
//! GA4GH-spec'd canonical form for VRS computed identifiers. It differs
//! from VEP-style minimal trim in repeat regions: where minimal trim
//! preserves the input's choice of position within the ambiguity, VOCA
//! expands the location to span the **entire** ambiguity region, giving
//! a representation that is unique per biological allele.
//!
//! Without VOCA, two pipelines describing the same biological insertion
//! / deletion in a homopolymer or microsatellite can hash to different
//! VRS IDs — defeating any cross-pipeline content-addressed index.
//!
//! # Algorithm
//!
//! 1. **Minimal trim:** strip shared prefix (left-first) then shared
//!    suffix to produce minimum-length ref/alt.
//! 2. **SNV / no-op / delins:** return the trimmed form (no rolling).
//!    SNVs and complex delins do not have ambiguous positions in the
//!    reference, so EXPAND collapses to the trimmed form.
//! 3. **Pure insertion / deletion:** roll right (extend `end` while the
//!    next reference base equals the indel's leading rotated base, with
//!    the indel rotating left each step) and then roll left (extend
//!    `start` while the previous reference base equals the indel's
//!    trailing rotated base, with the indel rotating right). This
//!    expands the location to cover the maximum interval over which the
//!    rotated indel remains an equivalent variant.
//! 4. **Reconstruct state:**
//!    - Insertion: `state = (rotated indel) ++ ref[start:end]`. The
//!      inserted bases anchor at the leftmost edge after rolling — this
//!      preserves variant equivalence under rotation.
//!    - Deletion: `state = ref[start:end][..len - indel_len]`. The
//!      deletion canonically removes the trailing `indel_len` bases of
//!      the expanded region; in a repeat all positions are equivalent.
//!
//! # Bounds
//!
//! Rolling runs to convergence — it stops only at a base mismatch or at
//! the chromosome edge. Defensive caps were considered but rejected: a
//! capped result would silently disagree with every other VRS-emitting
//! tool on the (rare but real) >1 kb tandem repeats, and the cost is
//! cheap (single-base mmap reads at ~5 ns each — even a 100 kb roll is
//! sub-millisecond).
//!
//! # FASTA reads
//!
//! Single-base reads via `FastaReader::fetch_base` are ~5 ns each
//! (memory-mapped). For typical indels (1-10 bp ambiguity) the total
//! cost is dominated by the SHA-512 hash, not the reads. We therefore
//! do not bother with batched / buffered reads.

use crate::error::VarEffectError;
use crate::fasta::FastaReader;

/// Result of VOCA normalization: a fully-justified VRS-canonical
/// allele expressed as 0-based interbase coordinates and the alt
/// sequence to substitute at that interval.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct NormalizedAllele {
    /// 0-based interbase start of the allele's location.
    pub start: u64,
    /// 0-based interbase end of the allele's location. `end == start`
    /// for pure insertions; otherwise `end > start`.
    pub end: u64,
    /// Alt bases to replace `ref[start:end]` with. Empty string for
    /// pure deletions in the trimmed form (rare after VOCA — most
    /// pure deletions in repeats expand to a non-trivial alt).
    pub alt: Vec<u8>,
}

/// Strip shared prefix (left-first, matches VEP `--minimal`) and
/// shared suffix from a (ref, alt) pair, returning (trimmed_ref,
/// trimmed_alt, prefix_strip_count). Caller adds `prefix_strip_count`
/// to the original 0-based position to get the trimmed start.
fn minimal_trim<'a>(ref_allele: &'a [u8], alt_allele: &'a [u8]) -> (&'a [u8], &'a [u8], u64) {
    let mut r = ref_allele;
    let mut a = alt_allele;
    let mut prefix: u64 = 0;
    while !r.is_empty() && !a.is_empty() && r[0] == a[0] {
        r = &r[1..];
        a = &a[1..];
        prefix += 1;
    }
    while !r.is_empty() && !a.is_empty() && r[r.len() - 1] == a[a.len() - 1] {
        r = &r[..r.len() - 1];
        a = &a[..a.len() - 1];
    }
    (r, a, prefix)
}

/// Compute the VOCA fully-justified normalized form of a variant.
///
/// Returns:
/// - `Ok(Some(_))` on success.
/// - `Ok(None)` for a no-op (`ref == alt`) — VRS does not assign an ID to
///   a non-variant — or when the chromosome is not present in the reader
///   (e.g. trimmed test fixtures).
/// - `Err(_)` for a real FASTA I/O error during rolling. Bounds checks
///   guarantee `fetch_base` is only called within the chromosome interior,
///   so any error here represents an underlying mmap / contig-table fault
///   rather than a normal stopping condition. Callers may choose to
///   downgrade `Err` to "no VRS ID available" — but the layer of
///   indirection lets them log/alert if they need to.
///
/// `pos` is the 0-based VCF position. `ref_allele` / `alt_allele` are
/// the original VCF allele bytes (uppercase ASCII).
pub(super) fn voca_normalize(
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    fasta: &FastaReader,
) -> Result<Option<NormalizedAllele>, VarEffectError> {
    let (trimmed_ref, trimmed_alt, prefix) = minimal_trim(ref_allele, alt_allele);

    // No-op (ref == alt). VRS does not assign an ID to a non-variant.
    if trimmed_ref.is_empty() && trimmed_alt.is_empty() {
        return Ok(None);
    }

    let mut start = pos + prefix;
    let mut end = start + trimmed_ref.len() as u64;

    // SNV: both lengths 1, ref != alt. No rolling, no FASTA reads.
    if trimmed_ref.len() == 1 && trimmed_alt.len() == 1 {
        return Ok(Some(NormalizedAllele {
            start,
            end,
            alt: trimmed_alt.to_vec(),
        }));
    }

    // Delins (both non-empty, lengths != 1). VRS EXPAND mode does not
    // roll complex delins; they retain their minimal-trim form.
    if !trimmed_ref.is_empty() && !trimmed_alt.is_empty() {
        return Ok(Some(NormalizedAllele {
            start,
            end,
            alt: trimmed_alt.to_vec(),
        }));
    }

    // Pure insertion or pure deletion — eligible for rolling.
    let is_insertion = trimmed_ref.is_empty();
    let mut indel: Vec<u8> = if is_insertion {
        trimmed_alt.to_vec()
    } else {
        trimmed_ref.to_vec()
    };

    // Chrom not in this reader (typical for test fixtures with a trimmed
    // contig set). Distinguished from a real I/O error: the variant
    // simply has no canonical form on this build.
    let Some(chrom_len) = fasta.chrom_length(chrom) else {
        return Ok(None);
    };

    // Roll right: extend `end` while the next ref base matches the
    // indel's leading base, rotating the indel left each step. Bounds
    // check before fetch so any fetch error is genuinely an I/O fault.
    while end < chrom_len {
        let next_base = fasta.fetch_base(chrom, end)?;
        if indel[0] != next_base {
            break;
        }
        end += 1;
        indel.rotate_left(1);
    }

    // Roll left: extend `start` while the previous ref base matches the
    // indel's trailing base, rotating right.
    while start > 0 {
        let prev_base = fasta.fetch_base(chrom, start - 1)?;
        let last_idx = indel.len() - 1;
        if indel[last_idx] != prev_base {
            break;
        }
        start -= 1;
        indel.rotate_right(1);
    }

    // `fetch_sequence` rejects zero-width ranges with an error. For a
    // pure insertion that didn't roll (start == end after rolling), the
    // expanded reference is just the empty slice.
    let expanded_ref = if start == end {
        Vec::new()
    } else {
        fasta.fetch_sequence(chrom, start, end)?
    };
    let alt = if is_insertion {
        // Inserted bases anchor at the leftmost edge of the expanded
        // region. After left-rolling, `indel` is rotated such that
        // prepending it to `expanded_ref` reconstructs the post-
        // insertion sequence at this location.
        let mut v = Vec::with_capacity(indel.len() + expanded_ref.len());
        v.extend_from_slice(&indel);
        v.extend_from_slice(&expanded_ref);
        v
    } else {
        // Deletion: trailing `indel_len` bases of the expanded region
        // are dropped. expanded_ref.len() >= indel.len() is invariant —
        // initial end-start equals indel.len() and rolling only grows
        // the interval — so the slice is always well-formed.
        debug_assert!(expanded_ref.len() >= indel.len());
        expanded_ref[..expanded_ref.len() - indel.len()].to_vec()
    };

    Ok(Some(NormalizedAllele { start, end, alt }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Assembly;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    /// Build a tempdir-backed synthetic FASTA with embedded repeat
    /// structures that exercise every VOCA branch. Layout (0-based):
    ///
    /// ```text
    /// >chr1   length 30
    /// pos  0..2:  GG          non-repeat prefix
    /// pos  2..7:  AAAAA       mono-A homopolymer (5 bases)
    /// pos  7..9:  CC          non-repeat spacer
    /// pos  9..15: ATATAT      di-AT repeat (3 units)
    /// pos 15..17: GG          non-repeat spacer
    /// pos 17..23: ATGATG      tri-ATG repeat (2 units)
    /// pos 23..30: TTTTTTT     mono-T tail (used for edge cases)
    ///
    /// >chr2   length 8
    /// pos  0..8:  ACGTACGT    no internal repeats; isolated bases
    /// ```
    fn synthetic_fasta() -> (TempDir, FastaReader) {
        let tmp = TempDir::new().unwrap();
        let bin = tmp.path().join("t.bin");
        let idx = tmp.path().join("t.bin.idx");
        let chr1: &[u8] = b"GGAAAAACCATATATGGATGATGTTTTTTT";
        let chr2: &[u8] = b"ACGTACGT";
        assert_eq!(chr1.len(), 30);
        assert_eq!(chr2.len(), 8);
        write_genome_binary(&[("chr1", chr1), ("chr2", chr2)], "test", &bin, &idx).unwrap();
        let reader = FastaReader::open_with_assembly(&bin, Assembly::GRCh38).unwrap();
        (tmp, reader)
    }

    /// Helper: unwrap a successful normalization. Tests that expect a
    /// hit chain `.unwrap().unwrap()` on `Result<Option<_>, _>`; this
    /// helper makes the intent explicit and the diagnostic clearer.
    fn ok_some(r: Result<Option<NormalizedAllele>, VarEffectError>) -> NormalizedAllele {
        r.expect("voca_normalize errored")
            .expect("voca_normalize returned None for an expected hit")
    }

    #[test]
    fn minimal_trim_strips_common_prefix_then_suffix() {
        // CTC -> CGC: trim leading C, then trailing C → T -> G.
        let (r, a, p) = minimal_trim(b"CTC", b"CGC");
        assert_eq!(r, b"T");
        assert_eq!(a, b"G");
        assert_eq!(p, 1);

        // Pure insertion (anchor base): C -> CAT
        let (r, a, p) = minimal_trim(b"C", b"CAT");
        assert_eq!(r, b"");
        assert_eq!(a, b"AT");
        assert_eq!(p, 1);

        // Pure deletion: CAT -> C
        let (r, a, p) = minimal_trim(b"CAT", b"C");
        assert_eq!(r, b"AT");
        assert_eq!(a, b"");
        assert_eq!(p, 1);

        // Identity: ref == alt
        let (r, a, p) = minimal_trim(b"ACGT", b"ACGT");
        assert_eq!(r, b"");
        assert_eq!(a, b"");
        assert_eq!(p, 4);

        // No common bases
        let (r, a, p) = minimal_trim(b"A", b"T");
        assert_eq!(r, b"A");
        assert_eq!(a, b"T");
        assert_eq!(p, 0);
    }

    #[test]
    fn no_op_returns_none() {
        let (_tmp, fasta) = synthetic_fasta();
        let r = voca_normalize("chr1", 0, b"GG", b"GG", &fasta).unwrap();
        assert!(r.is_none());
    }

    #[test]
    fn snv_does_not_roll() {
        // chr1:0 G>C — non-repeat, plain SNV.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 0, b"G", b"C", &fasta));
        assert_eq!(n.start, 0);
        assert_eq!(n.end, 1);
        assert_eq!(n.alt, b"C");
    }

    #[test]
    fn delins_keeps_minimal_form() {
        // chr1:0 GG>CT — delins (length-preserving but mixed).
        // VRS EXPAND mode does not roll delins.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 0, b"GG", b"CT", &fasta));
        assert_eq!(n.start, 0);
        assert_eq!(n.end, 2);
        assert_eq!(n.alt, b"CT");
    }

    #[test]
    fn ins_in_mono_homopolymer_expands_full_run() {
        // chr1 has AAAAA at [2,7). Insert one A — VOCA expands the
        // location to [2,7) and emits "AAAAAA" (6 bases).
        // Anchored input: ref="G", alt="GA" at pos=1. Trim of (G,GA)
        // leaves ref="" alt="A" with prefix=1, so trimmed pos = 2 — right
        // at the start of the AAAAA run.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 1, b"G", b"GA", &fasta));
        assert_eq!(n.start, 2);
        assert_eq!(n.end, 7);
        assert_eq!(n.alt, b"AAAAAA");
    }

    #[test]
    fn del_in_mono_homopolymer_expands_full_run() {
        // chr1 has AAAAA at [2,7). Delete one A.
        // Anchored input: ref="GA", alt="G" at pos=1.
        // Trim leaves ref="A" alt="" with prefix=1, pos=2.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 1, b"GA", b"G", &fasta));
        assert_eq!(n.start, 2);
        assert_eq!(n.end, 7);
        assert_eq!(n.alt, b"AAAA");
    }

    #[test]
    fn ins_in_di_repeat_expands_full_repeat_region() {
        // chr1 has ATATAT at [9,15). Insert "AT" between any unit
        // boundary — all positions in {9,11,13,15} produce the same
        // post-insertion sequence ATATATAT.
        // Anchored input at the inner boundary: pos=10, ref="T", alt="TAT".
        // Trim → ref="" alt="AT" prefix=1, trimmed pos=11.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 10, b"T", b"TAT", &fasta));
        assert_eq!(n.start, 9);
        assert_eq!(n.end, 15);
        // Alt is the 4-AT-unit form (8 bases).
        assert_eq!(n.alt, b"ATATATAT");
    }

    #[test]
    fn del_in_di_repeat_expands_full_repeat_region() {
        // chr1 has ATATAT at [9,15). Delete one "AT" unit.
        // Anchored input: pos=10, ref="TAT", alt="T".
        // Trim → ref="AT" alt="" prefix=1, trimmed pos=11.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 10, b"TAT", b"T", &fasta));
        assert_eq!(n.start, 9);
        assert_eq!(n.end, 15);
        // 6 bases - 2 (deleted unit) = 4 remaining = "ATAT".
        assert_eq!(n.alt, b"ATAT");
    }

    #[test]
    fn ins_in_tri_repeat_expands_full_repeat_region() {
        // chr1 has ATGATG at [17,23). Insert one "ATG" unit at the
        // upstream boundary (interbase 17).
        // Anchored input: pos=16, ref="G", alt="GATG".
        // Trim → ref="" alt="ATG" prefix=1, trimmed pos=17.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 16, b"G", b"GATG", &fasta));
        // After rolling right through the ATG repeat and rolling left
        // into the boundary G at pos 16, location spans [16,23) and alt
        // is the expanded GATGATGATG.
        assert_eq!(n.start, 16);
        assert_eq!(n.end, 23);
        assert_eq!(n.alt, b"GATGATGATG"); // 7-base expanded ref + 3-base indel
    }

    #[test]
    fn ins_outside_repeat_does_not_roll() {
        // Insert "C" between the GG at [0,2). No repeat to roll into.
        // Input: pos=0, ref="G", alt="GC". Trim → ref="" alt="C" prefix=1.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 0, b"G", b"GC", &fasta));
        assert_eq!(n.start, 1);
        assert_eq!(n.end, 1);
        // No expansion — alt is just the inserted base.
        assert_eq!(n.alt, b"C");
    }

    #[test]
    fn del_in_double_c_run_rolls_one_position() {
        // Delete the first C of CC at [7,9). The C-run is a 2-mer
        // mono-C repeat, so VOCA expands across both Cs and emits the
        // single-C deletion form.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 6, b"AC", b"A", &fasta));
        assert_eq!(n.start, 7);
        assert_eq!(n.end, 9);
        assert_eq!(n.alt, b"C");
    }

    #[test]
    fn del_isolated_base_does_not_roll() {
        // S0-5: an isolated base flanked by non-matching bases must NOT
        // roll. chr2 = "ACGTACGT" — each base is unique in its 3-mer
        // neighborhood. Deleting the C at chr2:1 (anchored: pos=0,
        // ref="AC", alt="A") leaves ref="C" alt="" with prefix=1,
        // trimmed pos=1. Roll right tests ref[2]='G' ≠ 'C' (stop). Roll
        // left tests ref[0]='A' ≠ 'C' (stop). Result: minimal-trim form.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr2", 0, b"AC", b"A", &fasta));
        assert_eq!(n.start, 1);
        assert_eq!(n.end, 2);
        assert_eq!(n.alt, b"");
    }

    #[test]
    fn roll_at_chrom_end_clamps_safely() {
        // chr1 ends with TTTTTTT (positions 23-30). Delete one T from
        // the homopolymer — rolling right will hit chrom end.
        // Input: pos=22, ref="GT", alt="G".
        // Trim → ref="T" alt="" prefix=1, trimmed pos=23.
        let (_tmp, fasta) = synthetic_fasta();
        let n = ok_some(voca_normalize("chr1", 22, b"GT", b"G", &fasta));
        // Rolls right to chrom_len=30, then attempts to roll left
        // through TTTTTT prefix until hitting non-T at pos 22 (G).
        assert_eq!(n.start, 23);
        assert_eq!(n.end, 30);
        assert_eq!(n.alt, b"TTTTTT"); // 7-1 = 6 Ts
    }

    #[test]
    fn missing_chrom_returns_ok_none_for_indel() {
        // Missing chrom is "no canonical form on this build", not an
        // I/O error — voca_normalize returns Ok(None), preserving the
        // distinction from real fetch errors.
        let (_tmp, fasta) = synthetic_fasta();
        let r = voca_normalize("chrUNKNOWN", 0, b"A", b"AT", &fasta).unwrap();
        assert!(r.is_none());
    }

    #[test]
    fn ins_in_long_homopolymer_runs_to_completion() {
        // S0-3: confirm the rolling no longer caps at 1000 bp. Build a
        // dedicated 1500-bp homopolymer fixture and insert one A; VOCA
        // must expand the location across the entire run, not stop at
        // an arbitrary cap.
        let tmp = TempDir::new().unwrap();
        let bin = tmp.path().join("t.bin");
        let idx = tmp.path().join("t.bin.idx");
        let mut chr_long = vec![b'G', b'G'];
        chr_long.extend(std::iter::repeat_n(b'A', 1500));
        chr_long.extend_from_slice(b"GG");
        let total_len = chr_long.len() as u64;
        write_genome_binary(&[("chr1", &chr_long)], "test", &bin, &idx).unwrap();
        let fasta = FastaReader::open_with_assembly(&bin, Assembly::GRCh38).unwrap();

        // Anchor at pos=1 (G), ref="G", alt="GA". Trim → trimmed pos=2,
        // start of the 1500-bp A run.
        let n = ok_some(voca_normalize("chr1", 1, b"G", b"GA", &fasta));
        assert_eq!(n.start, 2);
        assert_eq!(n.end, total_len - 2);
        assert_eq!(n.alt.len(), 1501); // 1500 ref + 1 inserted A
        assert!(n.alt.iter().all(|&b| b == b'A'));
    }
}
