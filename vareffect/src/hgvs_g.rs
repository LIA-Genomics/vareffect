//! HGVS genomic (`g.`) notation generation.
//!
//! Builds a canonical genomic-HGVS string (e.g. `NC_000017.11:g.7674220C>G`)
//! from a plus-strand VCF-style variant — a chromosome, a 0-based position,
//! and REF/ALT alleles — plus the reference genome for indel normalization.
//! This is the genomic counterpart of the [`crate::hgvs_c`] / [`crate::hgvs_p`]
//! formatters and the concordance target is Ensembl VEP's `--hgvsg` output.
//!
//! # Nomenclature rules implemented
//!
//! Per the HGVS recommendations (den Dunnen et al. 2016, Hum Mutat; HGVS
//! Nomenclature v21.1) on 1-based positions:
//!
//! - substitution — `g.{pos}{REF}>{ALT}` (single base → single base only)
//! - deletion — `g.{pos}del` / `g.{start}_{end}del` (deleted bases not written)
//! - duplication — `g.{pos}dup` / `g.{start}_{end}dup` (tandem 3'-flanking copy)
//! - insertion — `g.{pos}_{pos+1}ins{SEQ}`
//! - deletion-insertion — `g.{pos}delins{SEQ}` / `g.{start}_{end}delins{SEQ}`
//!
//! Indels and duplications are placed at the **most 3' (rightmost) genomic
//! position** (the "3' rule"); VCF left-aligns, so a naive conversion would be
//! non-canonical and disagree with VEP/ClinVar. An insertion whose 3'-shifted
//! sequence copies the immediately-preceding reference is emitted as a
//! `dup`, not an `ins`.
//!
//! # Safety
//!
//! Every classification/formatting inability returns `Ok(None)` — a missing
//! `g.` string is always preferred over a wrong one, since the output is a
//! clinical notation a user may copy verbatim. `Err` is reserved for genuine
//! reference-genome I/O faults.

use crate::chrom;
use crate::error::VarEffectError;
use crate::fasta::FastaReader;

/// Build canonical genomic-HGVS (`g.`) notation for a plus-strand variant.
///
/// # Arguments
///
/// * `chrom` — Chromosome name (UCSC `"chr17"`, or an `NC_*` accession). Bare
///   Ensembl names, patch/alt contigs, and unknown names yield `Ok(None)`.
/// * `pos` — 0-based genomic position of the first REF base (vareffect
///   convention, matching [`crate::VarEffect::annotate`]).
/// * `ref_allele` / `alt_allele` — VCF REF/ALT allele bytes. Case-insensitive;
///   non-`ACGTN` (symbolic / breakend) alleles yield `Ok(None)`.
/// * `fasta` — Reference genome, used for the 3'-shift and duplication test.
///
/// # Returns
///
/// * `Ok(Some(s))` — the canonical `g.` string.
/// * `Ok(None)` — the variant cannot be expressed as a canonical `g.` string
///   (unknown contig, symbolic allele, REF/genome mismatch, or a degenerate
///   contig-boundary case).
///
/// # Errors
///
/// [`VarEffectError`] on a genuine reference-genome read fault.
pub fn format_hgvs_g(
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    fasta: &FastaReader,
) -> Result<Option<String>, VarEffectError> {
    let ref_u = ref_allele.to_ascii_uppercase();
    let alt_u = alt_allele.to_ascii_uppercase();
    if ref_u.is_empty() || alt_u.is_empty() || !is_acgtn(&ref_u) || !is_acgtn(&alt_u) {
        return Ok(None);
    }

    let acc = match resolve_accession(chrom) {
        Some(acc) => acc,
        None => return Ok(None),
    };
    // Also confirms the contig is present, so the shift-loop `fetch_base`
    // calls below (all bounded by `chrom_len`) cannot raise a coordinate or
    // chromosome error.
    let chrom_len = match fasta.chrom_length(chrom) {
        Some(len) => len,
        None => return Ok(None),
    };

    // Fail closed if the stated REF disagrees with the reference genome — a
    // coordinate or reference-build mismatch must yield an absent notation,
    // never a wrong one. Verifying the full untrimmed REF here also covers the
    // insertion-anchor base: an insertion's trimmed ref is empty, so a
    // per-shape check could not catch a mismatched anchor.
    if !fasta.verify_ref(chrom, pos, &ref_u)? {
        return Ok(None);
    }

    let (start, r, a) = trim_to_minimal_edit(pos, &ref_u, &alt_u);

    match (r.as_slice(), a.as_slice()) {
        // Identical alleles — not a variant.
        ([], []) => Ok(None),
        // Substitution.
        ([rb], [ab]) => Ok(Some(format!(
            "{acc}:g.{}{}>{}",
            start + 1,
            char::from(*rb),
            char::from(*ab)
        ))),
        // Deletion.
        (del, []) => {
            let len = del.len() as u64;
            let s = shift_deletion(chrom, start, len, chrom_len, fasta)?;
            if len == 1 {
                Ok(Some(format!("{acc}:g.{}del", s + 1)))
            } else {
                Ok(Some(format!("{acc}:g.{}_{}del", s + 1, s + len)))
            }
        }
        // Insertion (or duplication after the 3'-shift).
        ([], ins) => format_insertion(&acc, chrom, start, ins, chrom_len, fasta),
        // Deletion-insertion. Not 3'-shifted (the replacement fixes the site).
        (del, ins) => {
            let len = del.len() as u64;
            let alt_str = ascii_to_string(ins);
            if len == 1 {
                Ok(Some(format!("{acc}:g.{}delins{}", start + 1, alt_str)))
            } else {
                Ok(Some(format!(
                    "{acc}:g.{}_{}delins{}",
                    start + 1,
                    start + len,
                    alt_str
                )))
            }
        }
    }
}

/// Resolve a chromosome name to its RefSeq `NC_*` accession, or `None` for a
/// non-primary contig (patch/alt/unknown) that has no canonical `g.` accession.
fn resolve_accession(chrom: &str) -> Option<String> {
    if chrom.starts_with("NC_") {
        return Some(chrom.to_string());
    }
    let acc = chrom::ucsc_to_refseq(chrom);
    // `ucsc_to_refseq` returns its input unchanged when the name is not a
    // primary GRCh38 chromosome (patch/alt contig, bare Ensembl name, or
    // anything unknown) — those have no canonical `NC_*` accession.
    if acc == chrom {
        None
    } else {
        Some(acc.to_string())
    }
}

/// Trim shared suffix then shared prefix of `ref`/`alt` to the minimal edit,
/// returning the 0-based start of the trimmed edit and the trimmed alleles.
fn trim_to_minimal_edit(pos: u64, ref_u: &[u8], alt_u: &[u8]) -> (u64, Vec<u8>, Vec<u8>) {
    let mut r_end = ref_u.len();
    let mut a_end = alt_u.len();
    while r_end > 0 && a_end > 0 && ref_u[r_end - 1] == alt_u[a_end - 1] {
        r_end -= 1;
        a_end -= 1;
    }
    let mut lead = 0;
    while lead < r_end && lead < a_end && ref_u[lead] == alt_u[lead] {
        lead += 1;
    }
    (
        pos + lead as u64,
        ref_u[lead..r_end].to_vec(),
        alt_u[lead..a_end].to_vec(),
    )
}

/// Slide a deletion of `len` bases at 0-based `start` to its most 3' position.
///
/// A deletion `[s, s+len)` is equivalent to `[s+1, s+1+len)` iff the base
/// leaving the 5' edge equals the base entering from the 3' side. The window
/// is clamped to `chrom_len`, so a deletion at the contig end simply does not
/// shift (never an out-of-range error). Returns the new 0-based start.
fn shift_deletion(
    chrom: &str,
    start: u64,
    len: u64,
    chrom_len: u64,
    fasta: &FastaReader,
) -> Result<u64, VarEffectError> {
    let mut s = start;
    while s + len < chrom_len {
        if fasta.fetch_base(chrom, s)? != fasta.fetch_base(chrom, s + len)? {
            break;
        }
        s += 1;
    }
    Ok(s)
}

/// Format an insertion at 0-based `start`, applying the 3'-shift and choosing
/// `dup` over `ins` when the shifted sequence copies the preceding reference.
fn format_insertion(
    acc: &str,
    chrom: &str,
    start: u64,
    ins: &[u8],
    chrom_len: u64,
    fasta: &FastaReader,
) -> Result<Option<String>, VarEffectError> {
    // 3'-shift: the insertion point slides right through a repeat while the
    // first inserted base equals the reference base being crossed; the
    // inserted sequence rotates by one at each step.
    let mut p = start;
    let mut seq = ins.to_vec();
    while p < chrom_len && fasta.fetch_base(chrom, p)? == seq[0] {
        seq.rotate_left(1);
        p += 1;
    }

    // Duplication: the inserted sequence copies the equal-length reference
    // immediately 5' of the (shifted) insertion point.
    let m = seq.len() as u64;
    if p >= m && fasta.fetch_sequence(chrom, p - m, p)? == seq {
        return if m == 1 {
            Ok(Some(format!("{acc}:g.{p}dup")))
        } else {
            Ok(Some(format!("{acc}:g.{}_{}dup", p - m + 1, p)))
        };
    }

    // Plain insertion needs both flanking bases: 5' flank 1-based `p` (so
    // `p >= 1`) and 3' flank 1-based `p+1` (so `p < chrom_len`). A degenerate
    // contig-boundary insertion has no valid flank pair.
    if p == 0 || p >= chrom_len {
        return Ok(None);
    }
    Ok(Some(format!(
        "{acc}:g.{}_{}ins{}",
        p,
        p + 1,
        ascii_to_string(&seq)
    )))
}

/// True if every byte is an uppercase IUPAC base this formatter can express.
fn is_acgtn(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// Convert validated `ACGTN` ASCII bytes to a `String`.
fn ascii_to_string(bytes: &[u8]) -> String {
    // Invariant: callers pass bytes already validated by `is_acgtn`, so this
    // is always valid ASCII/UTF-8.
    String::from_utf8(bytes.to_vec()).expect("ACGTN alleles are valid ASCII")
}

#[cfg(test)]
mod tests {
    //! CI-runnable exact-output tests over a tiny synthetic genome (no 3.1 GB
    //! FASTA). The `GGGG` run in `CONTIG` exercises the 3'-shift and the
    //! duplication-vs-insertion decision.

    use super::*;
    use crate::fasta::{FastaReader, write_genome_binary};
    use tempfile::TempDir;

    /// Synthetic contig "1" (queried as "chr1"):
    ///   1-based pos: 1 2 3 4 5 6 7 8 9 10
    ///   base:        A A T G G G G T A A
    const CONTIG: &[u8] = b"AATGGGGTAA";

    fn fasta_with(contigs: &[(&str, &[u8])]) -> (TempDir, FastaReader) {
        let tmp = TempDir::new().expect("tempdir");
        let bin = tmp.path().join("g.bin");
        let idx = tmp.path().join("g.bin.idx");
        write_genome_binary(contigs, "test", &bin, &idx).expect("write synthetic genome");
        let fasta = FastaReader::open(&bin).expect("open synthetic genome");
        (tmp, fasta)
    }

    fn fmt(chrom: &str, pos: u64, r: &[u8], a: &[u8]) -> Option<String> {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        format_hgvs_g(chrom, pos, r, a, &fasta).expect("format_hgvs_g should not error")
    }

    #[test]
    fn substitution() {
        // 0-based 2 (T) > A.
        assert_eq!(
            fmt("chr1", 2, b"T", b"A"),
            Some("NC_000001.11:g.3T>A".to_string())
        );
    }

    #[test]
    fn substitution_is_case_insensitive() {
        assert_eq!(
            fmt("chr1", 2, b"t", b"a"),
            Some("NC_000001.11:g.3T>A".to_string())
        );
    }

    #[test]
    fn substitution_ref_genome_mismatch_yields_none() {
        // Genome has T at 0-based 2, not C.
        assert_eq!(fmt("chr1", 2, b"C", b"A"), None);
    }

    #[test]
    fn single_base_deletion_shifts_3prime() {
        // Left-aligned deletion of one G in the GGGG run (T>del, VCF TG>T at
        // 0-based 2) must land on the 3'-most G at 1-based 7.
        assert_eq!(
            fmt("chr1", 2, b"TG", b"T"),
            Some("NC_000001.11:g.7del".to_string())
        );
    }

    #[test]
    fn range_deletion_shifts_3prime() {
        // Delete two G's from the run: 3'-most range is 1-based 6_7.
        assert_eq!(
            fmt("chr1", 2, b"TGG", b"T"),
            Some("NC_000001.11:g.6_7del".to_string())
        );
    }

    #[test]
    fn duplication_in_homopolymer() {
        // Left-aligned 1 bp insertion into the GGGG run (T>TG at 0-based 2)
        // 3'-shifts to a duplication of the last G at 1-based 7.
        assert_eq!(
            fmt("chr1", 2, b"T", b"TG"),
            Some("NC_000001.11:g.7dup".to_string())
        );
    }

    #[test]
    fn plain_insertion_not_a_duplication() {
        // Insert C after the T at 0-based 7 — no repeat, no dup.
        assert_eq!(
            fmt("chr1", 7, b"T", b"TC"),
            Some("NC_000001.11:g.8_9insC".to_string())
        );
    }

    #[test]
    fn deletion_insertion() {
        // TG (0-based 2_3) > CC — a delins, not 3'-shifted.
        assert_eq!(
            fmt("chr1", 2, b"TG", b"CC"),
            Some("NC_000001.11:g.3_4delinsCC".to_string())
        );
    }

    #[test]
    fn single_ref_base_delins() {
        // T (0-based 2) > CC — a 1-base-ref delins.
        assert_eq!(
            fmt("chr1", 2, b"T", b"CC"),
            Some("NC_000001.11:g.3delinsCC".to_string())
        );
    }

    #[test]
    fn multibase_tandem_duplication() {
        // Contig with a GTGTGT tandem run (0-based 2..8). A left-aligned
        // insertion of a "GT" copy (A>AGT at 0-based 1) 3'-shifts to a range
        // duplication of the final GT unit at 1-based 7_8.
        let (_tmp, fasta) = fasta_with(&[("1", b"AAGTGTGTCC")]);
        assert_eq!(
            format_hgvs_g("chr1", 1, b"A", b"AGT", &fasta).unwrap(),
            Some("NC_000001.11:g.7_8dup".to_string())
        );
    }

    #[test]
    fn contig_end_deletion_clamps_without_error() {
        // Delete an A from the trailing "AA" (VCF AA>A at 0-based 8) — the
        // 3'-shift reaches the final base (1-based 10) and stops at the contig
        // end rather than erroring.
        assert_eq!(
            fmt("chr1", 8, b"AA", b"A"),
            Some("NC_000001.11:g.10del".to_string())
        );
    }

    #[test]
    fn insertion_before_first_base_yields_none() {
        // Right-anchored insertion of a non-matching base before 0-based 0
        // (ref A, alt GA) cannot shift and has no 5' flank → no canonical g.
        assert_eq!(fmt("chr1", 0, b"A", b"GA"), None);
    }

    #[test]
    fn symbolic_allele_yields_none() {
        assert_eq!(fmt("chr1", 2, b"T", b"<DEL>"), None);
    }

    #[test]
    fn patch_contig_yields_none() {
        // Non-primary contig has no canonical NC_ accession.
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        assert_eq!(
            format_hgvs_g("chrUn_KI270302v1", 5, b"A", b"T", &fasta).unwrap(),
            None
        );
    }

    #[test]
    fn accession_resolution_covers_chrom_edges() {
        assert_eq!(resolve_accession("chr17").as_deref(), Some("NC_000017.11"));
        assert_eq!(resolve_accession("chrX").as_deref(), Some("NC_000023.11"));
        // chrY is .10, not .11.
        assert_eq!(resolve_accession("chrY").as_deref(), Some("NC_000024.10"));
        // Mitochondrion.
        assert_eq!(resolve_accession("chrM").as_deref(), Some("NC_012920.1"));
        // Already an accession — passthrough.
        assert_eq!(
            resolve_accession("NC_000017.11").as_deref(),
            Some("NC_000017.11")
        );
        // Bare Ensembl name and patch contigs have no canonical accession.
        assert_eq!(resolve_accession("17"), None);
        assert_eq!(resolve_accession("chr9_KN196479v1_fix"), None);
    }
}
