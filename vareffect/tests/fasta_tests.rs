//! Integration tests for [`vareffect::FastaReader`] against a real GRCh38
//! reference genome binary.
//!
//! These tests are `#[ignore]`-gated because they require a ~3.1 GB flat
//! binary genome on disk. Run them explicitly with:
//!
//! ```bash
//! FASTA_PATH=data/vareffect/GRCh38.bin cargo test -p vareffect -- --ignored
//! ```
//!
//! The expected reference base / length values come from the GRCh38.p14
//! primary assembly (NCBI RefSeq). They are not expected to change across
//! patch releases unless the entire assembly is rebased to a new patch
//! level — only a GRCh38 -> GRCh39 switch would alter them, and vareffect
//! targets GRCh38 exclusively for now.

use std::path::{Path, PathBuf};

use vareffect::{FastaReader, VarEffectError};

/// Helper: read `FASTA_PATH` from the environment and open a reader.
///
/// Panics with a clear message if the env var is unset or the path is
/// invalid — integration tests are explicitly opt-in, so this is fine.
fn open_reader() -> FastaReader {
    let path = std::env::var("FASTA_PATH").expect(
        "FASTA_PATH env var must point to a GRCh38 genome binary (.bin) \
         with its .bin.idx sidecar. Run `vareffect-cli setup` first, then \
         set FASTA_PATH=data/vareffect/GRCh38.bin.",
    );
    let path_buf = PathBuf::from(path);
    FastaReader::open(Path::new(&path_buf)).expect("opening the reference genome binary")
}

/// TP53 c.742C>T (p.Arg248Trp) lives at chr17:7674221 (1-based VCF) which is
/// chr17:7674220 in 0-based half-open coordinates. TP53 is on the minus
/// strand, so the coding-strand reference allele is `C` (first base of the
/// CGG Arg codon) — the FASTA stores the plus-strand complement, `G`.
///
/// This test locks in a well-known clinical position so future refactors
/// can't silently change the offset translation.
#[test]
#[ignore]
fn fetch_tp53_codon_248_first_base() {
    let reader = open_reader();
    let base = reader.fetch_base("chr17", 7674220).unwrap();
    assert_eq!(
        base, b'G',
        "expected plus-strand reference base G at chr17:7674220 (TP53 c.742C>T)",
    );
}

/// Sanity check: fetch a 10 bp window from BRCA2 (chr13) and confirm the
/// result is the expected length and contains only IUPAC DNA bases.
#[test]
#[ignore]
fn fetch_brca2_sequence() {
    let reader = open_reader();
    let seq = reader
        .fetch_sequence("chr13", 32_340_300, 32_340_310)
        .unwrap();
    assert_eq!(seq.len(), 10, "expected 10 bp window");
    assert!(
        seq.iter()
            .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')),
        "expected only IUPAC DNA bases, got {:?}",
        String::from_utf8_lossy(&seq),
    );
}

/// `verify_ref` should return `true` for the real base and `false` for any
/// other base at the same position.
#[test]
#[ignore]
fn verify_ref_allele() {
    let reader = open_reader();
    // Same TP53 c.742C>T position as above — plus strand is G.
    assert!(reader.verify_ref("chr17", 7674220, b"G").unwrap());
    assert!(!reader.verify_ref("chr17", 7674220, b"T").unwrap());
    // Case-insensitive compare.
    assert!(reader.verify_ref("chr17", 7674220, b"g").unwrap());
}

/// GRCh38 chr1 length is a well-known constant (248,956,422 bp). If this
/// changes, either the FASTA is the wrong assembly or the `.fai` was built
/// from a different file.
#[test]
#[ignore]
fn chrom_length_chr1() {
    let reader = open_reader();
    assert_eq!(reader.chrom_length("chr1"), Some(248_956_422));
}

/// Mitochondrial contig length is 16,569 bp on GRCh38 (NCBI's
/// `NC_012920.1` is the revised Cambridge reference sequence, identical
/// to Ensembl's `MT`). This test exercises the primary-chrom translation
/// path: against an NCBI FASTA, `chrM` maps to `NC_012920.1` via the
/// `ucsc_to_refseq` const table; against a legacy Ensembl FASTA, it maps
/// to the bare `MT` name.
#[test]
#[ignore]
fn chrom_length_chrm() {
    let reader = open_reader();
    assert_eq!(reader.chrom_length("chrM"), Some(16_569));
}

/// A bogus chromosome name should surface as `ChromNotFound`, not a panic
/// or a generic I/O error.
#[test]
#[ignore]
fn unknown_chrom_returns_error() {
    let reader = open_reader();
    let err = reader.fetch_base("chrZZ", 0).unwrap_err();
    assert!(
        matches!(err, VarEffectError::ChromNotFound { .. }),
        "expected ChromNotFound, got {err:?}",
    );
}
