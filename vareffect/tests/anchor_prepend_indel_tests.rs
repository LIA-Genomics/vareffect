//! Integration tests for [`vareffect::VarEffect::anchor_prepend_indel`]
//! against the real GRCh38 reference genome binary.
//!
//! These tests are `#[ignore]`-gated because they require the same ~3.1 GB
//! flat binary genome as [`fasta_tests`]. Run them explicitly with:
//!
//! ```bash
//! FASTA_PATH=data/vareffect/GRCh38.bin cargo test -p vareffect -- --ignored
//! ```
//!
//! All reference bases are derived from the GRCh38.p14 primary assembly
//! and cross-checked against the Ensembl VEP REST API (POLH frameshift,
//! BRCA1 insertion). If the reference genome changes, the expected anchor
//! bases here must be updated.

use std::path::{Path, PathBuf};

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect, VarEffectError};

/// Helper: build a minimal `VarEffect` with an empty transcript store and a
/// real FASTA reader. `anchor_prepend_indel` only touches the FASTA, so the
/// transcript store can safely be empty — the helper never queries it.
///
/// Panics with a clear message if `FASTA_PATH` is unset or the path is
/// invalid. Integration tests are explicitly opt-in, so this is fine.
fn open_var_effect() -> VarEffect {
    let path = std::env::var("FASTA_PATH").expect(
        "FASTA_PATH env var must point to a GRCh38 genome binary (.bin) \
         with its .bin.idx sidecar. Run `vareffect-cli setup` first, then \
         set FASTA_PATH=data/vareffect/GRCh38.bin.",
    );
    let path_buf = PathBuf::from(path);
    let fasta = FastaReader::open_with_assembly(Path::new(&path_buf), Assembly::GRCh38)
        .expect("opening the reference genome binary");
    let transcripts = TranscriptStore::from_transcripts(Assembly::GRCh38, Vec::new());
    VarEffect::builder()
        .with_handles(Assembly::GRCh38, transcripts, fasta)
        .expect("matching assemblies")
        .build()
        .expect("builder")
}

/// POLH deletion `NM_006772.2:c.1861_1862del` — VEP REST emits
/// `pos=33409450` with alleles `TG/-` on chr6. The expected VCF form is
/// `pos=33409449, ref="ATG", alt="A"` where the anchor base `A` lives at
/// 1-based `chr6:33409449` (0-based `33409448`).
#[test]
#[ignore]
fn anchor_prepend_deletion_polh() {
    let ve = open_var_effect();
    let (pos, ref_a, alt_a) = ve
        .anchor_prepend_indel(Assembly::GRCh38, "chr6", 33_409_450, "TG", "-")
        .unwrap()
        .expect("pure deletion should be normalized");
    assert_eq!(pos, 33_409_449);
    assert_eq!(ref_a, "ATG");
    assert_eq!(alt_a, "A");
}

/// BRCA1 dup `NM_007294.4:c.5266dupC` — VEP REST emits `pos=43057063`
/// with alleles `-/C` on chr17. The expected VCF form is
/// `pos=43057062, ref="T", alt="TC"` where the anchor base `T` lives at
/// 1-based `chr17:43057062` (0-based `43057061`).
#[test]
#[ignore]
fn anchor_prepend_insertion_brca1() {
    let ve = open_var_effect();
    let (pos, ref_a, alt_a) = ve
        .anchor_prepend_indel(Assembly::GRCh38, "chr17", 43_057_063, "-", "C")
        .unwrap()
        .expect("pure insertion should be normalized");
    assert_eq!(pos, 43_057_062);
    assert_eq!(ref_a, "T");
    assert_eq!(alt_a, "TC");
}

/// SNV (`TP53 c.742C>T`, plus-strand `G>A`) should pass through as
/// `Ok(None)` without touching the FASTA — the early return fires before
/// any base lookup.
#[test]
#[ignore]
fn anchor_prepend_snv_passthrough() {
    let ve = open_var_effect();
    assert!(
        ve.anchor_prepend_indel(Assembly::GRCh38, "chr17", 7_674_221, "G", "A")
            .unwrap()
            .is_none(),
    );
}

/// MNV / complex substitution (neither allele `"-"`) also passes through
/// as `Ok(None)` — the helper is strictly about anchor-prepending HGVS
/// placeholders, not about minimal-representation normalization.
#[test]
#[ignore]
fn anchor_prepend_complex_passthrough() {
    let ve = open_var_effect();
    assert!(
        ve.anchor_prepend_indel(Assembly::GRCh38, "chr17", 7_674_221, "GA", "TC")
            .unwrap()
            .is_none(),
    );
}

/// Guard: `pos_1based = 1` has no base immediately 5', so the helper
/// must surface `CoordinateOutOfRange` with the real `chrom_len` instead
/// of underflowing.
#[test]
#[ignore]
fn anchor_prepend_guard_pos_one() {
    let ve = open_var_effect();
    let err = ve
        .anchor_prepend_indel(Assembly::GRCh38, "chr1", 1, "TG", "-")
        .expect_err("pos=1 should fail out-of-range guard");
    match err {
        VarEffectError::CoordinateOutOfRange {
            chrom, chrom_len, ..
        } => {
            assert_eq!(chrom, "chr1");
            assert!(
                chrom_len > 0,
                "chrom_len should reflect the real chr1 length, got {chrom_len}",
            );
        }
        other => panic!("expected CoordinateOutOfRange, got {other:?}"),
    }
}

/// A bogus chromosome name must surface as `ChromNotFound`, whether via
/// the guard path (if `pos_1based < 2`) or via `fetch_base` (the common
/// indel case).
#[test]
#[ignore]
fn anchor_prepend_chrom_not_found() {
    let ve = open_var_effect();
    let err = ve
        .anchor_prepend_indel(Assembly::GRCh38, "chrZZ", 100, "TG", "-")
        .expect_err("bogus chrom should fail");
    assert!(
        matches!(err, VarEffectError::ChromNotFound { .. }),
        "expected ChromNotFound, got {err:?}",
    );
}
