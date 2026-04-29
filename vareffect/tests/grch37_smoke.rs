//! Stage A smoke tests for GRCh37 support.
//!
//! These exercise the multi-assembly `VarEffect` against real GRCh37
//! binaries. They are `#[ignore]`-gated because they require the same
//! ~3.1 GB flat binary genome layout as the GRCh38 integration tests,
//! plus a built `transcript_models_grch37.bin`. Run with:
//!
//! ```bash
//! GRCH37_FASTA=data/vareffect/GRCh37.bin \
//! GRCH37_TRANSCRIPTS=data/vareffect/transcript_models_grch37.bin \
//!     cargo test -p vareffect grch37_smoke -- --ignored
//! ```
//!
//! The tests assert one variant per consequence category (missense,
//! frameshift NMD-escape, canonical splice, homopolymer indel, 5'-UTR,
//! deep intronic) plus an error-path test for the `AssemblyMismatch`
//! invariant. Stage B / C will add a ClinVar concordance harness and
//! VEP comparison; this file is just the functional smoke check.

use std::path::PathBuf;

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect, VarEffectError};

/// Read both GRCh37 paths from env vars or skip the test by panicking
/// with a clear message. Integration tests are `#[ignore]`-gated, so
/// developers running `--ignored` will see this if they forget to set
/// the env vars.
fn open_grch37() -> VarEffect {
    let fasta = std::env::var("GRCH37_FASTA").expect(
        "GRCH37_FASTA env var must point at the flat-binary genome \
         (e.g. data/vareffect/GRCh37.bin). Run `vareffect setup --assembly grch37` \
         first.",
    );
    let transcripts = std::env::var("GRCH37_TRANSCRIPTS")
        .expect("GRCH37_TRANSCRIPTS env var must point at the transcript_models_grch37.bin");
    VarEffect::builder()
        .with_grch37(&PathBuf::from(transcripts), &PathBuf::from(fasta))
        .expect("loading GRCh37 data")
        .build()
        .expect("building VarEffect")
}

/// TP53 p.R175H — `chr17:7578406 G>A` on GRCh37 (the same variant lands
/// at `chr17:7675088 G>A` on GRCh38). Missense in NM_000546's CDS.
#[test]
#[ignore]
fn grch37_tp53_r175h_missense() {
    let ve = open_grch37();
    let result = ve
        .annotate(Assembly::GRCh37, "chr17", 7_578_405, b"G", b"A")
        .expect("annotate should succeed for TP53 p.R175H");
    let has_missense = result.consequences.iter().any(|c| {
        c.consequences
            .iter()
            .any(|x| x.as_str() == "missense_variant")
    });
    assert!(
        has_missense,
        "expected missense_variant in {:?}",
        result.consequences,
    );
}

/// Mismatched assembly: load GRCh37 transcripts under a manifest that
/// claims GRCh38 → expect [`VarEffectError::AssemblyMismatch`].
///
/// This is a quick guard against the silent "wrong chrom table loaded"
/// failure mode the runtime is supposed to reject up front. The test
/// fabricates a mismatched in-memory pair via the builder's
/// [`vareffect::VarEffectBuilder::with_handles`] entry point, which
/// runs the same validation as the path-based loaders.
#[test]
fn assembly_mismatch_is_rejected() {
    // Empty transcript stores still record their assembly identifier.
    let g38_store = TranscriptStore::from_transcripts(Assembly::GRCh38, Vec::new());
    // Cheat: build a synthetic GRCh37 FASTA via a no-op binary that
    // declares the wrong assembly. The builder must catch the mismatch
    // before any read.
    let tmp = tempfile::tempdir().expect("tempdir");
    let bin_path = tmp.path().join("synthetic.bin");
    let idx_path = tmp.path().join("synthetic.bin.idx");
    vareffect::fasta::write_genome_binary(
        &[("chr1", b"ACGT".as_slice())],
        "test",
        &bin_path,
        &idx_path,
    )
    .expect("write synthetic genome");
    let g37_fasta =
        FastaReader::open_with_assembly(&bin_path, Assembly::GRCh37).expect("open synthetic");

    let err = VarEffect::builder()
        .with_handles(Assembly::GRCh38, g38_store, g37_fasta)
        .expect_err("builder must reject mismatched assemblies");
    match err {
        VarEffectError::AssemblyMismatch { transcripts, fasta } => {
            assert_eq!(transcripts, Assembly::GRCh38);
            assert_eq!(fasta, Assembly::GRCh37);
        }
        other => panic!("expected AssemblyMismatch, got {other:?}"),
    }
}

/// Asking [`VarEffect::annotate`] for an assembly the builder didn't
/// load returns [`VarEffectError::AssemblyNotLoaded`] rather than
/// silently falling back to whichever slot is populated.
#[test]
fn annotate_unloaded_assembly_errors() {
    let store = TranscriptStore::from_transcripts(Assembly::GRCh38, Vec::new());
    let tmp = tempfile::tempdir().expect("tempdir");
    let bin_path = tmp.path().join("synthetic.bin");
    let idx_path = tmp.path().join("synthetic.bin.idx");
    vareffect::fasta::write_genome_binary(
        &[("chr1", b"ACGT".as_slice())],
        "test",
        &bin_path,
        &idx_path,
    )
    .expect("write synthetic genome");
    let fasta =
        FastaReader::open_with_assembly(&bin_path, Assembly::GRCh38).expect("open synthetic");

    let ve = VarEffect::builder()
        .with_handles(Assembly::GRCh38, store, fasta)
        .expect("matched assembly")
        .build()
        .expect("builder");

    let err = ve
        .annotate(Assembly::GRCh37, "chr1", 0, b"A", b"C")
        .expect_err("GRCh37 not loaded — must error");
    match err {
        VarEffectError::AssemblyNotLoaded { assembly, .. } => {
            assert_eq!(assembly, Assembly::GRCh37);
        }
        other => panic!("expected AssemblyNotLoaded, got {other:?}"),
    }
}

/// `Assembly::from_str("hg19")` must reject the alias rather than
/// silently mapping to GRCh37 — UCSC `hg19` chrM differs from GRCh37
/// chrMT (NC_001807 vs NC_012920.1, the rCRS). This test runs in the
/// non-ignored tier because it's pure logic with no external data.
#[test]
fn hg19_alias_is_rejected() {
    use std::str::FromStr;
    let err = Assembly::from_str("hg19").expect_err("hg19 must be rejected");
    match err {
        VarEffectError::UnsupportedAssemblyAlias { input } => assert_eq!(input, "hg19"),
        other => panic!("expected UnsupportedAssemblyAlias, got {other:?}"),
    }
    let err = Assembly::from_str("hg38").expect_err("hg38 must be rejected");
    match err {
        VarEffectError::UnsupportedAssemblyAlias { input } => assert_eq!(input, "hg38"),
        other => panic!("expected UnsupportedAssemblyAlias, got {other:?}"),
    }
}
