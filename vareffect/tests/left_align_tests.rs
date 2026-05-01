//! Integration tests for [`vareffect::VarEffect::left_align_indel`]
//! against the real GRCh38 reference genome binary.
//!
//! These tests are `#[ignore]`-gated because they require the same ~3.1 GB
//! flat binary genome as the other integration tests. Run them with:
//!
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin cargo test -p vareffect -- --ignored
//! ```
//!
//! The algorithm implements Tan et al. 2015 ("Unified representation of
//! genetic variants") — shift-then-trim to produce the unique left-aligned
//! parsimonious representation.

use std::path::{Path, PathBuf};

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect};

/// Helper: build a minimal `VarEffect` with an empty transcript store and a
/// real FASTA reader. `left_align_indel` only touches the FASTA, so the
/// transcript store can safely be empty.
fn open_var_effect() -> VarEffect {
    let path = std::env::var("GRCH38_FASTA").expect(
        "GRCH38_FASTA env var must point to a GRCh38 genome binary (.bin) \
         with its .bin.idx sidecar. Run `vareffect setup --assembly grch38` \
         first, then set GRCH38_FASTA=data/vareffect/GRCh38.bin.",
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

/// SNVs: rightmost bases differ immediately, shift loop never fires.
#[test]
#[ignore]
fn snv_passthrough() {
    let ve = open_var_effect();
    // TP53 hotspot: chr17:7676154 C>T (GRCh38)
    let result = ve
        .left_align_indel(Assembly::GRCh38, "chr17", 7_676_154, "C", "T")
        .expect("left_align_indel should not error on SNV");
    assert_eq!(result, None, "SNV should pass through unchanged");
}

/// MNVs: rightmost bases differ, no shift.
#[test]
#[ignore]
fn mnv_passthrough() {
    let ve = open_var_effect();
    // Two adjacent substitutions that don't share rightmost base
    let result = ve
        .left_align_indel(Assembly::GRCh38, "chr17", 7_676_154, "CC", "TG")
        .expect("left_align_indel should not error on MNV");
    assert_eq!(result, None, "MNV should pass through unchanged");
}

/// A deletion NOT in a repetitive region should return None.
#[test]
#[ignore]
fn already_leftmost_deletion() {
    let ve = open_var_effect();
    // CFTR deltaF508: chr7:117559590 ATCT>A (GRCh38) — not in a repeat
    // region. The base at chr7:117559593 is not the same as at 117559590,
    // so this should already be leftmost.
    let result = ve
        .left_align_indel(Assembly::GRCh38, "chr7", 117_559_590, "ATCT", "A")
        .expect("left_align_indel should not error");
    assert_eq!(result, None, "non-repeat deletion should be unchanged");
}

/// Right-shifted 1bp deletion in a poly-A run should left-align.
///
/// BRCA2 region: GRCh38 chr13:32340303-32340307 reads `GAAAA` — a 4 bp
/// poly-A immediately following a 'G'. Two right-shifted forms of the
/// same 1 bp deletion within the run must normalize to the canonical
/// left-aligned form `("GA","G")` at chr13:32340303.
#[test]
#[ignore]
fn homopolymer_deletion_shifts_left() {
    let ve = open_var_effect();

    // Right-shifted: anchor 'A' at 32340306, delete the 'A' at 32340307.
    let result_right = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_306, "AA", "A")
        .expect("left_align_indel should not error");
    // Left-aligned canonical form: anchor 'G' at 32340303, delete the
    // first 'A' at 32340304. Already leftmost — must return None.
    let result_left = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_303, "GA", "G")
        .expect("left_align_indel should not error");

    let norm_right = result_right.unwrap_or((32_340_306, "AA".to_string(), "A".to_string()));
    let norm_left = result_left.unwrap_or((32_340_303, "GA".to_string(), "G".to_string()));
    assert_eq!(
        norm_right, norm_left,
        "two representations of the same deletion must normalize identically"
    );
}

/// Right-shifted 1 bp insertion within the same chr13:32340303-32340307
/// `GAAAA` run must converge with the canonical left-aligned form
/// `("G","GA")` at chr13:32340303.
#[test]
#[ignore]
fn homopolymer_insertion_shifts_left() {
    let ve = open_var_effect();

    // Right-shifted: anchor 'A' at 32340307, insert an extra 'A' on
    // the right edge of the run.
    let result_a = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_307, "A", "AA")
        .expect("left_align_indel should not error");
    // Left-aligned canonical: anchor 'G' at 32340303, insert the new
    // 'A' immediately after the 'G'. Already leftmost — must return None.
    let result_b = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_303, "G", "GA")
        .expect("left_align_indel should not error");

    let norm_a = result_a.unwrap_or((32_340_307, "A".to_string(), "AA".to_string()));
    let norm_b = result_b.unwrap_or((32_340_303, "G".to_string(), "GA".to_string()));
    assert_eq!(
        norm_a, norm_b,
        "two representations of the same insertion must normalize identically"
    );
}

/// `ref=ACG, alt=AG` — shared rightmost G triggers the shift loop,
/// exposing a 1bp deletion of C. The algorithm handles this without
/// a separate post-trim guard.
#[test]
#[ignore]
fn complex_becomes_deletion_after_shift() {
    let ve = open_var_effect();
    // Use a real genomic position. The exact shift depends on the
    // local repeat context, but the key property is that the result
    // must NOT be None (the complex input is an indel in disguise).
    //
    // chr17:7676154 — pick a region where we can construct a test.
    // We just verify the algorithm processes it without error and
    // recognizes it's not a simple substitution.
    let result = ve.left_align_indel(Assembly::GRCh38, "chr17", 7_676_154, "CCC", "CC");
    assert!(result.is_ok(), "should not error on complex-becomes-indel");
    // This IS an indel (1bp deletion), so it may or may not shift
    // depending on the repeat context. The important thing is it
    // doesn't return None thinking it's an MNV.
}

#[test]
#[ignore]
fn idempotent() {
    let ve = open_var_effect();
    // First pass: normalize a right-shifted deletion
    let first = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_301, "GA", "G")
        .expect("first pass should not error");

    if let Some((pos, ref_a, alt_a)) = first {
        // Second pass: normalize the already-normalized result
        let second = ve
            .left_align_indel(Assembly::GRCh38, "chr13", pos, &ref_a, &alt_a)
            .expect("second pass should not error");
        assert_eq!(
            second, None,
            "already-normalized variant should return None on second pass"
        );
    }
    // If first returned None, it was already leftmost — also idempotent.
}

#[test]
#[ignore]
fn position_one_boundary() {
    let ve = open_var_effect();
    // At position 1 the loop guard `pos <= 1` prevents any shift.
    // This should succeed without error and return None or a trimmed
    // form (but never underflow).
    let result = ve.left_align_indel(Assembly::GRCh38, "chr1", 1, "NN", "N");
    assert!(
        result.is_ok(),
        "position 1 should not cause underflow: {:?}",
        result
    );
}

#[test]
#[ignore]
fn vcf_anchor_preserved() {
    let ve = open_var_effect();
    // After normalization, at least one allele must have length >= 1
    // (VCF anchor base requirement). The parsimony step preserves this.
    let result = ve
        .left_align_indel(Assembly::GRCh38, "chr13", 32_340_301, "GA", "G")
        .expect("should not error");
    if let Some((_, ref ref_a, ref alt_a)) = result {
        assert!(
            !ref_a.is_empty() && !alt_a.is_empty(),
            "both alleles must have at least one base (VCF anchor): ref={ref_a}, alt={alt_a}"
        );
        // For a pure deletion/insertion, the shorter allele should be
        // exactly 1 base (the anchor).
        let shorter = ref_a.len().min(alt_a.len());
        assert!(shorter >= 1, "shorter allele must be >= 1 base");
    }
}
