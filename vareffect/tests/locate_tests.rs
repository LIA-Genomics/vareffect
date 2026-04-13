//! Integration tests for [`vareffect::locate_variant`] against real
//! transcript models from the MessagePack store.
//!
//! These tests are `#[ignore]`-gated because they require the transcript
//! store file on disk. Run them explicitly with:
//!
//! ```bash
//! cargo test -p vareffect -- --ignored locate
//! ```
//!
//! Expected values are cross-checked against Ensembl VEP (release 115)
//! output for the same variants.

use std::path::Path;

use vareffect::{TranscriptStore, VariantLocation, locate_variant};

/// Load the transcript store from the workspace data directory.
fn load_store() -> TranscriptStore {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let workspace_root = manifest_dir
        .parent()
        .expect("Could not determine workspace root from CARGO_MANIFEST_DIR");
    let path = workspace_root.join("data/vareffect/transcript_models.bin");
    TranscriptStore::load_from_path(&path).unwrap_or_else(|e| {
        panic!(
            "Failed to load transcript store from {}: {}. \
             Run `cargo run -p vareffect-cli -- build` first.",
            path.display(),
            e,
        )
    })
}

/// TP53 c.742C>T (p.Arg248Trp) — a hotspot missense mutation.
///
/// TP53 NM_000546.6 is minus-strand on chr17.
/// CDS offset 741 (0-based) -> codon_number = 741/3 + 1 = 248,
/// codon_position = 741 % 3 = 0 (first base of codon 248).
#[test]
#[ignore]
fn tp53_codon_248_location() {
    let store = load_store();
    let (tx, idx) = store
        .get_by_accession("NM_000546.6")
        .expect("NM_000546.6 not found in store");

    let loc = locate_variant("chr17", 7_674_220, 7_674_221, tx, idx).unwrap();
    match loc {
        VariantLocation::CdsExon {
            cds_offset,
            codon_number,
            codon_position,
            exon_index,
            ..
        } => {
            assert_eq!(cds_offset, 741, "CDS offset should be 741 (c.742, 0-based)");
            assert_eq!(codon_number, 248, "Should be codon 248 (Arg248)");
            assert_eq!(codon_position, 0, "First base of codon 248");
            assert_eq!(exon_index, 6, "Should be exon index 6 (exon 7/11)");
        }
        other => panic!("Expected CdsExon for TP53 c.742C>T, got {:?}", other),
    }
}

/// SYNGAP1 NM_006772.3 — verify transcript structure.
#[test]
#[ignore]
fn syngap1_transcript_structure() {
    let store = load_store();
    let (tx, _) = store
        .get_by_accession("NM_006772.3")
        .expect("NM_006772.3 not found in store");

    assert_eq!(tx.exon_count, 19, "SYNGAP1 should have 19 exons");
    assert_eq!(tx.gene_symbol, "SYNGAP1");
    assert!(!tx.cds_segments.is_empty());
}

/// BRCA1 NM_007294.4 — splice donor at a known intron boundary.
#[test]
#[ignore]
fn brca1_splice_donor() {
    let store = load_store();
    let (tx, idx) = store
        .get_by_accession("NM_007294.4")
        .expect("NM_007294.4 not found in store");

    assert_eq!(tx.exon_count, 23, "BRCA1 should have 23 exons");

    let donor_pos = tx.exons[0].genomic_start - 1;
    let loc = locate_variant("chr17", donor_pos, donor_pos + 1, tx, idx).unwrap();
    match loc {
        VariantLocation::SpliceDonor {
            intron_index,
            offset,
        } => {
            assert_eq!(intron_index, 0, "Should be intron 0");
            assert_eq!(offset, 1, "Should be donor +1");
        }
        other => panic!(
            "Expected SpliceDonor at pos {} for BRCA1, got {:?}",
            donor_pos, other,
        ),
    }
}

/// BRCA2 NM_000059.4 — CDS offset at the start of coding.
#[test]
#[ignore]
fn brca2_cds_start() {
    let store = load_store();
    let (tx, idx) = store
        .get_by_accession("NM_000059.4")
        .expect("NM_000059.4 not found in store");

    assert_eq!(tx.exon_count, 27, "BRCA2 should have 27 exons");

    let first_cds_pos = tx.cds_segments[0].genomic_start;
    let loc = locate_variant("chr13", first_cds_pos, first_cds_pos + 1, tx, idx).unwrap();
    match loc {
        VariantLocation::CdsExon {
            cds_offset,
            codon_number,
            codon_position,
            ..
        } => {
            assert_eq!(cds_offset, 0, "First CDS base should be offset 0");
            assert_eq!(codon_number, 1, "First codon");
            assert_eq!(codon_position, 0, "First position of first codon");
        }
        other => panic!(
            "Expected CdsExon at first CDS base of BRCA2, got {:?}",
            other,
        ),
    }
}
