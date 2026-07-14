//! Concordance spot-check for the HGVS g. (genomic) reverse mapper.
//!
//! Validates that [`vareffect::VarEffect::resolve_hgvs_g`] produces correct
//! **left-aligned, VCF-canonical** 0-based genomic coordinates on the real
//! GRCh38 genome, across all six change types (substitution, deletion,
//! duplication, insertion, delins, inversion). `#[ignore]`-gated because it
//! requires the reference genome binary on disk.
//!
//! Run with:
//! ```bash
//! FASTA_PATH=/abs/path/to/data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect --test vep_concordance_hgvs_g_reverse -- --ignored
//! ```
//!
//! # Ground-truth provenance
//!
//! Expected coordinates are the **left-aligned parsimonious** (bcftools-canonical)
//! form — the same normalization `left_align_indel` implements, and what
//! biocommons `hgvs_to_vcf` (`shuffle_direction=5`) / gnomAD emit.
//!
//! - **substitution / delins** do not shift; the expected values are the direct
//!   0-based genomic form of the HGVS coordinates.
//! - **inversion** is a pure sequence operation: `ref = genome span`,
//!   `alt = reverse_complement(span)`, verified directly against the genome
//!   window (non-affix-sharing span, so no trim).
//! - **deletion / duplication** land in repeat contexts and were left-aligned by
//!   an **independent oracle** — `bcftools norm -f GRCh38.fa` — NOT by the code
//!   under test and NOT from SPDI/NCBI Variation Services (whose maximally-
//!   expanded form differs). Feeding the raw HGVS-stated forms to `bcftools norm`
//!   produced (1-based): CFTR `g.117559592_117559594del` to chr7 117559590
//!   `ATCT`>`A`; CFTR `g.117559595dup` to chr7 117559592 `C`>`CT`; BRCA1
//!   `g.43057065dup` to chr17 43057062 `T`>`TG`.
//!
//! The BRCA1 dup differs from `resolve_hgvs_c`'s raw 0-based 43057062 `G`>`GG`:
//! `resolve_hgvs_g` left-aligns, `resolve_hgvs_c` does not.

use std::path::Path;

use vareffect::{FastaReader, TranscriptStore, VarEffect, VarEffectError};

/// Open the reference genome from `FASTA_PATH` and pair it with an empty
/// transcript store — `resolve_hgvs_g` only touches the FASTA.
fn load_ve() -> VarEffect {
    let path = std::env::var("FASTA_PATH").expect(
        "FASTA_PATH env var must point to a GRCh38 genome binary (data/vareffect/GRCh38.bin)",
    );
    let fasta = FastaReader::open(Path::new(&path))
        .unwrap_or_else(|e| panic!("failed to open FASTA at {path}: {e}"));
    VarEffect::new(TranscriptStore::from_transcripts(Vec::new()), fasta)
}

/// Expected output for one reverse-mapped genomic variant.
struct ExpectedReverseG {
    label: &'static str,
    /// Full HGVS g. notation including accession.
    hgvs_input: &'static str,
    /// Expected UCSC chromosome (e.g. `"chr17"`).
    expected_chrom: &'static str,
    /// Expected **left-aligned** 0-based genomic position.
    expected_pos: u64,
    /// Expected reference allele (plus-strand, anchor-prepended for indels).
    expected_ref: &'static [u8],
    /// Expected alternate allele (plus-strand).
    expected_alt: &'static [u8],
    /// If `true`, also assert the round-trip `format_hgvs_g(resolve) == input`.
    round_trip: bool,
}

const VARIANTS: &[ExpectedReverseG] = &[
    // Substitution — direct, no shift.
    ExpectedReverseG {
        label: "#1 TP53 substitution",
        hgvs_input: "NC_000017.11:g.7674220C>G",
        expected_chrom: "chr17",
        expected_pos: 7_674_219,
        expected_ref: b"C",
        expected_alt: b"G",
        round_trip: true,
    },
    // Delins — direct, no shift (no shared affix).
    ExpectedReverseG {
        label: "#2 TP53 delins",
        hgvs_input: "NC_000017.11:g.7674220delinsAG",
        expected_chrom: "chr17",
        expected_pos: 7_674_219,
        expected_ref: b"C",
        expected_alt: b"AG",
        round_trip: false,
    },
    // Deletion in a repeat — left-aligned by bcftools norm.
    ExpectedReverseG {
        label: "#3 CFTR F508del",
        hgvs_input: "NC_000007.14:g.117559592_117559594del",
        expected_chrom: "chr7",
        expected_pos: 117_559_589,
        expected_ref: b"ATCT",
        expected_alt: b"A",
        round_trip: false,
    },
    // Duplication in a homopolymer — left-aligned by bcftools norm.
    ExpectedReverseG {
        label: "#4 CFTR homopolymer dup",
        hgvs_input: "NC_000007.14:g.117559595dup",
        expected_chrom: "chr7",
        expected_pos: 117_559_591,
        expected_ref: b"C",
        expected_alt: b"CT",
        round_trip: false,
    },
    // Clinical duplication (BRCA1 c.5266dup) — left-aligned by bcftools norm.
    ExpectedReverseG {
        label: "#5 BRCA1 c.5266dup",
        hgvs_input: "NC_000017.11:g.43057065dup",
        expected_chrom: "chr17",
        expected_pos: 43_057_061,
        expected_ref: b"T",
        expected_alt: b"TG",
        round_trip: false,
    },
    // Inversion — ref = genome span, alt = reverse_complement(span).
    ExpectedReverseG {
        label: "#6 chr17 inversion",
        hgvs_input: "NC_000017.11:g.7674220_7674225inv",
        expected_chrom: "chr17",
        expected_pos: 7_674_219,
        expected_ref: b"CGGTTC",
        expected_alt: b"GAACCG",
        round_trip: false,
    },
];

#[test]
#[ignore = "requires GRCh38 genome binary via FASTA_PATH"]
fn resolve_hgvs_g_matches_ground_truth() {
    let ve = load_ve();
    let mut failures = Vec::new();

    for exp in VARIANTS {
        let gv = match ve.resolve_hgvs_g(exp.hgvs_input) {
            Ok(gv) => gv,
            Err(e) => {
                failures.push(format!("[{}] {} -> ERROR {e}", exp.label, exp.hgvs_input));
                continue;
            }
        };
        if gv.chrom != exp.expected_chrom
            || gv.pos != exp.expected_pos
            || gv.ref_allele != exp.expected_ref
            || gv.alt_allele != exp.expected_alt
        {
            failures.push(format!(
                "[{}] {}\n    expected {} {} {}>{}\n    got      {} {} {}>{}",
                exp.label,
                exp.hgvs_input,
                exp.expected_chrom,
                exp.expected_pos,
                String::from_utf8_lossy(exp.expected_ref),
                String::from_utf8_lossy(exp.expected_alt),
                gv.chrom,
                gv.pos,
                String::from_utf8_lossy(&gv.ref_allele),
                String::from_utf8_lossy(&gv.alt_allele),
            ));
            continue;
        }

        if exp.round_trip {
            let formatted = ve
                .format_hgvs_g(&gv.chrom, gv.pos, &gv.ref_allele, &gv.alt_allele)
                .expect("format_hgvs_g should not error");
            if formatted.as_deref() != Some(exp.hgvs_input) {
                failures.push(format!(
                    "[{}] round-trip: format_hgvs_g gave {:?}, expected {:?}",
                    exp.label, formatted, exp.hgvs_input
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "concordance failures:\n{}",
        failures.join("\n")
    );
}

/// The reject-classes must fail closed on the real genome too (no coordinate is
/// ever guessed for identity / imprecise / non-chromosomal input).
#[test]
#[ignore = "requires GRCh38 genome binary via FASTA_PATH"]
fn resolve_hgvs_g_rejects_unsupported() {
    let ve = load_ve();
    for hgvs in [
        "NC_000017.11:g.7674220=",             // identity
        "NC_000017.11:g.(7674220_7674225)del", // imprecise
        "NG_017013.2:g.5000del",               // non-chromosomal reference
        "LRG_199:g.5000del",                   // LRG reference
    ] {
        assert!(
            matches!(
                ve.resolve_hgvs_g(hgvs),
                Err(VarEffectError::HgvsParseError(_))
            ),
            "expected HgvsParseError for {hgvs}"
        );
    }
}
