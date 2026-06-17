//! VEP concordance spot-check for structural-variant consequence prediction.
//!
//! Each case was validated against the Ensembl VEP REST API
//! (`POST /vep/human/region?refseq=1`, GRCh38, queried 2026-06-17) on MANE
//! Select transcripts. `#[ignore]`-gated because it requires the transcript
//! store + reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! cargo test -p vareffect --test vep_concordance_sv -- --ignored
//! ```
//!
//! # Concordance scope
//!
//! The assertion compares the per-transcript SO term **set** vareffect emits
//! against VEP's, after excluding three documented divergences (see
//! `VEP_DIVERGENCES.md`):
//! - precise codon consequences for SV breakpoints (`stop_lost`, `start_lost`,
//!   `stop_gained`, `missense_variant`, `frameshift_variant`, splice terms) —
//!   the SV path is geometry-only, no codon analysis;
//! - `NMD_transcript_variant` (modelled as the `predicts_nmd` flag elsewhere);
//! - regulatory-feature terms (`regulatory_region_*`, `TFBS_*`) — no regulatory
//!   store.
//!
//! All cases use intervals/points with multi-kb margins from region boundaries,
//! so the 1-based-VCF vs 0-based-vareffect convention (a ±1 shift) cannot change
//! which regions are hit.

use std::path::Path;

use vareffect::{SvKind, VarEffect};

/// The SV operation under test.
enum Op {
    /// Interval SV: `annotate_interval(start, end, kind)`.
    Interval(SvKind, u64, u64),
    /// Breakend breakpoint at a single 0-based position.
    Breakend(u64),
    /// Symbolic insertion at a single 0-based position.
    Insertion(u64),
}

/// One VEP-validated structural-variant concordance case.
struct Case {
    /// Human-readable label for assertion failures.
    label: &'static str,
    /// UCSC-style chromosome.
    chrom: &'static str,
    /// SV operation.
    op: Op,
    /// MANE Select transcript whose annotation is asserted.
    transcript: &'static str,
    /// Expected SO term set (sorted), per VEP minus the documented divergences.
    expected: &'static [&'static str],
}

/// Load the GRCh38 store + FASTA from the workspace data directory.
fn load_var_effect() -> VarEffect {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = manifest_dir.parent().expect("workspace root");
    VarEffect::open(
        &root.join("data/vareffect/transcript_models.bin"),
        &root.join("data/vareffect/GRCh38.bin"),
    )
    .expect("open VarEffect (run `vareffect-cli setup` to provision data/vareffect/)")
}

const CASES: &[Case] = &[
    // --- TP53 (minus strand, NM_000546.6) ---------------------------------
    Case {
        label: "TP53 inversion, 5' partial",
        chrom: "chr17",
        op: Op::Interval(SvKind::Inversion, 7_680_000, 7_700_000),
        transcript: "NM_000546.6",
        expected: &["5_prime_UTR_variant", "intron_variant"],
    },
    Case {
        label: "TP53 inversion, whole CDS span (no ablation: balanced)",
        chrom: "chr17",
        op: Op::Interval(SvKind::Inversion, 7_668_421, 7_687_550),
        transcript: "NM_000546.6",
        expected: &[
            "3_prime_UTR_variant",
            "5_prime_UTR_variant",
            "coding_sequence_variant",
            "intron_variant",
        ],
    },
    Case {
        label: "TP53 deletion, 5' partial -> feature_truncation",
        chrom: "chr17",
        op: Op::Interval(SvKind::Deletion, 7_680_000, 7_700_000),
        transcript: "NM_000546.6",
        expected: &[
            "5_prime_UTR_variant",
            "feature_truncation",
            "intron_variant",
        ],
    },
    Case {
        label: "TP53 deletion, whole transcript -> ablation only",
        chrom: "chr17",
        op: Op::Interval(SvKind::Deletion, 7_660_000, 7_690_000),
        transcript: "NM_000546.6",
        expected: &["transcript_ablation"],
    },
    Case {
        label: "TP53 duplication, partial extending out -> no headline",
        chrom: "chr17",
        op: Op::Interval(SvKind::Duplication, 7_680_000, 7_700_000),
        transcript: "NM_000546.6",
        expected: &["5_prime_UTR_variant", "intron_variant"],
    },
    Case {
        label: "TP53 duplication, intragenic exonic -> feature_elongation",
        chrom: "chr17",
        op: Op::Interval(SvKind::Duplication, 7_670_000, 7_672_000),
        transcript: "NM_000546.6",
        expected: &[
            "coding_sequence_variant",
            "feature_elongation",
            "intron_variant",
        ],
    },
    Case {
        label: "TP53 deletion, intragenic exonic -> feature_truncation",
        chrom: "chr17",
        op: Op::Interval(SvKind::Deletion, 7_670_000, 7_672_000),
        transcript: "NM_000546.6",
        expected: &[
            "coding_sequence_variant",
            "feature_truncation",
            "intron_variant",
        ],
    },
    Case {
        label: "TP53 breakend in intron -> feature_truncation",
        chrom: "chr17",
        op: Op::Breakend(7_674_000),
        transcript: "NM_000546.6",
        expected: &["feature_truncation", "intron_variant"],
    },
    // --- PTEN (plus strand, NM_000314.8, chr10) ---------------------------
    Case {
        label: "PTEN deletion, 5' partial -> feature_truncation",
        chrom: "chr10",
        op: Op::Interval(SvKind::Deletion, 87_863_000, 87_900_000),
        transcript: "NM_000314.8",
        expected: &[
            "5_prime_UTR_variant",
            "coding_sequence_variant",
            "feature_truncation",
            "intron_variant",
        ],
    },
    Case {
        label: "PTEN inversion, 5' partial (no headline)",
        chrom: "chr10",
        op: Op::Interval(SvKind::Inversion, 87_863_000, 87_900_000),
        transcript: "NM_000314.8",
        expected: &[
            "5_prime_UTR_variant",
            "coding_sequence_variant",
            "intron_variant",
        ],
    },
    Case {
        label: "PTEN duplication, intragenic 5'UTR exon -> feature_elongation",
        chrom: "chr10",
        op: Op::Interval(SvKind::Duplication, 87_863_700, 87_864_200),
        transcript: "NM_000314.8",
        expected: &["5_prime_UTR_variant", "feature_elongation"],
    },
    Case {
        label: "PTEN deletion, intronic-only -> no truncation",
        chrom: "chr10",
        op: Op::Interval(SvKind::Deletion, 87_899_000, 87_901_000),
        transcript: "NM_000314.8",
        expected: &["intron_variant"],
    },
    Case {
        label: "PTEN duplication, intronic-only -> no elongation",
        chrom: "chr10",
        op: Op::Interval(SvKind::Duplication, 87_899_000, 87_901_000),
        transcript: "NM_000314.8",
        expected: &["intron_variant"],
    },
    Case {
        label: "PTEN breakend in intron -> feature_truncation",
        chrom: "chr10",
        op: Op::Breakend(87_960_000),
        transcript: "NM_000314.8",
        expected: &["feature_truncation", "intron_variant"],
    },
    Case {
        label: "PTEN insertion in deep intron -> positional only",
        chrom: "chr10",
        op: Op::Insertion(87_900_000),
        transcript: "NM_000314.8",
        expected: &["intron_variant"],
    },
    // --- BRCA1 (minus strand, NM_007294.4) --------------------------------
    // VEP additionally reports `stop_lost` here (the deletion removes the stop
    // codon); the geometry-only SV path does not compute codon consequences —
    // a documented divergence, so `stop_lost` is excluded from `expected`.
    Case {
        label: "BRCA1 deletion, 3' partial -> feature_truncation",
        chrom: "chr17",
        op: Op::Interval(SvKind::Deletion, 43_043_000, 43_050_000),
        transcript: "NM_007294.4",
        expected: &[
            "3_prime_UTR_variant",
            "coding_sequence_variant",
            "feature_truncation",
            "intron_variant",
        ],
    },
];

#[test]
#[ignore = "requires data/vareffect/{transcript_models.bin,GRCh38.bin} on disk"]
fn vep_concordance_sv() {
    let ve = load_var_effect();

    for case in CASES {
        let results = match case.op {
            Op::Interval(kind, start, end) => ve
                .annotate_interval(case.chrom, start, end, kind)
                .unwrap_or_else(|e| panic!("{}: annotate_interval failed: {e}", case.label)),
            Op::Breakend(pos) => ve.annotate_breakend(case.chrom, pos),
            Op::Insertion(pos) => ve.annotate_sv_insertion(case.chrom, pos),
        };

        let result = results
            .iter()
            .find(|r| r.transcript == case.transcript)
            .unwrap_or_else(|| {
                panic!(
                    "{}: transcript {} not found in {} results",
                    case.label,
                    case.transcript,
                    results.len()
                )
            });

        let mut got: Vec<&str> = result.consequences.iter().map(|c| c.as_str()).collect();
        got.sort_unstable();

        assert_eq!(
            got, case.expected,
            "{}: vareffect terms diverged from VEP",
            case.label
        );
    }
}
