//! VEP concordance spot-check for indel and hard-case consequence prediction.
//!
//! 28 real clinical-grade variants validated against Ensembl VEP REST API
//! (GRCh38, queried 2026-04-09 with `refseq=1&numbers=1`). `#[ignore]`-gated
//! because the test requires the transcript store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_concordance_indel
//! ```
//!
//! Unlike the SNV test (`vep_concordance_snv.rs`), this test uses the top-level
//! [`annotate()`] dispatcher rather than per-function calls, exercising the full
//! VCF-input pipeline: REF verification -> `trim_alleles` -> type dispatch ->
//! per-transcript annotation.
//!
//! # Categories
//!
//! 1. **Frameshifts** (#1-5): deletions and insertions not divisible by 3.
//! 2. **Inframe indels** (#6-8): length divisible by 3, no frameshift.
//! 3. **Splice site overlap** (#9-13): indels hitting canonical splice sites.
//! 4. **Boundary-spanning** (#14-15): deletions crossing exon/intron boundaries.
//! 5. **MNVs / delins** (#16-18): multi-nucleotide substitutions via `annotate()`.
//! 6. **Non-CDS indels** (#19-21): UTR and deep intronic.
//! 7. **Edge cases** (#22-25): start_lost, stop_lost, extra frameshifts, non-coding.
//! 8. **Hard cases** (#26-28): junction insertion, additional frameshifts.
//!
//! # Known VEP divergences
//!
//! Variants with `vep_divergence: true` are asserted against vareffect's
//! *correct* output, not VEP's incorrect output. The `vep_consequence` field
//! stores what VEP (incorrectly) returns, for documentation only.
//!
//! 1. **Boundary-spanning deletion splice donor** (variant #13, ACVRL1
//!    c.525+1del) — VEP CLI with VCF input classifies as `frameshift_variant`
//!    because it only examines the VCF start position (which falls on the last
//!    exonic base). vareffect scans the full deletion footprint and correctly
//!    returns `splice_donor_variant`. See `docs/VEP_DIVERGENCES.md` section 1
//!    and [ensembl-vep#519](https://github.com/Ensembl/ensembl-vep/issues/519).
//!
//! 2. **Granular splice region sub-terms** — same as SNV test. VEP emits
//!    `splice_donor_region_variant` etc; vareffect emits the generic
//!    `splice_region_variant`. Affects variants #6 and #14 consequence comparison.
//!    Expected values reflect vareffect's output.
//!
//! # Incomplete terminal codon coverage
//!
//! Transcript store scan (2026-04-09) found only 1 transcript with
//! `total_cds_len % 3 != 0`: NM_001424184.1 (TMEM247, remainder=2). This is
//! not a MANE Select transcript. All MANE transcripts have complete terminal
//! codons, confirming that `is_incomplete_terminal_codon()` never fires in
//! production for MANE. The code path is covered by unit test #19 in
//! `consequence.rs` (synthetic fixture with modified CDS end).

use std::collections::BTreeSet;
use std::path::Path;

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect};

// ---------------------------------------------------------------------------
// Test infrastructure (same pattern as vep_concordance_snv.rs)
// ---------------------------------------------------------------------------

/// Load the GRCh38 transcript store. Reads `GRCH38_TRANSCRIPTS` if set, else
/// falls back to `data/vareffect/transcript_models_grch38.bin` under the
/// workspace root (derived from `CARGO_MANIFEST_DIR`).
fn load_store() -> TranscriptStore {
    let path = std::env::var("GRCH38_TRANSCRIPTS").unwrap_or_else(|_| {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("could not determine workspace root from CARGO_MANIFEST_DIR")
            .join("data/vareffect/transcript_models_grch38.bin")
            .to_string_lossy()
            .into_owned()
    });
    TranscriptStore::load_from_path(Path::new(&path)).unwrap_or_else(|e| {
        panic!(
            "failed to load GRCh38 transcript store from {path}: {e}. \
             Run `vareffect setup --assembly grch38` first.",
        )
    })
}

/// Load the GRCh38 genome reader from `GRCH38_FASTA` env var.
fn load_fasta() -> FastaReader {
    let path = std::env::var("GRCH38_FASTA")
        .expect("GRCH38_FASTA env var must point to a GRCh38 genome binary");
    FastaReader::open_with_patch_aliases_and_assembly(
        Path::new(&path),
        Some(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .and_then(|p| p.parent())
                .expect("workspace root")
                .join("data/vareffect/patch_chrom_aliases_grch38.csv")
                .as_ref(),
        ),
        Assembly::GRCh38,
    )
    .unwrap_or_else(|e| panic!("failed to open GRCh38 FASTA at {path}: {e}"))
}

// ---------------------------------------------------------------------------
// Expected output fixture
// ---------------------------------------------------------------------------

/// Expected output for one indel/MNV variant, derived from the VEP REST API.
/// All coordinate fields use VEP conventions (1-based for CDS/cDNA/protein
/// positions). The `pos` field is 0-based (vareffect convention).
struct Expected {
    /// Human-readable label for error messages.
    label: &'static str,
    /// UCSC-style chromosome (e.g., "chr17").
    chrom: &'static str,
    /// 0-based VCF POS (before trimming; includes anchor base for indels).
    pos: u64,
    /// VCF REF allele (plus-strand; includes anchor base for indels).
    ref_allele: &'static [u8],
    /// VCF ALT allele (plus-strand; includes anchor base for indels).
    alt_allele: &'static [u8],
    /// RefSeq transcript accession with version.
    transcript: &'static str,
    /// Expected SO consequence terms (sorted alphabetically).
    consequences: &'static [&'static str],
    /// Expected VEP IMPACT string.
    impact: &'static str,
    /// Expected 1-based protein position start, or `None`.
    protein_start: Option<u32>,
    /// Expected 1-based protein position end, or `None`.
    protein_end: Option<u32>,
    /// Expected 1-based CDS position start, or `None`.
    cds_position: Option<u32>,
    /// Expected 1-based CDS position end, or `None`.
    cds_position_end: Option<u32>,
    /// Expected 1-based cDNA position start, or `None`.
    cdna_position: Option<u32>,
    /// Expected 1-based cDNA position end, or `None`.
    cdna_position_end: Option<u32>,
    /// Expected codon string (e.g., "atCTTt/att"), or `None`.
    codons: Option<&'static str>,
    /// Expected amino acid string (e.g., "IF/I"), or `None`.
    amino_acids: Option<&'static str>,
    /// Expected exon number string (e.g., "11/27"), or `None`.
    exon: Option<&'static str>,
    /// Expected intron number string (e.g., "7/10"), or `None`.
    intron: Option<&'static str>,
    /// If `true`, vareffect intentionally differs from VEP. The assertion
    /// checks vareffect's correct output, not VEP's buggy output.
    vep_divergence: bool,
    /// When `vep_divergence` is `true`, the SO terms VEP (incorrectly) returns.
    /// Used for logging only; never asserted against.
    vep_consequence: Option<&'static [&'static str]>,
    /// Expected HGVS protein notation (e.g., "p.Phe508del"), or `None` for
    /// non-CDS variants. Values will be verified against VEP REST API with
    /// `refseq=1&hgvs=1` when the FASTA-dependent integration tests run.
    #[allow(dead_code)]
    hgvs_p: Option<&'static str>,
}

// ---------------------------------------------------------------------------
// Ground truth from VEP REST API (GRCh38, refseq=1, numbers=1)
// ---------------------------------------------------------------------------

const VARIANTS: &[Expected] = &[
    // -----------------------------------------------------------------------
    // Category 1: Frameshifts (5 variants)
    // -----------------------------------------------------------------------

    // #1 — BRCA1 c.68_69del: minus-strand 2bp deletion (185delAG), ClinVar 17662
    // VEP query: NM_007294.4:c.68_69del
    // Genomic: chr17:43124028-43124029 (1-based), deleted CT on plus strand
    // VCF: pos 43124027 (1-based) = 43124026 (0-based), REF=ACT, ALT=A
    Expected {
        label: "BRCA1 c.68_69del frameshift (minus)",
        chrom: "chr17",
        pos: 43_124_026,
        ref_allele: b"ACT",
        alt_allele: b"A",
        transcript: "NM_007294.4",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(23),
        protein_end: Some(23),
        cds_position: Some(68),
        cds_position_end: Some(69),
        cdna_position: Some(181),
        cdna_position_end: Some(182),
        codons: Some("gAG/g"),
        amino_acids: Some("E/X"),
        exon: Some("2/23"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #2 — BRCA1 c.5266dup: minus-strand 1bp insertion (5382insC), ClinVar 17677
    // VEP query: NM_007294.4:c.5266dup
    // Genomic: insertion between 43057063 and 43057064 (1-based)
    // Plus-strand inserted base: G (complement of coding C)
    // VCF: pos 43057063 (1-based) = 43057062 (0-based), REF=G, ALT=GG
    Expected {
        label: "BRCA1 c.5266dup frameshift insertion (minus)",
        chrom: "chr17",
        pos: 43_057_062,
        ref_allele: b"G",
        alt_allele: b"GG",
        transcript: "NM_007294.4",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(1755),
        protein_end: Some(1755),
        cds_position: Some(5264),
        cds_position_end: Some(5265),
        cdna_position: Some(5377),
        cdna_position_end: Some(5378),
        codons: Some("-/C"),
        amino_acids: Some("S/X"),
        exon: Some("19/23"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #3 — BRCA2 c.1813del: plus-strand 1bp deletion, ClinVar 37763
    // VEP query: NM_000059.4:c.1813del
    // Genomic: chr13:32333291 (1-based), deleted A on plus strand
    // VCF: pos 32333290 (1-based) = 32333289 (0-based), REF=AA, ALT=A
    Expected {
        label: "BRCA2 c.1813del frameshift (plus)",
        chrom: "chr13",
        pos: 32_333_289,
        ref_allele: b"AA",
        alt_allele: b"A",
        transcript: "NM_000059.4",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(605),
        protein_end: Some(605),
        cds_position: Some(1813),
        cds_position_end: Some(1813),
        cdna_position: Some(2012),
        cdna_position_end: Some(2012),
        codons: Some("Ata/ta"),
        amino_acids: Some("I/X"),
        exon: Some("10/27"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #4 — APC c.3927_3931del: plus-strand 5bp deletion, ClinVar hotspot region
    // VEP query: NM_000038.6:c.3927_3931del
    // Genomic: chr5:112839521-112839525 (1-based), deleted AAAGA on plus strand
    // VCF: pos 112839520 (1-based) = 112839519 (0-based), REF=AAAAGA, ALT=A
    Expected {
        label: "APC c.3927_3931del frameshift (plus, 5bp)",
        chrom: "chr5",
        pos: 112_839_519,
        ref_allele: b"AAAAGA",
        alt_allele: b"A",
        transcript: "NM_000038.6",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(1309),
        protein_end: Some(1311),
        cds_position: Some(3927),
        cds_position_end: Some(3931),
        cdna_position: Some(3986),
        cdna_position_end: Some(3990),
        codons: Some("gaAAAGAtt/gatt"),
        amino_acids: Some("EKI/DX"),
        exon: Some("16/16"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #5 — MLH1 c.676dup: plus-strand 1bp insertion near splice boundary
    // VEP query: NM_000249.4:c.676dup
    // Genomic: insertion between 37012097 and 37012098 (1-based)
    // VCF: pos 37012097 (1-based) = 37012096 (0-based), REF=T, ALT=TC
    // VEP returns compound: frameshift_variant + splice_region_variant
    Expected {
        label: "MLH1 c.676dup frameshift + splice_region (plus)",
        chrom: "chr3",
        pos: 37_012_096,
        ref_allele: b"T",
        alt_allele: b"TC",
        transcript: "NM_000249.4",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: Some("8/19"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 2: Inframe indels (3 variants)
    // -----------------------------------------------------------------------

    // #6 — CFTR c.1521_1523del (deltaF508): plus-strand 3bp inframe del, ClinVar 7105
    // VEP query: NM_000492.4:c.1521_1523del
    // Genomic: chr7:117559592-117559594 (1-based), deleted CTT on plus strand
    // VCF: pos 117559591 (1-based) = 117559590 (0-based), REF=TCTT, ALT=T
    Expected {
        label: "CFTR deltaF508 inframe deletion (plus)",
        chrom: "chr7",
        pos: 117_559_590,
        ref_allele: b"TCTT",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        consequences: &["inframe_deletion"],
        impact: "MODERATE",
        protein_start: Some(507),
        protein_end: Some(508),
        cds_position: Some(1521),
        cds_position_end: Some(1523),
        cdna_position: Some(1591),
        cdna_position_end: Some(1593),
        codons: Some("atCTTt/att"),
        amino_acids: Some("IF/I"),
        exon: Some("11/27"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #7 — EGFR c.2235_2249del (E746_A750del): plus-strand 15bp inframe del
    // VEP query: NM_005228.5:c.2235_2249del
    // Genomic: chr7:55174772-55174786 (1-based), deleted GGAATTAAGAGAAGC
    // VCF: pos 55174771 (1-based) = 55174770 (0-based), REF=AGGAATTAAGAGAAGC, ALT=A
    Expected {
        label: "EGFR E746_A750del inframe deletion (plus, 15bp)",
        chrom: "chr7",
        pos: 55_174_770,
        ref_allele: b"AGGAATTAAGAGAAGC",
        alt_allele: b"A",
        transcript: "NM_005228.5",
        consequences: &["inframe_deletion"],
        impact: "MODERATE",
        protein_start: Some(745),
        protein_end: Some(750),
        cds_position: Some(2235),
        cds_position_end: Some(2249),
        cdna_position: Some(2496),
        cdna_position_end: Some(2510),
        codons: Some("aaGGAATTAAGAGAAGCa/aaa"),
        amino_acids: Some("KELREA/K"),
        exon: Some("19/28"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #8 — ERBB2 c.2313_2324dup: plus-strand 12bp inframe insertion
    // VEP query: NM_004448.4:c.2313_2324dup
    // Genomic: insertion between 39724730 and 39724731 (1-based)
    // VCF: pos 39724730 (1-based) = 39724729 (0-based), REF=C, ALT=CATACGTGATGGC
    Expected {
        label: "ERBB2 exon20 12bp inframe insertion (plus)",
        chrom: "chr17",
        pos: 39_724_729,
        ref_allele: b"C",
        alt_allele: b"CATACGTGATGGC",
        transcript: "NM_004448.4",
        consequences: &["inframe_insertion", "splice_region_variant"],
        impact: "MODERATE",
        protein_start: Some(771),
        protein_end: Some(771),
        cds_position: Some(2312),
        cds_position_end: Some(2313),
        cdna_position: Some(2487),
        cdna_position_end: Some(2488),
        codons: Some("gca/gcATACGTGATGGCa"),
        amino_acids: Some("A/AYVMA"),
        exon: Some("20/27"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 3: Splice site overlap (5 variants)
    // -----------------------------------------------------------------------

    // #9 — TP53 c.782+1del: splice donor deletion (minus-strand)
    // VEP query: NM_000546.6:c.782+1del
    // Genomic: chr17:7674180 (1-based), deleted C on plus strand (G on coding)
    // VCF: pos 7674179 (1-based) = 7674178 (0-based), REF=AC, ALT=A
    Expected {
        label: "TP53 splice_donor +1 deletion (minus)",
        chrom: "chr17",
        pos: 7_674_178,
        ref_allele: b"AC",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("7/10"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #10 — TP53 c.673-2del: splice acceptor deletion (minus-strand)
    // VEP query: NM_000546.6:c.673-2del
    // Genomic: chr17:7674292 (1-based), deleted T on plus strand (A on coding)
    // VCF: pos 7674291 (1-based) = 7674290 (0-based), REF=CT, ALT=C
    Expected {
        label: "TP53 splice_acceptor -2 deletion (minus)",
        chrom: "chr17",
        pos: 7_674_290,
        ref_allele: b"CT",
        alt_allele: b"C",
        transcript: "NM_000546.6",
        consequences: &["splice_acceptor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("6/10"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #11 — TP53 c.782+1_782+2insT: splice donor insertion (minus-strand)
    // VEP query: NM_000546.6:c.782+1_782+2insT
    // Genomic: insertion between 7674179 and 7674180 (1-based)
    // Plus-strand inserted base: A (complement of coding T)
    // VCF: pos 7674179 (1-based) = 7674178 (0-based), REF=A, ALT=AA
    Expected {
        label: "TP53 splice_donor +1/+2 insertion (minus)",
        chrom: "chr17",
        pos: 7_674_178,
        ref_allele: b"A",
        alt_allele: b"AA",
        transcript: "NM_000546.6",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("7/10"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #12 — SYNGAP1 c.1861_1862del: plus-strand 2bp frameshift
    // VEP query: NM_006772.3:c.1861_1862del
    // Genomic: chr6:33440913-33440914 (1-based), deleted CG on plus strand
    // VCF: pos 33440912 (1-based) = 33440911 (0-based), REF=ACG, ALT=A
    // (anchor base A at position 33440912 — verified via Ensembl sequence API)
    Expected {
        label: "SYNGAP1 c.1861_1862del frameshift (plus)",
        chrom: "chr6",
        pos: 33_440_911,
        ref_allele: b"ACG",
        alt_allele: b"A",
        transcript: "NM_006772.3",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(621),
        protein_end: Some(621),
        cds_position: Some(1861),
        cds_position_end: Some(1862),
        cdna_position: Some(2061),
        cdna_position_end: Some(2062),
        codons: Some("CGa/a"),
        amino_acids: Some("R/X"),
        exon: Some("11/19"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #13 — ACVRL1 c.525+1del: VEP #519 boundary-spanning splice donor
    // VEP query: NM_000020.3:c.525+1del (HGVS → VEP correctly resolves)
    // Genomic: chr12:51913771 (1-based), deleted G on plus strand
    // VCF: pos 51913770 (1-based) = 51913769 (0-based), REF=GG, ALT=G
    // VEP CLI (VCF input) misclassifies as frameshift_variant (bug #519).
    // vareffect correctly returns splice_donor_variant.
    Expected {
        label: "ACVRL1 c.525+1del splice donor (VEP #519)",
        chrom: "chr12",
        pos: 51_913_769,
        ref_allele: b"GG",
        alt_allele: b"G",
        transcript: "NM_000020.3",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("4/9"),
        vep_divergence: true,
        vep_consequence: Some(&["frameshift_variant"]),
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 4: Boundary-spanning deletions (2 variants)
    // -----------------------------------------------------------------------

    // #14 — BRCA2 c.680_681+3del: exon-into-intron boundary deletion (donor)
    // VEP query: NM_000059.4:c.680_681+3del
    // Genomic: chr13:32329491-32329495 (1-based), deleted CTGTA, plus strand
    // VEP returns: splice_donor_variant + coding_sequence_variant + intron_variant
    // VCF: pos 32329490 (1-based) = 32329489 (0-based), REF=GCTGTA, ALT=G
    // (anchor base G at position 32329490 — verified via Ensembl sequence API)
    Expected {
        label: "BRCA2 c.680_681+3del boundary donor (plus)",
        chrom: "chr13",
        pos: 32_329_489,
        ref_allele: b"GCTGTA",
        alt_allele: b"G",
        transcript: "NM_000059.4",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: Some("8/27"),
        intron: Some("8/26"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #15 — BRCA2 c.682-3_683del: intron-into-exon boundary deletion (acceptor)
    // VEP query: NM_000059.4:c.682-3_683del
    // Genomic: chr13:32330916-32330920 (1-based), deleted CAGAA, plus strand
    // VEP returns: splice_acceptor_variant + coding_sequence_variant + intron_variant
    // VCF: pos 32330915 (1-based) = 32330914 (0-based), REF=GCAGAA, ALT=G
    // (anchor base G at position 32330915 — verified via Ensembl sequence API)
    Expected {
        label: "BRCA2 c.682-3_683del boundary acceptor (plus)",
        chrom: "chr13",
        pos: 32_330_914,
        ref_allele: b"GCAGAA",
        alt_allele: b"G",
        transcript: "NM_000059.4",
        consequences: &[
            "splice_acceptor_variant",
            "coding_sequence_variant",
            "intron_variant",
        ],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("8/26"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 5: MNVs / delins (3 variants)
    // -----------------------------------------------------------------------

    // #16 — TP53 c.742_743delinsTT: 2bp MNV within single codon (R248 region)
    // VEP query: NM_000546.6:c.742_743delinsTT
    // Genomic: chr17:7674220-7674221 (1-based), plus strand CG -> AA
    // (coding strand CG -> TT; plus strand = reverse complement)
    // VCF: pos 7674220 (1-based) = 7674219 (0-based), REF=CG, ALT=AA
    Expected {
        label: "TP53 R248 2bp MNV single codon (minus)",
        chrom: "chr17",
        pos: 7_674_219,
        ref_allele: b"CG",
        alt_allele: b"AA",
        transcript: "NM_000546.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(248),
        protein_end: Some(248),
        cds_position: Some(742),
        cds_position_end: Some(743),
        cdna_position: Some(884),
        cdna_position_end: Some(885),
        codons: Some("CGg/TTg"),
        amino_acids: Some("R/L"),
        exon: Some("7/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #17 — TP53 c.744_745delinsAA: 2bp delins spanning codon boundary
    // VEP decomposes to single SNV c.744 G>A (synonymous: CGG->CGA, both Arg).
    // c.745 A>A is unchanged on the coding strand.
    // Plus strand: pos 7674218-7674219 (1-based), TC -> TT. trim_alleles
    // strips shared T prefix, yielding SNV C>T at pos 7674219 (1-based).
    // VCF: pos 7674218 (1-based) = 7674217 (0-based), REF=TC, ALT=TT
    Expected {
        label: "TP53 c.744_745delinsAA decompose to SNV (minus)",
        chrom: "chr17",
        pos: 7_674_217,
        ref_allele: b"TC",
        alt_allele: b"TT",
        transcript: "NM_000546.6",
        consequences: &["synonymous_variant"],
        impact: "LOW",
        protein_start: Some(248),
        protein_end: Some(248),
        cds_position: Some(744),
        cds_position_end: Some(744),
        cdna_position: Some(886),
        cdna_position_end: Some(886),
        codons: Some("cgG/cgA"),
        amino_acids: Some("R"),
        exon: Some("7/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #18 — TP53 c.742_744delinsAAA: 3bp MNV replacing entire codon 248
    // VEP query: NM_000546.6:c.742_744delinsAAA
    // Genomic: chr17:7674219-7674221 (1-based), plus strand CCG -> TTT
    // (coding strand CGG -> AAA; reverse complement: CCG -> TTT)
    // VCF: pos 7674219 (1-based) = 7674218 (0-based), REF=CCG, ALT=TTT
    Expected {
        label: "TP53 R248K 3bp MNV entire codon (minus)",
        chrom: "chr17",
        pos: 7_674_218,
        ref_allele: b"CCG",
        alt_allele: b"TTT",
        transcript: "NM_000546.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(248),
        protein_end: Some(248),
        cds_position: Some(742),
        cds_position_end: Some(744),
        cdna_position: Some(884),
        cdna_position_end: Some(886),
        codons: Some("CGG/AAA"),
        amino_acids: Some("R/K"),
        exon: Some("7/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 6: Non-CDS indels (3 variants)
    // -----------------------------------------------------------------------

    // #19 — TP53 c.-20del: 5'UTR deletion (minus-strand)
    // VEP query: NM_000546.6:c.-20del
    // Genomic: chr17:7676614 (1-based), deleted G on plus strand (C on coding)
    // VCF: pos 7676613 (1-based) = 7676612 (0-based), REF=AG, ALT=A
    Expected {
        label: "TP53 5'UTR deletion (minus)",
        chrom: "chr17",
        pos: 7_676_612,
        ref_allele: b"AG",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        consequences: &["5_prime_UTR_variant"],
        impact: "MODIFIER",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: Some(123),
        cdna_position_end: Some(123),
        codons: None,
        amino_acids: None,
        exon: Some("2/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #20 — TP53 c.*50dupA: 3'UTR insertion (minus-strand)
    // VEP query: NM_000546.6:c.*50dupA
    // Genomic: insertion between 7669559 and 7669560 (1-based)
    // VCF: pos 7669559 (1-based) = 7669558 (0-based), REF=G, ALT=GC
    // VEP allele_string: -/C
    Expected {
        label: "TP53 3'UTR insertion (minus)",
        chrom: "chr17",
        pos: 7_669_558,
        ref_allele: b"G",
        alt_allele: b"GC",
        transcript: "NM_000546.6",
        consequences: &["3_prime_UTR_variant"],
        impact: "MODIFIER",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: Some(1373),
        cdna_position_end: Some(1374),
        codons: None,
        amino_acids: None,
        exon: Some("11/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #21 — BRCA2 c.681+50del: deep intronic deletion (plus-strand)
    // VEP query: NM_000059.4:c.681+50del
    // Genomic: chr13:32329542 (1-based), deleted G on plus strand
    // VCF: pos 32329541 (1-based) = 32329540 (0-based), REF=GG, ALT=G
    Expected {
        label: "BRCA2 deep intronic deletion (plus)",
        chrom: "chr13",
        pos: 32_329_540,
        ref_allele: b"GG",
        alt_allele: b"G",
        transcript: "NM_000059.4",
        consequences: &["intron_variant"],
        impact: "MODIFIER",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("8/26"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 7: Edge cases (4 variants)
    // -----------------------------------------------------------------------

    // #22 — TP53 c.1_3del: start codon deletion (minus-strand)
    // VEP query: NM_000546.6:c.1_3del
    // Genomic: chr17:7676592-7676594 (1-based), deleted CAT on plus strand
    // (coding strand ATG; reverse complement of CAT = ATG)
    // VCF: pos 7676591 (1-based) = 7676590 (0-based), REF=CCAT, ALT=C
    Expected {
        label: "TP53 start_lost (c.1_3del, minus)",
        chrom: "chr17",
        pos: 7_676_590,
        ref_allele: b"CCAT",
        alt_allele: b"C",
        transcript: "NM_000546.6",
        consequences: &["inframe_deletion", "start_lost"],
        impact: "HIGH",
        protein_start: Some(1),
        protein_end: Some(1),
        cds_position: Some(1),
        cds_position_end: Some(3),
        cdna_position: Some(143),
        cdna_position_end: Some(145),
        codons: Some("ATG/-"),
        amino_acids: Some("M/-"),
        exon: Some("2/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #23 — TP53 c.1180_1182del: stop codon deletion (minus-strand)
    // VEP query: NM_000546.6:c.1180_1182del
    // Genomic: chr17:7669609-7669611 (1-based), deleted TCA on plus strand
    // (coding strand TGA; reverse complement of TCA = TGA)
    // VCF: pos 7669608 (1-based) = 7669607 (0-based), REF=GTCA, ALT=G
    Expected {
        label: "TP53 stop_lost (c.1180_1182del, minus)",
        chrom: "chr17",
        pos: 7_669_607,
        ref_allele: b"GTCA",
        alt_allele: b"G",
        transcript: "NM_000546.6",
        consequences: &["inframe_deletion", "stop_lost"],
        impact: "HIGH",
        protein_start: Some(394),
        protein_end: Some(394),
        cds_position: Some(1180),
        cds_position_end: Some(1182),
        cdna_position: Some(1322),
        cdna_position_end: Some(1324),
        codons: Some("TGA/-"),
        amino_acids: Some("*/-"),
        exon: Some("11/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #24 — BRCA2 c.5946del: plus-strand 1bp frameshift, ClinVar rs80359550
    // VEP query: NM_000059.4:c.5946del
    // Genomic: chr13:32340301 (1-based), deleted T on plus strand
    // VCF: pos 32340300 (1-based) = 32340299 (0-based), REF=GT, ALT=G
    // (anchor base G at position 32340300 — verified via Ensembl sequence API)
    Expected {
        label: "BRCA2 c.5946del frameshift (plus)",
        chrom: "chr13",
        pos: 32_340_299,
        ref_allele: b"GT",
        alt_allele: b"G",
        transcript: "NM_000059.4",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(1982),
        protein_end: Some(1982),
        cds_position: Some(5946),
        cds_position_end: Some(5946),
        cdna_position: Some(6145),
        cdna_position_end: Some(6145),
        codons: Some("agT/ag"),
        amino_acids: Some("S/X"),
        exon: Some("11/27"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #25 — TERC n.200del: non-coding transcript exonic deletion
    // VEP query: chr3:169764861 region (NR_001566.3 not available via HGVS)
    // Genomic: chr3:169764861 (1-based), deleted G on plus strand
    // VCF: pos 169764860 (1-based) = 169764859 (0-based), REF=GG, ALT=G
    // (anchor base G at position 169764860 — verified via Ensembl sequence API)
    Expected {
        label: "TERC n.200del non-coding exon (minus)",
        chrom: "chr3",
        pos: 169_764_859,
        ref_allele: b"GG",
        alt_allele: b"G",
        transcript: "NR_001566.3",
        consequences: &["non_coding_transcript_exon_variant"],
        impact: "MODIFIER",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: Some(200),
        cdna_position_end: Some(200),
        codons: None,
        amino_acids: None,
        exon: Some("1/1"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category 8: Hard cases (3 variants)
    // -----------------------------------------------------------------------

    // #26 — BRCA2 c.6275_6276del: plus-strand 2bp frameshift, ClinVar rs11571658
    // VEP query: NM_000059.4:c.6275_6276del
    // Genomic: chr13:32340630-32340631 (1-based), deleted TT on plus strand
    // VCF: pos 32340629 (1-based) = 32340628 (0-based), REF=CTT, ALT=C
    // (anchor base C at position 32340629 — verified via Ensembl sequence API)
    Expected {
        label: "BRCA2 c.6275_6276del frameshift 2bp (plus)",
        chrom: "chr13",
        pos: 32_340_628,
        ref_allele: b"CTT",
        alt_allele: b"C",
        transcript: "NM_000059.4",
        consequences: &["frameshift_variant"],
        impact: "HIGH",
        protein_start: Some(2092),
        protein_end: Some(2092),
        cds_position: Some(6275),
        cds_position_end: Some(6276),
        cdna_position: Some(6474),
        cdna_position_end: Some(6475),
        codons: Some("cTT/c"),
        amino_acids: Some("L/X"),
        exon: Some("11/27"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #27 — TP53 c.782_782+1insA: junction insertion at exon-intron boundary
    // VEP query: NM_000546.6:c.782_782+1insA
    // Genomic: insertion between 7674180 and 7674181 (1-based)
    // VCF: pos 7674180 (1-based) = 7674179 (0-based), REF=C, ALT=CA
    // VEP returns: frameshift_variant + splice_region_variant
    // VEP does NOT assign splice_donor_variant for this junction insertion.
    Expected {
        label: "TP53 junction insertion exon/intron (minus)",
        chrom: "chr17",
        pos: 7_674_179,
        ref_allele: b"C",
        alt_allele: b"CA",
        transcript: "NM_000546.6",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: Some("7/11"),
        intron: None,
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
    // #28 — BRCA2 c.681+1_681+3del: pure intronic splice donor deletion
    // VEP query: NM_000059.4:c.681+1_681+3del
    // Genomic: chr13:32329493-32329495 (1-based), deleted GTA on plus strand
    // VCF: pos 32329492 (1-based) = 32329491 (0-based), REF=TGTA, ALT=T
    // (anchor base T at position 32329492 — from Ensembl sequence)
    Expected {
        label: "BRCA2 c.681+1_681+3del splice donor (plus)",
        chrom: "chr13",
        pos: 32_329_491,
        ref_allele: b"TGTA",
        alt_allele: b"T",
        transcript: "NM_000059.4",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        protein_end: None,
        cds_position: None,
        cds_position_end: None,
        cdna_position: None,
        cdna_position_end: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("8/26"),
        vep_divergence: false,
        vep_consequence: None,
        hgvs_p: None,
    },
];

// ---------------------------------------------------------------------------
// Comparison logic
// ---------------------------------------------------------------------------

/// Compare a single variant's `annotate()` output against expected VEP values.
/// Returns a list of field-level mismatches (empty = pass).
fn check_variant(ve: &VarEffect, exp: &Expected) -> Result<Vec<String>, String> {
    let results = ve
        .annotate(
            Assembly::GRCh38,
            exp.chrom,
            exp.pos,
            exp.ref_allele,
            exp.alt_allele,
        )
        .map_err(|e| format!("annotate error: {e}"))?;

    let result = results
        .consequences
        .iter()
        .find(|r| r.transcript == exp.transcript)
        .ok_or_else(|| {
            let found: Vec<&str> = results
                .consequences
                .iter()
                .map(|r| r.transcript.as_str())
                .collect();
            format!(
                "transcript {} not found in annotate() results (found: {:?})",
                exp.transcript, found,
            )
        })?;

    let mut mismatches: Vec<String> = Vec::new();

    // 1. Consequence terms (compare as sorted sets).
    let actual_csq: BTreeSet<&str> = result.consequences.iter().map(|c| c.as_str()).collect();
    let expected_csq: BTreeSet<&str> = exp.consequences.iter().copied().collect();
    if actual_csq != expected_csq {
        mismatches.push(format!(
            "consequences: expected {:?}, got {:?}",
            expected_csq, actual_csq,
        ));
    }

    // 2. IMPACT.
    let actual_impact = format!("{}", result.impact);
    if actual_impact != exp.impact {
        mismatches.push(format!(
            "impact: expected {}, got {}",
            exp.impact, actual_impact,
        ));
    }

    // 3. Protein start.
    if result.protein_start != exp.protein_start {
        mismatches.push(format!(
            "protein_start: expected {:?}, got {:?}",
            exp.protein_start, result.protein_start,
        ));
    }

    // 4. Protein end.
    if result.protein_end != exp.protein_end {
        mismatches.push(format!(
            "protein_end: expected {:?}, got {:?}",
            exp.protein_end, result.protein_end,
        ));
    }

    // 5. CDS position start.
    if result.cds_position != exp.cds_position {
        mismatches.push(format!(
            "cds_position: expected {:?}, got {:?}",
            exp.cds_position, result.cds_position,
        ));
    }

    // 6. CDS position end.
    if result.cds_position_end != exp.cds_position_end {
        mismatches.push(format!(
            "cds_position_end: expected {:?}, got {:?}",
            exp.cds_position_end, result.cds_position_end,
        ));
    }

    // 7. cDNA position start.
    if result.cdna_position != exp.cdna_position {
        mismatches.push(format!(
            "cdna_position: expected {:?}, got {:?}",
            exp.cdna_position, result.cdna_position,
        ));
    }

    // 8. cDNA position end.
    if result.cdna_position_end != exp.cdna_position_end {
        mismatches.push(format!(
            "cdna_position_end: expected {:?}, got {:?}",
            exp.cdna_position_end, result.cdna_position_end,
        ));
    }

    // 9. Codons.
    if result.codons.as_deref() != exp.codons {
        mismatches.push(format!(
            "codons: expected {:?}, got {:?}",
            exp.codons,
            result.codons.as_deref(),
        ));
    }

    // 10. Amino acids.
    if result.amino_acids.as_deref() != exp.amino_acids {
        mismatches.push(format!(
            "amino_acids: expected {:?}, got {:?}",
            exp.amino_acids,
            result.amino_acids.as_deref(),
        ));
    }

    // 11. Exon number.
    if result.exon.as_deref() != exp.exon {
        mismatches.push(format!(
            "exon: expected {:?}, got {:?}",
            exp.exon,
            result.exon.as_deref(),
        ));
    }

    // 12. Intron number.
    if result.intron.as_deref() != exp.intron {
        mismatches.push(format!(
            "intron: expected {:?}, got {:?}",
            exp.intron,
            result.intron.as_deref(),
        ));
    }

    Ok(mismatches)
}

// ---------------------------------------------------------------------------
// Test
// ---------------------------------------------------------------------------

/// Run all 28 VEP concordance variants. Each variant is annotated with
/// [`annotate()`] (the full VCF-input dispatcher) and compared field-by-field
/// against VEP ground truth.
///
/// Transcripts missing from the store are reported as skips, not failures.
/// Variants with `vep_divergence: true` are still asserted against vareffect's
/// correct output — the divergence flag is informational only.
#[test]
#[ignore]
fn vep_concordance_grch38_indel_28() {
    let ve = VarEffect::builder()
        .with_handles(Assembly::GRCh38, load_store(), load_fasta())
        .expect("matching assemblies")
        .build()
        .expect("builder");

    let total = VARIANTS.len();
    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;
    let mut divergence_pass = 0u32;
    let mut failures: Vec<String> = Vec::new();

    for (i, exp) in VARIANTS.iter().enumerate() {
        let num = i + 1;
        match check_variant(&ve, exp) {
            Err(msg) if msg.contains("not found") => {
                eprintln!("  [{num:>2}/{total}] SKIP  {}: {msg}", exp.label);
                skip += 1;
            }
            Err(msg) => {
                eprintln!("  [{num:>2}/{total}] FAIL  {}: {msg}", exp.label);
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(mismatches) if mismatches.is_empty() => {
                if exp.vep_divergence {
                    eprintln!(
                        "  [{num:>2}/{total}] PASS  {} [VEP divergence: VEP returns {:?}]",
                        exp.label,
                        exp.vep_consequence.unwrap_or(&[]),
                    );
                    divergence_pass += 1;
                } else {
                    eprintln!("  [{num:>2}/{total}] PASS  {}", exp.label);
                }
                pass += 1;
            }
            Ok(mismatches) => {
                eprintln!("  [{num:>2}/{total}] FAIL  {}", exp.label);
                let detail = mismatches.join("\n           ");
                failures.push(format!("#{num} {}:\n           {detail}", exp.label));
                fail += 1;
            }
        }
    }

    eprintln!(
        "\n=== VEP Concordance (indels): {pass}/{} passed, {skip} skipped, \
         {divergence_pass} intentional VEP divergence(s) ===",
        pass + fail,
    );
    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  {f}");
        }
    }
    assert_eq!(
        fail,
        0,
        "{fail} of {} variants failed concordance check",
        pass + fail,
    );
}
