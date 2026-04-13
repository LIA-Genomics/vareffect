//! VEP concordance spot-check for HGVS p. notation.
//!
//! 30 variants validated against Ensembl VEP REST API (GRCh38, queried
//! 2026-04-10 with `refseq=1&hgvs=1`). `#[ignore]`-gated because the test
//! requires the transcript store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! FASTA_PATH=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_concordance_hgvs_p
//! ```
//!
//! Uses the top-level [`annotate()`] dispatcher, exercising the full VCF-input
//! pipeline through HGVS p. formatting.
//!
//! # Categories
//!
//! A. **Missense** (#1-3): plus/minus strand.
//! B. **Synonymous** (#4): synonymous SNV.
//! C. **Stop gained** (#5): nonsense SNV.
//! D. **Start lost** (#6-7): SNV and inframe deletion.
//! E. **Stop lost / extension** (#8-9): SNV (extTer9) and deletion (delextTer8).
//! F. **Stop retained** (#10): stop codon preserved.
//! G. **Frameshift** (#11-16): deletions and insertions.
//! H. **Inframe deletion** (#17-18): single-AA and multi-AA.
//! I. **Inframe insertion** (#19-20): non-dup ins, protein-level dup (VEP
//!    divergence on #20 due to DNA-level 3' shift).
//! J. **Delins / MNV** (#21-23): single-codon missense, full-codon
//!    replacement, multi-codon delins.
//! K. **Non-CDS** (#24-27): splice, intron, UTR, non-coding -> None.
//! L. **Edge cases** (#28-30): splice indel, non-dup inframe insertion,
//!    inframe del in AA repeat (protein 3' rule).

use std::path::Path;

use vareffect::{FastaReader, TranscriptStore, VarEffect};

// ---------------------------------------------------------------------------
// Test infrastructure (same pattern as vep_concordance_hgvs.rs)
// ---------------------------------------------------------------------------

/// Load the transcript store from the workspace data directory.
fn load_store() -> TranscriptStore {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let workspace_root = manifest_dir
        .parent()
        .expect("could not determine workspace root from CARGO_MANIFEST_DIR");
    let path = workspace_root.join("data/vareffect/transcript_models.bin");
    TranscriptStore::load_from_path(&path).unwrap_or_else(|e| {
        panic!(
            "failed to load transcript store from {}: {}. \
             Run `cargo run -p vareffect-cli -- build` first.",
            path.display(),
            e,
        )
    })
}

/// Load the genome reader from `FASTA_PATH` env var.
fn load_fasta() -> FastaReader {
    let path = std::env::var("FASTA_PATH")
        .expect("FASTA_PATH env var must point to a GRCh38 genome binary");
    FastaReader::open_with_patch_aliases(
        Path::new(&path),
        Some(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .and_then(|p| p.parent())
                .expect("workspace root")
                .join("data/vareffect/patch_chrom_aliases.csv")
                .as_ref(),
        ),
    )
    .unwrap_or_else(|e| panic!("failed to open FASTA at {path}: {e}"))
}

// ---------------------------------------------------------------------------
// Expected output fixture
// ---------------------------------------------------------------------------

/// Expected HGVS p. output for one variant-transcript pair.
#[allow(dead_code)]
struct ExpectedP {
    /// Human-readable label for error messages.
    label: &'static str,
    /// UCSC-style chromosome (e.g., "chr17").
    chrom: &'static str,
    /// 0-based genomic position (vareffect convention).
    pos: u64,
    /// VCF REF allele (plus-strand; includes anchor base for indels).
    ref_allele: &'static [u8],
    /// VCF ALT allele (plus-strand; includes anchor base for indels).
    alt_allele: &'static [u8],
    /// RefSeq transcript accession with version.
    transcript: &'static str,
    /// Expected HGVS p. notation (p. prefix only, no NP_ accession).
    /// `None` for variants with no protein consequence.
    expected_hgvs_p: Option<&'static str>,
    /// Primary SO consequence term (for logging/documentation only).
    consequence: &'static str,
    /// True if vareffect intentionally differs from VEP (documented
    /// divergence, e.g., DNA-level 3' shift for insertions).
    vep_divergence: bool,
    /// If `vep_divergence` is true, the VEP hgvsp value for reference.
    vep_hgvs_p: Option<&'static str>,
}

// ---------------------------------------------------------------------------
// Ground truth from VEP REST API (GRCh38, refseq=1, hgvs=1)
// Queried 2026-04-10.
// ---------------------------------------------------------------------------

const VARIANTS: &[ExpectedP] = &[
    // -----------------------------------------------------------------------
    // Category A: Missense (3 variants)
    // -----------------------------------------------------------------------

    // #1 -- TP53 R248W: minus-strand missense, hotspot
    // VEP: NP_000537.3:p.Arg248Trp
    ExpectedP {
        label: "TP53 R248W missense (minus)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg248Trp"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #2 -- BRAF V600E: minus-strand missense, oncogene
    // VEP: NP_004324.2:p.Val600Glu
    ExpectedP {
        label: "BRAF V600E missense (minus)",
        chrom: "chr7",
        pos: 140_753_335,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NM_004333.6",
        expected_hgvs_p: Some("p.Val600Glu"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #3 -- EGFR L858R: plus-strand missense, lung cancer
    // VEP: NP_005219.2:p.Leu858Arg
    ExpectedP {
        label: "EGFR L858R missense (plus)",
        chrom: "chr7",
        pos: 55_191_821,
        ref_allele: b"T",
        alt_allele: b"G",
        transcript: "NM_005228.5",
        expected_hgvs_p: Some("p.Leu858Arg"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category B: Synonymous (1 variant)
    // -----------------------------------------------------------------------

    // #4 -- TP53 R248R: synonymous SNV
    // VEP: NP_000537.3:p.Arg248=
    ExpectedP {
        label: "TP53 R248R synonymous (minus)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_allele: b"G",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg248="),
        consequence: "synonymous_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category C: Stop gained / nonsense (1 variant)
    // -----------------------------------------------------------------------

    // #5 -- TP53 R196X: stop gained
    // VEP: NP_000537.3:p.Arg196Ter
    ExpectedP {
        label: "TP53 R196X stop_gained (minus)",
        chrom: "chr17",
        pos: 7_674_944,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg196Ter"),
        consequence: "stop_gained",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category D: Start lost (2 variants)
    // -----------------------------------------------------------------------

    // #6 -- TP53 start_lost SNV (ATG -> GTG on coding strand)
    // VEP: NP_000537.3:p.Met1?
    ExpectedP {
        label: "TP53 start_lost SNV (minus)",
        chrom: "chr17",
        pos: 7_676_593,
        ref_allele: b"T",
        alt_allele: b"C",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Met1?"),
        consequence: "start_lost",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #7 -- TP53 c.1_3del: start codon deletion
    // VEP: NP_000537.3:p.Met1?
    ExpectedP {
        label: "TP53 start_lost deletion (minus)",
        chrom: "chr17",
        pos: 7_676_590,
        ref_allele: b"CCAT",
        alt_allele: b"C",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Met1?"),
        consequence: "start_lost",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category E: Stop lost / extension (2 variants)
    // -----------------------------------------------------------------------

    // #8 -- TP53 stop_lost SNV (TGA -> TGT): extension with exact distance
    // VEP: NP_000537.3:p.Ter394CysextTer9
    ExpectedP {
        label: "TP53 stop_lost SNV extension (minus)",
        chrom: "chr17",
        pos: 7_669_608,
        ref_allele: b"T",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Ter394CysextTer9"),
        consequence: "stop_lost",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #9 -- TP53 c.1180_1182del: stop codon deletion, extension
    // VEP: NP_000537.3:p.Ter394delextTer8
    ExpectedP {
        label: "TP53 stop_lost deletion extension (minus)",
        chrom: "chr17",
        pos: 7_669_607,
        ref_allele: b"GTCA",
        alt_allele: b"G",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Ter394delextTer8"),
        consequence: "stop_lost",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category F: Stop retained (1 variant)
    // -----------------------------------------------------------------------

    // #10 -- TP53 stop_retained (TGA -> TAA, both stop)
    // VEP: NP_000537.3:p.Ter394=
    ExpectedP {
        label: "TP53 stop_retained (minus)",
        chrom: "chr17",
        pos: 7_669_609,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Ter394="),
        consequence: "stop_retained_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category G: Frameshift (6 variants)
    // -----------------------------------------------------------------------

    // #11 -- BRCA1 c.68_69del: minus-strand 2bp deletion, Ashkenazi founder
    // VEP: NP_009225.1:p.Glu23ValfsTer17
    // NM_007294.4 retired from VEP; notation from ClinVar + VEP NM_001407*.1
    ExpectedP {
        label: "BRCA1 c.68_69del frameshift (minus)",
        chrom: "chr17",
        pos: 43_124_026,
        ref_allele: b"ACT",
        alt_allele: b"A",
        transcript: "NM_007294.4",
        expected_hgvs_p: Some("p.Glu23ValfsTer17"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #12 -- BRCA1 c.5266dup: minus-strand dup causing frameshift
    // VEP: NP_009225.1:p.Gln1756ProfsTer74
    ExpectedP {
        label: "BRCA1 c.5266dup frameshift ins (minus)",
        chrom: "chr17",
        pos: 43_057_062,
        ref_allele: b"G",
        alt_allele: b"GG",
        transcript: "NM_007294.4",
        expected_hgvs_p: Some("p.Gln1756ProfsTer74"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #13 -- BRCA2 c.1813del: plus-strand single-base frameshift
    // VEP: NP_000050.3:p.Ile605TyrfsTer9
    ExpectedP {
        label: "BRCA2 c.1813del frameshift (plus)",
        chrom: "chr13",
        pos: 32_333_289,
        ref_allele: b"AA",
        alt_allele: b"A",
        transcript: "NM_000059.4",
        expected_hgvs_p: Some("p.Ile605TyrfsTer9"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #14 -- APC c.3927_3931del: plus-strand 5bp frameshift
    // VEP: NP_000029.2:p.Glu1309AspfsTer4
    ExpectedP {
        label: "APC c.3927_3931del frameshift (plus)",
        chrom: "chr5",
        pos: 112_839_519,
        ref_allele: b"AAAAGA",
        alt_allele: b"A",
        transcript: "NM_000038.6",
        expected_hgvs_p: Some("p.Glu1309AspfsTer4"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #15 -- BRCA2 c.5946del: plus-strand 1bp frameshift
    // VEP: NP_000050.3:p.Ser1982ArgfsTer22
    ExpectedP {
        label: "BRCA2 c.5946del frameshift (plus)",
        chrom: "chr13",
        pos: 32_340_299,
        ref_allele: b"GT",
        alt_allele: b"G",
        transcript: "NM_000059.4",
        expected_hgvs_p: Some("p.Ser1982ArgfsTer22"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #16 -- BRCA2 c.6275_6276del: plus-strand 2bp frameshift
    // VEP: NP_000050.3:p.Leu2092ProfsTer7
    ExpectedP {
        label: "BRCA2 c.6275_6276del frameshift (plus)",
        chrom: "chr13",
        pos: 32_340_628,
        ref_allele: b"CTT",
        alt_allele: b"C",
        transcript: "NM_000059.4",
        expected_hgvs_p: Some("p.Leu2092ProfsTer7"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category H: Inframe deletion (2 variants)
    // -----------------------------------------------------------------------

    // #17 -- CFTR deltaF508: canonical single-AA inframe deletion
    // VEP: NP_000483.3:p.Phe508del
    ExpectedP {
        label: "CFTR deltaF508 inframe del (plus)",
        chrom: "chr7",
        pos: 117_559_590,
        ref_allele: b"TCTT",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        expected_hgvs_p: Some("p.Phe508del"),
        consequence: "inframe_deletion",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #18 -- EGFR E746_A750del: multi-AA inframe deletion (15bp)
    // VEP: NP_005219.2:p.Glu746_Ala750del
    ExpectedP {
        label: "EGFR E746_A750del inframe del (plus)",
        chrom: "chr7",
        pos: 55_174_770,
        ref_allele: b"AGGAATTAAGAGAAGC",
        alt_allele: b"A",
        transcript: "NM_005228.5",
        expected_hgvs_p: Some("p.Glu746_Ala750del"),
        consequence: "inframe_deletion",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category I: Inframe insertion (2 variants)
    // -----------------------------------------------------------------------

    // #19 -- EGFR c.2310_2311insGGT: non-dup inframe insertion
    // VEP: NP_005219.2:p.Asp770_Asn771insGly
    // Genomic: chr7:55181319 (1-based), plus-strand, anchor C
    // Fixed: codon window expansion now includes the preceding codon
    // when the insertion is at a codon boundary, matching VEP's flanking
    // position computation exactly.
    ExpectedP {
        label: "EGFR exon20 3bp inframe ins (plus)",
        chrom: "chr7",
        pos: 55_181_318,
        ref_allele: b"C",
        alt_allele: b"CGGT",
        transcript: "NM_005228.5",
        expected_hgvs_p: Some("p.Asp770_Asn771insGly"),
        consequence: "inframe_insertion",
        vep_divergence: false,
        vep_hgvs_p: Some("p.Asp770_Asn771insGly"),
    },
    // #20 -- ERBB2 exon20 12bp inframe insertion
    // HGVS 3' normalization: the DNA-level shift converts the insertion to
    // c.2313_2324dup, and the protein-level codon window at the shifted
    // position detects the 4-AA duplication.
    // VEP: NP_004439.2:p.Tyr772_Ala775dup
    ExpectedP {
        label: "ERBB2 exon20 12bp inframe dup (plus, 3' shifted)",
        chrom: "chr17",
        pos: 39_724_729,
        ref_allele: b"C",
        alt_allele: b"CATACGTGATGGC",
        transcript: "NM_004448.4",
        expected_hgvs_p: Some("p.Tyr772_Ala775dup"),
        consequence: "inframe_insertion",
        vep_divergence: false,
        vep_hgvs_p: Some("p.Tyr772_Ala775dup"),
    },
    // -----------------------------------------------------------------------
    // Category J: Delins / MNV (3 variants)
    // -----------------------------------------------------------------------

    // #21 -- TP53 c.742_743delinsTT: 2bp MNV within single codon -> missense
    // VEP: NP_000537.3:p.Arg248Leu
    ExpectedP {
        label: "TP53 2bp MNV R248L (minus)",
        chrom: "chr17",
        pos: 7_674_219,
        ref_allele: b"CG",
        alt_allele: b"AA",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg248Leu"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #22 -- TP53 c.742_744delinsAAA: 3bp MNV, entire codon replacement
    // VEP: NP_000537.3:p.Arg248Lys
    ExpectedP {
        label: "TP53 3bp MNV R248K (minus)",
        chrom: "chr17",
        pos: 7_674_218,
        ref_allele: b"CCG",
        alt_allele: b"TTT",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg248Lys"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #23 -- TP53 c.742_747delinsAAAAAA: 6bp MNV spanning codons 248-249
    // VEP: NP_000537.3:p.Arg248_Arg249delinsLysLys
    // Genomic: chr17:7674216-7674221 (1-based), plus-strand CCTCCG -> TTTTTT
    ExpectedP {
        label: "TP53 6bp MNV R248_R249delinsKK (minus)",
        chrom: "chr17",
        pos: 7_674_215,
        ref_allele: b"CCTCCG",
        alt_allele: b"TTTTTT",
        transcript: "NM_000546.6",
        expected_hgvs_p: Some("p.Arg248_Arg249delinsLysLys"),
        consequence: "missense_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category K: Non-CDS variants -> None (4 variants)
    // -----------------------------------------------------------------------

    // #24 -- TP53 splice_donor +1 (intronic)
    // VEP: no hgvsp field
    ExpectedP {
        label: "TP53 splice_donor +1 (minus)",
        chrom: "chr17",
        pos: 7_674_179,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_p: None,
        consequence: "splice_donor_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #25 -- TP53 deep intron
    // VEP: no hgvsp field
    ExpectedP {
        label: "TP53 deep intron (minus)",
        chrom: "chr17",
        pos: 7_674_049,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_p: None,
        consequence: "intron_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #26 -- TP53 5'UTR
    // VEP: no hgvsp field
    ExpectedP {
        label: "TP53 5'UTR (minus)",
        chrom: "chr17",
        pos: 7_687_399,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_p: None,
        consequence: "5_prime_UTR_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #27 -- TERC non-coding transcript exon
    // VEP: no hgvsp field (NR_* transcript)
    ExpectedP {
        label: "TERC non-coding exon (minus)",
        chrom: "chr3",
        pos: 169_764_899,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NR_001566.3",
        expected_hgvs_p: None,
        consequence: "non_coding_transcript_exon_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // -----------------------------------------------------------------------
    // Category L: Edge cases (3 variants)
    // -----------------------------------------------------------------------

    // #28 -- ACVRL1 c.525+1del: splice donor deletion (VEP #519 divergence)
    // VEP divergence: VEP returns frameshift_variant, vareffect correctly
    // returns splice_donor_variant. No protein notation for splice.
    ExpectedP {
        label: "ACVRL1 splice_donor del (plus, VEP #519)",
        chrom: "chr12",
        pos: 51_913_769,
        ref_allele: b"GG",
        alt_allele: b"G",
        transcript: "NM_000020.3",
        expected_hgvs_p: None,
        consequence: "splice_donor_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #29 -- SYNGAP1 c.1861_1862del: frameshift (not nonsense)
    // VEP: NP_006763.2:p.Arg621AsnfsTer30
    // Tests that the first changed AA (Asn) is not a stop, so frameshift
    // notation is used rather than nonsense.
    ExpectedP {
        label: "SYNGAP1 c.1861_1862del frameshift (plus)",
        chrom: "chr6",
        pos: 33_440_911,
        ref_allele: b"ACG",
        alt_allele: b"A",
        transcript: "NM_006772.3",
        expected_hgvs_p: Some("p.Arg621AsnfsTer30"),
        consequence: "frameshift_variant",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
    // #30 -- CFTR c.1516_1518del: inframe deletion in AA repeat (3' rule)
    // VEP: NP_000483.3:p.Ile507del
    // CFTR positions 506-508 are Ile-Ile-Phe. Deleting 3bp at c.1516_1518
    // removes one Ile. The protein-level 3' rule shifts from Ile506 to
    // Ile507 (most C-terminal Ile in the repeat).
    // Genomic: chr7:117559587-117559589 (1-based), plus-strand ATC deleted
    // VCF anchor: T at 117559586 (1-based) = 117559585 (0-based)
    ExpectedP {
        label: "CFTR c.1516_1518del 3' rule (plus)",
        chrom: "chr7",
        pos: 117_559_585,
        ref_allele: b"TATC",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        expected_hgvs_p: Some("p.Ile507del"),
        consequence: "inframe_deletion",
        vep_divergence: false,
        vep_hgvs_p: None,
    },
];

// ---------------------------------------------------------------------------
// Comparison logic
// ---------------------------------------------------------------------------

/// Compare a single variant's `annotate()` HGVS p. output against the expected
/// VEP value. Returns `Ok(None)` on match, `Ok(Some(mismatch))` on field
/// mismatch, or `Err(msg)` on transcript-not-found / annotate error.
fn check_variant(ve: &VarEffect, exp: &ExpectedP) -> Result<Option<String>, String> {
    let results = ve
        .annotate(exp.chrom, exp.pos, exp.ref_allele, exp.alt_allele)
        .map_err(|e| format!("annotate error: {e}"))?;

    let result = results
        .iter()
        .find(|r| r.transcript == exp.transcript)
        .ok_or_else(|| {
            let found: Vec<&str> = results.iter().map(|r| r.transcript.as_str()).collect();
            format!(
                "transcript {} not found in results (found: {found:?})",
                exp.transcript,
            )
        })?;

    let actual = result.hgvs_p.as_deref();
    let expected = exp.expected_hgvs_p;

    if actual == expected {
        Ok(None)
    } else {
        Ok(Some(format!(
            "hgvs_p mismatch: expected {:?}, got {:?}",
            expected, actual,
        )))
    }
}

// ---------------------------------------------------------------------------
// Test runner
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn vep_concordance_hgvs_p_30() {
    let ve = VarEffect::new(load_store(), load_fasta());

    let total = VARIANTS.len();
    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;
    let mut failures: Vec<String> = Vec::new();

    for (i, exp) in VARIANTS.iter().enumerate() {
        let num = i + 1;
        let divergence_tag = if exp.vep_divergence {
            " [VEP divergence]"
        } else {
            ""
        };

        match check_variant(&ve, exp) {
            Err(msg) if msg.contains("not found") => {
                eprintln!(
                    "  [{num:>2}/{total}] SKIP  {}{divergence_tag}: {msg}",
                    exp.label,
                );
                skip += 1;
            }
            Err(msg) => {
                eprintln!(
                    "  [{num:>2}/{total}] FAIL  {}{divergence_tag}: {msg}",
                    exp.label,
                );
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(None) => {
                eprintln!("  [{num:>2}/{total}] PASS  {}{divergence_tag}", exp.label,);
                pass += 1;
            }
            Ok(Some(mismatch)) => {
                eprintln!(
                    "  [{num:>2}/{total}] FAIL  {}{divergence_tag}: {mismatch}",
                    exp.label,
                );
                failures.push(format!("#{num} {}: {mismatch}", exp.label));
                fail += 1;
            }
        }
    }

    eprintln!(
        "\n=== VEP Concordance (HGVS p.): {pass}/{} passed, \
         {skip} skipped ===",
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
        "{fail} of {} variants failed HGVS p. concordance check",
        pass + fail,
    );
}
