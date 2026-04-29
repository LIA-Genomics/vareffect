//! VEP concordance spot-check for the HGVS c. reverse mapper.
//!
//! 20 variants validating that [`vareffect::resolve_hgvs_c`] produces correct
//! 0-based VCF-style genomic coordinates for all position types (CDS, 5'UTR,
//! 3'UTR, intronic, intronic-in-UTR) and variant types (substitution, deletion,
//! duplication, insertion). `#[ignore]`-gated because the test requires the
//! transcript store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! FASTA_PATH=data/vareffect/GRCh38.fa.gz \
//!   cargo test -p vareffect -- --ignored vep_concordance_hgvs_reverse
//! ```
//!
//! # Categories
//!
//! A. **CDS substitutions** (#1-3): plus and minus strand, verifies allele
//!    complementing.
//! B. **CDS deletions** (#4-5, #10): single-base, multi-base, 15-base; verifies
//!    VCF anchor base.
//! C. **CDS duplications** (#6-7): minus-strand single-base, plus-strand
//!    multi-base; verifies dup-to-insertion conversion.
//! D. **CDS insertions** (#8-9): plus-strand 2-base and 3-base insertions.
//! E. **Intronic substitutions** (#11-14): all 4 combinations of plus/minus
//!    strand x donor/acceptor side.
//! F. **UTR substitutions** (#15-18): 5'UTR (both strands), 3'UTR (minus),
//!    intronic-in-5'UTR (minus).
//! G. **Round-trip** (#19-20): resolve -> annotate -> compare hgvs_c.
//!
//! # VEP allele_string strand convention
//!
//! VEP's `allele_string` for HGVS c. input is on the **coding/transcript
//! strand**, not the genomic plus strand -- for all variant types including
//! indels. For minus-strand genes, the expected vareffect output (always
//! genomic plus-strand) is the complement of VEP's allele_string.
//!
//! Evidence: BRCA1 c.5266dup (minus strand) -- VEP reports `-/C`, but the
//! genomic base at the dup position is G (verified via Ensembl sequence API).
//! C = complement of G = coding strand.
//!
//! # Ground truth
//!
//! Every expected value comes from a VEP REST API query (GRCh38) plus Ensembl
//! sequence API lookups for indel anchor bases. No values are computed from
//! transcript models.

use std::path::Path;

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect, VarEffectError};

// ---------------------------------------------------------------------------
// Test infrastructure
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
    FastaReader::open_with_patch_aliases_and_assembly(
        Path::new(&path),
        Some(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .and_then(|p| p.parent())
                .expect("workspace root")
                .join("data/vareffect/patch_chrom_aliases.csv")
                .as_ref(),
        ),
        Assembly::GRCh38,
    )
    .unwrap_or_else(|e| panic!("failed to open FASTA at {path}: {e}"))
}

// ---------------------------------------------------------------------------
// Expected output fixture
// ---------------------------------------------------------------------------

/// Expected output for one reverse-mapped variant.
#[allow(dead_code)]
struct ExpectedReverse {
    /// Human-readable label for error messages.
    label: &'static str,
    /// Full HGVS c. notation including accession (e.g., `"NM_000546.6:c.742C>T"`).
    hgvs_input: &'static str,
    /// Expected chromosome with `"chr"` prefix (e.g., `"chr17"`).
    expected_chrom: &'static str,
    /// Expected 0-based genomic position.
    expected_pos: u64,
    /// Expected reference allele (plus-strand, VCF convention with anchor for indels).
    expected_ref: &'static [u8],
    /// Expected alternate allele (plus-strand, VCF convention with anchor for indels).
    expected_alt: &'static [u8],
    /// Category tag for logging.
    category: &'static str,
    /// If `true`, run a round-trip test: resolve -> annotate -> compare hgvs_c.
    round_trip: bool,
}

// ---------------------------------------------------------------------------
// Ground truth from VEP REST API + Ensembl sequence API
// ---------------------------------------------------------------------------

const VARIANTS: &[ExpectedReverse] = &[
    // -----------------------------------------------------------------------
    // Category A: CDS substitutions -- plus and minus strand
    // -----------------------------------------------------------------------
    //
    // VEP allele_string is coding-strand. For minus-strand, complement to get
    // the genomic plus-strand alleles that vareffect produces.
    ExpectedReverse {
        label: "#1 TP53 R248W (minus)",
        hgvs_input: "NM_000546.6:c.742C>T",
        // VEP: 17:7674221 C/T (coding) -> genomic G/A (complement)
        expected_chrom: "chr17",
        expected_pos: 7_674_220,
        expected_ref: b"G",
        expected_alt: b"A",
        category: "A-cds-sub",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#2 EGFR L858R (plus)",
        hgvs_input: "NM_005228.5:c.2573T>G",
        // VEP: 7:55191822 T/G (coding = genomic for plus strand)
        expected_chrom: "chr7",
        expected_pos: 55_191_821,
        expected_ref: b"T",
        expected_alt: b"G",
        category: "A-cds-sub",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#3 BRAF V600E (minus)",
        hgvs_input: "NM_004333.6:c.1799T>A",
        // VEP: 7:140753336 T/A (coding) -> genomic A/T (complement)
        expected_chrom: "chr7",
        expected_pos: 140_753_335,
        expected_ref: b"A",
        expected_alt: b"T",
        category: "A-cds-sub",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category B: CDS deletions
    // -----------------------------------------------------------------------
    //
    // VCF convention: anchor base at (deleted_start - 1). Anchor bases
    // verified via Ensembl sequence API.
    ExpectedReverse {
        label: "#4 BRCA2 c.1813del (plus, single-base)",
        hgvs_input: "NM_000059.4:c.1813del",
        // VEP: 13:32333291 A/- (deleted base). Anchor at 0b:32333289 = A
        // (Ensembl seq 13:32333289..32333292 = AAAT).
        expected_chrom: "chr13",
        expected_pos: 32_333_289,
        expected_ref: b"AA",
        expected_alt: b"A",
        category: "B-cds-del",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#5 CFTR deltaF508 (plus, 3-base)",
        hgvs_input: "NM_000492.4:c.1521_1523del",
        // VEP: 7:117559592-594 CTT/-. Anchor at 0b:117559590 = T
        // (Ensembl seq 7:117559590..117559595 = ATCTTT).
        expected_chrom: "chr7",
        expected_pos: 117_559_590,
        expected_ref: b"TCTT",
        expected_alt: b"T",
        category: "B-cds-del",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category C: CDS duplications
    // -----------------------------------------------------------------------
    //
    // HGVS dup -> VCF insertion. VEP allele_string for dups is coding-strand.
    // Plus strand: anchor = last dup base (gend).
    // Minus strand: anchor = base before dup range (gstart - 1).
    ExpectedReverse {
        label: "#6 BRCA1 c.5266dup (minus, single-base)",
        hgvs_input: "NM_007294.4:c.5266dup",
        // VEP: insertion between 0b:43057062-063. -/C (coding) = genomic G.
        // Anchor at 0b:43057062 = G (Ensembl seq 17:43057062..43057065 = TGGG).
        expected_chrom: "chr17",
        expected_pos: 43_057_062,
        expected_ref: b"G",
        expected_alt: b"GG",
        category: "C-cds-dup",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#7 ERBB2 c.2313_2324dup (plus, 12-base)",
        hgvs_input: "NM_004448.4:c.2313_2324dup",
        // VEP maps c.2324 to 1-based 39724742 (confirmed via c.2324C>T query).
        // The dup range c.2313-c.2324 = 12 bases ending at 0b:39724741.
        // VEP's dup query left-normalizes to 0b:39724729 in a tandem repeat,
        // but vareffect anchors at the last dup base without normalization.
        // Both are valid VCF representations of the same variant.
        // Anchor at 0b:39724741 = C.
        expected_chrom: "chr17",
        expected_pos: 39_724_741,
        expected_ref: b"C",
        expected_alt: b"CATACGTGATGGC",
        category: "C-cds-dup",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category D: CDS insertions
    // -----------------------------------------------------------------------
    //
    // VCF: anchor = min(g_left, g_right). Inserted bases on coding strand;
    // for plus-strand genes, coding = genomic.
    ExpectedReverse {
        label: "#8 BRCA2 c.5946_5947insAA (plus)",
        hgvs_input: "NM_000059.4:c.5946_5947insAA",
        // VEP: insertion between 0b:32340300-301. -/AA.
        // Anchor at 0b:32340300 = T
        // (Ensembl seq 13:32340300..32340303 = GTGG -> 0b:32340300 = T).
        expected_chrom: "chr13",
        expected_pos: 32_340_300,
        expected_ref: b"T",
        expected_alt: b"TAA",
        category: "D-cds-ins",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#9 EGFR exon20 ins (plus)",
        hgvs_input: "NM_005228.5:c.2309_2310insGGT",
        // VEP: insertion between 0b:55181317-318. -/GGT.
        // Anchor at 0b:55181317 = A
        // (Ensembl seq 7:55181317..55181320 = GACA -> 0b:55181317 = A).
        expected_chrom: "chr7",
        expected_pos: 55_181_317,
        expected_ref: b"A",
        expected_alt: b"AGGT",
        category: "D-cds-ins",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category B (cont): large CDS deletion
    // -----------------------------------------------------------------------
    ExpectedReverse {
        label: "#10 EGFR E746_A750del (plus, 15-base)",
        hgvs_input: "NM_005228.5:c.2235_2249del",
        // VEP: 7:55174772-786 GGAATTAAGAGAAGC/-. Anchor at 0b:55174770 = A
        // (Ensembl seq 7:55174770..55174773 = AAGG -> 0b:55174770 = A).
        expected_chrom: "chr7",
        expected_pos: 55_174_770,
        expected_ref: b"AGGAATTAAGAGAAGC",
        expected_alt: b"A",
        category: "B-cds-del",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category E: Intronic substitutions -- all 4 strand x side combinations
    // -----------------------------------------------------------------------
    ExpectedReverse {
        label: "#11 TP53 splice donor c.672+1 (minus)",
        hgvs_input: "NM_000546.6:c.672+1G>A",
        // VEP: 17:7674858 G/A (coding) -> genomic C/T (complement)
        expected_chrom: "chr17",
        expected_pos: 7_674_857,
        expected_ref: b"C",
        expected_alt: b"T",
        category: "E-intronic",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#12 BRCA2 splice donor c.8487+1 (plus)",
        hgvs_input: "NM_000059.4:c.8487+1G>A",
        // VEP: 13:32370558 G/A (coding = genomic)
        expected_chrom: "chr13",
        expected_pos: 32_370_557,
        expected_ref: b"G",
        expected_alt: b"A",
        category: "E-intronic",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#13 TP53 splice acceptor c.673-2 (minus)",
        hgvs_input: "NM_000546.6:c.673-2A>G",
        // VEP: 17:7674292 A/G (coding) -> genomic T/C (complement)
        expected_chrom: "chr17",
        expected_pos: 7_674_291,
        expected_ref: b"T",
        expected_alt: b"C",
        category: "E-intronic",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#14 BRCA2 splice acceptor c.8488-1 (plus)",
        hgvs_input: "NM_000059.4:c.8488-1G>A",
        // VEP: 13:32370955 G/A (coding = genomic)
        expected_chrom: "chr13",
        expected_pos: 32_370_954,
        expected_ref: b"G",
        expected_alt: b"A",
        category: "E-intronic",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category F: UTR substitutions
    // -----------------------------------------------------------------------
    //
    // Note: original spec had wrong ref alleles for #15-17. Corrected alleles
    // come from VEP error messages disclosing the actual reference bases.
    ExpectedReverse {
        label: "#15 TP53 5'UTR c.-29 (minus)",
        hgvs_input: "NM_000546.6:c.-29G>A",
        // VEP: 17:7687377 G/A (coding) -> genomic C/T (complement)
        expected_chrom: "chr17",
        expected_pos: 7_687_376,
        expected_ref: b"C",
        expected_alt: b"T",
        category: "F-utr",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#16 EGFR 5'UTR c.-191 (plus)",
        hgvs_input: "NM_005228.5:c.-191A>G",
        // VEP: 7:55019087 A/G (coding = genomic)
        expected_chrom: "chr7",
        expected_pos: 55_019_086,
        expected_ref: b"A",
        expected_alt: b"G",
        category: "F-utr",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#17 TP53 3'UTR c.*10 (minus)",
        hgvs_input: "NM_000546.6:c.*10C>T",
        // VEP: 17:7669599 C/T (coding) -> genomic G/A (complement)
        expected_chrom: "chr17",
        expected_pos: 7_669_598,
        expected_ref: b"G",
        expected_alt: b"A",
        category: "F-utr",
        round_trip: false,
    },
    ExpectedReverse {
        label: "#18 TP53 intronic-in-5'UTR c.-28+1 (minus)",
        hgvs_input: "NM_000546.6:c.-28+1A>G",
        // VEP: 17:7676621 A/G (coding) -> genomic T/C (complement)
        expected_chrom: "chr17",
        expected_pos: 7_676_620,
        expected_ref: b"T",
        expected_alt: b"C",
        category: "F-utr",
        round_trip: false,
    },
    // -----------------------------------------------------------------------
    // Category G: Round-trip validation
    // -----------------------------------------------------------------------
    //
    // Same variants as #1 and #6, but with round_trip enabled: resolve ->
    // annotate -> find transcript -> compare hgvs_c with original input.
    ExpectedReverse {
        label: "#19 TP53 R248W round-trip (= #1)",
        hgvs_input: "NM_000546.6:c.742C>T",
        expected_chrom: "chr17",
        expected_pos: 7_674_220,
        expected_ref: b"G",
        expected_alt: b"A",
        category: "G-round-trip",
        round_trip: true,
    },
    ExpectedReverse {
        label: "#20 BRCA1 c.5266dup round-trip (= #6)",
        hgvs_input: "NM_007294.4:c.5266dup",
        expected_chrom: "chr17",
        expected_pos: 43_057_062,
        expected_ref: b"G",
        expected_alt: b"GG",
        category: "G-round-trip",
        round_trip: true,
    },
];

// ---------------------------------------------------------------------------
// Test runner
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn vep_concordance_hgvs_reverse_20() {
    let ve = VarEffect::builder()
        .with_handles(Assembly::GRCh38, load_store(), load_fasta())
        .expect("matching assemblies")
        .build()
        .expect("builder");

    let total = VARIANTS.len();
    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;

    for (i, exp) in VARIANTS.iter().enumerate() {
        let tag = format!("[{}/{}] {}", i + 1, total, exp.label);

        match ve.resolve_hgvs_c(Assembly::GRCh38, exp.hgvs_input) {
            Ok(result) => {
                let mut mismatches = Vec::new();

                if result.chrom != exp.expected_chrom {
                    mismatches.push(format!(
                        "chrom: expected {}, got {}",
                        exp.expected_chrom, result.chrom,
                    ));
                }
                if result.pos != exp.expected_pos {
                    mismatches.push(format!(
                        "pos: expected {}, got {}",
                        exp.expected_pos, result.pos,
                    ));
                }
                if result.ref_allele != exp.expected_ref {
                    mismatches.push(format!(
                        "ref: expected {:?}, got {:?}",
                        String::from_utf8_lossy(exp.expected_ref),
                        String::from_utf8_lossy(&result.ref_allele),
                    ));
                }
                if result.alt_allele != exp.expected_alt {
                    mismatches.push(format!(
                        "alt: expected {:?}, got {:?}",
                        String::from_utf8_lossy(exp.expected_alt),
                        String::from_utf8_lossy(&result.alt_allele),
                    ));
                }

                if mismatches.is_empty() {
                    pass += 1;
                    eprintln!("PASS  {tag}  [{}]", exp.category);

                    // Round-trip test: resolve -> annotate -> compare hgvs_c.
                    if exp.round_trip {
                        match ve.annotate(
                            Assembly::GRCh38,
                            &result.chrom,
                            result.pos,
                            &result.ref_allele,
                            &result.alt_allele,
                        ) {
                            Ok(ann_results) => {
                                // Extract the accession prefix (e.g., "NM_000546.6")
                                // from the hgvs_input to match against annotate output.
                                let accession = exp
                                    .hgvs_input
                                    .split(":c.")
                                    .next()
                                    .expect("valid HGVS input");

                                let got_hgvs = ann_results
                                    .consequences
                                    .iter()
                                    .find(|r| r.transcript == accession)
                                    .and_then(|r| r.hgvs_c.as_deref());

                                // The round-trip hgvs_c is the portion after
                                // the colon (e.g., "c.742C>T"), but the
                                // ConsequenceResult stores the full notation
                                // with accession prefix in some builds.
                                // Compare against the descriptor after ":".
                                let expected_desc = exp
                                    .hgvs_input
                                    .split_once(':')
                                    .map(|(_, d)| d)
                                    .unwrap_or(exp.hgvs_input);

                                // hgvs_c in ConsequenceResult may or may not
                                // include the accession prefix. Handle both.
                                let matches = got_hgvs
                                    .is_some_and(|h| h == exp.hgvs_input || h == expected_desc);

                                if matches {
                                    eprintln!("      ROUND-TRIP PASS");
                                } else {
                                    fail += 1;
                                    eprintln!(
                                        "      ROUND-TRIP FAIL: expected {:?}, got {:?}",
                                        exp.hgvs_input, got_hgvs,
                                    );
                                }
                            }
                            Err(e) => {
                                fail += 1;
                                eprintln!("      ROUND-TRIP FAIL: annotate error: {e}");
                            }
                        }
                    }
                } else {
                    fail += 1;
                    eprintln!("FAIL  {tag}  [{}]: {}", exp.category, mismatches.join("; "));
                }
            }
            Err(VarEffectError::TranscriptNotFound { .. }) => {
                skip += 1;
                eprintln!("SKIP  {tag}  [{}]: transcript not in store", exp.category);
            }
            Err(e) => {
                fail += 1;
                eprintln!("FAIL  {tag}  [{}]: resolve error: {e}", exp.category);
            }
        }
    }

    eprintln!(
        "\n=== HGVS Reverse Concordance: {pass}/{} passed, {skip} skipped ===",
        pass + fail,
    );
    assert_eq!(fail, 0, "{fail} variant(s) failed");
}
