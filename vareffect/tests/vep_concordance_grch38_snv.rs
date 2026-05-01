//! VEP concordance spot-check for SNV consequence prediction.
//!
//! 20 real clinical-grade SNVs validated against Ensembl VEP REST API (GRCh38,
//! queried 2026-04-09 with `refseq=1&numbers=1`). `#[ignore]`-gated because
//! the test requires the transcript store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_concordance
//! ```
//!
//! # Known VEP divergences
//!
//! These are intentional differences between vareffect and VEP. Each is
//! documented here rather than in the test fixture so the expected values
//! reflect what vareffect *should* produce.
//!
//! 1. **Granular splice region terms** — VEP (release 113+) emits
//!    `splice_polypyrimidine_tract_variant` (acceptor side),
//!    `splice_donor_region_variant` (donor +3/+4/+6), and
//!    `splice_donor_5th_base_variant` (donor +5) alongside the generic
//!    `splice_region_variant`. Vareffect emits only the core
//!    `splice_region_variant` + `intron_variant` pair. These granular terms
//!    are VEP-specific extensions and not required by SO or ClinGen.
//!
//! 2. **`non_coding_transcript_variant`** — VEP may emit this term for
//!    intronic positions in non-coding transcripts. Vareffect currently does
//!    not; it uses the standard `intron_variant` term alone.

use std::collections::BTreeSet;
use std::path::Path;

use vareffect::{Assembly, FastaReader, TranscriptStore, annotate_snv};

/// Load the GRCh38 transcript store.
///
/// Reads the path from `GRCH38_TRANSCRIPTS` (matching the GRCh37 sibling's
/// `GRCH37_TRANSCRIPTS` convention) or falls back to the workspace default
/// `data/vareffect/transcript_models_grch38.bin`. Panics with a clear
/// message if the store is missing.
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

/// Load the GRCh38 genome reader from the `GRCH38_FASTA` environment variable.
fn load_fasta() -> FastaReader {
    let path = std::env::var("GRCH38_FASTA").expect(
        "GRCH38_FASTA env var must point to a GRCh38 genome binary (.bin) \
         with its .bin.idx sidecar. Run `vareffect setup --assembly grch38` \
         first, then set GRCH38_FASTA=data/vareffect/GRCh38.bin.",
    );
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

/// Expected output for one variant-transcript pair, derived from the VEP REST
/// API. All coordinate fields use VEP's conventions (1-based for CDS/cDNA/
/// protein positions). The `pos` field is 0-based (vareffect convention).
struct Expected {
    /// Human-readable label for error messages.
    label: &'static str,
    /// UCSC-style chromosome (e.g., "chr17").
    chrom: &'static str,
    /// 0-based genomic position (vareffect convention).
    pos: u64,
    /// Plus-strand reference base.
    ref_base: u8,
    /// Plus-strand alternate base.
    alt_base: u8,
    /// RefSeq transcript accession with version.
    transcript: &'static str,
    /// Expected SO consequence terms (sorted alphabetically).
    consequences: &'static [&'static str],
    /// Expected VEP IMPACT string.
    impact: &'static str,
    /// Expected 1-based protein position, or `None`.
    protein_start: Option<u32>,
    /// Expected 1-based CDS position, or `None`.
    cds_position: Option<u32>,
    /// Expected 1-based cDNA position, or `None`.
    cdna_position: Option<u32>,
    /// Expected codon string (e.g., "Cgg/Tgg"), or `None`.
    codons: Option<&'static str>,
    /// Expected amino acid string (e.g., "R/W"), or `None`.
    amino_acids: Option<&'static str>,
    /// Expected exon number string (e.g., "7/11"), or `None`.
    exon: Option<&'static str>,
    /// Expected intron number string (e.g., "7/10"), or `None`.
    intron: Option<&'static str>,
    /// Expected HGVS protein notation (e.g., "p.Arg248Trp"), or `None` for
    /// non-CDS variants. Values verified against VEP REST API (2026-04-10).
    hgvs_p: Option<&'static str>,
}

const VARIANTS: &[Expected] = &[
    // #1 — TP53 R248W: minus-strand missense, hotspot
    // VEP query: 17:7674221/A
    Expected {
        label: "TP53 R248W missense (minus)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_base: b'G',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(248),
        cds_position: Some(742),
        cdna_position: Some(884),
        codons: Some("Cgg/Tgg"),
        amino_acids: Some("R/W"),
        exon: Some("7/11"),
        intron: None,
        hgvs_p: Some("p.Arg248Trp"),
    },
    // #2 — TP53 R248R: synonymous at same codon (different alt allele)
    // VEP query: 17:7674221/T
    Expected {
        label: "TP53 R248R synonymous (minus)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_base: b'G',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["synonymous_variant"],
        impact: "LOW",
        protein_start: Some(248),
        cds_position: Some(742),
        cdna_position: Some(884),
        codons: Some("Cgg/Agg"),
        amino_acids: Some("R"),
        exon: Some("7/11"),
        intron: None,
        hgvs_p: Some("p.Arg248="),
    },
    // #3 — BRAF V600E: minus-strand missense, famous oncogene
    // VEP query: 7:140753336/T
    Expected {
        label: "BRAF V600E missense (minus)",
        chrom: "chr7",
        pos: 140_753_335,
        ref_base: b'A',
        alt_base: b'T',
        transcript: "NM_004333.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(600),
        cds_position: Some(1799),
        cdna_position: Some(2025),
        codons: Some("gTg/gAg"),
        amino_acids: Some("V/E"),
        exon: Some("15/18"),
        intron: None,
        hgvs_p: Some("p.Val600Glu"),
    },
    // #4 — CFTR I507F: plus-strand missense
    // VEP query: 7:117559590/T
    Expected {
        label: "CFTR I507F missense (plus)",
        chrom: "chr7",
        pos: 117_559_589,
        ref_base: b'A',
        alt_base: b'T',
        transcript: "NM_000492.4",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(507),
        cds_position: Some(1519),
        cdna_position: Some(1589),
        codons: Some("Atc/Ttc"),
        amino_acids: Some("I/F"),
        exon: Some("11/27"),
        intron: None,
        hgvs_p: Some("p.Ile507Phe"),
    },
    // #5 — TP53 R196X: stop gained (CGA → TGA)
    // VEP query: 17:7674945/A
    Expected {
        label: "TP53 R196X stop_gained (minus)",
        chrom: "chr17",
        pos: 7_674_944,
        ref_base: b'G',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["stop_gained"],
        impact: "HIGH",
        protein_start: Some(196),
        cds_position: Some(586),
        cdna_position: Some(728),
        codons: Some("Cga/Tga"),
        amino_acids: Some("R/*"),
        exon: Some("6/11"),
        intron: None,
        hgvs_p: Some("p.Arg196Ter"),
    },
    // #6 — TP53 start lost (ATG → GTG on coding strand)
    // VEP query: 17:7676594/C
    Expected {
        label: "TP53 start_lost (minus)",
        chrom: "chr17",
        pos: 7_676_593,
        ref_base: b'T',
        alt_base: b'C',
        transcript: "NM_000546.6",
        consequences: &["start_lost"],
        impact: "HIGH",
        protein_start: Some(1),
        cds_position: Some(1),
        cdna_position: Some(143),
        codons: Some("Atg/Gtg"),
        amino_acids: Some("M/V"),
        exon: Some("2/11"),
        intron: None,
        hgvs_p: Some("p.Met1?"),
    },
    // #7 — TP53 stop lost (TGA → TGT on coding strand)
    // VEP query: 17:7669609/A
    Expected {
        label: "TP53 stop_lost (minus)",
        chrom: "chr17",
        pos: 7_669_608,
        ref_base: b'T',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["stop_lost"],
        impact: "HIGH",
        protein_start: Some(394),
        cds_position: Some(1182),
        cdna_position: Some(1324),
        codons: Some("tgA/tgT"),
        amino_acids: Some("*/C"),
        exon: Some("11/11"),
        intron: None,
        // VEP: p.Ter394CysextTer9 — the actual extension distance is
        // computed by a 3'UTR stop-scan.
        hgvs_p: Some("p.Ter394CysextTer9"),
    },
    // #8 — TP53 stop retained (TGA → TAA on coding strand, both stop)
    // VEP query: 17:7669610/T
    Expected {
        label: "TP53 stop_retained (minus)",
        chrom: "chr17",
        pos: 7_669_609,
        ref_base: b'C',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["stop_retained_variant"],
        impact: "LOW",
        protein_start: Some(394),
        cds_position: Some(1181),
        cdna_position: Some(1323),
        codons: Some("tGa/tAa"),
        amino_acids: Some("*"),
        exon: Some("11/11"),
        intron: None,
        hgvs_p: Some("p.Ter394="),
    },
    // #9 — EGFR L858R: plus-strand missense, famous lung cancer variant
    // VEP query: 7:55191822/G
    Expected {
        label: "EGFR L858R missense (plus)",
        chrom: "chr7",
        pos: 55_191_821,
        ref_base: b'T',
        alt_base: b'G',
        transcript: "NM_005228.5",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(858),
        cds_position: Some(2573),
        cdna_position: Some(2834),
        codons: Some("cTg/cGg"),
        amino_acids: Some("L/R"),
        exon: Some("21/28"),
        intron: None,
        hgvs_p: Some("p.Leu858Arg"),
    },
    // #10 — TP53 S261R: split codon across exon 7/8 boundary (CDS 781-782
    //   in exon 7, CDS 783 in exon 8). Tests cross-exon codon extraction.
    //   Also triggers splice_region_variant (first base of exon 8).
    // VEP query: 17:7673837/T
    Expected {
        label: "TP53 S261R split codon (minus)",
        chrom: "chr17",
        pos: 7_673_836,
        ref_base: b'A',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["missense_variant", "splice_region_variant"],
        impact: "MODERATE",
        protein_start: Some(261),
        cds_position: Some(783),
        cdna_position: Some(925),
        codons: Some("agT/agA"),
        amino_acids: Some("S/R"),
        exon: Some("8/11"),
        intron: None,
        hgvs_p: Some("p.Ser261Arg"),
    },
    // #11 — Splice donor +1 (intronic, exon 7 donor side)
    // VEP query: 17:7674180/T
    Expected {
        label: "TP53 splice_donor +1 (minus)",
        chrom: "chr17",
        pos: 7_674_179,
        ref_base: b'C',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["splice_donor_variant"],
        impact: "HIGH",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("7/10"),
        hgvs_p: None,
    },
    // #12 — Splice acceptor -1 (intronic, exon 7 acceptor side from intron 6)
    // VEP query: 17:7674291/T
    Expected {
        label: "TP53 splice_acceptor -1 (minus)",
        chrom: "chr17",
        pos: 7_674_290,
        ref_base: b'C',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["splice_acceptor_variant"],
        impact: "HIGH",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("6/10"),
        hgvs_p: None,
    },
    // #13 — Splice region + intron (acceptor -3 from exon 7)
    // VEP returns: splice_region_variant + splice_polypyrimidine_tract_variant
    //   + intron_variant. vareffect intentionally omits the granular
    //   splice_polypyrimidine_tract_variant term (see known divergence #1).
    // VEP query: 17:7674293/G
    Expected {
        label: "TP53 splice_region + intron (minus)",
        chrom: "chr17",
        pos: 7_674_292,
        ref_base: b'A',
        alt_base: b'G',
        transcript: "NM_000546.6",
        consequences: &["intron_variant", "splice_region_variant"],
        impact: "LOW",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("6/10"),
        hgvs_p: None,
    },
    // #14 — Exonic splice region + synonymous (3 bases from exon 7 boundary)
    // VEP query: 17:7674183/A
    Expected {
        label: "TP53 exonic splice_region + synonymous (minus)",
        chrom: "chr17",
        pos: 7_674_182,
        ref_base: b'G',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["splice_region_variant", "synonymous_variant"],
        impact: "LOW",
        protein_start: Some(260),
        cds_position: Some(780),
        cdna_position: Some(922),
        codons: Some("tcC/tcT"),
        amino_acids: Some("S"),
        exon: Some("7/11"),
        intron: None,
        hgvs_p: Some("p.Ser260="),
    },
    // #15 — Deep intronic (>8 bp from nearest exon boundary)
    // VEP query: 17:7674050/A
    Expected {
        label: "TP53 deep intron (minus)",
        chrom: "chr17",
        pos: 7_674_049,
        ref_base: b'G',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["intron_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: Some("7/10"),
        hgvs_p: None,
    },
    // #16 — 5' UTR (exon 1, before CDS start)
    // VEP query: 17:7687400/A
    Expected {
        label: "TP53 5_prime_UTR (minus)",
        chrom: "chr17",
        pos: 7_687_399,
        ref_base: b'G',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["5_prime_UTR_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: Some(91),
        codons: None,
        amino_acids: None,
        exon: Some("1/11"),
        intron: None,
        hgvs_p: None,
    },
    // #17 — 3' UTR (exon 11, after stop codon)
    // VEP query: 17:7669580/A
    Expected {
        label: "TP53 3_prime_UTR (minus)",
        chrom: "chr17",
        pos: 7_669_579,
        ref_base: b'T',
        alt_base: b'A',
        transcript: "NM_000546.6",
        consequences: &["3_prime_UTR_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: Some(1353),
        codons: None,
        amino_acids: None,
        exon: Some("11/11"),
        intron: None,
        hgvs_p: None,
    },
    // #18 — Upstream gene variant (>5' of transcript, within 5 kb)
    // VEP query: 17:7688000/G
    Expected {
        label: "TP53 upstream (minus)",
        chrom: "chr17",
        pos: 7_687_999,
        ref_base: b'A',
        alt_base: b'G',
        transcript: "NM_000546.6",
        consequences: &["upstream_gene_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: None,
        hgvs_p: None,
    },
    // #19 — Downstream gene variant (>3' of transcript, within 5 kb)
    // VEP query: 17:7668000/T
    Expected {
        label: "TP53 downstream (minus)",
        chrom: "chr17",
        pos: 7_667_999,
        ref_base: b'G',
        alt_base: b'T',
        transcript: "NM_000546.6",
        consequences: &["downstream_gene_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: None,
        codons: None,
        amino_acids: None,
        exon: None,
        intron: None,
        hgvs_p: None,
    },
    // #20 — Non-coding transcript exon variant (TERC, telomerase RNA)
    // VEP query: 3:169764900/T
    // NOTE: This variant requires NR_001566.3 in the transcript store.
    //   If the store was built without NR_* transcripts, this test will
    //   report "transcript not found" and be counted as a skip, not a fail.
    Expected {
        label: "TERC non_coding_transcript_exon (minus)",
        chrom: "chr3",
        pos: 169_764_899,
        ref_base: b'A',
        alt_base: b'T',
        transcript: "NR_001566.3",
        consequences: &["non_coding_transcript_exon_variant"],
        impact: "MODIFIER",
        protein_start: None,
        cds_position: None,
        cdna_position: Some(161),
        codons: None,
        amino_acids: None,
        exon: Some("1/1"),
        intron: None,
        hgvs_p: None,
    },
];

/// Compare a single variant's `annotate_snv` output against the expected VEP
/// values. Returns a list of field-level mismatches (empty = pass).
fn check_variant(
    store: &TranscriptStore,
    fasta: &FastaReader,
    exp: &Expected,
) -> Result<Vec<String>, String> {
    let (tx, idx) = store
        .get_by_accession(exp.transcript)
        .ok_or_else(|| format!("transcript {} not found in store", exp.transcript))?;

    let result = annotate_snv(
        exp.chrom,
        exp.pos,
        exp.ref_base,
        exp.alt_base,
        tx,
        idx,
        fasta,
    )
    .map_err(|e| format!("annotate_snv error: {e}"))?;

    let mut mismatches: Vec<String> = Vec::new();

    // Consequence terms (compare as sorted sets).
    let actual_csq: BTreeSet<&str> = result.consequences.iter().map(|c| c.as_str()).collect();
    let expected_csq: BTreeSet<&str> = exp.consequences.iter().copied().collect();
    if actual_csq != expected_csq {
        mismatches.push(format!(
            "consequences: expected {:?}, got {:?}",
            expected_csq, actual_csq,
        ));
    }

    // IMPACT.
    let actual_impact = format!("{}", result.impact);
    if actual_impact != exp.impact {
        mismatches.push(format!(
            "impact: expected {}, got {}",
            exp.impact, actual_impact,
        ));
    }

    // Protein position.
    if result.protein_start != exp.protein_start {
        mismatches.push(format!(
            "protein_start: expected {:?}, got {:?}",
            exp.protein_start, result.protein_start,
        ));
    }

    // CDS position.
    if result.cds_position != exp.cds_position {
        mismatches.push(format!(
            "cds_position: expected {:?}, got {:?}",
            exp.cds_position, result.cds_position,
        ));
    }

    // cDNA position.
    if result.cdna_position != exp.cdna_position {
        mismatches.push(format!(
            "cdna_position: expected {:?}, got {:?}",
            exp.cdna_position, result.cdna_position,
        ));
    }

    // Codons.
    if result.codons.as_deref() != exp.codons {
        mismatches.push(format!(
            "codons: expected {:?}, got {:?}",
            exp.codons,
            result.codons.as_deref(),
        ));
    }

    // Amino acids.
    if result.amino_acids.as_deref() != exp.amino_acids {
        mismatches.push(format!(
            "amino_acids: expected {:?}, got {:?}",
            exp.amino_acids,
            result.amino_acids.as_deref(),
        ));
    }

    // Exon number.
    if result.exon.as_deref() != exp.exon {
        mismatches.push(format!(
            "exon: expected {:?}, got {:?}",
            exp.exon,
            result.exon.as_deref(),
        ));
    }

    // Intron number.
    if result.intron.as_deref() != exp.intron {
        mismatches.push(format!(
            "intron: expected {:?}, got {:?}",
            exp.intron,
            result.intron.as_deref(),
        ));
    }

    // HGVS protein notation.
    if result.hgvs_p.as_deref() != exp.hgvs_p {
        mismatches.push(format!(
            "hgvs_p: expected {:?}, got {:?}",
            exp.hgvs_p,
            result.hgvs_p.as_deref(),
        ));
    }

    Ok(mismatches)
}

/// Run all 20 VEP concordance variants. Each variant is annotated with
/// [`annotate_snv`] and compared field-by-field against VEP ground truth.
///
/// Transcripts that are missing from the store are reported as skips, not
/// failures — the TERC non-coding transcript (NR_001566.3) may not be
/// present depending on the store build configuration.
#[test]
#[ignore]
fn vep_concordance_grch38_snv_20() {
    let store = load_store();
    let fasta = load_fasta();

    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;
    let mut failures: Vec<String> = Vec::new();

    for (i, exp) in VARIANTS.iter().enumerate() {
        let num = i + 1;
        match check_variant(&store, &fasta, exp) {
            Err(msg) if msg.contains("not found in store") => {
                eprintln!("  [{num:>2}/20] SKIP  {}: {msg}", exp.label);
                skip += 1;
            }
            Err(msg) => {
                eprintln!("  [{num:>2}/20] FAIL  {}: {msg}", exp.label);
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(mismatches) if mismatches.is_empty() => {
                eprintln!("  [{num:>2}/20] PASS  {}", exp.label);
                pass += 1;
            }
            Ok(mismatches) => {
                eprintln!("  [{num:>2}/20] FAIL  {}", exp.label);
                let detail = mismatches.join("\n           ");
                failures.push(format!("#{num} {}:\n           {detail}", exp.label));
                fail += 1;
            }
        }
    }

    eprintln!(
        "\n=== VEP Concordance: {pass}/{} passed, {skip} skipped ===",
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
