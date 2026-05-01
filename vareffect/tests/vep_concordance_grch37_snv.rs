//! VEP concordance spot-check for SNV consequence prediction on GRCh37.
//!
//! Mirrors `vep_concordance_snv.rs` (GRCh38) — same self-contained file
//! pattern, same `Expected`-struct + `VARIANTS` array + per-variant
//! comparison loop. Variants validated against Ensembl VEP REST API
//! responses captured at `https://grch37.rest.ensembl.org` (release
//! 115/116, queried 2026-04-30 with `refseq=1&hgvs=1&shift_hgvs=1&numbers=1`).
//!
//! Run with:
//! ```bash
//! GRCH37_FASTA=/abs/path/to/data/vareffect/GRCh37.bin \
//! GRCH37_TRANSCRIPTS=/abs/path/to/data/vareffect/transcript_models_grch37.bin \
//!   cargo test -p vareffect --release -- --ignored vep_concordance_grch37_snv
//! ```
//!
//! # Transcript-version tolerance
//!
//! Ensembl VEP's GRCh37 REST cache is updated independently of NCBI's
//! GRCh37.p13 RefSeq snapshot, so VEP returns versioned accessions like
//! `NM_000546.6` while our store carries `NM_000546.5`. The lookup path
//! tries strict accession-with-version first and falls back to
//! accession-without-version on miss. HGVS values are compared on the
//! bare `c.`/`p.` suffix (accession prefix stripped) for the same reason.
//!
//! # Known intentional divergences
//!
//! Same as the GRCh38 spot-check (`vep_concordance_snv.rs:13-29`):
//! granular splice region terms and `non_coding_transcript_variant` are
//! VEP-specific extensions that vareffect deliberately does not emit.
//!
//! # Adding fixtures
//!
//! Same workflow as the GRCh38 spot-checks: query the Ensembl VEP REST
//! API for each variant (using `https://grch37.rest.ensembl.org/vep/human/region/...
//! ?refseq=1&hgvs=1&shift_hgvs=1&numbers=1`), paste the response into
//! the fixture comment, hand-write the `Expected { … }` block, and
//! re-run the test. **Important**: VEP REST ignores the input REF
//! allele and infers from its own reference, so confirm the plus-strand
//! REF against ClinVar / UCSC before committing — a mismatched REF hits
//! a `RefMismatch` error from `annotate`.

use std::collections::BTreeSet;
use std::path::PathBuf;

use vareffect::{Assembly, ConsequenceResult, VarEffect};

/// Build a `VarEffect` instance with the GRCh37 store + reference loaded
/// from env vars.
fn open_grch37() -> VarEffect {
    let fasta = std::env::var("GRCH37_FASTA").expect(
        "GRCH37_FASTA env var must point at the flat-binary genome \
         (e.g. data/vareffect/GRCh37.bin). Run `vareffect setup --assembly grch37` first.",
    );
    let transcripts = std::env::var("GRCH37_TRANSCRIPTS")
        .expect("GRCH37_TRANSCRIPTS env var must point at the transcript_models_grch37.bin");
    VarEffect::builder()
        .with_grch37(&PathBuf::from(transcripts), &PathBuf::from(fasta))
        .expect("loading GRCh37")
        .build()
        .expect("building VarEffect")
}

/// Strip the `accession:` prefix from a full HGVS string. VEP returns
/// `NM_000546.6:c.524G>T`; vareffect returns `NM_000546.5:c.524G>T`. We
/// compare only the bare `c.`/`p.` suffix.
fn strip_hgvs_accession(full: &str) -> &str {
    full.split_once(':').map(|(_, rest)| rest).unwrap_or(full)
}

/// Find the per-transcript consequence row matching `expected_accession`.
/// Tries strict match (versioned) first, then falls back to
/// accession-without-version. Returns the matched row and a flag
/// indicating whether the fallback was used.
fn find_result_for_transcript<'a>(
    consequences: &'a [ConsequenceResult],
    expected_accession: &str,
) -> Option<(&'a ConsequenceResult, bool)> {
    if let Some(r) = consequences
        .iter()
        .find(|r| r.transcript == expected_accession)
    {
        return Some((r, false));
    }
    let base = expected_accession.split('.').next().unwrap_or("");
    consequences
        .iter()
        .find(|r| r.transcript.split('.').next().unwrap_or("") == base)
        .map(|r| (r, true))
}

/// Expected output for one variant–transcript pair, derived from the VEP
/// REST API. Coordinate fields use VEP's conventions (1-based for
/// CDS/cDNA/protein); `pos` is 0-based per vareffect's annotate API.
struct Expected {
    label: &'static str,
    chrom: &'static str,
    /// 0-based genomic position.
    pos: u64,
    ref_allele: &'static [u8],
    alt_allele: &'static [u8],
    /// RefSeq accession with VEP's reported version. The lookup falls
    /// back to accession-without-version if the store carries a different
    /// version.
    transcript: &'static str,
    consequences: &'static [&'static str],
    impact: &'static str,
    protein_start: Option<u32>,
    cds_position: Option<u32>,
    cdna_position: Option<u32>,
    codons: Option<&'static str>,
    amino_acids: Option<&'static str>,
    exon: Option<&'static str>,
    intron: Option<&'static str>,
    /// Bare HGVS `c.` notation (accession prefix stripped).
    hgvs_c: Option<&'static str>,
    /// Bare HGVS `p.` notation (accession prefix stripped). `None` for
    /// non-CDS variants.
    hgvs_p: Option<&'static str>,
}

const VARIANTS: &[Expected] = &[
    // #1 — TP53 R175L: minus-strand missense at the c.524 position.
    // Plus-strand reference at chr17:7578406 is C; the variant
    // `chr17:7578406 C>A` corresponds to coding c.524G>T on the minus
    // strand, CGC→CTC, p.Arg175Leu. (Distinct from the R175H hotspot,
    // which would be chr17:7578406 C>T → c.524G>A.)
    // VEP query: https://grch37.rest.ensembl.org/vep/human/region/17:7578406-7578406/A?refseq=1&hgvs=1&shift_hgvs=1&numbers=1
    Expected {
        label: "TP53 R175L missense (minus)",
        chrom: "chr17",
        pos: 7_578_405,
        ref_allele: b"C",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(175),
        cds_position: Some(524),
        cdna_position: Some(666),
        codons: Some("cGc/cTc"),
        amino_acids: Some("R/L"),
        exon: Some("5/11"),
        intron: None,
        hgvs_c: Some("c.524G>T"),
        hgvs_p: Some("p.Arg175Leu"),
    },
    // #2 — BRAF V600E: minus-strand missense, the canonical activating
    // BRAF mutation in melanoma / colorectal / thyroid cancers.
    // VEP query: https://grch37.rest.ensembl.org/vep/human/region/7:140453136-140453136/T?refseq=1&hgvs=1&shift_hgvs=1&numbers=1
    Expected {
        label: "BRAF V600E missense (minus)",
        chrom: "chr7",
        pos: 140_453_135,
        ref_allele: b"A",
        alt_allele: b"T",
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
        hgvs_c: Some("c.1799T>A"),
        hgvs_p: Some("p.Val600Glu"),
    },
    // #3 — TP53 R273S: minus-strand missense at c.817C>A. Plus-strand
    // reference at chr17:7577121 is G; the variant `chr17:7577121 G>T`
    // corresponds to coding c.817C>A on the minus strand, CGT→AGT,
    // p.Arg273Ser.
    // VEP query: https://grch37.rest.ensembl.org/vep/human/region/17:7577121-7577121/T?refseq=1&hgvs=1&shift_hgvs=1&numbers=1
    Expected {
        label: "TP53 R273S missense (minus)",
        chrom: "chr17",
        pos: 7_577_120,
        ref_allele: b"G",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        consequences: &["missense_variant"],
        impact: "MODERATE",
        protein_start: Some(273),
        cds_position: Some(817),
        cdna_position: Some(959),
        codons: Some("Cgt/Agt"),
        amino_acids: Some("R/S"),
        exon: Some("8/11"),
        intron: None,
        hgvs_c: Some("c.817C>A"),
        hgvs_p: Some("p.Arg273Ser"),
    },
];

fn check_variant(ve: &VarEffect, exp: &Expected) -> Result<Vec<String>, String> {
    let result = ve
        .annotate(
            Assembly::GRCh37,
            exp.chrom,
            exp.pos,
            exp.ref_allele,
            exp.alt_allele,
        )
        .map_err(|e| format!("annotate error: {e}"))?;

    let (row, version_fallback) = find_result_for_transcript(&result.consequences, exp.transcript)
        .ok_or_else(|| format!("transcript {} not found in store", exp.transcript))?;

    if version_fallback {
        eprintln!(
            "      version fallback: expected {}, store has {}",
            exp.transcript, row.transcript,
        );
    }

    let mut mismatches: Vec<String> = Vec::new();

    let actual_csq: BTreeSet<&str> = row.consequences.iter().map(|c| c.as_str()).collect();
    let expected_csq: BTreeSet<&str> = exp.consequences.iter().copied().collect();
    if actual_csq != expected_csq {
        mismatches.push(format!(
            "consequences: expected {expected_csq:?}, got {actual_csq:?}",
        ));
    }

    let actual_impact = format!("{}", row.impact);
    if actual_impact != exp.impact {
        mismatches.push(format!(
            "impact: expected {}, got {}",
            exp.impact, actual_impact,
        ));
    }
    if row.protein_start != exp.protein_start {
        mismatches.push(format!(
            "protein_start: expected {:?}, got {:?}",
            exp.protein_start, row.protein_start,
        ));
    }
    if row.cds_position != exp.cds_position {
        mismatches.push(format!(
            "cds_position: expected {:?}, got {:?}",
            exp.cds_position, row.cds_position,
        ));
    }
    if row.cdna_position != exp.cdna_position {
        mismatches.push(format!(
            "cdna_position: expected {:?}, got {:?}",
            exp.cdna_position, row.cdna_position,
        ));
    }
    if row.codons.as_deref() != exp.codons {
        mismatches.push(format!(
            "codons: expected {:?}, got {:?}",
            exp.codons,
            row.codons.as_deref(),
        ));
    }
    if row.amino_acids.as_deref() != exp.amino_acids {
        mismatches.push(format!(
            "amino_acids: expected {:?}, got {:?}",
            exp.amino_acids,
            row.amino_acids.as_deref(),
        ));
    }
    if row.exon.as_deref() != exp.exon {
        mismatches.push(format!(
            "exon: expected {:?}, got {:?}",
            exp.exon,
            row.exon.as_deref(),
        ));
    }
    if row.intron.as_deref() != exp.intron {
        mismatches.push(format!(
            "intron: expected {:?}, got {:?}",
            exp.intron,
            row.intron.as_deref(),
        ));
    }
    let actual_c = row.hgvs_c.as_deref().map(strip_hgvs_accession);
    if actual_c != exp.hgvs_c {
        mismatches.push(format!(
            "hgvs_c: expected {:?}, got {actual_c:?}",
            exp.hgvs_c,
        ));
    }
    let actual_p = row.hgvs_p.as_deref().map(strip_hgvs_accession);
    if actual_p != exp.hgvs_p {
        mismatches.push(format!(
            "hgvs_p: expected {:?}, got {actual_p:?}",
            exp.hgvs_p,
        ));
    }

    Ok(mismatches)
}

#[test]
#[ignore]
fn vep_concordance_grch37_snv() {
    let ve = open_grch37();
    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;
    let mut failures: Vec<String> = Vec::new();
    let total = VARIANTS.len();

    for (i, exp) in VARIANTS.iter().enumerate() {
        let num = i + 1;
        match check_variant(&ve, exp) {
            Err(msg) if msg.contains("not found in store") => {
                eprintln!("  [{num:>2}/{total}] SKIP  {}: {msg}", exp.label);
                skip += 1;
            }
            Err(msg) => {
                eprintln!("  [{num:>2}/{total}] FAIL  {}: {msg}", exp.label);
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(mismatches) if mismatches.is_empty() => {
                eprintln!("  [{num:>2}/{total}] PASS  {}", exp.label);
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

    eprintln!("\n=== VEP Concordance (GRCh37 SNV): {pass}/{total} passed, {skip} skipped ===");
    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  {f}");
        }
    }
    assert_eq!(
        fail, 0,
        "{fail} of {total} variants failed concordance check"
    );
}
