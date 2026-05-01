//! VEP concordance spot-check for HGVS c. notation.
//!
//! 20 variants validated against Ensembl VEP REST API (GRCh38, queried
//! 2026-04-09 with `refseq=1&hgvs=1`). `#[ignore]`-gated because the test
//! requires the transcript store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_concordance_hgvs
//! ```
//!
//! Uses the top-level [`annotate()`] dispatcher, exercising the full VCF-input
//! pipeline (REF verification -> `trim_alleles` -> type dispatch -> HGVS c.
//! formatting via `hgvs_c.rs`).
//!
//! # Categories
//!
//! A. **CDS substitutions** (#1-4): plus-strand and minus-strand missense.
//! B. **UTR substitutions** (#5-6): 5'UTR (c.-N) and 3'UTR (c.*N).
//! C. **Intronic substitutions** (#7-10): splice donor/acceptor, deep intron,
//!    5'UTR intronic (c.-N-M notation).
//! D. **CDS deletions** (#11-12): single-base and multi-base.
//! E. **Duplications** (#13-14): single-base minus-strand, multi-base
//!    plus-strand.
//! F. **Insertion (non-dup)** (#15): multi-base insertion that is NOT a dup.
//! G. **Delins / MNV** (#16-17): 2bp and 3bp multi-nucleotide variants.
//! H. **Frameshift deletion** (#18): minus-strand 2bp deletion.
//! I. **Non-coding transcript** (#19): n. prefix on NR_* transcript.
//! J. **Boundary / edge** (#20): intronic deletion near splice donor.
//!
//! # BRCA1 NM_007294.4 note
//!
//! VEP no longer carries NM_007294.4 (superseded by NM_001407*.1 isoforms).
//! Variants #13 and #18 use the universally agreed-upon clinical HGVS from
//! ClinVar / HGMD. Genomic coordinates are verified in the existing indel
//! concordance test.

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

/// Expected HGVS c./n. output for one variant-transcript pair.
struct Expected {
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
    /// Full HGVS c./n. string from VEP (e.g., "NM_000546.6:c.742C>T").
    expected_hgvs_c: &'static str,
}

// ---------------------------------------------------------------------------
// Ground truth from VEP REST API (GRCh38, refseq=1, hgvs=1)
// ---------------------------------------------------------------------------

const VARIANTS: &[Expected] = &[
    // -----------------------------------------------------------------------
    // Category A: CDS substitutions (4 variants)
    // -----------------------------------------------------------------------

    // #1 -- TP53 R248W: minus-strand missense
    // VEP query: 17:7674221:7674221/A
    Expected {
        label: "TP53 R248W c.742C>T (minus)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.742C>T",
    },
    // #2 -- BRAF V600E: minus-strand missense
    // VEP query: 7:140753336:140753336/T
    Expected {
        label: "BRAF V600E c.1799T>A (minus)",
        chrom: "chr7",
        pos: 140_753_335,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NM_004333.6",
        expected_hgvs_c: "NM_004333.6:c.1799T>A",
    },
    // #3 -- EGFR L858R: plus-strand missense
    // VEP query: 7:55191822:55191822/G
    Expected {
        label: "EGFR L858R c.2573T>G (plus)",
        chrom: "chr7",
        pos: 55_191_821,
        ref_allele: b"T",
        alt_allele: b"G",
        transcript: "NM_005228.5",
        expected_hgvs_c: "NM_005228.5:c.2573T>G",
    },
    // #4 -- CFTR I507F: plus-strand missense
    // VEP query: 7:117559590:117559590/T
    Expected {
        label: "CFTR I507F c.1519A>T (plus)",
        chrom: "chr7",
        pos: 117_559_589,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        expected_hgvs_c: "NM_000492.4:c.1519A>T",
    },
    // -----------------------------------------------------------------------
    // Category B: UTR substitutions (2 variants)
    // -----------------------------------------------------------------------

    // #5 -- TP53 5'UTR: minus-strand, c.-N notation (negative from CDS start)
    // VEP query: 17:7687400:7687400/A
    Expected {
        label: "TP53 5'UTR c.-52C>T (minus)",
        chrom: "chr17",
        pos: 7_687_399,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.-52C>T",
    },
    // #6 -- TP53 3'UTR: minus-strand, c.*N notation (after stop codon)
    // VEP query: 17:7669580:7669580/A
    Expected {
        label: "TP53 3'UTR c.*29A>T (minus)",
        chrom: "chr17",
        pos: 7_669_579,
        ref_allele: b"T",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.*29A>T",
    },
    // -----------------------------------------------------------------------
    // Category C: Intronic substitutions (4 variants)
    // -----------------------------------------------------------------------

    // #7 -- TP53 splice donor +1: c.N+1 notation
    // VEP query: 17:7674180:7674180/T
    Expected {
        label: "TP53 splice donor c.782+1G>A (minus)",
        chrom: "chr17",
        pos: 7_674_179,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.782+1G>A",
    },
    // #8 -- TP53 splice acceptor -1: c.N-1 notation
    // VEP query: 17:7674291:7674291/T
    Expected {
        label: "TP53 splice acceptor c.673-1G>A (minus)",
        chrom: "chr17",
        pos: 7_674_290,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.673-1G>A",
    },
    // #9 -- TP53 deep intron: >8bp from nearest exon, c.N+M notation
    // VEP query: 17:7674050:7674050/A
    Expected {
        label: "TP53 deep intron c.782+131C>T (minus)",
        chrom: "chr17",
        pos: 7_674_049,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.782+131C>T",
    },
    // #10 -- TP53 5'UTR intron: intron 1 lies entirely in 5'UTR (exon 1 is
    //   non-coding, CDS starts in exon 2). Tests c.-N-M notation where the
    //   anchor exon boundary is in the 5'UTR and the offset is into the intron.
    // VEP query: 17:7676650:7676650/T (given_ref=A, used_ref=A)
    Expected {
        label: "TP53 5'UTR intron c.-28-28T>A (minus)",
        chrom: "chr17",
        pos: 7_676_649,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.-28-28T>A",
    },
    // -----------------------------------------------------------------------
    // Category D: CDS deletions (2 variants)
    // -----------------------------------------------------------------------

    // #11 -- BRCA2 c.1813del: plus-strand single-base deletion (frameshift)
    // VEP query: NM_000059.4:c.1813del
    Expected {
        label: "BRCA2 c.1813del single-base (plus)",
        chrom: "chr13",
        pos: 32_333_289,
        ref_allele: b"AA",
        alt_allele: b"A",
        transcript: "NM_000059.4",
        expected_hgvs_c: "NM_000059.4:c.1813del",
    },
    // #12 -- CFTR deltaF508: plus-strand 3bp inframe deletion
    // VEP query: NM_000492.4:c.1521_1523del
    Expected {
        label: "CFTR deltaF508 c.1521_1523del (plus)",
        chrom: "chr7",
        pos: 117_559_590,
        ref_allele: b"TCTT",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        expected_hgvs_c: "NM_000492.4:c.1521_1523del",
    },
    // -----------------------------------------------------------------------
    // Category E: Duplications (2 variants)
    // -----------------------------------------------------------------------

    // #13 -- BRCA1 c.5266dup (5382insC): minus-strand single-base dup.
    //   CRITICAL: Ashkenazi founder variant. Tests `is_duplication()` on
    //   minus-strand (5' = higher genomic coords). VCF inserts G on plus
    //   strand; on coding strand the duplicated base is C in a tandem CC.
    //   NM_007294.4 is no longer in VEP RefSeq — notation from ClinVar.
    //
    //   HGVS 3' normalization: both c.5265 and c.5266 are C on the coding
    //   strand. The dup is shifted to c.5266 (the 3'-most position),
    //   matching VEP's `--shift_hgvs` default.
    Expected {
        label: "BRCA1 c.5266dup single-base (minus, 3' shifted)",
        chrom: "chr17",
        pos: 43_057_062,
        ref_allele: b"G",
        alt_allele: b"GG",
        transcript: "NM_007294.4",
        expected_hgvs_c: "NM_007294.4:c.5266dup",
    },
    // #14 -- ERBB2 exon 20 12bp insertion: plus-strand.
    //   The 12bp inserted sequence matches the 3' flanking reference at
    //   [ins_pos, ins_pos+12).
    //
    //   HGVS 3' normalization: insertion point shifts +12bp, at which point
    //   the 5' flanking matches the inserted sequence and it becomes a dup.
    //   Matches VEP's `--shift_hgvs` default.
    Expected {
        label: "ERBB2 c.2313_2324dup 12bp (plus, 3' shifted)",
        chrom: "chr17",
        pos: 39_724_729,
        ref_allele: b"C",
        alt_allele: b"CATACGTGATGGC",
        transcript: "NM_004448.4",
        expected_hgvs_c: "NM_004448.4:c.2313_2324dup",
    },
    // -----------------------------------------------------------------------
    // Category F: Insertion (non-dup) (1 variant)
    // -----------------------------------------------------------------------

    // #15 -- BRCA2 c.5946_5947insAA: plus-strand 2bp insertion that is NOT
    //   a duplication. The 5' flanking 2 bases are "GT" (from indel test #24
    //   REF at that region), which does not match "AA".
    // VEP query: NM_000059.4:c.5946_5947insAA
    Expected {
        label: "BRCA2 c.5946_5947insAA non-dup (plus)",
        chrom: "chr13",
        pos: 32_340_300,
        ref_allele: b"T",
        alt_allele: b"TAA",
        transcript: "NM_000059.4",
        expected_hgvs_c: "NM_000059.4:c.5946_5947insAA",
    },
    // -----------------------------------------------------------------------
    // Category G: Delins / MNV (2 variants)
    // -----------------------------------------------------------------------

    // #16 -- TP53 c.742_743delinsTT: minus-strand 2bp MNV
    // VEP query: NM_000546.6:c.742_743delinsTT
    Expected {
        label: "TP53 c.742_743delinsTT MNV (minus)",
        chrom: "chr17",
        pos: 7_674_219,
        ref_allele: b"CG",
        alt_allele: b"AA",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.742_743delinsTT",
    },
    // #17 -- TP53 c.742_744delinsAAA: minus-strand 3bp MNV (full codon)
    // VEP query: NM_000546.6:c.742_744delinsAAA
    Expected {
        label: "TP53 c.742_744delinsAAA MNV (minus)",
        chrom: "chr17",
        pos: 7_674_218,
        ref_allele: b"CCG",
        alt_allele: b"TTT",
        transcript: "NM_000546.6",
        expected_hgvs_c: "NM_000546.6:c.742_744delinsAAA",
    },
    // -----------------------------------------------------------------------
    // Category H: Frameshift deletion (1 variant)
    // -----------------------------------------------------------------------

    // #18 -- BRCA1 c.68_69del: minus-strand 2bp frameshift deletion.
    //   NM_007294.4 is no longer in VEP RefSeq — notation from ClinVar.
    // VEP query: NM_007294.4:c.68_69del (returns NM_001407*.1:c.68_69del)
    Expected {
        label: "BRCA1 c.68_69del frameshift (minus)",
        chrom: "chr17",
        pos: 43_124_026,
        ref_allele: b"ACT",
        alt_allele: b"A",
        transcript: "NM_007294.4",
        expected_hgvs_c: "NM_007294.4:c.68_69del",
    },
    // -----------------------------------------------------------------------
    // Category I: Non-coding transcript (1 variant)
    // -----------------------------------------------------------------------

    // #19 -- TERC n.161T>A: minus-strand non-coding (telomerase RNA).
    //   Tests `n.` prefix instead of `c.` for NR_* transcripts.
    //   NR_001566.3 may not be in the transcript store — counted as skip.
    // VEP query: 3:169764900:169764900/T
    Expected {
        label: "TERC n.161T>A non-coding (minus)",
        chrom: "chr3",
        pos: 169_764_899,
        ref_allele: b"A",
        alt_allele: b"T",
        transcript: "NR_001566.3",
        expected_hgvs_c: "NR_001566.3:n.161T>A",
    },
    // -----------------------------------------------------------------------
    // Category J: Boundary / edge case (1 variant)
    // -----------------------------------------------------------------------

    // #20 -- BRCA2 c.681+1_681+3del: plus-strand pure intronic deletion
    //   at splice donor site. Tests intronic HGVS notation on both endpoints
    //   of a multi-base deletion: c.N+start_N+end format.
    // VEP query: NM_000059.4:c.681+1_681+3del
    Expected {
        label: "BRCA2 c.681+1_681+3del splice donor (plus)",
        chrom: "chr13",
        pos: 32_329_491,
        ref_allele: b"TGTA",
        alt_allele: b"T",
        transcript: "NM_000059.4",
        expected_hgvs_c: "NM_000059.4:c.681+1_681+3del",
    },
];

// ---------------------------------------------------------------------------
// Comparison logic
// ---------------------------------------------------------------------------

/// Compare a single variant's `annotate()` HGVS c. output against the expected
/// VEP value. Returns `Ok(None)` on match, `Ok(Some(mismatch))` on field
/// mismatch, or `Err(msg)` on transcript-not-found / annotate error.
fn check_variant(ve: &VarEffect, exp: &Expected) -> Result<Option<String>, String> {
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
                "transcript {} not found in results (found: {found:?})",
                exp.transcript,
            )
        })?;

    let actual = result.hgvs_c.as_deref();
    let expected = Some(exp.expected_hgvs_c);

    if actual == expected {
        Ok(None)
    } else {
        Ok(Some(format!(
            "hgvs_c mismatch: expected {:?}, got {:?}",
            expected, actual,
        )))
    }
}

// ---------------------------------------------------------------------------
// Test runner
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn vep_concordance_grch38_hgvs_20() {
    let ve = VarEffect::builder()
        .with_handles(Assembly::GRCh38, load_store(), load_fasta())
        .expect("matching assemblies")
        .build()
        .expect("builder");

    let total = VARIANTS.len();
    let mut pass = 0u32;
    let mut fail = 0u32;
    let mut skip = 0u32;
    let mut failures: Vec<String> = Vec::new();

    for (i, exp) in VARIANTS.iter().enumerate() {
        let num = i + 1;
        match check_variant(&ve, exp) {
            Err(msg) if msg.contains("not found") => {
                eprintln!("  [{num:>2}/{total}] SKIP  {}: {msg}", exp.label,);
                skip += 1;
            }
            Err(msg) => {
                eprintln!("  [{num:>2}/{total}] FAIL  {}: {msg}", exp.label,);
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(None) => {
                eprintln!("  [{num:>2}/{total}] PASS  {}", exp.label,);
                pass += 1;
            }
            Ok(Some(mismatch)) => {
                eprintln!("  [{num:>2}/{total}] FAIL  {}: {mismatch}", exp.label,);
                failures.push(format!("#{num} {}: {mismatch}", exp.label,));
                fail += 1;
            }
        }
    }

    eprintln!(
        "\n=== VEP Concordance (HGVS c.): {pass}/{} passed, \
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
        "{fail} of {} variants failed HGVS c. concordance check",
        pass + fail,
    );
}
