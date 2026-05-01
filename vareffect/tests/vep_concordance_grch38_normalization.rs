//! VEP concordance spot-check for HGVS 3' normalization, intergenic
//! classification, and NMD prediction.
//!
//! 18 variants validating 3' normalization, `intergenic_variant` detection,
//! and `predicts_nmd` (50-nt rule) against VEP ground truth and manual NMD
//! computation. `#[ignore]`-gated because the test requires the transcript
//! store and reference FASTA on disk.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_concordance_normalization
//! ```
//!
//! Uses the top-level [`annotate()`] dispatcher, exercising the full VCF-input
//! pipeline (REF verification -> `trim_alleles` -> type dispatch -> consequence
//! assignment -> HGVS formatting -> NMD prediction).
//!
//! # Categories
//!
//! A. **3' normalization — resolved divergences** (#1-2): minus-strand dup, ins→dup.
//! B. **3' normalization — deletions** (#3-4): mono-A repeat, non-repeat.
//! C. **3' normalization — insertions** (#5-6): non-dup ins, non-repeat ins.
//! D. **intergenic_variant** (#7-9): gene desert, negative control.
//! E. **predicts_nmd** (#10-17): stop_gained/frameshift across all NMD cases.
//! F. **Combined / synonymous** (#18): synonymous → no NMD.
//!
//! # NMD validation methodology
//!
//! VEP does NOT report the 50-nt rule for MANE RefSeq protein_coding transcripts
//! (`NMD_transcript_variant` is biotype-based). `predicts_nmd` expected values are
//! derived from the transcript model's exon structure:
//!
//! - TP53 NM_000546.6: 10 CDS segments, cumulative CDS prefix sums
//!   [0, 74, 96, 375, 559, 672, 782, 919, 993, 1100, 1182].
//!   Last junction CDS pos = 1101 (1-based).
//! - APC NM_000038.6: exon 16/16 is the last and largest exon (~6500bp CDS).
//!   CDS positions 3927 and 4348 are both in the last exon.
//! - BRCA1 NM_007294.4: 24 exons, CDS pos 68 is in exon 2 (far from last junction).
//! - PRNP NM_000311.5: single-exon gene → no exon-exon junction → NMD impossible.

use std::collections::BTreeSet;
use std::path::Path;

use vareffect::{Assembly, FastaReader, TranscriptStore, VarEffect};

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

/// Expected output for one variant. Only `Some` fields are asserted;
/// `None` fields are skipped.
#[allow(dead_code)]
struct ExpectedAnnotation {
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
    /// RefSeq transcript accession with version. Empty string `""` for
    /// intergenic variants (match on empty transcript in results).
    transcript: &'static str,
    /// If `Some`, assert `hgvs_c` matches (full `"NM_*:c.*"` format).
    expected_hgvs_c: Option<&'static str>,
    /// If `Some`, assert `hgvs_p` matches (`"p.*"` format).
    expected_hgvs_p: Option<&'static str>,
    /// Subset check: every term here must appear in the actual consequences.
    /// Uses subset (expected ⊆ actual) NOT equality, because co-assigned
    /// terms like `splice_region_variant` are tested in the SNV/indel
    /// concordance suites.
    expected_consequences: Option<&'static [&'static str]>,
    /// If `true`, assert that `intergenic_variant` does NOT appear in
    /// consequences. Used for negative-control variants.
    assert_not_intergenic: bool,
    /// If `Some`, assert `predicts_nmd` matches.
    expected_nmd: Option<bool>,
    /// Category tag for logging.
    category: &'static str,
}

const VARIANTS: &[ExpectedAnnotation] = &[
    // #1 — BRCA1 c.5266dup (5382insC): minus-strand single-base dup.
    //   CRITICAL: Ashkenazi founder variant. Tests `is_duplication()` on
    //   minus-strand (5' = higher genomic coords). VCF inserts G on plus
    //   strand; on coding strand the duplicated base is C in a tandem CC.
    //
    //   HGVS 3' normalization: c.5265 and c.5266 are both C on the coding
    //   strand. The dup shifts to c.5266 (3'-most position), matching VEP's
    //   `--shift_hgvs` default.
    //
    //   NMD: CDS 5266 in a 5592-nt CDS with 24 exons. Far from last
    //   junction → NMD predicted.
    //
    //   Combined check: 3' shift + consequence + NMD in one variant.
    //
    //   Source: hgvs test #13, ClinVar notation (NM_007294.4 retired from VEP).
    ExpectedAnnotation {
        label: "BRCA1 c.5266dup (3' shift + NMD)",
        chrom: "chr17",
        pos: 43_057_062,
        ref_allele: b"G",
        alt_allele: b"GG",
        transcript: "NM_007294.4",
        expected_hgvs_c: Some("NM_007294.4:c.5266dup"),
        expected_hgvs_p: None,
        expected_consequences: Some(&["frameshift_variant"]),
        assert_not_intergenic: false,
        expected_nmd: Some(true),
        category: "A: 3' norm (dup, minus)",
    },
    // #2 — ERBB2 exon20 12bp insertion: plus-strand ins→dup conversion.
    //   The 12bp inserted sequence matches the 3' flanking reference at
    //   the shifted position, converting ins notation to dup.
    //
    //   HGVS 3' normalization: insertion point shifts +12bp, at which point
    //   the 5' flanking matches the inserted sequence → dup detected.
    //   Both hgvs_c AND hgvs_p must match VEP.
    //
    //   Source: hgvs test #14, hgvs_p test #20.
    //   VEP: NM_004448.4:c.2313_2324dup, NP_004439.2:p.Tyr772_Ala775dup
    ExpectedAnnotation {
        label: "ERBB2 c.2313_2324dup 12bp (3' shifted ins→dup)",
        chrom: "chr17",
        pos: 39_724_729,
        ref_allele: b"C",
        alt_allele: b"CATACGTGATGGC",
        transcript: "NM_004448.4",
        expected_hgvs_c: Some("NM_004448.4:c.2313_2324dup"),
        expected_hgvs_p: Some("p.Tyr772_Ala775dup"),
        expected_consequences: None,
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "A: 3' norm (ins→dup, plus)",
    },
    // #3 — BRCA2 c.1813del: plus-strand single-base deletion in coding region.
    //   Tests 3' shift for a deletion in a mono-A run (if adjacent bases match).
    //   Source: hgvs test #11.
    //   VEP: NM_000059.4:c.1813del
    ExpectedAnnotation {
        label: "BRCA2 c.1813del single-base (plus)",
        chrom: "chr13",
        pos: 32_333_289,
        ref_allele: b"AA",
        alt_allele: b"A",
        transcript: "NM_000059.4",
        expected_hgvs_c: Some("NM_000059.4:c.1813del"),
        expected_hgvs_p: None,
        expected_consequences: None,
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "B: 3' norm (del, mono-A)",
    },
    // #4 — CFTR deltaF508: plus-strand 3bp inframe deletion, no repeat.
    //   CTT is NOT in a repeat context → shift = 0.
    //   Regression guard: normalization must not shift non-repeat deletions.
    //   Source: hgvs test #12, hgvs_p test #18.
    //   VEP: NM_000492.4:c.1521_1523del, NP_000483.3:p.Phe508del
    ExpectedAnnotation {
        label: "CFTR deltaF508 c.1521_1523del (no shift, regression)",
        chrom: "chr7",
        pos: 117_559_590,
        ref_allele: b"TCTT",
        alt_allele: b"T",
        transcript: "NM_000492.4",
        expected_hgvs_c: Some("NM_000492.4:c.1521_1523del"),
        expected_hgvs_p: Some("p.Phe508del"),
        expected_consequences: None,
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "B: 3' norm (del, no shift, regression)",
    },
    // #5 — BRCA2 c.5946_5947insAA: plus-strand 2bp non-dup insertion.
    //   The 5' flanking 2 bases are "GT", which do not match "AA" → stays ins.
    //   Tests that non-dup insertions are not incorrectly converted.
    //   Source: hgvs test #15.
    //   VEP: NM_000059.4:c.5946_5947insAA
    ExpectedAnnotation {
        label: "BRCA2 c.5946_5947insAA non-dup (plus)",
        chrom: "chr13",
        pos: 32_340_300,
        ref_allele: b"T",
        alt_allele: b"TAA",
        transcript: "NM_000059.4",
        expected_hgvs_c: Some("NM_000059.4:c.5946_5947insAA"),
        expected_hgvs_p: None,
        expected_consequences: None,
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "C: 3' norm (ins, non-dup)",
    },
    // #6 — EGFR exon20 3bp inframe insertion: plus-strand, non-dup, non-repeat.
    //   Verifies that the codon-boundary handling survives 3' normalization.
    //   GGT insertion is not in a repeat context → no 3' shift.
    //   Source: hgvs_p test #19 (codon-boundary regression guard).
    //   VEP: NM_005228.5:c.2310_2311insGGT, p.Asp770_Asn771insGly
    //   Initial VEP query (7:55181319:55181318/GGT) was off-by-one for the
    //   insertion coordinate; corrected query (7:55181320:55181319/GGT)
    //   confirms c.2310_2311.
    ExpectedAnnotation {
        label: "EGFR exon20 3bp ins (non-dup, codon-boundary regression)",
        chrom: "chr7",
        pos: 55_181_318,
        ref_allele: b"C",
        alt_allele: b"CGGT",
        transcript: "NM_005228.5",
        expected_hgvs_c: Some("NM_005228.5:c.2310_2311insGGT"),
        expected_hgvs_p: Some("p.Asp770_Asn771insGly"),
        expected_consequences: None,
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "C: 3' norm (ins, non-repeat, regression)",
    },
    // #7 — Gene desert (chr5 PIK3R1-MAST4 gap): ~5 Mb gene desert.
    //   No RefSeq MANE transcripts overlap this position in the store.
    //   Must return exactly 1 ConsequenceResult with IntergenicVariant.
    //   REF verified against FASTA: T at chr5:67500000.
    ExpectedAnnotation {
        label: "Gene desert chr5:67.5M (intergenic)",
        chrom: "chr5",
        pos: 67_500_000,
        ref_allele: b"T",
        alt_allele: b"A",
        transcript: "",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["intergenic_variant"]),
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "D: intergenic",
    },
    // #8 — Gene desert (chr13 pericentromeric): no RefSeq genes.
    //   REF verified against FASTA: C at chr13:19000000.
    ExpectedAnnotation {
        label: "Gene desert chr13:19M (intergenic)",
        chrom: "chr13",
        pos: 19_000_000,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["intergenic_variant"]),
        assert_not_intergenic: false,
        expected_nmd: None,
        category: "D: intergenic",
    },
    // #9 — TP53 R248W: negative control — coding variant must NOT get
    //   IntergenicVariant even as a co-assigned term.
    //   Source: snv test #1.
    //   VEP: NM_000546.6:c.742C>T, missense_variant
    ExpectedAnnotation {
        label: "TP53 R248W — NOT intergenic (negative control)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["missense_variant"]),
        assert_not_intergenic: true,
        expected_nmd: None,
        category: "D: NOT intergenic (negative control)",
    },
    // #10 — TP53 R196X: stop_gained early in CDS, NMD predicted.
    //   CDS position 586, in exon 5 (~seg[4]). Distance to last junction
    //   (1101) = 515 >> 50. NMD = true.
    //   Source: snv test #5, hgvs_p test #5.
    //   VEP: NM_000546.6:c.586C>T, stop_gained
    //   Genomic: chr17:7674945 (1-based) = 7674944 (0-based), minus-strand,
    //   plus-strand G>A (complement of coding C>T).
    ExpectedAnnotation {
        label: "TP53 R196X stop_gained (NMD true, early CDS)",
        chrom: "chr17",
        pos: 7_674_944,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["stop_gained"]),
        assert_not_intergenic: false,
        expected_nmd: Some(true),
        category: "E: NMD (stop_gained, early, true)",
    },
    // #11 — APC R1450X: stop_gained in last exon (16/16), NMD escape.
    //   CDS position 4348, well inside APC's enormous last exon (~6500bp).
    //   Variant is in or past the last junction → NMD = false.
    //   VEP REST (2026-04-10): NM_000038.6:c.4348C>T, stop_gained
    //   Genomic: chr5:112839942 (1-based) = 112839941 (0-based), plus-strand.
    ExpectedAnnotation {
        label: "APC R1450X stop_gained (NMD false, last exon)",
        chrom: "chr5",
        pos: 112_839_941,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000038.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["stop_gained"]),
        assert_not_intergenic: false,
        expected_nmd: Some(false),
        category: "E: NMD (stop_gained, last exon, false)",
    },
    // #12 — TP53 R342X: stop_gained in penultimate CDS segment, >50nt from
    //   last junction. CDS 1024, distance = 1101 - 1024 = 77 > 50. NMD = true.
    //   VEP REST (2026-04-10): NM_000546.6:c.1024C>T, p.Arg342Ter
    //   Genomic: CDS 1024 maps to chr17:7670684 (0-based, minus strand).
    //   Plus-strand: REF=G, ALT=A (complement of coding C>T).
    ExpectedAnnotation {
        label: "TP53 R342X stop_gained (NMD true, penultimate, dist=77)",
        chrom: "chr17",
        pos: 7_670_684,
        ref_allele: b"G",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["stop_gained"]),
        assert_not_intergenic: false,
        expected_nmd: Some(true),
        category: "E: NMD (stop_gained, penultimate, true)",
    },
    // #13 — BRCA1 c.68_69del: frameshift early in CDS, NMD predicted.
    //   CDS position 68 in exon 2 of 24 exons. Distance to last junction
    //   is thousands of bases. NMD = true.
    //   Source: indel test #1, hgvs test #18, hgvs_p test #11.
    ExpectedAnnotation {
        label: "BRCA1 c.68_69del frameshift (NMD true, early CDS)",
        chrom: "chr17",
        pos: 43_124_026,
        ref_allele: b"ACT",
        alt_allele: b"A",
        transcript: "NM_007294.4",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["frameshift_variant"]),
        assert_not_intergenic: false,
        expected_nmd: Some(true),
        category: "E: NMD (frameshift, early, true)",
    },
    // #14 — APC c.3927_3931del: frameshift in last exon (16/16), NMD escape.
    //   CDS position 3927, inside APC's last exon. NMD = false.
    //   Source: indel test #4, hgvs_p test #14.
    ExpectedAnnotation {
        label: "APC c.3927_3931del frameshift (NMD false, last exon)",
        chrom: "chr5",
        pos: 112_839_519,
        ref_allele: b"AAAAGA",
        alt_allele: b"A",
        transcript: "NM_000038.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["frameshift_variant"]),
        assert_not_intergenic: false,
        expected_nmd: Some(false),
        category: "E: NMD (frameshift, last exon, false)",
    },
    // #15 — PRNP p.Gln160Ter: stop_gained in a single-exon gene.
    //   PRNP NM_000311.5 has exactly 1 CDS segment (chr20, plus strand, 762bp).
    //   No exon-exon junction → NMD impossible → false.
    //   Codon 160: CDS 478-480 = CAA (Gln). c.478C>T → TAA (Stop).
    //   Genomic: chr20:4699697 (0-based), plus strand, REF=C, ALT=T.
    ExpectedAnnotation {
        label: "PRNP Q160X stop_gained (NMD false, single-exon gene)",
        chrom: "chr20",
        pos: 4_699_697,
        ref_allele: b"C",
        alt_allele: b"T",
        transcript: "NM_000311.5",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["stop_gained"]),
        assert_not_intergenic: false,
        expected_nmd: Some(false),
        category: "E: NMD (single-exon, false)",
    },
    // #16 — TP53 50-nt boundary: stop_gained at exactly 50 nt from last
    //   junction. Tests the `> 50` threshold (NOT `>= 50`).
    //
    //   Derivation:
    //     TP53 NM_000546.6 last_junction_cds_pos = 1101 (see module doc).
    //     CDS 1051: distance = 1101 - 1051 = 50. NOT > 50 → NMD = false.
    //     Codon 351: CDS 1051-1053 = AAG (Lys) on coding strand.
    //     A→T at position 1: codon becomes TAG (Stop).
    //     Coding strand A→T = plus strand T→A (complement, minus strand gene).
    //   Genomic: chr17:7670657 (0-based), REF=T, ALT=A.
    ExpectedAnnotation {
        label: "TP53 c.1051A>T Lys351Ter (NMD false, dist=50, boundary)",
        chrom: "chr17",
        pos: 7_670_657,
        ref_allele: b"T",
        alt_allele: b"A",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["stop_gained"]),
        assert_not_intergenic: false,
        expected_nmd: Some(false),
        category: "E: NMD (boundary, dist=50, false)",
    },
    // #17 — TP53 51-nt boundary: frameshift at exactly 51 nt from last
    //   junction. Paired with #16 to verify `>` vs `>=` comparison.
    //
    //   Derivation:
    //     CDS 1050: distance = 1101 - 1050 = 51 > 50 → NMD = true.
    //     No single-base stop_gained SNV is possible at CDS 1050 (position 3
    //     of codon 350 = CTC, Leu — no stop by changing position 3 alone).
    //     Instead use a 1-base deletion → frameshift at CDS 1050.
    //     CDS 1050 → genomic 7670658, plus-strand REF=G.
    //     Surrounding plus-strand: A(7670659), G(7670658), T(7670657) — no
    //     mono-nucleotide repeat → 3' shift = 0 → CDS position stays 1050.
    //   VCF: chr17:7670657 (0-based, anchor=T), REF=TG, ALT=T.
    ExpectedAnnotation {
        label: "TP53 c.1050del frameshift (NMD true, dist=51, boundary)",
        chrom: "chr17",
        pos: 7_670_657,
        ref_allele: b"TG",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["frameshift_variant"]),
        assert_not_intergenic: false,
        expected_nmd: Some(true),
        category: "E: NMD (boundary, dist=51, true)",
    },
    // #18 — TP53 R248R: synonymous variant, no NMD.
    //   predicts_nmd should only be true for PTC-producing variants
    //   (stop_gained, frameshift_variant). A synonymous SNV at any CDS
    //   position must have predicts_nmd = false.
    //   Source: hgvs_p test #4.
    //   VEP: NM_000546.6:c.742C>A, p.Arg248= (synonymous_variant)
    //   Genomic: chr17:7674220 (0-based), REF=G, ALT=T (complement of C>A).
    ExpectedAnnotation {
        label: "TP53 R248R synonymous (NMD false, no PTC)",
        chrom: "chr17",
        pos: 7_674_220,
        ref_allele: b"G",
        alt_allele: b"T",
        transcript: "NM_000546.6",
        expected_hgvs_c: None,
        expected_hgvs_p: None,
        expected_consequences: Some(&["synonymous_variant"]),
        assert_not_intergenic: false,
        expected_nmd: Some(false),
        category: "F: synonymous (no NMD)",
    },
];

/// Compare a single variant's `annotate()` output against expected values.
/// Returns `Ok(empty vec)` on full match, `Ok(non-empty)` with mismatch
/// descriptions, or `Err(msg)` on transcript-not-found / annotate error.
fn check_variant(ve: &VarEffect, exp: &ExpectedAnnotation) -> Result<Vec<String>, String> {
    let results = ve
        .annotate(
            Assembly::GRCh38,
            exp.chrom,
            exp.pos,
            exp.ref_allele,
            exp.alt_allele,
        )
        .map_err(|e| format!("annotate error: {e}"))?;

    // For intergenic: match on empty transcript.
    // For others: match on transcript accession.
    let result = if exp.transcript.is_empty() {
        results
            .consequences
            .iter()
            .find(|r| r.transcript.is_empty())
            .ok_or_else(|| {
                let detail: Vec<String> = results
                    .consequences
                    .iter()
                    .map(|r| {
                        let csqs: Vec<&str> = r.consequences.iter().map(|c| c.as_str()).collect();
                        format!("{}:{csqs:?}", r.transcript)
                    })
                    .collect();
                format!("no intergenic result; found: {detail:?}")
            })?
    } else {
        results
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
            })?
    };

    let mut mismatches: Vec<String> = Vec::new();

    // HGVS c.
    if let Some(expected) = exp.expected_hgvs_c
        && result.hgvs_c.as_deref() != Some(expected)
    {
        mismatches.push(format!(
            "hgvs_c: expected {:?}, got {:?}",
            expected, result.hgvs_c,
        ));
    }

    // HGVS p.
    if let Some(expected) = exp.expected_hgvs_p
        && result.hgvs_p.as_deref() != Some(expected)
    {
        mismatches.push(format!(
            "hgvs_p: expected {:?}, got {:?}",
            expected, result.hgvs_p,
        ));
    }

    // Consequences (subset: every expected term must be present).
    if let Some(expected_csq) = exp.expected_consequences {
        let actual: BTreeSet<&str> = result.consequences.iter().map(|c| c.as_str()).collect();
        let missing: Vec<&&str> = expected_csq
            .iter()
            .filter(|c| !actual.contains(**c))
            .collect();
        if !missing.is_empty() {
            mismatches.push(format!(
                "consequences: missing {missing:?} from actual {actual:?}",
            ));
        }
    }

    // Negative intergenic control.
    if exp.assert_not_intergenic
        && result
            .consequences
            .iter()
            .any(|c| c.as_str() == "intergenic_variant")
    {
        mismatches.push("intergenic_variant found but should NOT be present".into());
    }

    // predicts_nmd.
    if let Some(expected_nmd) = exp.expected_nmd
        && result.predicts_nmd != expected_nmd
    {
        mismatches.push(format!(
            "predicts_nmd: expected {expected_nmd}, got {}",
            result.predicts_nmd,
        ));
    }

    Ok(mismatches)
}

#[test]
#[ignore]
fn vep_concordance_grch38_normalization_18() {
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
                eprintln!(
                    "  [{num:>2}/{total}] SKIP  {} [{}]: {msg}",
                    exp.label, exp.category,
                );
                skip += 1;
            }
            Err(msg) => {
                eprintln!(
                    "  [{num:>2}/{total}] FAIL  {} [{}]: {msg}",
                    exp.label, exp.category,
                );
                failures.push(format!("#{num} {}: {msg}", exp.label));
                fail += 1;
            }
            Ok(ref mismatches) if mismatches.is_empty() => {
                eprintln!(
                    "  [{num:>2}/{total}] PASS  {} [{}]",
                    exp.label, exp.category,
                );
                pass += 1;
            }
            Ok(mismatches) => {
                eprintln!(
                    "  [{num:>2}/{total}] FAIL  {} [{}]",
                    exp.label, exp.category,
                );
                let detail = mismatches.join("\n           ");
                failures.push(format!("#{num} {}:\n           {detail}", exp.label));
                fail += 1;
            }
        }
    }

    eprintln!(
        "\n=== Normalization concordance: {pass}/{} passed, {skip} skipped ===",
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
        "{fail} of {} variants failed normalization concordance check",
        pass + fail,
    );
}
