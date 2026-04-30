//! Stage C — GRCh37 ClinVar self-concordance harness.
//!
//! This test reads the dual-coordinate ClinVar pairing TSV produced by
//! `vareffect/scripts/generate_clinvar_pairs.py` and verifies that
//! `vareffect.annotate(GRCh37, …)` and `vareffect.annotate(GRCh38, …)`
//! agree on `hgvs_c`, `hgvs_p`, and the consequence-term set for every
//! transcript accession that exists in both assemblies' stores. It is
//! `#[ignore]`-gated because it requires the GRCh38 + GRCh37 binaries
//! and the committed pairing TSV.
//!
//! Unlike `vep_large_concordance.rs` (which compares vareffect against
//! VEP's REST output on GRCh38 alone), this harness is *self-concordance*:
//! both halves are produced by vareffect, so divergence here points
//! squarely at GRCh37-specific logic — chrom mapping, RefSeq Select
//! handling, divergent transcripts, the chrM codon table.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//! GRCH37_FASTA=data/vareffect/GRCh37.bin \
//! GRCH37_TRANSCRIPTS=data/vareffect/transcript_models_grch37.bin \
//!   cargo test -p vareffect --release -- --ignored grch37_clinvar_concordance --nocapture
//! ```
//!
//! # Output files
//!
//! * Report — printed to stderr (`--nocapture` to see it live)
//! * Mismatch dump — `vareffect/tests/data/grch37_clinvar_mismatches.log`
//!   (gitignored, overwritten every run)
//!
//! # Exclusion logic
//!
//! Per `VEP_DIVERGENCES.md`, ~5 % of NCBI RefSeq Select GRCh37 transcripts
//! carry `genome_transcript_divergent: true` (transcript sequence diverges
//! from the GRCh37 reference at one or more positions). HGVS positions
//! against a divergent transcript may not round-trip to the same genomic
//! position they would against the reference — so for those accessions,
//! disagreement between vareffect-37 and vareffect-38 is *expected by
//! construction* and counting them in the concordance metric would mask
//! real bugs.
//!
//! At setup we build a `HashSet` of divergent accessions from the GRCh37
//! store and assert its size lands between 4 % and 6 % of the store
//! (NCBI's published ~5 % figure ± 1 percentage point). If the GFF3
//! `Note=` parser regex stops matching, the set silently shrinks toward
//! zero; this regression target catches that.
//!
//! chrM transcripts get *no* special exclusion — both ClinVar VCFs use
//! `NC_012920.1` (the rCRS), and the chrom module already routes them
//! identically across builds. The ~37 chrM transcripts that Stage B's
//! UCSC cross-validation couldn't compare (UCSC `hg19` chrM is the 1981
//! Anderson reference, NCBI is the rCRS) are validated here.

use std::any::Any;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;
use vareffect::{AnnotationResult, Assembly, ConsequenceResult, VarEffect};

// ---------------------------------------------------------------------------
// Paths and constants
// ---------------------------------------------------------------------------

/// Hard lower bound on each per-metric concordance rate. Below this the
/// test fails. Set to 0.99 because both sides are vareffect — divergence
/// here is a real bug, not a VEP/vareffect interpretive disagreement.
const METRIC_THRESHOLD: f64 = 0.99;

/// Acceptable range for the divergent-accessions fraction of the GRCh37
/// store. NCBI publishes ~5 % for RefSeq Select on GRCh37; we tolerate
/// ±1 percentage point of drift across weekly RefSeq refreshes.
const DIVERGENT_FRACTION_LO: f64 = 0.04;
const DIVERGENT_FRACTION_HI: f64 = 0.06;

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("could not determine workspace root from CARGO_MANIFEST_DIR")
        .to_path_buf()
}

fn pairs_tsv_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("clinvar_grch37_grch38_pairs.tsv")
}

fn mismatches_log_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("grch37_clinvar_mismatches.log")
}

// ---------------------------------------------------------------------------
// VarEffect bring-up
// ---------------------------------------------------------------------------

fn open_var_effect() -> VarEffect {
    let g37_fasta = std::env::var("GRCH37_FASTA").expect(
        "GRCH37_FASTA env var must point at the flat-binary genome \
         (e.g. data/vareffect/GRCh37.bin). Run `vareffect setup --assembly grch37` first.",
    );
    let g37_transcripts = std::env::var("GRCH37_TRANSCRIPTS")
        .expect("GRCH37_TRANSCRIPTS env var must point at the transcript_models_grch37.bin");
    let g38_fasta =
        std::env::var("GRCH38_FASTA").expect("GRCH38_FASTA env var must point at GRCh38.bin");

    let root = workspace_root();
    let g38_transcripts = root.join("data/vareffect/transcript_models_grch38.bin");
    let g38_aliases = root.join("data/vareffect/patch_chrom_aliases_grch38.csv");

    VarEffect::builder()
        .with_grch37(&PathBuf::from(g37_transcripts), &PathBuf::from(&g37_fasta))
        .expect("loading GRCh37")
        .with_grch38_and_patch_aliases(&g38_transcripts, &PathBuf::from(&g38_fasta), &g38_aliases)
        .expect("loading GRCh38")
        .build()
        .expect("building VarEffect")
}

/// Build the divergent-accession set from the GRCh37 store already loaded
/// inside the [`VarEffect`] handle. Asserts the fraction is within the
/// regression band (catches a silent parser regression where the GFF3
/// `Note=` patterns stop matching).
fn divergent_accessions_with_assertion(ve: &VarEffect) -> HashSet<String> {
    let store = ve
        .transcripts(Assembly::GRCh37)
        .expect("GRCh37 store must be loaded into VarEffect by open_var_effect");
    let total = store.transcripts().len();
    let set: HashSet<String> = store
        .transcripts()
        .iter()
        .filter(|t| t.genome_transcript_divergent)
        .map(|t| t.accession.clone())
        .collect();
    let fraction = set.len() as f64 / total.max(1) as f64;
    eprintln!(
        "  divergent transcripts: {} / {} ({:.2}%)",
        set.len(),
        total,
        fraction * 100.0,
    );
    assert!(
        (DIVERGENT_FRACTION_LO..=DIVERGENT_FRACTION_HI).contains(&fraction),
        "divergent-accession fraction {:.2}% is outside [{:.0}%, {:.0}%] -- \
         the GFF3 `Note=` parser may have stopped matching, or NCBI changed \
         the wording (see vareffect-cli/src/builders/transcript_models/gff3_attrs.rs).",
        fraction * 100.0,
        DIVERGENT_FRACTION_LO * 100.0,
        DIVERGENT_FRACTION_HI * 100.0,
    );
    set
}

// ---------------------------------------------------------------------------
// TSV parsing
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct Pair {
    variation_id: String,
    chrom_38: String,
    /// 1-based VCF position.
    pos_38: u64,
    chrom_37: String,
    /// 1-based VCF position.
    pos_37: u64,
    ref_allele: String,
    alt_allele: String,
}

fn parse_tsv(path: &Path) -> Vec<Pair> {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("opening pairing TSV at {}: {e}", path.display()));
    let reader = BufReader::new(file);
    let mut out = Vec::new();

    for (lineno, line_result) in reader.lines().enumerate() {
        let line = line_result.unwrap_or_else(|e| panic!("reading TSV line {lineno}: {e}"));
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("variation_id\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            panic!(
                "TSV line {lineno} has {} fields, expected 8: {line}",
                fields.len()
            );
        }
        let pos_38: u64 = fields[2]
            .parse()
            .unwrap_or_else(|_| panic!("TSV line {lineno}: invalid pos_38 {:?}", fields[2]));
        let pos_37: u64 = fields[4]
            .parse()
            .unwrap_or_else(|_| panic!("TSV line {lineno}: invalid pos_37 {:?}", fields[4]));
        // VCF positions are 1-based; pos=0 is invalid in the VCF spec.
        // Catch malformed input here rather than masking it with a
        // saturating subtraction at the call site.
        assert!(
            pos_38 >= 1,
            "TSV line {lineno}: pos_38 must be 1-based (>=1), got {pos_38}"
        );
        assert!(
            pos_37 >= 1,
            "TSV line {lineno}: pos_37 must be 1-based (>=1), got {pos_37}"
        );
        out.push(Pair {
            variation_id: fields[0].to_owned(),
            chrom_38: fields[1].to_owned(),
            pos_38,
            chrom_37: fields[3].to_owned(),
            pos_37,
            ref_allele: fields[5].to_owned(),
            alt_allele: fields[6].to_owned(),
        });
    }
    out
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

/// Compare a single (accession, GRCh37 row, GRCh38 row) triplet on the
/// three concordance dimensions. Returns `None` when everything matches
/// — the per-row log only records divergences.
struct TranscriptCompare {
    accession: String,
    consequences_match: bool,
    hgvs_c_match: bool,
    hgvs_p_match: bool,
    /// Compact mismatch detail; empty when all three matched.
    detail: Option<String>,
}

/// Strip `accession:` prefix from an HGVS string. vareffect's `hgvs_c`
/// always carries the accession; `hgvs_p` does not. Stripping makes the
/// comparison insensitive to GRCh37 vs GRCh38 transcript-version drift
/// when the bare c./p. suffix is identical.
fn strip_hgvs_accession(full: &str) -> &str {
    full.split_once(':').map(|(_, rest)| rest).unwrap_or(full)
}

fn compare_transcript(
    accession: &str,
    r37: &ConsequenceResult,
    r38: &ConsequenceResult,
) -> TranscriptCompare {
    let csq37: BTreeSet<String> = r37
        .consequences
        .iter()
        .map(|c| c.as_str().to_owned())
        .collect();
    let csq38: BTreeSet<String> = r38
        .consequences
        .iter()
        .map(|c| c.as_str().to_owned())
        .collect();
    let consequences_match = csq37 == csq38;

    let c37 = r37.hgvs_c.as_deref().map(strip_hgvs_accession);
    let c38 = r38.hgvs_c.as_deref().map(strip_hgvs_accession);
    let hgvs_c_match = c37 == c38;

    let p37 = r37.hgvs_p.as_deref().map(strip_hgvs_accession);
    let p38 = r38.hgvs_p.as_deref().map(strip_hgvs_accession);
    let hgvs_p_match = p37 == p38;

    let detail = if consequences_match && hgvs_c_match && hgvs_p_match {
        None
    } else {
        let mut parts = Vec::new();
        if !consequences_match {
            parts.push(format!("csq: g37={csq37:?} g38={csq38:?}"));
        }
        if !hgvs_c_match {
            parts.push(format!("hgvs_c: g37={c37:?} g38={c38:?}"));
        }
        if !hgvs_p_match {
            parts.push(format!("hgvs_p: g37={p37:?} g38={p38:?}"));
        }
        Some(parts.join(" | "))
    };

    TranscriptCompare {
        accession: accession.to_owned(),
        consequences_match,
        hgvs_c_match,
        hgvs_p_match,
        detail,
    }
}

// ---------------------------------------------------------------------------
// Per-row outcome and reporting structs
// ---------------------------------------------------------------------------

enum RowResult {
    /// `VarEffect::annotate` returned `Err` on either side.
    AnnotationError {
        lineno: usize,
        side: &'static str,
        message: String,
    },
    /// Annotation panicked on either side (hard bug).
    Panic {
        lineno: usize,
        side: &'static str,
        message: String,
    },
    /// Annotation succeeded; carries per-shared-transcript comparisons
    /// and the count of divergent accessions skipped.
    Completed {
        lineno: usize,
        compares: Vec<TranscriptCompare>,
        divergent_skipped: usize,
        no_overlap_grch37: bool,
        no_overlap_grch38: bool,
    },
}

#[derive(Default)]
struct Stats {
    total_tsv_rows: usize,
    annotation_errors: usize,
    panics: usize,
    /// Variants where one side or the other returned no transcripts at
    /// all (likely intergenic or chrom alias miss).
    no_overlap_either_side: usize,
    /// Variants where neither side shared a transcript accession with
    /// the other. Tracked separately from the per-transcript stats.
    no_shared_transcripts: usize,
    /// Variants where every shared accession was on the divergent list.
    only_divergent_shared: usize,
    /// Sum of (shared accessions skipped due to divergent flag) across
    /// all rows. Informational.
    divergent_skipped: usize,
    transcripts_compared: usize,
    consequences_match: usize,
    hgvs_c_match: usize,
    hgvs_p_match: usize,
    /// Top mismatch dimensions for the summary table.
    mismatch_patterns: BTreeMap<String, usize>,
}

fn pct(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        (num as f64 / den as f64) * 100.0
    }
}

fn mismatch_pattern_key(c: &TranscriptCompare) -> String {
    let mut dims = Vec::new();
    if !c.consequences_match {
        dims.push("csq");
    }
    if !c.hgvs_c_match {
        dims.push("hgvs_c");
    }
    if !c.hgvs_p_match {
        dims.push("hgvs_p");
    }
    dims.join("+")
}

fn panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_owned()
    } else {
        "<non-string panic payload>".to_owned()
    }
}

// ---------------------------------------------------------------------------
// Annotation (parallel)
// ---------------------------------------------------------------------------

fn annotate_one_side(
    ve: &VarEffect,
    assembly: Assembly,
    chrom: &str,
    pos_one_based: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Result<AnnotationResult, vareffect::VarEffectError> {
    // VCF is 1-based; vareffect annotate is 0-based. The TSV preserves
    // VCF convention so we subtract here at the call boundary.
    // `parse_tsv` already asserts `pos >= 1` so the subtraction can't
    // underflow.
    let pos = pos_one_based - 1;
    ve.annotate(assembly, chrom, pos, ref_allele, alt_allele)
}

fn process_pair(
    lineno: usize,
    pair: &Pair,
    ve: &VarEffect,
    divergent: &HashSet<String>,
) -> RowResult {
    let ref_b = pair.ref_allele.as_bytes();
    let alt_b = pair.alt_allele.as_bytes();

    let g37 = panic::catch_unwind(AssertUnwindSafe(|| {
        annotate_one_side(
            ve,
            Assembly::GRCh37,
            &pair.chrom_37,
            pair.pos_37,
            ref_b,
            alt_b,
        )
    }));
    let g37 = match g37 {
        Ok(Ok(r)) => r,
        Ok(Err(e)) => {
            return RowResult::AnnotationError {
                lineno,
                side: "GRCh37",
                message: e.to_string(),
            };
        }
        Err(payload) => {
            return RowResult::Panic {
                lineno,
                side: "GRCh37",
                message: panic_message(payload.as_ref()),
            };
        }
    };

    let g38 = panic::catch_unwind(AssertUnwindSafe(|| {
        annotate_one_side(
            ve,
            Assembly::GRCh38,
            &pair.chrom_38,
            pair.pos_38,
            ref_b,
            alt_b,
        )
    }));
    let g38 = match g38 {
        Ok(Ok(r)) => r,
        Ok(Err(e)) => {
            return RowResult::AnnotationError {
                lineno,
                side: "GRCh38",
                message: e.to_string(),
            };
        }
        Err(payload) => {
            return RowResult::Panic {
                lineno,
                side: "GRCh38",
                message: panic_message(payload.as_ref()),
            };
        }
    };

    let no_overlap_grch37 = g37.consequences.is_empty();
    let no_overlap_grch38 = g38.consequences.is_empty();

    let g37_by_acc: HashMap<&str, &ConsequenceResult> = g37
        .consequences
        .iter()
        .map(|c| (c.transcript.as_str(), c))
        .collect();
    let g38_by_acc: HashMap<&str, &ConsequenceResult> = g38
        .consequences
        .iter()
        .map(|c| (c.transcript.as_str(), c))
        .collect();

    let mut compares = Vec::new();
    let mut divergent_skipped = 0usize;
    for (acc, r37) in &g37_by_acc {
        let Some(r38) = g38_by_acc.get(acc) else {
            continue;
        };
        if divergent.contains(*acc) {
            divergent_skipped += 1;
            continue;
        }
        compares.push(compare_transcript(acc, r37, r38));
    }

    RowResult::Completed {
        lineno,
        compares,
        divergent_skipped,
        no_overlap_grch37,
        no_overlap_grch38,
    }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

fn print_report(stats: &Stats, elapsed_secs: f64, num_threads: usize) {
    let tested = stats.transcripts_compared;
    eprintln!();
    eprintln!("============================================================");
    eprintln!("GRCh37 ClinVar SELF-CONCORDANCE REPORT");
    eprintln!("============================================================");
    eprintln!(
        "Total TSV rows:                   {:>6}",
        stats.total_tsv_rows
    );
    eprintln!(
        "Skipped (annotation error):       {:>6}",
        stats.annotation_errors
    );
    eprintln!(
        "Skipped (PANIC):                  {:>6}  <- hard bugs, see mismatches.log",
        stats.panics,
    );
    eprintln!(
        "Skipped (no overlap on either):   {:>6}",
        stats.no_overlap_either_side,
    );
    eprintln!(
        "Skipped (no shared transcripts):  {:>6}",
        stats.no_shared_transcripts,
    );
    eprintln!(
        "Skipped (only divergent shared):  {:>6}",
        stats.only_divergent_shared,
    );
    eprintln!(
        "Divergent transcripts skipped:    {:>6}  (informational, not counted)",
        stats.divergent_skipped,
    );
    eprintln!("Transcripts compared:             {tested:>6}");
    eprintln!();
    eprintln!(
        "Consequence concordance: {:>6}/{tested} ({:.2}%)",
        stats.consequences_match,
        pct(stats.consequences_match, tested),
    );
    eprintln!(
        "HGVS c. concordance:     {:>6}/{tested} ({:.2}%)",
        stats.hgvs_c_match,
        pct(stats.hgvs_c_match, tested),
    );
    eprintln!(
        "HGVS p. concordance:     {:>6}/{tested} ({:.2}%)",
        stats.hgvs_p_match,
        pct(stats.hgvs_p_match, tested),
    );
    eprintln!();

    let mut patterns: Vec<(&String, &usize)> = stats.mismatch_patterns.iter().collect();
    patterns.sort_by(|a, b| b.1.cmp(a.1));
    if !patterns.is_empty() {
        eprintln!("TOP MISMATCH DIMENSIONS:");
        for (pattern, count) in patterns.iter().take(20) {
            eprintln!("  {count:>5}x  {pattern}");
        }
        eprintln!();
    }

    eprintln!("Threads:    {num_threads}");
    if tested > 0 {
        eprintln!(
            "Throughput: {:.0} comparisons/second",
            tested as f64 / elapsed_secs
        );
    }
    eprintln!("Elapsed:    {elapsed_secs:.1}s");
    eprintln!("============================================================");
}

fn write_mismatches_log(
    path: &Path,
    mismatches: &[(usize, &Pair, &TranscriptCompare)],
    annotation_errors: &[(usize, &'static str, &str, &Pair)],
    panics: &[(usize, &'static str, &str, &Pair)],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = File::create(path)?;
    writeln!(
        f,
        "# GRCh37 ClinVar self-concordance — mismatch log ({} mismatches, {} errors, {} panics)",
        mismatches.len(),
        annotation_errors.len(),
        panics.len(),
    )?;

    if !panics.is_empty() {
        writeln!(f, "\n# === PANICS ===")?;
        writeln!(
            f,
            "# lineno\tside\tvariation_id\tchrom:pos(g37/g38)\tref>alt\tmessage"
        )?;
        for (lineno, side, msg, p) in panics {
            writeln!(
                f,
                "{lineno}\t{side}\t{}\t{}:{}/{}:{}\t{}>{}\t{msg}",
                p.variation_id,
                p.chrom_37,
                p.pos_37,
                p.chrom_38,
                p.pos_38,
                p.ref_allele,
                p.alt_allele,
            )?;
        }
    }

    if !annotation_errors.is_empty() {
        writeln!(f, "\n# === ANNOTATION ERRORS ===")?;
        writeln!(
            f,
            "# lineno\tside\tvariation_id\tchrom:pos(g37/g38)\tref>alt\tmessage"
        )?;
        for (lineno, side, msg, p) in annotation_errors {
            writeln!(
                f,
                "{lineno}\t{side}\t{}\t{}:{}/{}:{}\t{}>{}\t{msg}",
                p.variation_id,
                p.chrom_37,
                p.pos_37,
                p.chrom_38,
                p.pos_38,
                p.ref_allele,
                p.alt_allele,
            )?;
        }
    }

    writeln!(f, "\n# === MISMATCHES ===")?;
    writeln!(f, "# lineno\tvariation_id\ttranscript\tdetail")?;
    for (lineno, p, c) in mismatches {
        writeln!(
            f,
            "{lineno}\t{}\t{}\t{}",
            p.variation_id,
            c.accession,
            c.detail.as_deref().unwrap_or(""),
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Test entry point
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn grch37_clinvar_concordance() {
    let tsv = pairs_tsv_path();
    if !tsv.exists() {
        eprintln!(
            "SKIP: pairing TSV not found at {} -- run \
             `vareffect/scripts/generate_clinvar_pairs.py` first.",
            tsv.display()
        );
        return;
    }

    eprintln!("Loading both assemblies (this can take a few seconds)...");
    let load_start = Instant::now();
    let ve = open_var_effect();
    let divergent = divergent_accessions_with_assertion(&ve);
    eprintln!("  loaded in {:.1}s", load_start.elapsed().as_secs_f64());

    let pairs = parse_tsv(&tsv);
    eprintln!("Parsed {} pairs from {}", pairs.len(), tsv.display());

    let num_threads = std::env::var("CONCORDANCE_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|n| *n > 0)
        .unwrap_or_else(rayon::current_num_threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("rayon pool");

    let start = Instant::now();
    let results: Vec<RowResult> = pool.install(|| {
        pairs
            .par_iter()
            .enumerate()
            .map(|(i, pair)| process_pair(i + 1, pair, &ve, &divergent))
            .collect()
    });
    let elapsed = start.elapsed().as_secs_f64();

    // Sequential reduction.
    let mut stats = Stats {
        total_tsv_rows: pairs.len(),
        ..Stats::default()
    };
    let mut mismatches: Vec<(usize, &Pair, TranscriptCompare)> = Vec::new();
    let mut annotation_errors: Vec<(usize, &'static str, String, &Pair)> = Vec::new();
    let mut panics: Vec<(usize, &'static str, String, &Pair)> = Vec::new();

    for result in results {
        match result {
            RowResult::AnnotationError {
                lineno,
                side,
                message,
            } => {
                stats.annotation_errors += 1;
                annotation_errors.push((lineno, side, message, &pairs[lineno - 1]));
            }
            RowResult::Panic {
                lineno,
                side,
                message,
            } => {
                stats.panics += 1;
                panics.push((lineno, side, message, &pairs[lineno - 1]));
            }
            RowResult::Completed {
                lineno,
                compares,
                divergent_skipped,
                no_overlap_grch37,
                no_overlap_grch38,
            } => {
                stats.divergent_skipped += divergent_skipped;
                if no_overlap_grch37 || no_overlap_grch38 {
                    stats.no_overlap_either_side += 1;
                    continue;
                }
                if compares.is_empty() {
                    if divergent_skipped > 0 {
                        stats.only_divergent_shared += 1;
                    } else {
                        stats.no_shared_transcripts += 1;
                    }
                    continue;
                }
                for c in compares {
                    stats.transcripts_compared += 1;
                    if c.consequences_match {
                        stats.consequences_match += 1;
                    }
                    if c.hgvs_c_match {
                        stats.hgvs_c_match += 1;
                    }
                    if c.hgvs_p_match {
                        stats.hgvs_p_match += 1;
                    }
                    if c.detail.is_some() {
                        let key = mismatch_pattern_key(&c);
                        *stats.mismatch_patterns.entry(key).or_default() += 1;
                        mismatches.push((lineno, &pairs[lineno - 1], c));
                    }
                }
            }
        }
    }

    print_report(&stats, elapsed, num_threads);

    let mismatches_view: Vec<(usize, &Pair, &TranscriptCompare)> =
        mismatches.iter().map(|(l, p, c)| (*l, *p, c)).collect();
    let annotation_errors_view: Vec<(usize, &'static str, &str, &Pair)> = annotation_errors
        .iter()
        .map(|(l, s, m, p)| (*l, *s, m.as_str(), *p))
        .collect();
    let panics_view: Vec<(usize, &'static str, &str, &Pair)> = panics
        .iter()
        .map(|(l, s, m, p)| (*l, *s, m.as_str(), *p))
        .collect();
    let log_path = mismatches_log_path();
    if let Err(e) = write_mismatches_log(
        &log_path,
        &mismatches_view,
        &annotation_errors_view,
        &panics_view,
    ) {
        eprintln!(
            "WARN: failed to write mismatch log to {}: {e}",
            log_path.display()
        );
    } else {
        eprintln!("Mismatch log: {}", log_path.display());
    }

    // Hard panics are always fatal.
    assert_eq!(
        stats.panics, 0,
        "annotation panicked on {} rows",
        stats.panics
    );

    // Per-metric thresholds.
    let csq_pct = stats.consequences_match as f64 / stats.transcripts_compared.max(1) as f64;
    let c_pct = stats.hgvs_c_match as f64 / stats.transcripts_compared.max(1) as f64;
    let p_pct = stats.hgvs_p_match as f64 / stats.transcripts_compared.max(1) as f64;
    assert!(
        csq_pct >= METRIC_THRESHOLD,
        "consequence concordance {:.4} below threshold {METRIC_THRESHOLD}",
        csq_pct,
    );
    assert!(
        c_pct >= METRIC_THRESHOLD,
        "hgvs_c concordance {:.4} below threshold {METRIC_THRESHOLD}",
        c_pct,
    );
    assert!(
        p_pct >= METRIC_THRESHOLD,
        "hgvs_p concordance {:.4} below threshold {METRIC_THRESHOLD}",
        p_pct,
    );
}
