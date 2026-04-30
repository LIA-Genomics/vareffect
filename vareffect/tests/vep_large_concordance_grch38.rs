//! Large-scale VEP concordance harness.
//!
//! This test reads a pre-computed TSV of VEP ground truth (see
//! `vareffect/scripts/generate_vep_ground_truth.py`), runs vareffect on
//! every row, compares per-field results, classifies divergences, and prints a
//! concordance report. It is `#[ignore]`-gated because it requires the
//! transcript store, reference FASTA, and the committed ground-truth TSV.
//!
//! Unlike the hand-picked spot-check tests (`vep_concordance_*.rs`), this
//! harness is statistical: it asserts a ≥95% consequence-concordance *rate*
//! across ~50k real ClinVar variants rather than per-variant equality. The
//! threshold is deliberately softer than the ≥99% long-term goal so the
//! harness can land before every divergence pattern has been triaged; fixes
//! are scheduled as follow-up work and the threshold is tightened once the
//! baseline is characterised.
//!
//! Run with:
//! ```bash
//! GRCH38_FASTA=data/vareffect/GRCh38.bin \
//!   cargo test -p vareffect -- --ignored vep_large_concordance --nocapture
//! ```
//!
//! # Output files
//!
//! * Report — printed to stderr (`--nocapture` to see it live)
//! * Mismatch dump — `crates/vareffect/tests/data/vep_large_concordance_mismatches.log`
//!   (gitignored, overwritten every run)
//!
//! # Divergence taxonomy
//!
//! Before comparison, the test *normalises* VEP's consequence set into
//! vareffect's vocabulary:
//!
//! * **Granular splice terms** — VEP emits `splice_donor_5th_base_variant`,
//!   `splice_donor_region_variant`, `splice_polypyrimidine_tract_variant`;
//!   vareffect collapses these into `splice_region_variant`. Each such term
//!   is folded into `splice_region_variant` on the VEP side.
//! * **`NMD_transcript_variant`** — VEP emits this SO term for NMD-biotype
//!   transcripts; vareffect tracks the 50-nt rule via the
//!   [`ConsequenceResult::predicts_nmd`] bool instead of a consequence term.
//!   It is filtered out of VEP's set.
//!
//! After normalisation, a variant passes if
//! `normalized_vep_consequences.is_subset(&vareffect_consequences)`. Any
//! strict subset mismatch is logged as a real divergence.
//!
//! Separately, the test classifies two *non*-mismatch buckets that don't
//! count against the concordance rate:
//!
//! * **Transcript not found** — VEP's chosen RefSeq transcript is absent from
//!   vareffect's store.
//! * **Transcript version mismatch** — same accession, different version
//!   (e.g. VEP has `NM_000546.6`, vareffect has `NM_000546.5`).

use std::any::Any;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;
use vareffect::{Assembly, ConsequenceResult, FastaReader, TranscriptStore, VarEffect};

// ---------------------------------------------------------------------------
// Paths and constants
// ---------------------------------------------------------------------------

/// Hard lower bound on consequence concordance. Below this the test fails.
/// Tighten to 0.99 once divergences are triaged.
const CONSEQUENCE_THRESHOLD: f64 = 0.95;

/// Granular VEP splice sub-terms that vareffect collapses into
/// `splice_region_variant`.
const GRANULAR_SPLICE_TERMS: &[&str] = &[
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
];

/// VEP-only consequence terms that vareffect deliberately does not emit.
/// Filtered out of VEP's set before comparison.
const VEP_ONLY_TERMS: &[&str] = &["NMD_transcript_variant"];

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("could not determine workspace root from CARGO_MANIFEST_DIR")
        .to_path_buf()
}

fn ground_truth_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("vep_ground_truth.tsv")
}

fn mismatches_log_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("vep_large_concordance_mismatches.log")
}

fn store_accessions_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("store_accessions.txt")
}

// ---------------------------------------------------------------------------
// Test infrastructure (copied from vep_concordance_indel.rs for consistency)
// ---------------------------------------------------------------------------

fn load_store() -> TranscriptStore {
    let path = std::env::var("GRCH38_TRANSCRIPTS").unwrap_or_else(|_| {
        workspace_root()
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

fn load_fasta() -> FastaReader {
    let path = std::env::var("GRCH38_FASTA")
        .expect("GRCH38_FASTA env var must point to a GRCh38 genome binary");
    let aliases = workspace_root().join("data/vareffect/patch_chrom_aliases_grch38.csv");
    FastaReader::open_with_patch_aliases_and_assembly(
        Path::new(&path),
        Some(aliases.as_ref()),
        Assembly::GRCh38,
    )
    .unwrap_or_else(|e| panic!("failed to open GRCh38 FASTA at {path}: {e}"))
}

// ---------------------------------------------------------------------------
// TSV parsing
// ---------------------------------------------------------------------------

/// One row of ground truth from the Python generator.
///
/// Lines that begin with `#` are treated as comments and skipped. Field
/// ordering must match the generator's `TSV_COLUMNS` constant exactly.
#[derive(Debug)]
struct VepGroundTruth {
    chrom: String,
    /// 1-based (VCF convention).
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    /// Full RefSeq accession with version, e.g. `"NM_000546.6"`. An empty
    /// string means VEP returned no RefSeq transcript for this variant — the
    /// row is skipped rather than counted as a miss.
    transcript_id: String,
    /// Unnormalised VEP SO terms, pipe-split as they appeared in the TSV.
    consequence: Vec<String>,
    /// `"HIGH"` / `"MODERATE"` / `"LOW"` / `"MODIFIER"`.
    impact: String,
    /// Full VEP HGVS with accession prefix (e.g. `"NM_000546.6:c.742C>T"`).
    hgvsc_full: Option<String>,
    /// Full VEP HGVS protein with accession prefix.
    hgvsp_full: Option<String>,
    protein_start: Option<u32>,
}

fn parse_tsv(path: &Path) -> Vec<VepGroundTruth> {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("failed to open TSV at {}: {}", path.display(), e));
    let reader = BufReader::new(file);

    let mut out = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|e| panic!("IO error reading TSV line {lineno}: {e}"));
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        // Header row starts with "chrom" — skip it.
        if line.starts_with("chrom\t") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 17 {
            panic!(
                "TSV line {lineno} has {} fields, expected 17: {line}",
                fields.len()
            );
        }

        let opt = |s: &str| -> Option<String> {
            if s == "-" || s.is_empty() {
                None
            } else {
                Some(s.to_owned())
            }
        };
        let opt_u32 = |s: &str| -> Option<u32> { opt(s).and_then(|v| v.parse().ok()) };

        let pos: u64 = fields[1]
            .parse()
            .unwrap_or_else(|_| panic!("TSV line {lineno}: invalid pos: {}", fields[1]));

        let consequence: Vec<String> = if fields[5] == "-" {
            Vec::new()
        } else {
            fields[5].split('|').map(str::to_owned).collect()
        };

        out.push(VepGroundTruth {
            chrom: fields[0].to_owned(),
            pos,
            ref_allele: fields[2].to_owned(),
            alt_allele: fields[3].to_owned(),
            transcript_id: fields[4].to_owned(),
            consequence,
            impact: fields[6].to_owned(),
            hgvsc_full: opt(fields[7]),
            hgvsp_full: opt(fields[8]),
            protein_start: opt_u32(fields[9]),
            // fields[10..=16] — protein_end, amino_acids, codons, exon, intron,
            // strand, biotype — are present in the TSV for future passes but
            // aren't consumed by today's report, so they're skipped to keep
            // the struct free of unused-field warnings.
        });
    }
    out
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

/// Normalise VEP's consequence set into vareffect's vocabulary.
///
/// * Drop `VEP_ONLY_TERMS` (e.g. `NMD_transcript_variant`).
/// * Fold `GRANULAR_SPLICE_TERMS` into `splice_region_variant`.
fn normalize_vep_consequences(terms: &[String]) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for term in terms {
        if VEP_ONLY_TERMS.contains(&term.as_str()) {
            continue;
        }
        if GRANULAR_SPLICE_TERMS.contains(&term.as_str()) {
            out.insert("splice_region_variant".to_owned());
        } else {
            out.insert(term.clone());
        }
    }
    out
}

/// Outcome classification for a single variant row.
///
/// `VarEffectError` is not represented here because annotation errors are
/// handled in the test loop itself, before `compare_variant` is ever called.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Outcome {
    /// VEP's chosen transcript isn't in vareffect's store.
    TranscriptNotFound,
    /// Same accession, different version — counted separately from real misses.
    TranscriptVersionMismatch,
    /// Compared cleanly; fields below indicate whether each dimension passed.
    Compared(Comparison),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Comparison {
    /// `true` if strict equality of raw consequence sets (no normalisation).
    consequence_exact: bool,
    /// `true` after normalisation (NMD filter + granular splice fold).
    consequence_match: bool,
    /// Whether the normalisation changed anything — useful for the "known
    /// divergence" bucket count.
    consequence_normalized: bool,
    hgvsc_match: bool,
    hgvsp_match: bool,
    impact_match: bool,
    protein_start_match: bool,
    /// Mismatch description for the per-row log; empty if everything matched.
    detail: Option<String>,
}

/// Strip an optional `accession:` prefix from an HGVS string.
///
/// VEP's REST response always prefixes with the accession (e.g.
/// `NM_000546.6:c.742C>T`, `NP_000537.3:p.Arg248Trp`). vareffect's output
/// is **asymmetric**: `hgvs_c` also carries the transcript accession, but
/// `hgvs_p` returns bare `p....` notation. Stripping from both sides makes
/// the comparison insensitive to that asymmetry — the bare c./p. suffix is
/// what we actually want to diff.
fn strip_hgvs_accession(full: &str) -> &str {
    full.split_once(':').map(|(_, rest)| rest).unwrap_or(full)
}

/// Extract a human-readable message from a caught-panic payload.
///
/// `std::panic::catch_unwind` returns a `Box<dyn Any + Send>` which is
/// typically either a `String` (from `panic!("...")`) or a `&'static str`
/// (from `panic!("literal")`). Anything else falls back to a generic label.
fn panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_owned()
    } else {
        "<non-string panic payload>".to_owned()
    }
}

fn compare_variant(vep: &VepGroundTruth, ve_results: &[ConsequenceResult]) -> Outcome {
    // Find the transcript row matching VEP's pick.
    let ve = ve_results
        .iter()
        .find(|r| r.transcript == vep.transcript_id);
    let ve = match ve {
        Some(r) => r,
        None => {
            // Fall back to accession-without-version match for the version
            // mismatch bucket.
            let vep_accession = vep.transcript_id.split('.').next().unwrap_or("");
            if ve_results
                .iter()
                .any(|r| r.transcript.split('.').next().unwrap_or("") == vep_accession)
            {
                return Outcome::TranscriptVersionMismatch;
            }
            return Outcome::TranscriptNotFound;
        }
    };

    // --- Consequence comparison -------------------------------------------------
    let vep_raw: BTreeSet<String> = vep.consequence.iter().cloned().collect();
    let vep_normalized = normalize_vep_consequences(&vep.consequence);
    let ve_csq: BTreeSet<String> = ve
        .consequences
        .iter()
        .map(|c| c.as_str().to_owned())
        .collect();

    let consequence_exact = vep_raw == ve_csq;
    let consequence_match = vep_normalized.is_subset(&ve_csq);
    let consequence_normalized = vep_raw != vep_normalized;

    // --- HGVS c/p comparison ----------------------------------------------------
    // Strip any `accession:` prefix from both sides so the comparison sees
    // only the c./p. suffix. vareffect's `hgvs_c` carries the accession,
    // `hgvs_p` does not — both asymmetries dissolve here.
    let vep_c = vep.hgvsc_full.as_deref().map(strip_hgvs_accession);
    let ve_c = ve.hgvs_c.as_deref().map(strip_hgvs_accession);
    let hgvsc_match = vep_c == ve_c;

    let vep_p = vep.hgvsp_full.as_deref().map(strip_hgvs_accession);
    let ve_p = ve.hgvs_p.as_deref().map(strip_hgvs_accession);
    let hgvsp_match = vep_p == ve_p;

    // --- Impact comparison ------------------------------------------------------
    let impact_match = format!("{}", ve.impact) == vep.impact;

    // --- Protein position -------------------------------------------------------
    let protein_start_match = vep.protein_start == ve.protein_start;

    // Build a compact detail string for the log if anything diverged.
    let detail =
        if consequence_match && hgvsc_match && hgvsp_match && impact_match && protein_start_match {
            None
        } else {
            let mut parts = Vec::new();
            if !consequence_match {
                parts.push(format!("csq: vep={:?} ve={:?}", vep_normalized, ve_csq));
            }
            if !hgvsc_match {
                parts.push(format!("hgvsc: vep={:?} ve={:?}", vep_c, ve_c));
            }
            if !hgvsp_match {
                parts.push(format!("hgvsp: vep={:?} ve={:?}", vep_p, ve_p));
            }
            if !impact_match {
                parts.push(format!("impact: vep={} ve={}", vep.impact, ve.impact));
            }
            if !protein_start_match {
                parts.push(format!(
                    "protein_start: vep={:?} ve={:?}",
                    vep.protein_start, ve.protein_start
                ));
            }
            Some(parts.join(" | "))
        };

    Outcome::Compared(Comparison {
        consequence_exact,
        consequence_match,
        consequence_normalized,
        hgvsc_match,
        hgvsp_match,
        impact_match,
        protein_start_match,
        detail,
    })
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

#[derive(Default)]
struct Stats {
    total_tsv_rows: usize,
    skipped_no_refseq: usize,
    skipped_vareffect_error: usize,
    /// Panics caught from `VarEffect::annotate`. Tracked separately from
    /// `Result::Err` returns because an uncaught panic signals a hard bug
    /// (out-of-bounds slice, integer overflow, etc.) and is always worth
    /// surfacing in the report.
    vareffect_panic: usize,
    transcript_not_found: usize,
    transcript_version_mismatch: usize,
    compared: usize,
    consequence_exact: usize,
    consequence_match: usize,
    consequence_normalized: usize,
    hgvsc_match: usize,
    hgvsc_total_with_vep: usize,
    hgvsp_match: usize,
    hgvsp_total_with_vep: usize,
    impact_match: usize,
    protein_start_match: usize,
    protein_start_total_with_vep: usize,
    /// Top mismatch patterns for the summary.
    mismatch_patterns: BTreeMap<String, usize>,
}

fn pct(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        (num as f64 / den as f64) * 100.0
    }
}

fn write_mismatches_log(
    path: &Path,
    mismatches: &[(usize, &VepGroundTruth, &Comparison)],
    panics: &[(usize, &VepGroundTruth, &str)],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = File::create(path)?;
    writeln!(
        f,
        "# Large-scale concordance — mismatch log ({} mismatches, {} panics)",
        mismatches.len(),
        panics.len(),
    )?;

    if !panics.is_empty() {
        writeln!(f)?;
        writeln!(f, "# === PANICS ===")?;
        writeln!(f, "# lineno\tchrom:pos\tref>alt\ttranscript\tpanic_message")?;
        for (lineno, vep, message) in panics {
            writeln!(
                f,
                "{lineno}\t{}:{}\t{}>{}\t{}\t{}",
                vep.chrom, vep.pos, vep.ref_allele, vep.alt_allele, vep.transcript_id, message,
            )?;
        }
    }

    writeln!(f)?;
    writeln!(f, "# === MISMATCHES ===")?;
    writeln!(f, "# lineno\tchrom:pos\tref>alt\ttranscript\tdetail")?;
    for (lineno, vep, cmp) in mismatches {
        writeln!(
            f,
            "{lineno}\t{}:{}\t{}>{}\t{}\t{}",
            vep.chrom,
            vep.pos,
            vep.ref_allele,
            vep.alt_allele,
            vep.transcript_id,
            cmp.detail.as_deref().unwrap_or(""),
        )?;
    }
    Ok(())
}

fn print_report(stats: &Stats, elapsed_secs: f64, num_threads: usize) {
    let tested = stats.compared;
    eprintln!();
    eprintln!("============================================================");
    eprintln!("LARGE-SCALE VEP CONCORDANCE REPORT");
    eprintln!("============================================================");
    eprintln!("Total TSV rows:               {:>6}", stats.total_tsv_rows);
    eprintln!(
        "Skipped (no RefSeq tc):       {:>6}",
        stats.skipped_no_refseq
    );
    eprintln!(
        "Skipped (vareffect error):    {:>6}",
        stats.skipped_vareffect_error
    );
    eprintln!(
        "Skipped (vareffect PANIC):    {:>6}  <- hard bugs, see mismatches.log",
        stats.vareffect_panic
    );
    eprintln!(
        "Transcript not found:         {:>6}",
        stats.transcript_not_found
    );
    eprintln!(
        "Transcript version mismatch:  {:>6}",
        stats.transcript_version_mismatch
    );
    eprintln!("Compared:                     {tested:>6}");
    eprintln!();
    eprintln!(
        "Consequence concordance:   {:>5}/{tested} ({:.2}%)",
        stats.consequence_match,
        pct(stats.consequence_match, tested),
    );
    eprintln!(
        "  Exact match:             {:>5}       ({:.2}%)",
        stats.consequence_exact,
        pct(stats.consequence_exact, tested),
    );
    eprintln!(
        "  Normalised (splice/NMD): {:>5}       ({:.2}%)",
        stats.consequence_normalized,
        pct(stats.consequence_normalized, tested),
    );
    eprintln!(
        "  Real mismatch:           {:>5}       ({:.2}%)",
        tested.saturating_sub(stats.consequence_match),
        pct(tested.saturating_sub(stats.consequence_match), tested),
    );
    eprintln!();
    eprintln!(
        "HGVS c. concordance:       {:>5}/{} ({:.2}%)",
        stats.hgvsc_match,
        stats.hgvsc_total_with_vep,
        pct(stats.hgvsc_match, stats.hgvsc_total_with_vep),
    );
    eprintln!(
        "HGVS p. concordance:       {:>5}/{} ({:.2}%)",
        stats.hgvsp_match,
        stats.hgvsp_total_with_vep,
        pct(stats.hgvsp_match, stats.hgvsp_total_with_vep),
    );
    eprintln!(
        "Impact concordance:        {:>5}/{tested} ({:.2}%)",
        stats.impact_match,
        pct(stats.impact_match, tested),
    );
    eprintln!(
        "Protein start concordance: {:>5}/{} ({:.2}%)",
        stats.protein_start_match,
        stats.protein_start_total_with_vep,
        pct(
            stats.protein_start_match,
            stats.protein_start_total_with_vep
        ),
    );
    eprintln!();

    let mut patterns: Vec<(&String, &usize)> = stats.mismatch_patterns.iter().collect();
    patterns.sort_by(|a, b| b.1.cmp(a.1));
    if !patterns.is_empty() {
        eprintln!("TOP MISMATCH PATTERNS (up to 20):");
        for (pattern, count) in patterns.iter().take(20) {
            eprintln!("  {count:>5}x  {pattern}");
        }
        eprintln!();
    }

    eprintln!("Threads:    {num_threads}");
    if tested > 0 {
        let throughput = tested as f64 / elapsed_secs;
        eprintln!("Throughput: {throughput:.0} variants/second");
    }
    eprintln!("Elapsed:    {elapsed_secs:.1}s");
    eprintln!("============================================================");
}

/// Summarise a detail string into a pattern key for the top-N table.
///
/// Strips quoted values and keeps only the field labels + raw SO term names,
/// so "csq: vep={\"missense_variant\"} ve={\"missense_variant\", \"splice_region_variant\"}"
/// collapses alongside every other case of vareffect adding splice_region to a
/// missense call.
fn mismatch_pattern_key(cmp: &Comparison) -> String {
    let mut dims = Vec::new();
    if !cmp.consequence_match {
        dims.push("csq");
    }
    if !cmp.hgvsc_match {
        dims.push("hgvsc");
    }
    if !cmp.hgvsp_match {
        dims.push("hgvsp");
    }
    if !cmp.impact_match {
        dims.push("impact");
    }
    if !cmp.protein_start_match {
        dims.push("protein_start");
    }
    dims.join("+")
}

// ---------------------------------------------------------------------------
// Parallel annotation result
// ---------------------------------------------------------------------------

/// Per-row result produced during the parallel annotation phase.
///
/// Carries all data needed for sequential stats accumulation so the parallel
/// phase is pure map (no shared mutable state).
enum RowResult {
    /// VEP returned no RefSeq transcript for this variant.
    SkippedNoRefseq,
    /// `VarEffect::annotate` returned `Err`.
    VareffectError,
    /// `VarEffect::annotate` panicked (hard bug).
    VareffectPanic { lineno: usize, message: String },
    /// Annotation succeeded; `outcome` carries the comparison result.
    Completed { lineno: usize, outcome: Outcome },
}

// ---------------------------------------------------------------------------
// Test entry point
// ---------------------------------------------------------------------------

#[test]
#[ignore]
fn vep_large_concordance_grch38() {
    let tsv_path = ground_truth_path();
    if !tsv_path.exists() {
        eprintln!(
            "SKIP: VEP ground-truth TSV not found at {}",
            tsv_path.display()
        );
        eprintln!(
            "      Run: uv run vareffect/scripts/generate_vep_ground_truth.py \\
             --assembly grch38 \\
             --input data/clinvar/clinvar_grch38.vcf.gz \\
             --output {} \\
             --sample-size 50000",
            tsv_path.display()
        );
        return;
    }

    let ground_truth = parse_tsv(&tsv_path);
    eprintln!(
        "Loaded {} ground-truth rows from {}",
        ground_truth.len(),
        tsv_path.display()
    );

    // Thread count: default to 1 (single-threaded) unless CONCORDANCE_THREADS
    // is set. This preserves existing behavior; opt in to parallelism with e.g.
    // `CONCORDANCE_THREADS=8`.
    let num_threads = std::env::var("CONCORDANCE_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1)
        .max(1);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("failed to build rayon thread pool");

    eprintln!("Using {num_threads} thread(s) for annotation");

    let store = load_store();
    let fasta = load_fasta();
    let ve = VarEffect::builder()
        .with_handles(Assembly::GRCh38, store, fasta)
        .expect("matching assemblies")
        .build()
        .expect("builder");

    // Silence the default panic hook for the duration of the scan. Without
    // this, every caught panic from `VarEffect::annotate` would print a full
    // stack trace to stderr, drowning out the concordance report. The hook
    // is restored before the final assert so any panic inside the report
    // itself still prints normally.
    let prev_hook = panic::take_hook();
    panic::set_hook(Box::new(|_info| {}));

    let start = Instant::now();

    // Phase 1 — Parallel annotation + comparison.
    //
    // Each row is independently annotated and compared. The VarEffect
    // instance is read-only (`Send + Sync` via Arc<Mmap>, Arc<[T]>),
    // so `&ve` is safely shared across rayon worker threads.
    let row_results: Vec<RowResult> = pool.install(|| {
        ground_truth
            .par_iter()
            .enumerate()
            .map(|(lineno, vep)| {
                if vep.transcript_id.is_empty() {
                    return RowResult::SkippedNoRefseq;
                }

                let ve_chrom = format!("chr{}", vep.chrom);
                let ve_pos = vep.pos.saturating_sub(1);

                let annotate_result = panic::catch_unwind(AssertUnwindSafe(|| {
                    ve.annotate(
                        Assembly::GRCh38,
                        &ve_chrom,
                        ve_pos,
                        vep.ref_allele.as_bytes(),
                        vep.alt_allele.as_bytes(),
                    )
                }));

                let results = match annotate_result {
                    Ok(Ok(r)) => r,
                    Ok(Err(_)) => return RowResult::VareffectError,
                    Err(panic_payload) => {
                        let message = panic_message(panic_payload.as_ref());
                        return RowResult::VareffectPanic { lineno, message };
                    }
                };

                RowResult::Completed {
                    lineno,
                    outcome: compare_variant(vep, &results.consequences),
                }
            })
            .collect()
    });

    let elapsed = start.elapsed().as_secs_f64();

    // Restore the default panic hook now that the scan is done.
    panic::set_hook(prev_hook);

    // Phase 2 — Sequential stats accumulation.
    let mut stats = Stats {
        total_tsv_rows: ground_truth.len(),
        ..Stats::default()
    };
    let mut mismatch_rows: Vec<(usize, &VepGroundTruth, Comparison)> = Vec::new();
    let mut panic_rows: Vec<(usize, &VepGroundTruth, String)> = Vec::new();

    for (idx, result) in row_results.into_iter().enumerate() {
        match result {
            RowResult::SkippedNoRefseq => {
                stats.skipped_no_refseq += 1;
            }
            RowResult::VareffectError => {
                stats.skipped_vareffect_error += 1;
            }
            RowResult::VareffectPanic { lineno, message } => {
                stats.vareffect_panic += 1;
                if stats.vareffect_panic <= 10 {
                    let vep = &ground_truth[idx];
                    eprintln!(
                        "PANIC chr{}:{} {}>{} ({}): {}",
                        vep.chrom,
                        vep.pos.saturating_sub(1),
                        vep.ref_allele,
                        vep.alt_allele,
                        vep.transcript_id,
                        message,
                    );
                }
                panic_rows.push((lineno, &ground_truth[idx], message));
            }
            RowResult::Completed { lineno, outcome } => match outcome {
                Outcome::TranscriptNotFound => {
                    stats.transcript_not_found += 1;
                }
                Outcome::TranscriptVersionMismatch => {
                    stats.transcript_version_mismatch += 1;
                }
                Outcome::Compared(cmp) => {
                    let vep = &ground_truth[idx];
                    stats.compared += 1;
                    if cmp.consequence_exact {
                        stats.consequence_exact += 1;
                    }
                    if cmp.consequence_match {
                        stats.consequence_match += 1;
                    }
                    if cmp.consequence_normalized {
                        stats.consequence_normalized += 1;
                    }
                    if vep.hgvsc_full.is_some() {
                        stats.hgvsc_total_with_vep += 1;
                        if cmp.hgvsc_match {
                            stats.hgvsc_match += 1;
                        }
                    }
                    if vep.hgvsp_full.is_some() {
                        stats.hgvsp_total_with_vep += 1;
                        if cmp.hgvsp_match {
                            stats.hgvsp_match += 1;
                        }
                    }
                    if cmp.impact_match {
                        stats.impact_match += 1;
                    }
                    if vep.protein_start.is_some() {
                        stats.protein_start_total_with_vep += 1;
                        if cmp.protein_start_match {
                            stats.protein_start_match += 1;
                        }
                    }
                    if cmp.detail.is_some() {
                        let key = mismatch_pattern_key(&cmp);
                        *stats.mismatch_patterns.entry(key).or_default() += 1;
                        mismatch_rows.push((lineno, vep, cmp));
                    }
                }
            },
        }
    }

    // Dump the per-row mismatch and panic log to disk for follow-up triage.
    let mismatch_refs: Vec<(usize, &VepGroundTruth, &Comparison)> = mismatch_rows
        .iter()
        .map(|(lineno, vep, cmp)| (*lineno, *vep, cmp))
        .collect();
    let panic_refs: Vec<(usize, &VepGroundTruth, &str)> = panic_rows
        .iter()
        .map(|(lineno, vep, msg)| (*lineno, *vep, msg.as_str()))
        .collect();
    if let Err(e) = write_mismatches_log(&mismatches_log_path(), &mismatch_refs, &panic_refs) {
        eprintln!("WARN: failed to write mismatch log: {e}");
    }

    print_report(&stats, elapsed, num_threads);

    // Hard assertion — deliberately softer than the long-term ≥99% goal so
    // the harness lands before every divergence is triaged.
    assert!(
        stats.compared > 0,
        "no variants compared; is the TSV empty or corrupt?"
    );
    let rate = stats.consequence_match as f64 / stats.compared as f64;
    assert!(
        rate >= CONSEQUENCE_THRESHOLD,
        "consequence concordance {:.2}% below {:.0}% threshold ({} / {})",
        rate * 100.0,
        CONSEQUENCE_THRESHOLD * 100.0,
        stats.consequence_match,
        stats.compared,
    );
}

/// One-shot helper that dumps every transcript accession in vareffect's
/// store to `tests/data/store_accessions.txt` (one per line). The file is
/// consumed by the Python ground-truth generator to restrict VEP's transcript
/// pick to accessions vareffect actually has — otherwise ~72% of VEP's best
/// picks land on transcripts absent from the store and the comparison
/// corpus shrinks to ~27% of the TSV rows.
///
/// Run after rebuilding `data/vareffect/transcript_models_grch38.bin`:
/// ```bash
/// cargo test -p vareffect -- --ignored dump_store_accessions --nocapture
/// ```
#[test]
#[ignore]
fn dump_store_accessions() {
    let store = load_store();
    let out_path = store_accessions_path();
    if let Some(parent) = out_path.parent() {
        std::fs::create_dir_all(parent).expect("create tests/data dir");
    }
    let mut f = File::create(&out_path)
        .unwrap_or_else(|e| panic!("failed to create {}: {}", out_path.display(), e));

    // Collect and sort so the output is deterministic across runs. The
    // store's internal order is already accession-sorted at build time, but
    // guarding against future changes is cheap.
    let mut accessions: Vec<&str> = store
        .transcripts()
        .iter()
        .map(|t| t.accession.as_str())
        .collect();
    accessions.sort_unstable();

    for accession in &accessions {
        writeln!(f, "{accession}").expect("write accession");
    }

    eprintln!(
        "Dumped {} accessions to {}",
        accessions.len(),
        out_path.display()
    );
}
