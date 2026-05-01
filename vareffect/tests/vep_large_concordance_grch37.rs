//! Large-scale VEP concordance harness — GRCh37.
//!
//! Counterpart to `vep_large_concordance.rs` (GRCh38). Reads a pre-computed
//! TSV of VEP REST API ground truth from `https://grch37.rest.ensembl.org`,
//! runs vareffect on every row, compares per-field, and prints a
//! concordance report. Statistical: asserts a ≥ 95 % consequence-concordance
//! rate across the full TSV, matching the GRCh38 sibling's threshold (the
//! ≥ 99 % long-term goal kicks in once divergences are triaged).
//!
//! `#[ignore]`-gated because it requires the GRCh37 transcript store, the
//! GRCh37 reference FASTA, and the committed ground-truth TSV.
//!
//! Run with:
//! ```bash
//! GRCH37_FASTA=/abs/path/to/data/vareffect/GRCh37.bin \
//! GRCH37_TRANSCRIPTS=/abs/path/to/data/vareffect/transcript_models_grch37.bin \
//!   cargo test -p vareffect --release -- --ignored vep_large_concordance_grch37 --nocapture
//! ```
//!
//! # Output files
//!
//! * Report — printed to stderr (`--nocapture` to see it live)
//! * Mismatch dump — `vareffect/tests/data/vep_large_concordance_grch37_mismatches.log`
//!   (gitignored, overwritten every run)
//!
//! # Generating the TSV
//!
//! ```bash
//! uv run vareffect/scripts/generate_vep_ground_truth.py \
//!   --assembly grch37 \
//!   --input  data/clinvar/clinvar_grch37.vcf.gz \
//!   --output vareffect/tests/data/vep_ground_truth_grch37.tsv \
//!   --sample-size 50000
//! ```
//!
//! # Divergence taxonomy
//!
//! Same normalization rules as the GRCh38 harness:
//! * **Granular splice terms** folded into `splice_region_variant`
//!   (VEP's `splice_donor_5th_base_variant`, `splice_donor_region_variant`,
//!   `splice_polypyrimidine_tract_variant`).
//! * **`NMD_transcript_variant`** filtered out (vareffect tracks NMD via
//!   the `predicts_nmd` flag instead).

use std::any::Any;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;
use vareffect::{Assembly, ConsequenceResult, FastaReader, TranscriptStore, VarEffect};

const CONSEQUENCE_THRESHOLD: f64 = 0.95;

const GRANULAR_SPLICE_TERMS: &[&str] = &[
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
];

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
        .join("vep_ground_truth_grch37.tsv")
}

fn mismatches_log_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("vep_large_concordance_grch37_mismatches.log")
}

fn load_store() -> TranscriptStore {
    let path = std::env::var("GRCH37_TRANSCRIPTS").unwrap_or_else(|_| {
        workspace_root()
            .join("data/vareffect/transcript_models_grch37.bin")
            .to_string_lossy()
            .into_owned()
    });
    TranscriptStore::load_from_path(Path::new(&path)).unwrap_or_else(|e| {
        panic!(
            "failed to load GRCh37 transcript store from {path}: {e}. \
             Run `vareffect setup --assembly grch37` first.",
        )
    })
}

fn load_fasta() -> FastaReader {
    let path = std::env::var("GRCH37_FASTA")
        .expect("GRCH37_FASTA env var must point to a GRCh37 genome binary (.bin)");
    let aliases = workspace_root().join("data/vareffect/patch_chrom_aliases_grch37.csv");
    FastaReader::open_with_patch_aliases_and_assembly(
        Path::new(&path),
        Some(aliases.as_ref()),
        Assembly::GRCh37,
    )
    .unwrap_or_else(|e| panic!("failed to open GRCh37 FASTA at {path}: {e}"))
}

#[derive(Debug)]
struct VepGroundTruth {
    chrom: String,
    /// 1-based VCF position.
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    transcript_id: String,
    consequence: Vec<String>,
    impact: String,
    hgvsc_full: Option<String>,
    hgvsp_full: Option<String>,
    protein_start: Option<u32>,
}

fn parse_tsv(path: &Path) -> Vec<VepGroundTruth> {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("failed to open TSV at {}: {e}", path.display()));
    let reader = BufReader::new(file);

    let mut out = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|e| panic!("IO error reading TSV line {lineno}: {e}"));
        if line.starts_with('#') || line.is_empty() || line.starts_with("chrom\t") {
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
        let consequence: Vec<String> = if fields[5] == "-" || fields[5].is_empty() {
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
        });
    }
    out
}

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

#[derive(Debug, Clone, PartialEq, Eq)]
enum Outcome {
    TranscriptNotFound,
    TranscriptVersionMismatch,
    Compared(Comparison),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Comparison {
    consequence_exact: bool,
    consequence_match: bool,
    consequence_normalized: bool,
    hgvsc_match: bool,
    hgvsp_match: bool,
    impact_match: bool,
    protein_start_match: bool,
    detail: Option<String>,
}

fn strip_hgvs_accession(full: &str) -> &str {
    full.split_once(':').map(|(_, rest)| rest).unwrap_or(full)
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

fn compare_variant(vep: &VepGroundTruth, ve_results: &[ConsequenceResult]) -> Outcome {
    let ve = ve_results
        .iter()
        .find(|r| r.transcript == vep.transcript_id);
    let ve = match ve {
        Some(r) => r,
        None => {
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

    let vep_c = vep.hgvsc_full.as_deref().map(strip_hgvs_accession);
    let ve_c = ve.hgvs_c.as_deref().map(strip_hgvs_accession);
    let hgvsc_match = vep_c == ve_c;

    let vep_p = vep.hgvsp_full.as_deref().map(strip_hgvs_accession);
    let ve_p = ve.hgvs_p.as_deref().map(strip_hgvs_accession);
    let hgvsp_match = vep_p == ve_p;

    let impact_match = format!("{}", ve.impact) == vep.impact;
    let protein_start_match = vep.protein_start == ve.protein_start;

    let detail =
        if consequence_match && hgvsc_match && hgvsp_match && impact_match && protein_start_match {
            None
        } else {
            let mut parts = Vec::new();
            if !consequence_match {
                parts.push(format!("csq: vep={vep_normalized:?} ve={ve_csq:?}"));
            }
            if !hgvsc_match {
                parts.push(format!("hgvsc: vep={vep_c:?} ve={ve_c:?}"));
            }
            if !hgvsp_match {
                parts.push(format!("hgvsp: vep={vep_p:?} ve={ve_p:?}"));
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

#[derive(Default)]
struct Stats {
    total_tsv_rows: usize,
    skipped_no_refseq: usize,
    skipped_vareffect_error: usize,
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
    mismatch_patterns: BTreeMap<String, usize>,
}

fn pct(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        (num as f64 / den as f64) * 100.0
    }
}

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
        "# Large-scale concordance (GRCh37) — mismatch log ({} mismatches, {} panics)",
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
    eprintln!("LARGE-SCALE VEP CONCORDANCE REPORT (GRCh37)");
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

enum RowResult {
    SkippedNoRefseq,
    VareffectError,
    VareffectPanic { lineno: usize, message: String },
    Completed { lineno: usize, outcome: Outcome },
}

#[test]
#[ignore]
fn vep_large_concordance_grch37() {
    let tsv_path = ground_truth_path();
    if !tsv_path.exists() {
        eprintln!(
            "SKIP: GRCh37 VEP ground-truth TSV not found at {}",
            tsv_path.display()
        );
        eprintln!(
            "      Run: uv run vareffect/scripts/generate_vep_ground_truth.py \\
                  --assembly grch37 \\
                  --input data/clinvar/clinvar_grch37.vcf.gz \\
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
        .with_handles(Assembly::GRCh37, store, fasta)
        .expect("matching assemblies")
        .build()
        .expect("builder");

    let prev_hook = panic::take_hook();
    panic::set_hook(Box::new(|_info| {}));

    let start = Instant::now();
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
                        Assembly::GRCh37,
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
    panic::set_hook(prev_hook);

    let mut stats = Stats {
        total_tsv_rows: ground_truth.len(),
        ..Stats::default()
    };
    let mut mismatch_rows: Vec<(usize, &VepGroundTruth, Comparison)> = Vec::new();
    let mut panic_rows: Vec<(usize, &VepGroundTruth, String)> = Vec::new();

    for (idx, result) in row_results.into_iter().enumerate() {
        match result {
            RowResult::SkippedNoRefseq => stats.skipped_no_refseq += 1,
            RowResult::VareffectError => stats.skipped_vareffect_error += 1,
            RowResult::VareffectPanic { lineno, message } => {
                stats.vareffect_panic += 1;
                if stats.vareffect_panic <= 10 {
                    let vep = &ground_truth[idx];
                    eprintln!(
                        "PANIC chr{}:{} {}>{} ({}): {message}",
                        vep.chrom,
                        vep.pos.saturating_sub(1),
                        vep.ref_allele,
                        vep.alt_allele,
                        vep.transcript_id,
                    );
                }
                panic_rows.push((lineno, &ground_truth[idx], message));
            }
            RowResult::Completed { lineno, outcome } => match outcome {
                Outcome::TranscriptNotFound => stats.transcript_not_found += 1,
                Outcome::TranscriptVersionMismatch => stats.transcript_version_mismatch += 1,
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
