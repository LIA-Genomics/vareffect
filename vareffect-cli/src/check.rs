//! Handler for the `vareffect check` subcommand.
//!
//! Validates that the config file exists and parses, the data directory
//! is present, and the genome binary + transcript model store are in
//! place and readable. Reports a pass/fail summary and returns a
//! non-zero exit code when any check fails.

use std::fs::{self, File};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Args;
use console::style;

use crate::common::format_size;
use crate::config;

/// Arguments for the `check` subcommand.
#[derive(Debug, Args)]
pub struct CheckArgs {
    /// Path to config file. Auto-discovered if omitted.
    #[arg(long)]
    pub config: Option<PathBuf>,

    /// Only print failures.
    #[arg(long)]
    pub quiet: bool,
}

/// Status of a single validation check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CheckStatus {
    Pass,
    Fail,
}

/// Result of a single validation check.
#[derive(Debug)]
struct CheckResult {
    status: CheckStatus,
    detail: String,
}

impl CheckResult {
    fn pass(detail: impl Into<String>) -> Self {
        Self {
            status: CheckStatus::Pass,
            detail: detail.into(),
        }
    }

    fn fail(detail: impl Into<String>) -> Self {
        Self {
            status: CheckStatus::Fail,
            detail: detail.into(),
        }
    }
}

/// Print a single check result to stderr.
fn print_result(result: &CheckResult, quiet: bool) {
    if quiet && result.status == CheckStatus::Pass {
        return;
    }
    match result.status {
        CheckStatus::Pass => {
            eprintln!("  {}  {}", style("OK").green().bold(), result.detail,);
        }
        CheckStatus::Fail => {
            eprintln!("  {}  {}", style("FAIL").red().bold(), result.detail,);
        }
    }
}

/// Check whether a file exists and is readable (by opening it).
///
/// Returns `Ok(())` if the file can be opened for reading, or an
/// error describing why it cannot.
fn check_file_readable(path: &Path) -> Result<()> {
    File::open(path).with_context(|| format!("{}", path.display()))?;
    Ok(())
}

/// Read the `record_count` field from a `.manifest.json` sidecar.
///
/// Returns `None` if the manifest does not exist or cannot be parsed.
fn read_manifest_record_count(bin_path: &Path) -> Option<u64> {
    let manifest_path = bin_path.with_extension("manifest.json");
    let content = fs::read_to_string(&manifest_path).ok()?;
    let json: serde_json::Value = serde_json::from_str(&content).ok()?;
    json.get("record_count")?.as_u64()
}

/// Run the `vareffect check` subcommand.
///
/// # Returns
///
/// `Ok(true)` if all checks pass, `Ok(false)` if any check fails.
/// `Err(...)` only for unexpected IO errors during the check itself.
pub fn run(args: &CheckArgs) -> Result<bool> {
    let mut results: Vec<CheckResult> = Vec::new();

    // C1: Config file exists and parses.
    let config_path = config::find_config(args.config.as_deref());
    let parsed_config = match &config_path {
        Ok(path) => match fs::read_to_string(path) {
            Ok(content) => match toml::from_str::<config::VareffectConfig>(&content) {
                Ok(cfg) => {
                    results.push(CheckResult::pass(format!(
                        "Config loaded from {}",
                        path.display(),
                    )));
                    Some((cfg, path.clone()))
                }
                Err(e) => {
                    results.push(CheckResult::fail(format!(
                        "Config parse error in {}: {e}",
                        path.display(),
                    )));
                    None
                }
            },
            Err(e) => {
                results.push(CheckResult::fail(format!(
                    "Config not readable at {}: {e}",
                    path.display(),
                )));
                None
            }
        },
        Err(_) => {
            results.push(CheckResult::fail("Config not found. Run: vareffect init"));
            None
        }
    };

    // C2-C5 depend on a parsed config.
    if let Some((cfg, _config_path)) = &parsed_config {
        let output_dir = PathBuf::from(&cfg.vareffect.output_dir);

        // C2: output_dir exists and is readable.
        match fs::read_dir(&output_dir) {
            Ok(_) => {
                results.push(CheckResult::pass(format!(
                    "Data directory: {}",
                    output_dir.display(),
                )));
            }
            Err(_) => {
                results.push(CheckResult::fail(format!(
                    "Data directory missing: {}",
                    output_dir.display(),
                )));
            }
        }

        // C3..C5: Per-assembly artifacts. Iterate over every entry the
        // config defines so the check command works for builds that
        // include only one assembly.
        let entries: Vec<(&'static str, &crate::config::AssemblyEntry)> = [
            cfg.vareffect.grch38.as_ref().map(|e| ("GRCh38", e)),
            cfg.vareffect.grch37.as_ref().map(|e| ("GRCh37", e)),
        ]
        .into_iter()
        .flatten()
        .collect();
        if entries.is_empty() {
            results.push(CheckResult::fail(
                "Config defines no assembly sub-table; \
                 add [vareffect.grch38] or [vareffect.grch37] then run `vareffect setup`",
            ));
        }
        for (label, entry) in &entries {
            let genome_path = output_dir.join(&entry.fasta_filename);
            match check_file_readable(&genome_path) {
                Ok(()) => {
                    let size = fs::metadata(&genome_path).map(|m| m.len()).unwrap_or(0);
                    results.push(CheckResult::pass(format!(
                        "{label} genome: {} ({})",
                        genome_path.display(),
                        format_size(size),
                    )));
                }
                Err(_) => {
                    results.push(CheckResult::fail(format!(
                        "{label} genome binary not found: {}",
                        genome_path.display(),
                    )));
                }
            }

            let idx_path = genome_path.with_extension("bin.idx");
            match check_file_readable(&idx_path) {
                Ok(()) => {
                    results.push(CheckResult::pass(format!(
                        "{label} genome index: {}",
                        idx_path.display(),
                    )));
                }
                Err(_) => {
                    results.push(CheckResult::fail(format!(
                        "{label} genome index not found: {}",
                        idx_path.display(),
                    )));
                }
            }

            let transcript_path = output_dir.join(&entry.transcript_models_filename);
            match check_file_readable(&transcript_path) {
                Ok(()) => {
                    let count_str = read_manifest_record_count(&transcript_path)
                        .map(|n| format!(" ({n} transcripts)"))
                        .unwrap_or_default();
                    results.push(CheckResult::pass(format!(
                        "{label} transcript DB: {}{count_str}",
                        transcript_path.display(),
                    )));
                }
                Err(_) => {
                    results.push(CheckResult::fail(format!(
                        "{label} transcript DB not found: {}",
                        transcript_path.display(),
                    )));
                }
            }
        }
    } else {
        // Config failed to load; mark all dependent checks as failed.
        let dependent = [
            "Data directory",
            "Genome binary",
            "Genome index",
            "Transcript DB",
        ];
        for name in &dependent {
            results.push(CheckResult::fail(format!(
                "{name}: skipped (no config). Run: vareffect init",
            )));
        }
    }

    // Print results.
    eprintln!();
    for result in &results {
        print_result(result, args.quiet);
    }

    let passed = results
        .iter()
        .filter(|r| r.status == CheckStatus::Pass)
        .count();
    let failed = results.len() - passed;

    eprintln!();
    if failed == 0 {
        eprintln!("{}/{} checks passed.", passed, results.len(),);
    } else {
        eprintln!(
            "{}/{} checks passed. {} {} found -- see above.",
            passed,
            results.len(),
            failed,
            if failed == 1 { "issue" } else { "issues" },
        );
    }

    Ok(failed == 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::fs;

    use tempfile::TempDir;

    /// Write a minimal valid config to a temp directory and return
    /// the path to the config file.
    fn write_test_config(dir: &Path, output_dir: &str) -> PathBuf {
        let config_path = dir.join("vareffect_build.toml");
        let output_dir = output_dir.replace('\\', "/");
        let content = format!(
            r#"
[vareffect]
output_dir = "{output_dir}"

[vareffect.grch38]
fasta_url = "https://example.com/GRCh38.fna.gz"
fasta_filename = "GRCh38.bin"
fasta_build = "GRCh38"
gff_url = "https://example.com/mane.gff.gz"
transcript_models_filename = "transcript_models_grch38.bin"
transcript_models_basename = "transcript_models_grch38"
transcript_models_version = "1.5"
transcript_source = "mane"
assembly_report_url = "https://example.com/report.txt"
assembly_report_input = "report.txt"
patch_aliases_filename = "patch_chrom_aliases_grch38.csv"
"#,
        );
        fs::write(&config_path, content).unwrap();
        config_path
    }

    #[test]
    fn check_passes_for_config_and_data_dir() {
        let tmp = TempDir::new().unwrap();
        let data_dir = tmp.path().join("data");
        fs::create_dir_all(&data_dir).unwrap();
        let config_path = write_test_config(tmp.path(), &data_dir.display().to_string());

        let args = CheckArgs {
            config: Some(config_path),
            quiet: false,
        };
        let result = run(&args).expect("check should not error");
        // Config and data dir pass, but genome/transcript files are
        // missing so overall result is false.
        assert!(!result);
    }

    #[test]
    fn check_fails_for_missing_config() {
        let tmp = TempDir::new().unwrap();
        let config_path = tmp.path().join("nonexistent.toml");

        let args = CheckArgs {
            config: Some(config_path),
            quiet: false,
        };
        let result = run(&args).expect("check should not error");
        assert!(!result);
    }

    #[test]
    fn check_passes_when_all_files_present() {
        let tmp = TempDir::new().unwrap();
        let data_dir = tmp.path().join("data");
        fs::create_dir_all(&data_dir).unwrap();

        // Create dummy data files matching the per-assembly filenames in
        // the test config (GRCh38 only).
        fs::write(data_dir.join("GRCh38.bin"), b"genome").unwrap();
        fs::write(data_dir.join("GRCh38.bin.idx"), b"index").unwrap();
        fs::write(
            data_dir.join("transcript_models_grch38.bin"),
            b"transcripts",
        )
        .unwrap();

        let config_path = write_test_config(tmp.path(), &data_dir.display().to_string());

        let args = CheckArgs {
            config: Some(config_path),
            quiet: false,
        };
        let result = run(&args).expect("check should not error");
        assert!(result);
    }

    #[test]
    fn manifest_record_count_read() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("transcript_models.bin");
        fs::write(&bin_path, b"data").unwrap();

        let manifest = serde_json::json!({
            "store_name": "transcript_models",
            "record_count": 1234,
        });
        fs::write(
            tmp.path().join("transcript_models.manifest.json"),
            serde_json::to_string(&manifest).unwrap(),
        )
        .unwrap();

        assert_eq!(read_manifest_record_count(&bin_path), Some(1234));
    }

    #[test]
    fn manifest_record_count_missing_manifest() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("transcript_models.bin");
        assert_eq!(read_manifest_record_count(&bin_path), None);
    }
}
