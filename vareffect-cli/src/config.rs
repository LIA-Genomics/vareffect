//! Configuration discovery and TOML parsing for `vareffect-cli`.
//!
//! The config file (`vareffect_build.toml`) defines URLs, filenames, and
//! output directories for the `setup` subcommand. The parser uses
//! `#[serde(default)]` on optional sections and silently ignores unknown
//! keys, so a single config file can be shared with other build tools
//! that add their own top-level sections.

use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use serde::Deserialize;

/// Config filename auto-discovered by [`find_config`].
const CONFIG_FILENAME: &str = "vareffect_build.toml";

/// Minimal view of the config TOML used by the `setup` subcommand.
///
/// `serde` silently ignores any other top-level sections, so a single
/// `vareffect_build.toml` can be shared with other build tools that add
/// their own sections without forcing this parser to know about them.
#[derive(Debug, Deserialize)]
pub struct VareffectConfig {
    /// Default directory paths (optional `[paths]` table).
    #[serde(default)]
    pub paths: PathsConfig,
    /// Required `[vareffect]` section with all URLs and filenames.
    pub vareffect: VareffectSection,
    /// Optional `[mane]` section -- only the `input` field is read, to
    /// override the default GFF3 cache filename when a downstream pipeline
    /// needs to share a pre-fetched copy.
    mane: Option<ManeSection>,
}

impl VareffectConfig {
    /// GFF3 cache filename: from `[mane].input` if present, else the
    /// default filename used by the `setup` subcommand.
    pub fn gff_filename(&self) -> &str {
        self.mane
            .as_ref()
            .map(|m| m.input.as_str())
            .unwrap_or(DEFAULT_MANE_GFF_FILENAME)
    }
}

/// Default filename for the cached MANE GFF3 download when the optional
/// `[mane]` section is absent.
const DEFAULT_MANE_GFF_FILENAME: &str = "mane_grch38.gff.gz";

/// Shared `[paths]` table. Only `raw_dir` is used by `vareffect-cli`.
#[derive(Debug, Default, Deserialize)]
pub struct PathsConfig {
    /// Directory for downloaded raw source files (e.g. `data/raw`).
    pub raw_dir: Option<String>,
}

/// Optional `[mane]` section -- only `input` is consumed; other fields
/// (`version`, `url`, ...) are silently ignored by serde.
#[derive(Debug, Deserialize)]
struct ManeSection {
    input: String,
}

/// The `[vareffect]` TOML section.
///
/// Field names match the TOML keys verbatim. All fields are required --
/// there is no meaningful default for a URL or filename, and defaulting
/// would silently mask a misconfigured config.
#[derive(Debug, Deserialize)]
pub struct VareffectSection {
    // Reference FASTA
    /// NCBI FTP URL for the gzipped reference FASTA.
    pub fasta_url: String,
    /// Output filename for the flat binary genome (e.g. `"GRCh38.bin"`).
    pub fasta_filename: String,
    /// Assembly build name (e.g. `"GRCh38"`). Sanity-checked against
    /// `fasta_url` to prevent assembly mismatches.
    pub fasta_build: String,

    // Transcript models
    /// URL to the MANE GFF3 release file.
    pub gff_url: String,
    /// URL to the MANE summary TSV file.
    pub summary_url: String,
    /// Cache filename for the summary TSV in `raw_dir`.
    pub summary_input: String,
    /// MANE release version string (e.g. `"1.5"`).
    pub transcript_models_version: String,

    // Patch-contig aliases (NCBI assembly report -> patch_chrom_aliases.csv)
    /// URL to the NCBI assembly report.
    pub assembly_report_url: String,
    /// Cache filename for the assembly report in `raw_dir`.
    pub assembly_report_input: String,

    // Output
    /// Directory for the final runtime files (`GRCh38.fa.gz`, `.fai`,
    /// `transcript_models.bin`, etc.).
    pub output_dir: String,
}

/// Locate the config file for `vareffect-cli`.
///
/// Resolution order:
/// 1. Explicit `--config` flag
/// 2. `VAREFFECT_BUILD_CONFIG` environment variable
/// 3. Well-known paths: `crates/vareffect-cli/vareffect_build.toml` (repo
///    root), `/etc/vareffect/vareffect_build.toml` (system),
///    `~/.config/vareffect/vareffect_build.toml` (user)
///
/// # Errors
///
/// Returns an error if no config file is found at any location.
pub fn find_config(explicit: Option<&Path>) -> Result<PathBuf> {
    // 1. Explicit --config flag.
    if let Some(p) = explicit {
        return Ok(p.to_path_buf());
    }

    // 2. Environment variable.
    if let Ok(p) = std::env::var("VAREFFECT_BUILD_CONFIG") {
        let path = PathBuf::from(p);
        if path.exists() {
            return Ok(path);
        }
        bail!(
            "VAREFFECT_BUILD_CONFIG is set to {} but the file does not exist",
            path.display()
        );
    }

    // 3. Well-known paths.
    let candidates: Vec<PathBuf> = [
        Some(PathBuf::from(format!("vareffect-cli/{CONFIG_FILENAME}"))),
        Some(PathBuf::from(format!("/etc/vareffect/{CONFIG_FILENAME}"))),
        dirs::config_dir().map(|d| d.join(format!("vareffect/{CONFIG_FILENAME}"))),
    ]
    .into_iter()
    .flatten()
    .collect();

    for path in &candidates {
        if path.exists() {
            return Ok(path.clone());
        }
    }

    bail!(
        "no config found -- pass --config, set VAREFFECT_BUILD_CONFIG, \
         or place {} in one of:\n{}",
        CONFIG_FILENAME,
        candidates
            .iter()
            .map(|p| format!("  - {}", p.display()))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_minimal_vareffect_config() {
        let toml = r#"
[vareffect]
fasta_url = "https://example.com/genome.fna.gz"
fasta_filename = "GRCh38.bin"
fasta_build = "GRCh38"
gff_url = "https://example.com/mane.gff.gz"
summary_url = "https://example.com/summary.txt.gz"
summary_input = "summary.tsv.gz"
transcript_models_version = "1.5"
assembly_report_url = "https://example.com/report.txt"
assembly_report_input = "report.txt"
output_dir = "data/vareffect"
"#;
        let config: VareffectConfig = toml::from_str(toml).unwrap();
        assert_eq!(config.vareffect.fasta_build, "GRCh38");
        assert_eq!(config.gff_filename(), DEFAULT_MANE_GFF_FILENAME);
        assert!(config.paths.raw_dir.is_none());
    }

    #[test]
    fn parse_config_with_mane_section_overrides_gff_filename() {
        let toml = r#"
[paths]
raw_dir = "data/raw"

[mane]
input = "custom_mane.gff.gz"
version = "1.3"

[vareffect]
fasta_url = "https://example.com/genome.fna.gz"
fasta_filename = "GRCh38.bin"
fasta_build = "GRCh38"
gff_url = "https://example.com/mane.gff.gz"
summary_url = "https://example.com/summary.txt.gz"
summary_input = "summary.tsv.gz"
transcript_models_version = "1.5"
assembly_report_url = "https://example.com/report.txt"
assembly_report_input = "report.txt"
output_dir = "data/vareffect"
"#;
        let config: VareffectConfig = toml::from_str(toml).unwrap();
        assert_eq!(config.gff_filename(), "custom_mane.gff.gz");
        assert_eq!(config.paths.raw_dir.as_deref(), Some("data/raw"));
    }

    #[test]
    fn parse_ignores_unknown_sections() {
        // Compatibility: extra top-level sections are silently skipped so
        // a shared config file with other tools' sections still parses.
        let toml = r#"
[constraint]
input = "gnomad.tsv"
version = "4.1"

[vareffect]
fasta_url = "https://example.com/genome.fna.gz"
fasta_filename = "GRCh38.bin"
fasta_build = "GRCh38"
gff_url = "https://example.com/mane.gff.gz"
summary_url = "https://example.com/summary.txt.gz"
summary_input = "summary.tsv.gz"
transcript_models_version = "1.5"
assembly_report_url = "https://example.com/report.txt"
assembly_report_input = "report.txt"
output_dir = "data/vareffect"
"#;
        let config: VareffectConfig = toml::from_str(toml).unwrap();
        assert_eq!(config.vareffect.fasta_build, "GRCh38");
    }
}
