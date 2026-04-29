//! Configuration discovery and TOML parsing for `vareffect-cli`.
//!
//! The config file (`vareffect_build.toml`) defines URLs, filenames, and
//! output directories for the `setup` subcommand. The schema is now
//! per-assembly: `[vareffect.grch38]` and `[vareffect.grch37]` are
//! parallel sub-tables under `[vareffect]`. Either or both may be
//! present; `vareffect setup --assembly <name>` operates on the matching
//! sub-table.
//!
//! The parser uses `#[serde(default)]` on optional sections and silently
//! ignores unknown keys, so a single config file can be shared with
//! other build tools that add their own top-level sections.

use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use serde::Deserialize;

use vareffect::Assembly;

/// Config filename auto-discovered by [`find_config`].
const CONFIG_FILENAME: &str = "vareffect_build.toml";

/// Top-level view of `vareffect_build.toml`.
#[derive(Debug, Deserialize)]
pub struct VareffectConfig {
    /// Default directory paths (optional `[paths]` table).
    #[serde(default)]
    pub paths: PathsConfig,
    /// Required `[vareffect]` section with shared output directory and
    /// per-assembly sub-tables.
    pub vareffect: VareffectSection,
    /// Optional `[mane]` section -- only the `input` field is read, to
    /// override the default GFF3 cache filename when a downstream pipeline
    /// needs to share a pre-fetched copy.
    mane: Option<ManeSection>,
}

impl VareffectConfig {
    /// GFF3 cache filename: from `[mane].input` if present, else the
    /// default filename used by the `setup` subcommand for GRCh38.
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

/// The `[vareffect]` TOML section. Holds the shared output directory and
/// up to two per-assembly sub-tables.
#[derive(Debug, Deserialize)]
pub struct VareffectSection {
    /// Directory for the final runtime files. Shared across both
    /// assemblies; per-assembly disambiguation comes from the
    /// per-section `*_filename` fields.
    pub output_dir: String,

    /// GRCh38 assembly entry. Optional so a config that only builds
    /// GRCh37 still parses.
    #[serde(default)]
    pub grch38: Option<AssemblyEntry>,

    /// GRCh37 assembly entry. Optional so a config that only builds
    /// GRCh38 still parses.
    #[serde(default)]
    pub grch37: Option<AssemblyEntry>,
}

impl VareffectSection {
    /// Look up the entry for a given assembly. Returns `None` if the
    /// caller asked for an assembly the config doesn't define.
    pub fn entry(&self, assembly: Assembly) -> Option<&AssemblyEntry> {
        match assembly {
            Assembly::GRCh38 => self.grch38.as_ref(),
            Assembly::GRCh37 => self.grch37.as_ref(),
        }
    }
}

/// Per-assembly `[vareffect.<build>]` sub-table.
///
/// Field names match the TOML keys verbatim. Only `cross_validation_url`
/// / `cross_validation_input` are optional — Stage A omits them on the
/// GRCh37 entry by design (no authoritative second source identified
/// yet; see plan §6).
#[derive(Debug, Deserialize)]
pub struct AssemblyEntry {
    // Reference FASTA
    /// NCBI FTP URL for the gzipped reference FASTA.
    pub fasta_url: String,
    /// Output filename for the flat binary genome (e.g. `"GRCh38.bin"`).
    pub fasta_filename: String,
    /// Assembly build name (e.g. `"GRCh38"`). Sanity-checked against
    /// `fasta_url` to prevent assembly mismatches.
    pub fasta_build: String,

    // Transcript models
    /// URL to the GFF3 release file (MANE for GRCh38, NCBI RefSeq for GRCh37).
    pub gff_url: String,
    /// URL to the cross-validation source (MANE summary TSV / future
    /// GRCh37 second source). `None` skips cross-validation with a build
    /// log warning.
    pub cross_validation_url: Option<String>,
    /// Cache filename for the cross-validation source in `raw_dir`.
    pub cross_validation_input: Option<String>,
    /// Output filename for the per-assembly transcript models bin.
    pub transcript_models_filename: String,
    /// File-stem (no extension) for the transcript models build, used
    /// by [`crate::common::serialize_and_finalize`] to derive the
    /// `.bin` / `.bin.sha256` / `.manifest.json` sibling filenames.
    pub transcript_models_basename: String,
    /// Upstream data version string (e.g. `"1.5"` for MANE,
    /// `"GRCh37.p13"` for the NCBI RefSeq GRCh37 release).
    pub transcript_models_version: String,
    /// Transcript source for the GFF3 builder. Selects which `tag=`
    /// values to admit and which tier each maps to.
    pub transcript_source: TranscriptSource,

    // Patch-contig aliases (NCBI assembly report -> patch_chrom_aliases.csv)
    /// URL to the NCBI assembly report.
    pub assembly_report_url: String,
    /// Cache filename for the assembly report in `raw_dir`.
    pub assembly_report_input: String,
    /// Output filename for the per-assembly patch_chrom_aliases CSV
    /// (e.g. `"patch_chrom_aliases_grch38.csv"`).
    pub patch_aliases_filename: String,
}

/// Selects which `tag=` values the GFF3 transcript builder admits.
///
/// The user-facing config exposes this as a string enum
/// (`transcript_source = "mane"` / `"refseq_select"`); the parser
/// dispatches to the corresponding admit-tag slice in
/// `crate::builders::transcript_models`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
pub enum TranscriptSource {
    /// `tag=MANE Select` and `tag=MANE Plus Clinical`. GRCh38-only.
    #[serde(rename = "mane")]
    Mane,
    /// `tag=RefSeq Select`. Used by NCBI's GRCh37.p13 GFF3 (and works
    /// on later GRCh38 releases too, but MANE is preferred there).
    #[serde(rename = "refseq_select")]
    RefSeqSelect,
}

/// Locate the config file for `vareffect-cli`.
///
/// Resolution order:
/// 1. Explicit `--config` flag
/// 2. `VAREFFECT_BUILD_CONFIG` environment variable
/// 3. Well-known paths: `vareffect-cli/vareffect_build.toml` (repo
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

    const FULL_CONFIG: &str = r#"
[paths]
raw_dir = "data/raw"

[vareffect]
output_dir = "data/vareffect"

[vareffect.grch38]
fasta_url = "https://example.com/GRCh38.fna.gz"
fasta_filename = "GRCh38.bin"
fasta_build = "GRCh38"
gff_url = "https://example.com/mane.gff.gz"
cross_validation_url = "https://example.com/mane.summary.tsv.gz"
cross_validation_input = "mane.summary.tsv.gz"
transcript_models_filename = "transcript_models_grch38.bin"
transcript_models_basename = "transcript_models_grch38"
transcript_models_version = "1.5"
transcript_source = "mane"
assembly_report_url = "https://example.com/grch38.report.txt"
assembly_report_input = "grch38.report.txt"
patch_aliases_filename = "patch_chrom_aliases_grch38.csv"

[vareffect.grch37]
fasta_url = "https://example.com/GRCh37.fna.gz"
fasta_filename = "GRCh37.bin"
fasta_build = "GRCh37"
gff_url = "https://example.com/grch37.gff.gz"
transcript_models_filename = "transcript_models_grch37.bin"
transcript_models_basename = "transcript_models_grch37"
transcript_models_version = "GRCh37.p13"
transcript_source = "refseq_select"
assembly_report_url = "https://example.com/grch37.report.txt"
assembly_report_input = "grch37.report.txt"
patch_aliases_filename = "patch_chrom_aliases_grch37.csv"
"#;

    #[test]
    fn parse_full_config_with_both_assemblies() {
        let config: VareffectConfig = toml::from_str(FULL_CONFIG).unwrap();
        let g38 = config.vareffect.grch38.as_ref().expect("grch38 entry");
        let g37 = config.vareffect.grch37.as_ref().expect("grch37 entry");
        assert_eq!(g38.fasta_build, "GRCh38");
        assert_eq!(g37.fasta_build, "GRCh37");
        assert_eq!(g38.transcript_source, TranscriptSource::Mane);
        assert_eq!(g37.transcript_source, TranscriptSource::RefSeqSelect);
        assert!(g37.cross_validation_url.is_none());
        assert_eq!(config.vareffect.output_dir, "data/vareffect");
    }

    #[test]
    fn parse_grch38_only_config() {
        // GRCh37 sub-table absent; build should still parse.
        let toml = r#"
[vareffect]
output_dir = "data/vareffect"

[vareffect.grch38]
fasta_url = "https://example.com/genome.fna.gz"
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
"#;
        let config: VareffectConfig = toml::from_str(toml).unwrap();
        assert!(config.vareffect.grch38.is_some());
        assert!(config.vareffect.grch37.is_none());
    }

    #[test]
    fn entry_lookup_returns_correct_sub_table() {
        let config: VareffectConfig = toml::from_str(FULL_CONFIG).unwrap();
        assert_eq!(
            config
                .vareffect
                .entry(Assembly::GRCh38)
                .unwrap()
                .fasta_build,
            "GRCh38"
        );
        assert_eq!(
            config
                .vareffect
                .entry(Assembly::GRCh37)
                .unwrap()
                .fasta_build,
            "GRCh37"
        );
    }

    #[test]
    fn parse_ignores_unknown_sections() {
        let toml = r#"
[constraint]
input = "gnomad.tsv"
version = "4.1"

[vareffect]
output_dir = "data/vareffect"

[vareffect.grch38]
fasta_url = "https://example.com/genome.fna.gz"
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
"#;
        let config: VareffectConfig = toml::from_str(toml).unwrap();
        assert_eq!(
            config.vareffect.grch38.as_ref().unwrap().fasta_build,
            "GRCh38"
        );
    }
}
