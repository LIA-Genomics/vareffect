//! Configuration discovery and TOML parsing for `vareffect-cli`.
//!
//! The config file (`vareffect_build.toml`) defines URLs, filenames, and
//! output directories for the `setup` subcommand, and the transcript source
//! that `setup` builds and `check` validates. The parser uses
//! `#[serde(default)]` on optional sections and silently ignores unknown
//! top-level sections, so a single config file can be shared with other
//! build tools that add their own. Inside `[vareffect]` and its tables an
//! unknown key is an error, so a misspelled `[vareffect.transcripts]` cannot
//! silently fall back to MANE.

use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow, bail};
use serde::Deserialize;
use vareffect::TranscriptSource;

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

    /// The transcript source selected by `[vareffect.transcripts] source`
    /// (MANE when the table is absent).
    pub fn source(&self) -> SourceName {
        self.vareffect.transcripts.source
    }

    /// Release of the selected source: `transcript_models_version` for MANE,
    /// the source table's `release` for RefSeq and GENCODE.
    ///
    /// # Errors
    ///
    /// Returns an error if the selected source is RefSeq or GENCODE and its
    /// `[vareffect.refseq]` / `[vareffect.gencode]` table is missing.
    pub fn source_release(&self) -> Result<&str> {
        match self.source() {
            SourceName::Mane => Ok(&self.vareffect.transcript_models_version),
            SourceName::Refseq | SourceName::Gencode => Ok(&self.pinned_source()?.release),
        }
    }

    /// Pinned download of the selected RefSeq or GENCODE source.
    ///
    /// # Errors
    ///
    /// Returns an error if the selected source is MANE (configured by the
    /// `[vareffect]` `gff_url` and `summary_url` keys instead) or if its
    /// `[vareffect.refseq]` / `[vareffect.gencode]` table is missing.
    pub fn pinned_source(&self) -> Result<&PinnedSource> {
        let (table, pinned) = match self.source() {
            SourceName::Mane => bail!(
                "the MANE source is configured by [vareffect] gff_url and summary_url, \
                 not a pinned source table"
            ),
            SourceName::Refseq => ("refseq", self.vareffect.refseq.as_ref()),
            SourceName::Gencode => ("gencode", self.vareffect.gencode.as_ref()),
        };
        pinned.ok_or_else(|| {
            anyhow!(
                "[vareffect.transcripts] source = \"{table}\" needs a [vareffect.{table}] \
                 table with url, release and sha256"
            )
        })
    }
}

/// Transcript source, as named in `[vareffect.transcripts] source` and by
/// `build-transcripts --source`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum SourceName {
    /// MANE Select and MANE Plus Clinical (NCBI MANE GFF3).
    #[default]
    Mane,
    /// Curated RefSeq `NM_` / `NR_` on the primary assembly.
    Refseq,
    /// GENCODE comprehensive.
    Gencode,
}

impl From<SourceName> for TranscriptSource {
    fn from(source: SourceName) -> Self {
        match source {
            SourceName::Mane => Self::Mane,
            SourceName::Refseq => Self::RefSeq,
            SourceName::Gencode => Self::Gencode,
        }
    }
}

/// Optional `[vareffect.transcripts]` table.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptsSection {
    /// Source that `setup` builds and `check` validates; `mane` when absent.
    /// An unknown value is a parse error.
    #[serde(default)]
    pub source: SourceName,
}

/// A pinned source GFF3: the `[vareffect.refseq]` or `[vareffect.gencode]`
/// table. `setup` refuses a file whose SHA-256 or header release differs.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PinnedSource {
    /// Download URL of the gzipped GFF3.
    pub url: String,
    /// Source release, recorded in the store and checked against the file
    /// header (e.g. `"RS_2025_08"` for RefSeq, `"50"` for GENCODE).
    pub release: String,
    /// Expected hex SHA-256 of the downloaded file.
    pub sha256: String,
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
/// Field names match the TOML keys verbatim. The URL and filename fields
/// are required -- there is no meaningful default for a URL or filename, and
/// defaulting would silently mask a misconfigured config. The transcript
/// source tables are optional, so a config written before they existed
/// still parses and selects MANE. Unknown keys are refused.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
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

    // Transcript source
    /// Optional `[vareffect.transcripts]` table selecting the source.
    #[serde(default)]
    pub transcripts: TranscriptsSection,
    /// Optional `[vareffect.refseq]` table; required when the source is
    /// `refseq`.
    pub refseq: Option<PinnedSource>,
    /// Optional `[vareffect.gencode]` table; required when the source is
    /// `gencode`.
    pub gencode: Option<PinnedSource>,
}

/// Locate the config file for `vareffect-cli`.
///
/// Resolution order:
/// 1. Explicit `--config` flag
/// 2. `VAREFFECT_BUILD_CONFIG` environment variable
/// 3. Well-known paths: `vareffect-cli/vareffect_build.toml` (repo root),
///    `/etc/vareffect/vareffect_build.toml` (system),
///    `vareffect/vareffect_build.toml` in the user config directory
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

    /// The required `[vareffect]` keys, followed by `extra` TOML.
    fn config_with(extra: &str) -> String {
        format!(
            r#"
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
{extra}"#
        )
    }

    const PINNED_TABLES: &str = r#"
[vareffect.refseq]
url = "https://example.com/refseq.gff.gz"
release = "RS_2025_08"
sha256 = "aa"

[vareffect.gencode]
url = "https://example.com/gencode.gff3.gz"
release = "50"
sha256 = "bb"
"#;

    #[test]
    fn source_defaults_to_mane_with_its_release() {
        let config: VareffectConfig = toml::from_str(&config_with("")).unwrap();
        assert_eq!(config.source(), SourceName::Mane);
        assert_eq!(config.source_release().unwrap(), "1.5");
        assert!(config.pinned_source().is_err());
    }

    #[test]
    fn source_selects_the_pinned_table() {
        for (name, source, release) in [
            ("refseq", SourceName::Refseq, "RS_2025_08"),
            ("gencode", SourceName::Gencode, "50"),
        ] {
            let text = config_with(&format!(
                "\n[vareffect.transcripts]\nsource = \"{name}\"\n{PINNED_TABLES}"
            ));
            let config: VareffectConfig = toml::from_str(&text).unwrap();
            assert_eq!(config.source(), source);
            assert_eq!(config.source_release().unwrap(), release);
            assert_eq!(config.pinned_source().unwrap().release, release);
        }
    }

    #[test]
    fn source_without_its_table_is_an_error() {
        let text = config_with("\n[vareffect.transcripts]\nsource = \"gencode\"\n");
        let config: VareffectConfig = toml::from_str(&text).unwrap();
        let error = config.source_release().expect_err("missing table");
        assert!(error.to_string().contains("[vareffect.gencode]"), "{error}");
    }

    #[test]
    fn unknown_source_or_key_is_a_parse_error() {
        let unknown_source = config_with("\n[vareffect.transcripts]\nsource = \"ensembl\"\n");
        assert!(toml::from_str::<VareffectConfig>(&unknown_source).is_err());
        let unknown_key = config_with(
            "\n[vareffect.refseq]\nurl = \"u\"\nrelease = \"r\"\nsha256 = \"s\"\nsha = \"s\"\n",
        );
        assert!(toml::from_str::<VareffectConfig>(&unknown_key).is_err());
        let misspelled_table = config_with("\n[vareffect.transcript]\nsource = \"gencode\"\n");
        assert!(toml::from_str::<VareffectConfig>(&misspelled_table).is_err());
    }

    #[test]
    fn source_name_converts_to_transcript_source() {
        assert_eq!(
            TranscriptSource::from(SourceName::Mane),
            TranscriptSource::Mane
        );
        assert_eq!(
            TranscriptSource::from(SourceName::Refseq),
            TranscriptSource::RefSeq
        );
        assert_eq!(
            TranscriptSource::from(SourceName::Gencode),
            TranscriptSource::Gencode
        );
    }
}
