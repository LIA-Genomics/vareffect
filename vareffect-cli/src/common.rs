//! Shared build utilities for store construction.
//!
//! Every store build produces three sidecar files:
//! - `{name}.bin` -- the store file (for the transcript store: a 12-byte
//!   versioned header followed by MessagePack, see
//!   [`vareffect::encode_transcript_store`])
//! - `{name}.bin.sha256` -- hex-encoded SHA-256 of the whole .bin file
//! - `{name}.manifest.json` -- build metadata (version, record count, etc.)

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use chrono::Utc;
use sha2::{Digest, Sha256};
use vareffect::transcript::{TRANSCRIPT_STORE_FORMAT_VERSION, TRANSCRIPT_STORE_HEADER_LEN};
use vareffect::{TranscriptModel, TranscriptSource};

/// Metadata returned after a successful store build.
#[derive(Debug)]
pub struct BuildOutput {
    /// Path to the written store file (`.bin`).
    #[allow(dead_code)]
    pub store_path: PathBuf,
    /// Hex-encoded SHA-256 of the store file contents.
    pub sha256: String,
    /// Number of records serialized.
    pub record_count: u64,
    /// Source data version string (included for provenance).
    #[allow(dead_code)]
    pub source_version: String,
    /// SHA-256 of the raw source file, if recorded during build.
    #[allow(dead_code)]
    pub source_file_sha256: Option<String>,
}

/// Provenance and format fields written into a store's manifest.
#[derive(Debug, Clone, Default)]
pub struct ManifestFields {
    /// Number of records in the store.
    pub record_count: u64,
    /// Upstream data version string (e.g. `"1.5"` for MANE v1.5).
    pub source_version: String,
    /// Original input filename.
    pub source_file: String,
    /// Hex SHA-256 of the raw input file, when known.
    pub source_file_sha256: Option<String>,
    /// Store file format version, for stores that carry a versioned header.
    pub format_version: Option<u32>,
    /// Hex SHA-256 of the payload after the header. It links builds whose
    /// files differ only by the header: across the move from unversioned
    /// stores to format version 1 the payload was untouched. Across later
    /// format versions the link key is (`format_version`, `payload_sha256`),
    /// since a version bump may change meaning without changing bytes.
    pub payload_sha256: Option<String>,
}

/// Write the .bin file, compute SHA-256, and create sidecar files.
///
/// # Arguments
///
/// * `name` -- Store name used as the file stem (e.g. "transcript_models").
/// * `output_dir` -- Directory to write files into (created if absent).
/// * `data` -- Complete store file contents.
/// * `fields` -- Manifest provenance and format fields.
///
/// # Returns
///
/// [`BuildOutput`] with the path, SHA-256, and record count.
///
/// # Errors
///
/// Returns an error if directory creation or file I/O fails.
pub fn finalize_store(
    name: &str,
    output_dir: &Path,
    data: &[u8],
    fields: &ManifestFields,
) -> Result<BuildOutput> {
    std::fs::create_dir_all(output_dir)
        .with_context(|| format!("creating output directory {}", output_dir.display()))?;

    // 1. Write .bin
    let bin_path = output_dir.join(format!("{name}.bin"));
    std::fs::write(&bin_path, data).with_context(|| format!("writing {}", bin_path.display()))?;

    // 2. Compute + write .sha256
    let sha256 = hex::encode(Sha256::digest(data));
    let sha_path = output_dir.join(format!("{name}.bin.sha256"));
    std::fs::write(&sha_path, &sha256)
        .with_context(|| format!("writing {}", sha_path.display()))?;

    // 3. Write manifest
    let mut manifest = serde_json::json!({
        "store_name": name,
        "source_version": fields.source_version,
        "source_file": fields.source_file,
        "record_count": fields.record_count,
        "bin_size_bytes": data.len(),
        "sha256": sha256,
        "built_at": Utc::now().to_rfc3339(),
        "built_by": format!("vareffect-cli {}", env!("CARGO_PKG_VERSION")),
    });
    if let Some(src_sha) = &fields.source_file_sha256 {
        manifest["source_file_sha256"] = serde_json::Value::String(src_sha.clone());
    }
    if let Some(format_version) = fields.format_version {
        manifest["format_version"] = serde_json::Value::from(format_version);
    }
    if let Some(payload_sha) = &fields.payload_sha256 {
        manifest["payload_sha256"] = serde_json::Value::String(payload_sha.clone());
    }
    let manifest_path = output_dir.join(format!("{name}.manifest.json"));
    std::fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)
        .with_context(|| format!("writing {}", manifest_path.display()))?;

    Ok(BuildOutput {
        store_path: bin_path,
        sha256,
        record_count: fields.record_count,
        source_version: fields.source_version.clone(),
        source_file_sha256: fields.source_file_sha256.clone(),
    })
}

/// Format a byte count as a human-readable size string (e.g. "4.8 MB").
///
/// Includes a GB tier for the ~3 GB decompressed reference FASTA produced
/// by the `setup` subcommand.
pub fn format_size(bytes: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = KB * 1024;
    const GB: u64 = MB * 1024;
    if bytes < KB {
        format!("{bytes} B")
    } else if bytes < MB {
        format!("{:.1} KB", bytes as f64 / KB as f64)
    } else if bytes < GB {
        format!("{:.1} MB", bytes as f64 / MB as f64)
    } else {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    }
}

/// Stream a file through SHA-256 without holding it in memory.
///
/// # Arguments
///
/// * `path` -- File to hash.
///
/// # Returns
///
/// The hex-encoded digest.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or read.
fn sha256_file(path: &Path) -> Result<String> {
    use std::io::Read;

    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening {} to checksum it", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let read = match file.read(&mut buffer) {
            Ok(read) => read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => {
                return Err(error).with_context(|| format!("checksumming {}", path.display()));
            }
        };
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(hex::encode(hasher.finalize()))
}

/// Encode transcript models as a versioned store and write its sidecar files.
///
/// Encodes through [`vareffect::encode_transcript_store`], so the builder and
/// the library's reader share one format implementation. The store payload
/// records `source` and `version`; the manifest records the format version,
/// the payload SHA-256 and the source file's SHA-256.
///
/// # Arguments
///
/// * `source` -- Annotation source; names the file stem (see [`store_stem`]).
/// * `transcripts` -- The models to store.
/// * `input` -- Original source file path (name and SHA-256 recorded for
///   provenance).
/// * `output_dir` -- Directory to write files into.
/// * `version` -- Upstream data version string, recorded as the store's
///   source release.
///
/// # Returns
///
/// Tuple of ([`BuildOutput`], store file size in bytes).
///
/// # Errors
///
/// Returns an error if `transcripts` is empty, encoding fails, the source has
/// no file stem, or the source file or output files cannot be read or
/// written.
pub(crate) fn write_transcript_store(
    source: TranscriptSource,
    transcripts: &[TranscriptModel],
    input: &Path,
    output_dir: &Path,
    version: &str,
) -> Result<(BuildOutput, usize)> {
    let name = store_stem(source)?;
    let data = vareffect::encode_transcript_store(Some(source), Some(version), transcripts)
        .with_context(|| format!("encoding {name} store"))?;
    let source_sha256 = sha256_file(input)?;
    let fields = ManifestFields {
        record_count: transcripts.len() as u64,
        source_version: version.to_string(),
        source_file: input
            .file_name()
            .map(|f| f.to_string_lossy().to_string())
            .unwrap_or_default(),
        source_file_sha256: Some(source_sha256),
        format_version: Some(TRANSCRIPT_STORE_FORMAT_VERSION),
        payload_sha256: Some(hex::encode(Sha256::digest(
            &data[TRANSCRIPT_STORE_HEADER_LEN..],
        ))),
    };
    let output = finalize_store(name, output_dir, &data, &fields)?;
    Ok((output, data.len()))
}

/// File stem of the transcript store built from `source`.
///
/// # Errors
///
/// Returns an error for a source this CLI has no store name for (the
/// library's `TranscriptSource` is `#[non_exhaustive]`).
pub(crate) fn store_stem(source: TranscriptSource) -> Result<&'static str> {
    match source {
        TranscriptSource::Mane => Ok("transcript_models"),
        TranscriptSource::RefSeq => Ok("transcript_models_refseq"),
        TranscriptSource::Gencode => Ok("transcript_models_gencode"),
        other => bail!("no store file name for transcript source {other}"),
    }
}

/// Print the summary table header for build output.
pub fn print_build_header() {
    eprintln!(
        "  {:<24} {:>8}    {:>10}    {:>15}    {}",
        console::style("Store").bold(),
        console::style("Records").bold(),
        console::style("Size").bold(),
        console::style("SHA256").bold(),
        console::style("Time").bold(),
    );
    eprintln!("  {}", console::style("-".repeat(78)).dim());
}

/// Print a build summary line for one store.
pub fn print_build_summary(name: &str, output: &BuildOutput, data_size: usize, elapsed: Instant) {
    let secs = elapsed.elapsed().as_secs();
    let time_str = if secs < 1 {
        "<1s".to_string()
    } else {
        format!("{secs}s")
    };
    eprintln!(
        "  {:<24} {:>8}    {:>10}    {}...    {}",
        console::style(name).cyan(),
        output.record_count,
        format_size(data_size as u64),
        console::style(&output.sha256[..12]).dim(),
        console::style(time_str).dim(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn store_stems_are_distinct_per_source() {
        assert_eq!(
            store_stem(TranscriptSource::Mane).expect("mane"),
            "transcript_models"
        );
        assert_eq!(
            store_stem(TranscriptSource::RefSeq).expect("refseq"),
            "transcript_models_refseq"
        );
        assert_eq!(
            store_stem(TranscriptSource::Gencode).expect("gencode"),
            "transcript_models_gencode"
        );
    }
}
