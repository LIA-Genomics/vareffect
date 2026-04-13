//! Shared build utilities for store construction.
//!
//! Every store build produces three sidecar files:
//! - `{name}.bin` -- MessagePack serialized data
//! - `{name}.bin.sha256` -- hex-encoded SHA-256 of the .bin file
//! - `{name}.manifest.json` -- build metadata (version, record count, etc.)

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use chrono::Utc;
use serde::Serialize;
use sha2::{Digest, Sha256};

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

/// Write the .bin file, compute SHA-256, and create sidecar files.
///
/// # Arguments
///
/// * `name` -- Store name used as the file stem (e.g. "transcript_models").
/// * `output_dir` -- Directory to write files into (created if absent).
/// * `data` -- Pre-serialized MessagePack bytes (`rmp_serde::to_vec_named()`).
/// * `record_count` -- Number of records for the manifest.
/// * `source_version` -- Upstream data version string.
/// * `source_file` -- Original input filename for provenance.
/// * `source_file_sha256` -- SHA-256 of the raw input file. Pass `None` for
///   individual builder invocations where diffing is not needed.
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
    record_count: u64,
    source_version: &str,
    source_file: &str,
    source_file_sha256: Option<&str>,
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
        "source_version": source_version,
        "source_file": source_file,
        "record_count": record_count,
        "bin_size_bytes": data.len(),
        "sha256": sha256,
        "built_at": Utc::now().to_rfc3339(),
        "built_by": format!("vareffect-cli {}", env!("CARGO_PKG_VERSION")),
    });
    if let Some(src_sha) = source_file_sha256 {
        manifest["source_file_sha256"] = serde_json::Value::String(src_sha.to_string());
    }
    let manifest_path = output_dir.join(format!("{name}.manifest.json"));
    std::fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)
        .with_context(|| format!("writing {}", manifest_path.display()))?;

    Ok(BuildOutput {
        store_path: bin_path,
        sha256,
        record_count,
        source_version: source_version.to_string(),
        source_file_sha256: source_file_sha256.map(String::from),
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

/// Serialize a store to MessagePack and write sidecar files.
///
/// Combines `rmp_serde::to_vec_named`, source filename extraction, and
/// [`finalize_store`] into a single call.
///
/// # Arguments
///
/// * `name` -- Store name used as the file stem (e.g. `"transcript_models"`).
/// * `store` -- Any `Serialize`-able store struct.
/// * `record_count` -- Number of records for the manifest.
/// * `input` -- Original source file path (filename extracted for provenance).
/// * `output_dir` -- Directory to write files into.
/// * `version` -- Upstream data version string.
///
/// # Returns
///
/// Tuple of ([`BuildOutput`], serialized byte count).
///
/// # Errors
///
/// Returns an error if serialization or file I/O fails.
pub(crate) fn serialize_and_finalize(
    name: &str,
    store: &impl Serialize,
    record_count: u64,
    input: &Path,
    output_dir: &Path,
    version: &str,
) -> Result<(BuildOutput, usize)> {
    let data =
        rmp_serde::to_vec_named(store).with_context(|| format!("serializing {name} store"))?;
    let data_len = data.len();

    let source_file = input
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_default();

    let output = finalize_store(
        name,
        output_dir,
        &data,
        record_count,
        version,
        &source_file,
        None,
    )?;
    Ok((output, data_len))
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
