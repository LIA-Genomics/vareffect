//! Reference genome downloader + flat binary builder for the `vareffect-cli
//! setup` command.
//!
//! Produces `GRCh38.bin` + `GRCh38.bin.idx` files in `data/vareffect/` from
//! the NCBI GRCh38 primary assembly FASTA.
//!
//! # Flow
//!
//! 1. If `{output_dir}/{filename}` and its `.idx` sidecar already exist,
//!    return [`GenomeOutcome::AlreadyPresent`].
//! 2. Otherwise, download the gzipped FASTA to `{raw_dir}/...`.
//! 3. Stream-decompress to a transient `{raw_dir}/GRCh38.staging.fa`.
//! 4. Convert the staging FASTA to flat binary: parse lines, uppercase,
//!    validate `{A,C,G,T,N}`, write `.bin` + `.bin.idx`.
//! 5. Delete the staging file.
//!
//! # Disk requirement
//!
//! Step 3 creates a ~3.2 GB transient file alongside the cached download,
//! and step 4 writes the ~3.1 GB `.bin` file, so `raw_dir` + `output_dir`
//! must have roughly 7 GB free during setup.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::thread::sleep;
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use vareffect::fasta::{ContigEntry, GENOME_BIN_INDEX_VERSION, GenomeBinIndex};

/// Number of download attempts before giving up.
const DOWNLOAD_RETRY_ATTEMPTS: u32 = 3;

/// Initial backoff before the first retry. Doubles on every subsequent attempt.
const DOWNLOAD_RETRY_BASE_DELAY: Duration = Duration::from_millis(250);

/// Paths to the sidecar file that sits next to the flat binary genome.
#[derive(Debug, Clone)]
pub struct GenomeSidecars {
    /// Path to the MessagePack `.bin.idx` index.
    pub idx: PathBuf,
}

/// Result of an [`ensure_genome`] call.
#[derive(Debug)]
pub enum GenomeOutcome {
    /// `.bin` and `.bin.idx` were already on disk.
    AlreadyPresent {
        /// Size of the `.bin` file in bytes.
        size_bytes: u64,
        /// Paths to the sidecars for logging.
        sidecars: GenomeSidecars,
    },
    /// `.bin` and/or `.bin.idx` were missing and had to be rebuilt.
    Built {
        /// Size of the final `.bin` in bytes.
        size_bytes: u64,
        /// Number of contigs in the freshly-built index.
        contigs: usize,
        /// Wall-clock duration for the full pipeline.
        elapsed: Duration,
        /// Paths to the freshly-written sidecars.
        sidecars: GenomeSidecars,
    },
}

/// Ensure the flat binary genome and its `.idx` sidecar exist in `output_dir`.
///
/// Idempotent: the first call does the work, subsequent calls return
/// [`GenomeOutcome::AlreadyPresent`].
///
/// # Arguments
///
/// * `client` -- Shared blocking HTTP client.
/// * `url` -- URL of the gzipped primary-assembly FASTA.
/// * `output_dir` -- Directory for the final `.bin` and `.bin.idx` files.
/// * `filename` -- Output filename for the flat binary (e.g. `"GRCh38.bin"`).
/// * `build` -- Assembly build label (e.g. `"GRCh38"`).
/// * `raw_dir` -- Directory for the cached `.gz` download and staging file.
pub fn ensure_genome(
    client: &reqwest::blocking::Client,
    url: &str,
    output_dir: &Path,
    filename: &str,
    build: &str,
    raw_dir: &Path,
) -> Result<GenomeOutcome> {
    let start = Instant::now();

    let final_bin = output_dir.join(filename);
    let sidecars = GenomeSidecars {
        idx: append_extension(&final_bin, "idx"),
    };

    // Fast path: both files already present.
    if final_bin.exists() && sidecars.idx.exists() {
        let size_bytes = std::fs::metadata(&final_bin)
            .with_context(|| format!("stat {}", final_bin.display()))?
            .len();
        return Ok(GenomeOutcome::AlreadyPresent {
            size_bytes,
            sidecars,
        });
    }

    std::fs::create_dir_all(output_dir)
        .with_context(|| format!("creating {}", output_dir.display()))?;
    std::fs::create_dir_all(raw_dir).with_context(|| format!("creating {}", raw_dir.display()))?;

    // Download the source .gz if not cached.
    let gz_basename = url
        .rsplit('/')
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or("reference.fa.gz");
    let gz_cache = raw_dir.join(gz_basename);

    if !gz_cache.exists() {
        download_to_path(client, url, &gz_cache, "reference FASTA")
            .with_context(|| format!("downloading {url}"))?;
    } else {
        eprintln!(
            "  {}  reference FASTA .gz cached at {}",
            console::style("SKIP").dim(),
            gz_cache.display(),
        );
    }

    let contigs = prepare_flat_binary(&gz_cache, &final_bin, &sidecars, build, raw_dir)
        .with_context(|| format!("preparing flat binary genome at {}", final_bin.display()))?;

    let size_bytes = std::fs::metadata(&final_bin)
        .with_context(|| format!("stat {}", final_bin.display()))?
        .len();

    Ok(GenomeOutcome::Built {
        size_bytes,
        contigs,
        elapsed: start.elapsed(),
        sidecars,
    })
}

/// Stream-download a URL to a local path with a progress bar.
///
/// Retries up to [`DOWNLOAD_RETRY_ATTEMPTS`] times with exponential backoff
/// on transport or truncation errors, and verifies the byte count matches
/// the server-advertised `Content-Length` before the final rename.
pub(crate) fn download_to_path(
    client: &reqwest::blocking::Client,
    url: &str,
    dest: &Path,
    display_name: &str,
) -> Result<u64> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating {}", parent.display()))?;
    }

    let tmp_path = dest.with_extension("tmp");

    let mut last_err: Option<anyhow::Error> = None;
    for attempt in 0..DOWNLOAD_RETRY_ATTEMPTS {
        match try_download_once(client, url, &tmp_path, display_name) {
            Ok(total) => {
                std::fs::rename(&tmp_path, dest).with_context(|| {
                    format!("renaming {} -> {}", tmp_path.display(), dest.display())
                })?;
                return Ok(total);
            }
            Err(err) => {
                let _ = std::fs::remove_file(&tmp_path);
                tracing::warn!(
                    attempt = attempt + 1,
                    max_attempts = DOWNLOAD_RETRY_ATTEMPTS,
                    download = %display_name,
                    "download attempt failed: {err:#}",
                );
                last_err = Some(err);
                if attempt + 1 < DOWNLOAD_RETRY_ATTEMPTS {
                    let backoff = DOWNLOAD_RETRY_BASE_DELAY * (1 << attempt);
                    sleep(backoff);
                }
            }
        }
    }
    Err(last_err
        .unwrap_or_else(|| anyhow::anyhow!("download failed for unknown reasons"))
        .context(format!(
            "{display_name}: exhausted {DOWNLOAD_RETRY_ATTEMPTS} attempts"
        )))
}

/// Single-attempt stream download with a progress bar and content-length
/// verification.
fn try_download_once(
    client: &reqwest::blocking::Client,
    url: &str,
    tmp_path: &Path,
    display_name: &str,
) -> Result<u64> {
    let resp = client
        .get(url)
        .send()
        .with_context(|| format!("HTTP GET {url}"))?;
    if !resp.status().is_success() {
        bail!("HTTP {} for {url}", resp.status());
    }

    let content_length = resp.content_length();

    let pb = if let Some(len) = content_length {
        let pb = ProgressBar::new(len);
        pb.set_style(
            ProgressStyle::with_template(
                "  DL    {msg} [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({percent}%)",
            )
            .expect("valid template")
            .progress_chars("=> "),
        );
        pb
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("  DL    {msg} {bytes} {elapsed}")
                .expect("valid template"),
        );
        pb
    };
    pb.set_message(display_name.to_string());

    let mut reader = resp;
    let mut file =
        File::create(tmp_path).with_context(|| format!("creating {}", tmp_path.display()))?;
    let mut buf = vec![0u8; 64 * 1024];
    let mut total: u64 = 0;

    loop {
        let n = reader
            .read(&mut buf)
            .with_context(|| format!("reading HTTP body for {display_name}"))?;
        if n == 0 {
            break;
        }
        file.write_all(&buf[..n])
            .with_context(|| format!("writing {}", tmp_path.display()))?;
        total += n as u64;
        pb.set_position(total);
    }

    file.sync_all()
        .with_context(|| format!("sync {}", tmp_path.display()))?;
    drop(file);

    pb.finish_and_clear();

    if let Some(expected) = content_length
        && total != expected
    {
        bail!("truncated download for {display_name}: expected {expected} bytes, wrote {total}");
    }

    Ok(total)
}

/// Full pipeline: gzip source -> staging .fa -> flat binary .bin + .bin.idx.
///
/// Returns the number of contigs in the freshly-built index.
fn prepare_flat_binary(
    source_gz: &Path,
    final_bin: &Path,
    sidecars: &GenomeSidecars,
    build: &str,
    raw_dir: &Path,
) -> Result<usize> {
    let staging_fa = raw_dir.join("GRCh38.staging.fa");

    eprintln!(
        "  {}  decompressing {} -> staging {}",
        console::style("GUNZIP").cyan(),
        source_gz.display(),
        staging_fa.display(),
    );
    decompress_to_staging(source_gz, &staging_fa)
        .with_context(|| format!("decompressing {}", source_gz.display()))?;

    eprintln!(
        "  {}  converting staging FASTA -> flat binary {}",
        console::style("BIN").cyan(),
        final_bin.display(),
    );
    let contigs = convert_fasta_to_flat_binary(&staging_fa, final_bin, &sidecars.idx, build)
        .with_context(|| {
            format!(
                "converting {} -> {}",
                staging_fa.display(),
                final_bin.display()
            )
        })?;

    // Clean up the ~3.2 GB transient staging file.
    if let Err(err) = std::fs::remove_file(&staging_fa) {
        tracing::warn!(
            path = %staging_fa.display(),
            error = %err,
            "failed to delete staging FASTA; you can remove it manually",
        );
    }

    Ok(contigs)
}

/// Stream-decompress a gzipped file into a plain-text destination using
/// the `.tmp` -> rename atomic pattern.
fn decompress_to_staging(gz_path: &Path, dest: &Path) -> Result<()> {
    let tmp_path = dest.with_extension("tmp");

    let gz_file = File::open(gz_path).with_context(|| format!("opening {}", gz_path.display()))?;
    let mut decoder = MultiGzDecoder::new(gz_file);
    let mut out =
        File::create(&tmp_path).with_context(|| format!("creating {}", tmp_path.display()))?;

    let mut buf = vec![0u8; 64 * 1024];
    loop {
        let n = decoder
            .read(&mut buf)
            .with_context(|| format!("reading {}", gz_path.display()))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])
            .with_context(|| format!("writing {}", tmp_path.display()))?;
    }

    out.sync_all()
        .with_context(|| format!("sync {}", tmp_path.display()))?;
    drop(out);

    std::fs::rename(&tmp_path, dest)
        .with_context(|| format!("renaming {} -> {}", tmp_path.display(), dest.display()))?;
    Ok(())
}

/// Convert a plain FASTA file to a flat binary genome (.bin + .bin.idx).
///
/// Single streaming pass: reads lines, detects `>` headers, uppercases and
/// validates sequence bytes (standard IUPAC nucleotide alphabet), writes
/// contiguous bytes to the .bin file, then serializes the MessagePack index.
///
/// Returns the number of contigs written.
fn convert_fasta_to_flat_binary(
    fasta_path: &Path,
    bin_path: &Path,
    idx_path: &Path,
    build: &str,
) -> Result<usize> {
    let fasta_file =
        File::open(fasta_path).with_context(|| format!("opening {}", fasta_path.display()))?;
    let reader = BufReader::new(fasta_file);

    let tmp_bin = bin_path.with_extension("tmp");
    let bin_file =
        File::create(&tmp_bin).with_context(|| format!("creating {}", tmp_bin.display()))?;
    let mut writer = BufWriter::new(bin_file);

    let mut entries: Vec<ContigEntry> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_length: u64 = 0;
    let mut total_bytes: u64 = 0;

    for line_result in reader.lines() {
        let line = line_result.with_context(|| format!("reading {}", fasta_path.display()))?;

        if let Some(name) = line.strip_prefix('>') {
            // Finalize previous contig.
            if let Some(prev_name) = current_name.take() {
                entries.push(ContigEntry {
                    name: prev_name,
                    offset: current_offset,
                    length: current_length,
                });
            }
            // Start new contig. The header line may contain a description
            // after the first whitespace — take only the name part.
            let contig_name = name.split_whitespace().next().unwrap_or(name).to_string();
            current_name = Some(contig_name);
            current_offset = total_bytes;
            current_length = 0;
        } else {
            // Sequence line: uppercase, validate, write as a single chunk.
            let line_bytes = line.into_bytes();
            let mut buf = line_bytes;
            buf.make_ascii_uppercase();
            for (i, &b) in buf.iter().enumerate() {
                if !vareffect::fasta::is_iupac_nucleotide(b) {
                    bail!(
                        "non-IUPAC byte 0x{:02X} ('{}') in {} at byte offset {}",
                        b,
                        b as char,
                        fasta_path.display(),
                        total_bytes + i as u64,
                    );
                }
            }
            let n = buf.len() as u64;
            writer
                .write_all(&buf)
                .with_context(|| format!("writing {}", tmp_bin.display()))?;
            total_bytes += n;
            current_length += n;
        }
    }

    // Finalize the last contig.
    if let Some(last_name) = current_name.take() {
        entries.push(ContigEntry {
            name: last_name,
            offset: current_offset,
            length: current_length,
        });
    }

    let contigs = entries.len();

    writer
        .flush()
        .with_context(|| format!("flushing {}", tmp_bin.display()))?;
    writer
        .into_inner()
        .map_err(|e| anyhow::anyhow!("flush error: {}", e.into_error()))?
        .sync_all()
        .with_context(|| format!("sync {}", tmp_bin.display()))?;

    // Atomic rename.
    std::fs::rename(&tmp_bin, bin_path)
        .with_context(|| format!("renaming {} -> {}", tmp_bin.display(), bin_path.display()))?;

    // Write the MessagePack index.
    let index = GenomeBinIndex {
        version: GENOME_BIN_INDEX_VERSION,
        build: build.to_string(),
        expected_size: total_bytes,
        contigs: entries,
    };
    let idx_bytes = rmp_serde::to_vec(&index).context("serializing genome binary index")?;

    let tmp_idx = idx_path.with_extension("tmp");
    std::fs::write(&tmp_idx, &idx_bytes)
        .with_context(|| format!("writing {}", tmp_idx.display()))?;
    std::fs::rename(&tmp_idx, idx_path)
        .with_context(|| format!("renaming {} -> {}", tmp_idx.display(), idx_path.display()))?;

    Ok(contigs)
}

/// Append a dot + extension to a path without replacing existing extensions.
fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut os = path.as_os_str().to_os_string();
    os.push(".");
    os.push(ext);
    PathBuf::from(os)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vareffect::FastaReader;
    use vareffect::fasta::GenomeBinIndex;

    /// `append_extension` must *append*, not replace. This is the property
    /// the rest of the builder relies on when deriving `.bin.idx` from the
    /// main binary path.
    #[test]
    fn append_extension_appends_suffix() {
        assert_eq!(
            append_extension(Path::new("/tmp/GRCh38.bin"), "idx"),
            PathBuf::from("/tmp/GRCh38.bin.idx"),
        );
        assert_eq!(
            append_extension(Path::new("/tmp/ref"), "idx"),
            PathBuf::from("/tmp/ref.idx"),
        );
    }

    /// Write a tiny synthetic FASTA, convert it via `convert_fasta_to_flat_binary`,
    /// and verify the result is readable by the runtime `FastaReader`.
    #[test]
    fn convert_fasta_to_flat_binary_round_trip() {
        let tmp = tempfile::tempdir().expect("tempdir");

        // Synthetic multi-contig FASTA (UCSC-style names, 80-col wrap).
        let fasta_path = tmp.path().join("test.fa");
        let fasta = b">chr1 Homo sapiens chromosome 1\nACGTACGT\nNNNNNNNN\n>chrM\nGGGGCCCC\n";
        std::fs::write(&fasta_path, fasta).unwrap();

        let bin_path = tmp.path().join("test.bin");
        let idx_path = tmp.path().join("test.bin.idx");

        let contigs =
            convert_fasta_to_flat_binary(&fasta_path, &bin_path, &idx_path, "test").unwrap();

        assert_eq!(contigs, 2, "expected two contigs");
        assert!(bin_path.exists());
        assert!(idx_path.exists());

        // Open via the runtime reader and verify sequence content.
        let reader =
            FastaReader::open_with_assembly(&bin_path, vareffect::Assembly::GRCh38).unwrap();
        assert_eq!(reader.fetch_sequence("chr1", 0, 4).unwrap(), b"ACGT");
        assert_eq!(reader.fetch_sequence("chr1", 8, 16).unwrap(), b"NNNNNNNN");
        assert_eq!(reader.chrom_length("chr1"), Some(16));
        assert_eq!(reader.fetch_sequence("chrM", 0, 4).unwrap(), b"GGGG");
        assert_eq!(reader.chrom_length("chrM"), Some(8));
    }

    /// A FASTA with a non-IUPAC byte (e.g. `X`) must be rejected at convert
    /// time. IUPAC ambiguity codes (R, Y, M, etc.) are accepted since the
    /// NCBI GRCh38.p14 assembly contains them.
    #[test]
    fn convert_fasta_to_flat_binary_rejects_non_iupac() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let fasta_path = tmp.path().join("bad.fa");
        std::fs::write(&fasta_path, b">chr1\nACGTXN\n").unwrap();

        let bin_path = tmp.path().join("bad.bin");
        let idx_path = tmp.path().join("bad.bin.idx");

        let err = convert_fasta_to_flat_binary(&fasta_path, &bin_path, &idx_path, "test");
        assert!(err.is_err(), "non-IUPAC byte must cause an error");
        let msg = format!("{:#}", err.unwrap_err());
        assert!(msg.contains("non-IUPAC"), "error message: {msg}");
    }

    /// IUPAC ambiguity codes (R, Y, M, etc.) must be accepted since the
    /// NCBI GRCh38.p14 assembly uses them in some patch-scaffold regions.
    #[test]
    fn convert_fasta_to_flat_binary_accepts_iupac_ambiguity() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let fasta_path = tmp.path().join("iupac.fa");
        std::fs::write(&fasta_path, b">chr1\nACGTNRYSWKMBDHV\n").unwrap();

        let bin_path = tmp.path().join("iupac.bin");
        let idx_path = tmp.path().join("iupac.bin.idx");

        let contigs =
            convert_fasta_to_flat_binary(&fasta_path, &bin_path, &idx_path, "test").unwrap();
        assert_eq!(contigs, 1);

        let reader =
            FastaReader::open_with_assembly(&bin_path, vareffect::Assembly::GRCh38).unwrap();
        assert_eq!(
            reader.fetch_sequence("chr1", 0, 15).unwrap(),
            b"ACGTNRYSWKMBDHV",
        );
    }

    /// The index written by the CLI builder must use the same format version
    /// constant as the runtime reader.
    #[test]
    fn convert_fasta_to_flat_binary_uses_shared_version() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let fasta_path = tmp.path().join("ver.fa");
        std::fs::write(&fasta_path, b">chr1\nACGT\n").unwrap();

        let bin_path = tmp.path().join("ver.bin");
        let idx_path = tmp.path().join("ver.bin.idx");

        convert_fasta_to_flat_binary(&fasta_path, &bin_path, &idx_path, "test").unwrap();

        let idx_bytes = std::fs::read(&idx_path).unwrap();
        let index: GenomeBinIndex = rmp_serde::from_slice(&idx_bytes).unwrap();
        assert_eq!(
            index.version, GENOME_BIN_INDEX_VERSION,
            "CLI builder must use the shared GENOME_BIN_INDEX_VERSION constant",
        );
    }
}
