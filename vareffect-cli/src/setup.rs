//! Orchestrator for the `vareffect-cli setup` subcommand.
//!
//! Responsible for provisioning every runtime input the `vareffect` crate
//! needs in a single idempotent invocation:
//!
//! 1. Download + decompress + convert the GRCh38 reference FASTA to flat binary.
//! 2. Download + parse the NCBI assembly report into patch-chrom aliases.
//! 3. Download + parse + validate the GFF3 of the configured transcript
//!    source (`[vareffect.transcripts] source`): MANE ->
//!    `transcript_models.bin`, RefSeq -> `transcript_models_refseq.bin`,
//!    GENCODE -> `transcript_models_gencode.bin`. RefSeq and GENCODE files
//!    are verified against their pinned SHA-256, release and genome build.
//!
//! Outputs live under the `output_dir` configured in `vareffect_build.toml`
//! (default `data/vareffect/`). Source archives (`.gz` FASTA, the source
//! GFF3, the MANE summary TSV) land in `raw_dir` (default `data/raw/`) as a
//! cache so reruns skip the download step.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use console::style;

use vareffect::TranscriptSource;

use crate::builders;
use crate::builders::reference_genome::{self, GenomeOutcome};
use crate::common::{self, format_size};
use crate::config::{self, SourceName};
use crate::source_file;

/// Run the `vareffect-cli setup` subcommand.
///
/// # Arguments
///
/// * `config_path` -- Explicit `--config` path, or `None` to auto-discover.
/// * `fasta_only` -- Skip the transcript-model build step.
/// * `models_only` -- Skip the reference FASTA step.
/// * `output_override` -- When `Some`, overrides `[vareffect].output_dir`
///   from the config file so runtime files land in the given directory.
///
/// `fasta_only` and `models_only` are mutually exclusive at the CLI layer
/// (clap enforces this via `conflicts_with`).
pub fn run(
    config_path: Option<&Path>,
    fasta_only: bool,
    models_only: bool,
    output_override: Option<&Path>,
) -> Result<()> {
    let total_start = Instant::now();

    // 1. Locate and parse the config.
    let config_file = config::find_config(config_path)?;
    let config_text = std::fs::read_to_string(&config_file)
        .with_context(|| format!("reading config {}", config_file.display()))?;
    let config: config::VareffectConfig = toml::from_str(&config_text)
        .with_context(|| format!("parsing [vareffect] section in {}", config_file.display()))?;

    // 2. Resolve directories. `output_dir` comes from the [vareffect]
    //    section (it's vareffect-specific). `raw_dir` comes from the shared
    //    [paths] block, falling back to `data/raw` if unset.
    let output_dir = output_override
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(&config.vareffect.output_dir));
    let raw_dir = config
        .paths
        .raw_dir
        .as_deref()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("data/raw"));

    std::fs::create_dir_all(&output_dir)
        .with_context(|| format!("creating {}", output_dir.display()))?;
    std::fs::create_dir_all(&raw_dir).with_context(|| format!("creating {}", raw_dir.display()))?;

    // Sanity-check that `fasta_url` matches the declared `fasta_build`.
    // A MANE / reference-FASTA assembly mismatch would silently corrupt
    // every coordinate lookup downstream; better to fail loudly here than
    // to debug an off-by-millions variant call months later.
    if !config
        .vareffect
        .fasta_url
        .contains(&config.vareffect.fasta_build)
    {
        bail!(
            "config mismatch: [vareffect].fasta_url={:?} does not contain the declared \
             fasta_build={:?} -- refusing to download a FASTA that may be the wrong assembly",
            config.vareffect.fasta_url,
            config.vareffect.fasta_build,
        );
    }

    eprintln!();
    eprintln!(
        "{}",
        style("Setting up vareffect runtime data").bold().cyan()
    );
    eprintln!("  output_dir  = {}", output_dir.display());
    eprintln!("  raw_dir     = {}", raw_dir.display());
    eprintln!("  fasta_build = {}", config.vareffect.fasta_build);
    if !fasta_only {
        eprintln!(
            "  transcripts = {} {}",
            TranscriptSource::from(config.source()),
            config.source_release()?,
        );
    }
    eprintln!();

    // 3. Build a single HTTP client for all downloads. 1800 s body timeout
    //    covers the ~3 GB gzipped reference FASTA on a slow link; 30 s
    //    connect keeps misconfigured hosts from hanging the build.
    let http_client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(1800))
        .connect_timeout(Duration::from_secs(30))
        .user_agent("vareffect-cli/0.1.0")
        .build()
        .context("building HTTP client")?;

    // -- Step 1/3: Flat binary reference genome -------------------------
    if !models_only {
        eprintln!(
            "{}",
            style("[1/3] Flat binary reference genome (NCBI)").bold()
        );
        let outcome = reference_genome::ensure_genome(
            &http_client,
            &config.vareffect.fasta_url,
            &output_dir,
            &config.vareffect.fasta_filename,
            &config.vareffect.fasta_build,
            &raw_dir,
        )
        .context("preparing reference genome binary")?;
        report_genome_outcome(&outcome, &output_dir.join(&config.vareffect.fasta_filename));

        // Warn (non-destructive) if a legacy BGZF FASTA is still present.
        let legacy_bgzf = output_dir.join("GRCh38.fa.gz");
        if legacy_bgzf.exists() {
            tracing::warn!(
                path = %legacy_bgzf.display(),
                "legacy BGZF FASTA detected; delete it (and .fai/.gzi sidecars) \
                 manually to complete the flat binary migration",
            );
        }
        let legacy_plain = output_dir.join("GRCh38.fa");
        if legacy_plain.exists() {
            tracing::warn!(
                path = %legacy_plain.display(),
                "legacy plain FASTA detected; delete it manually",
            );
        }
        eprintln!();
    } else {
        eprintln!(
            "{}",
            style("[1/3] Flat binary genome (skipped via --models-only)").dim()
        );
        eprintln!();
    }

    // -- Step 2/3: Patch-chrom aliases (feeds transcript cross-validation)
    let aliases_path = if !fasta_only {
        eprintln!("{}", style("[2/3] Patch-chrom aliases").bold());

        let assembly_report_path = raw_dir.join(&config.vareffect.assembly_report_input);
        ensure_cached(
            &http_client,
            &config.vareffect.assembly_report_url,
            &assembly_report_path,
            "NCBI assembly report",
        )?;

        let (csv_path, row_count) = builders::patch_chrom_aliases::build_from_assembly_report(
            &assembly_report_path,
            &output_dir,
        )
        .context("building patch_chrom_aliases.csv")?;

        eprintln!(
            "  {}  wrote {} ({} aliases)",
            style("OK").green(),
            style(
                csv_path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(builders::patch_chrom_aliases::OUTPUT_FILENAME)
            )
            .cyan(),
            row_count,
        );
        eprintln!();
        Some(csv_path)
    } else {
        eprintln!(
            "{}",
            style("[2/3] Patch-chrom aliases (skipped via --fasta-only)").dim()
        );
        eprintln!();
        None
    };

    // -- Step 3/3: Transcript models ------------------------------------
    if !fasta_only {
        eprintln!("{}", style("[3/3] Transcript models").bold());

        let source = config.source();
        let release = config.source_release()?;
        let (gff_path, summary_path) = match source {
            SourceName::Mane => {
                // GFF3 cache filename from the optional `[mane]` section, so
                // a downstream pipeline can pre-stage the ~80 MB MANE GFF3
                // and have the `setup` command reuse it instead of
                // re-downloading.
                let gff_path = raw_dir.join(config.gff_filename());
                ensure_cached(
                    &http_client,
                    &config.vareffect.gff_url,
                    &gff_path,
                    "MANE GFF3",
                )?;
                let summary_path = raw_dir.join(&config.vareffect.summary_input);
                ensure_cached(
                    &http_client,
                    &config.vareffect.summary_url,
                    &summary_path,
                    "MANE summary TSV",
                )?;
                (gff_path, Some(summary_path))
            }
            SourceName::Refseq | SourceName::Gencode => {
                let pinned = config.pinned_source()?;
                let fasta_build = &config.vareffect.fasta_build;
                let gff_path = raw_dir.join(source_file::cache_filename(
                    source,
                    fasta_build,
                    &pinned.release,
                )?);
                let display_name = format!("{} GFF3", TranscriptSource::from(source));
                ensure_cached(&http_client, &pinned.url, &gff_path, &display_name)?;
                // Cached files are verified too: a stale or swapped cache
                // must not build a store labelled with the pinned release.
                let header = source_file::verify(&gff_path, source, pinned, fasta_build)?;
                eprintln!(
                    "  {}  verified {} (SHA-256, release {}, {})",
                    style("OK").green(),
                    display_name,
                    header.release,
                    header.genome_build,
                );
                (gff_path, None)
            }
        };

        // Always rebuild the transcript model store: the parse is cheap
        // (seconds) and the cross-validation is cheap too, so an
        // always-fresh build is simpler than trying to decide "is this .bin
        // up to date?".
        let build_start = Instant::now();
        // `aliases_path` is `Some` here by construction -- Step 2/3 only
        // returns `None` under `fasta_only`, which short-circuits this block.
        let aliases_ref = aliases_path
            .as_deref()
            .expect("patch aliases must be built before transcript models in vareffect-cli setup");
        let transcript_source = TranscriptSource::from(source);
        let (out, size_bytes) =
            builders::transcript_models::build(&builders::transcript_models::BuildRequest {
                source: transcript_source,
                input: &gff_path,
                summary_input: summary_path.as_deref(),
                patch_aliases_input: Some(aliases_ref),
                output_dir: &output_dir,
                release,
            })
            .with_context(|| format!("building the {transcript_source} transcript store"))?;
        let store_name = format!("{}.bin", common::store_stem(transcript_source)?);

        eprintln!(
            "  {}  wrote {} ({} records, {} in {}s)",
            style("OK").green(),
            style(store_name).cyan(),
            out.record_count,
            format_size(size_bytes as u64),
            build_start.elapsed().as_secs(),
        );
        eprintln!();
    } else {
        eprintln!(
            "{}",
            style("[3/3] Transcript models (skipped via --fasta-only)").dim()
        );
        eprintln!();
    }

    eprintln!(
        "{}  vareffect-cli setup complete in {}s",
        style("OK").green().bold(),
        total_start.elapsed().as_secs(),
    );

    Ok(())
}

/// Download `url` to `dest` if `dest` does not already exist, then return
/// silently. Used for the GFF3 and summary TSV cache entries.
fn ensure_cached(
    client: &reqwest::blocking::Client,
    url: &str,
    dest: &Path,
    display_name: &str,
) -> Result<()> {
    if dest.exists() {
        eprintln!(
            "  {}  {} cached at {}",
            style("SKIP").dim(),
            display_name,
            dest.display(),
        );
        return Ok(());
    }
    let bytes = reference_genome::download_to_path(client, url, dest, display_name)?;
    eprintln!(
        "  {}  downloaded {} ({}) -> {}",
        style("OK").green(),
        display_name,
        format_size(bytes),
        dest.display(),
    );
    Ok(())
}

/// Pretty-print an [`GenomeOutcome`].
fn report_genome_outcome(outcome: &GenomeOutcome, bin_path: &Path) {
    match outcome {
        GenomeOutcome::AlreadyPresent {
            size_bytes,
            sidecars,
        } => {
            eprintln!(
                "  {}  {} already present ({})",
                style("SKIP").dim(),
                bin_path.display(),
                format_size(*size_bytes),
            );
            log_genome_sidecar_sizes(sidecars);
        }
        GenomeOutcome::Built {
            size_bytes,
            contigs,
            elapsed,
            sidecars,
        } => {
            eprintln!(
                "  {}  built {} ({}, {} contigs) in {}s",
                style("OK").green(),
                bin_path.display(),
                format_size(*size_bytes),
                contigs,
                elapsed.as_secs(),
            );
            log_genome_sidecar_sizes(sidecars);
        }
    }
}

/// Log `.bin.idx` size.
fn log_genome_sidecar_sizes(sidecars: &reference_genome::GenomeSidecars) {
    let size = std::fs::metadata(&sidecars.idx)
        .map(|m| m.len())
        .unwrap_or(0);
    eprintln!(
        "    {}  {} ({})",
        style("idx").dim(),
        sidecars.idx.display(),
        format_size(size),
    );
}
