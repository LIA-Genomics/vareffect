//! Orchestrator for the `vareffect-cli setup` subcommand.
//!
//! Responsible for provisioning every runtime input the `vareffect` crate
//! needs, in a single idempotent invocation, for one or both supported
//! assemblies.
//!
//! Per assembly:
//! 1. Download + decompress + convert the reference FASTA to flat binary.
//! 2. Download + parse the NCBI assembly report into per-assembly
//!    patch-chrom aliases.
//! 3. Download + parse + (when supported) cross-validate the GFF3 →
//!    `transcript_models_<assembly>.bin`.
//!
//! Output artifacts live under the `output_dir` configured in
//! `vareffect_build.toml` (default `data/vareffect/`). Source archives
//! land in `raw_dir` (default `data/raw/`) as a cache so reruns skip the
//! download step.

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use console::style;

use crate::builders;
use crate::builders::reference_genome::{self, GenomeOutcome};
use crate::common::format_size;
use crate::config::{self, AssemblyEntry, TranscriptSource};
use vareffect::Assembly;

/// Selection passed in by the CLI: build a single assembly, or every
/// assembly defined in the config file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyFilter {
    /// Build the named assembly only. Errors out if the config does not
    /// define a sub-table for it.
    Single(Assembly),
    /// Build every assembly defined in the config (`grch38`, `grch37`,
    /// or both).
    All,
}

/// Run the `vareffect-cli setup` subcommand.
///
/// `fasta_only` and `models_only` are mutually exclusive at the CLI
/// layer (clap enforces this via `conflicts_with`).
pub fn run(
    config_path: Option<&Path>,
    fasta_only: bool,
    models_only: bool,
    output_override: Option<&Path>,
    assemblies: AssemblyFilter,
) -> Result<()> {
    let total_start = Instant::now();

    let config_file = config::find_config(config_path)?;
    let config_text = std::fs::read_to_string(&config_file)
        .with_context(|| format!("reading config {}", config_file.display()))?;
    let config: config::VareffectConfig = toml::from_str(&config_text)
        .with_context(|| format!("parsing {}", config_file.display()))?;

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

    let http_client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(1800))
        .connect_timeout(Duration::from_secs(30))
        .user_agent(concat!("vareffect-cli/", env!("CARGO_PKG_VERSION")))
        .build()
        .context("building HTTP client")?;

    // Resolve the list of (assembly, entry) pairs to build.
    let targets = resolve_targets(&config.vareffect, assemblies)?;

    eprintln!();
    eprintln!(
        "{}  ({} assembly target{})",
        style("Setting up vareffect runtime data").bold().cyan(),
        targets.len(),
        if targets.len() == 1 { "" } else { "s" },
    );
    eprintln!("  output_dir = {}", output_dir.display());
    eprintln!("  raw_dir    = {}", raw_dir.display());
    for (assembly, _) in &targets {
        eprintln!("  build      = {assembly}");
    }
    eprintln!();

    for (assembly, entry) in &targets {
        run_one_assembly(
            &http_client,
            *assembly,
            entry,
            &config,
            &output_dir,
            &raw_dir,
            fasta_only,
            models_only,
        )?;
    }

    eprintln!(
        "{}  vareffect-cli setup complete in {}s",
        style("OK").green().bold(),
        total_start.elapsed().as_secs(),
    );

    Ok(())
}

/// Resolve the `(Assembly, &AssemblyEntry)` list the CLI asked for,
/// erroring out for entries the config doesn't define.
fn resolve_targets(
    section: &config::VareffectSection,
    filter: AssemblyFilter,
) -> Result<Vec<(Assembly, &AssemblyEntry)>> {
    let mut out: Vec<(Assembly, &AssemblyEntry)> = Vec::new();
    match filter {
        AssemblyFilter::Single(asm) => {
            let entry = section.entry(asm).with_context(|| {
                format!(
                    "config has no [vareffect.{}] sub-table; \
                     define one or pass `--assembly` for a different build",
                    asm.as_str().to_ascii_lowercase()
                )
            })?;
            out.push((asm, entry));
        }
        AssemblyFilter::All => {
            if let Some(entry) = &section.grch38 {
                out.push((Assembly::GRCh38, entry));
            }
            if let Some(entry) = &section.grch37 {
                out.push((Assembly::GRCh37, entry));
            }
            if out.is_empty() {
                bail!(
                    "config has no [vareffect.grch38] or [vareffect.grch37] \
                     sub-table; nothing to build"
                );
            }
        }
    }
    Ok(out)
}

/// Build all artifacts for a single assembly. Mirrors the previous
/// monolithic flow but takes the per-assembly entry as an argument so
/// the orchestrator can fan out across multiple builds.
#[allow(clippy::too_many_arguments)]
fn run_one_assembly(
    http_client: &reqwest::blocking::Client,
    assembly: Assembly,
    entry: &AssemblyEntry,
    config: &config::VareffectConfig,
    output_dir: &Path,
    raw_dir: &Path,
    fasta_only: bool,
    models_only: bool,
) -> Result<()> {
    eprintln!(
        "{}",
        style(format!("=== Building {assembly} ===")).bold().cyan()
    );

    // Sanity-check that `fasta_url` matches the declared `fasta_build`.
    if !entry.fasta_url.contains(&entry.fasta_build) {
        bail!(
            "config mismatch for {assembly}: fasta_url={:?} does not contain \
             fasta_build={:?} -- refusing to download a FASTA that may be \
             the wrong assembly",
            entry.fasta_url,
            entry.fasta_build,
        );
    }

    // -- Step 1/3: Flat binary reference genome -------------------------
    if !models_only {
        eprintln!(
            "{}",
            style("[1/3] Flat binary reference genome (NCBI)").bold()
        );
        let outcome = reference_genome::ensure_genome(
            http_client,
            &entry.fasta_url,
            output_dir,
            &entry.fasta_filename,
            &entry.fasta_build,
            raw_dir,
        )
        .with_context(|| format!("preparing reference genome binary for {assembly}"))?;
        report_genome_outcome(&outcome, &output_dir.join(&entry.fasta_filename));
        eprintln!();
    } else {
        eprintln!(
            "{}",
            style("[1/3] Flat binary genome (skipped via --models-only)").dim()
        );
        eprintln!();
    }

    // -- Step 2/3: Patch-chrom aliases ---------------------------------
    let aliases_path = if !fasta_only {
        eprintln!("{}", style("[2/3] Patch-chrom aliases").bold());

        let assembly_report_path = raw_dir.join(&entry.assembly_report_input);
        ensure_cached(
            http_client,
            &entry.assembly_report_url,
            &assembly_report_path,
            "NCBI assembly report",
        )?;

        let (csv_path, row_count) = builders::patch_chrom_aliases::build_from_assembly_report(
            &assembly_report_path,
            output_dir,
            &entry.patch_aliases_filename,
        )
        .with_context(|| format!("building {}", entry.patch_aliases_filename))?;

        eprintln!(
            "  {}  wrote {} ({} aliases)",
            style("OK").green(),
            style(&entry.patch_aliases_filename).cyan(),
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

        // GFF3 cache filename. For GRCh38 use the optional `[mane].input`
        // override; for GRCh37 derive a default name from the assembly.
        let gff_filename = match assembly {
            Assembly::GRCh38 => config.gff_filename().to_string(),
            Assembly::GRCh37 => "ncbi_grch37.gff.gz".to_string(),
        };
        let gff_path = raw_dir.join(&gff_filename);
        ensure_cached(http_client, &entry.gff_url, &gff_path, "RefSeq GFF3")?;

        let cross_validation_path = if let Some(input) = &entry.cross_validation_input {
            let path = raw_dir.join(input);
            if let Some(url) = &entry.cross_validation_url {
                ensure_cached(http_client, url, &path, "cross-validation source")?;
            }
            Some(path)
        } else {
            None
        };

        let admit_tags = match entry.transcript_source {
            TranscriptSource::Mane => builders::transcript_models::ADMIT_TAGS_MANE,
            TranscriptSource::RefSeqSelect => builders::transcript_models::ADMIT_TAGS_REFSEQ_SELECT,
        };

        let build_start = Instant::now();
        let aliases_ref = aliases_path.as_deref();
        let (out, size_bytes) = builders::transcript_models::build(
            &gff_path,
            cross_validation_path.as_deref(),
            aliases_ref,
            output_dir,
            &entry.transcript_models_version,
            assembly,
            admit_tags,
            &entry.transcript_models_basename,
        )
        .with_context(|| format!("building transcript models for {assembly}"))?;

        eprintln!(
            "  {}  wrote {} ({} records, {} in {}s)",
            style("OK").green(),
            style(&entry.transcript_models_filename).cyan(),
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

/// Pretty-print a [`GenomeOutcome`].
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
