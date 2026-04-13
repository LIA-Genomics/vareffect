//! vareffect-cli -- build tooling for vareffect runtime data.
//!
//! Downloads reference FASTA, builds transcript model stores, and provisions
//! everything the `vareffect` crate needs at runtime. Ships as a standalone
//! binary with zero `lia-*` dependencies so `cargo install vareffect-cli`
//! works for external users.

mod annotate;
mod builders;
mod check;
mod common;
mod config;
mod csq;
mod init;
mod setup;
mod vcf;

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};

use crate::common::{print_build_header, print_build_summary};

/// Build tooling for vareffect runtime data.
///
/// Downloads and indexes a GRCh38 reference FASTA, builds transcript model
/// stores from MANE GFF3 releases, and provisions everything the `vareffect`
/// crate needs at runtime under `data/vareffect/`.
#[derive(Parser)]
#[command(name = "vareffect", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

/// Available subcommands.
#[derive(Subcommand)]
enum Command {
    /// Initialize a configuration file at a discoverable location.
    ///
    /// Writes a default `vareffect_build.toml` so that `setup`, `check`,
    /// and other subcommands can find it without an explicit `--config`
    /// flag. Run this once after installing.
    Init(init::InitArgs),

    /// Validate configuration and reference data files.
    ///
    /// Checks that the config file exists and parses, the data directory
    /// is present, and the genome binary + transcript model store are in
    /// place and readable. Exits with code 1 if any check fails.
    Check(check::CheckArgs),

    /// Set up all runtime data for vareffect in one command.
    ///
    /// Downloads the GRCh38 reference FASTA, decompresses and indexes it
    /// (produces `GRCh38.fa.gz` + `.fai` + `.gzi`), downloads the NCBI
    /// assembly report and builds patch-chrom aliases, then downloads the
    /// MANE GFF3 + summary TSV and builds `transcript_models.bin`.
    ///
    /// Idempotent: the reference FASTA is skipped if already present;
    /// transcript models are always rebuilt. Source archives land in
    /// `data/raw/` as a cache so repeated runs skip the download step.
    Setup {
        /// Path to `vareffect_build.toml`. Auto-discovered if omitted.
        #[arg(long)]
        config: Option<PathBuf>,
        /// Override the output directory from `vareffect_build.toml`.
        /// When set, runtime files (genome binary, transcript models, etc.)
        /// are written here instead of `[vareffect].output_dir`.
        #[arg(long)]
        output: Option<PathBuf>,
        /// Only download + index the reference FASTA; skip the transcript
        /// model build.
        #[arg(long, conflicts_with = "models_only")]
        fasta_only: bool,
        /// Only build the transcript model store; skip the reference FASTA
        /// download + index.
        #[arg(long, conflicts_with = "fasta_only")]
        models_only: bool,
    },

    /// Annotate a VCF file with consequence predictions.
    ///
    /// Reads a VCF (plain text or gzip-compressed), annotates each variant
    /// against the vareffect transcript model and reference genome, and writes
    /// an annotated VCF with a CSQ INFO field matching VEP's `--vcf` output
    /// format.
    Annotate {
        /// Path to the input VCF file (`.vcf` or `.vcf.gz`).
        #[arg(long)]
        input: PathBuf,
        /// Path to the output VCF file (`.vcf` or `.vcf.gz`).
        #[arg(long)]
        output: PathBuf,
        /// Path to `GRCh38.bin` (flat binary genome).
        #[arg(long)]
        fasta: PathBuf,
        /// Path to `transcript_models.bin`.
        #[arg(long)]
        transcripts: PathBuf,
        /// Number of threads for parallel annotation.
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Path to `patch_chrom_aliases.csv` (optional).
        #[arg(long)]
        patch_aliases: Option<PathBuf>,
    },

    /// Build transcript model store from a MANE GFF3 file.
    ///
    /// Produces `transcript_models.bin` containing full
    /// `Vec<vareffect::TranscriptModel>` with exon structure, CDS bounds,
    /// protein accessions, and MANE tier. When `--summary-input` is
    /// provided, every transcript is cross-validated against the MANE
    /// summary TSV and the build fails on mismatches.
    ///
    /// Prefer `setup` for production builds -- it downloads the source
    /// files automatically. This standalone subcommand is kept for
    /// iterative development.
    BuildTranscripts {
        /// Path to `MANE.GRCh38.vX.X.refseq_genomic.gff.gz`.
        #[arg(long)]
        input: PathBuf,
        /// Optional path to `MANE.GRCh38.vX.X.summary.txt.gz` for
        /// GFF3-vs-summary cross-validation.
        #[arg(long)]
        summary_input: Option<PathBuf>,
        /// Path to `patch_chrom_aliases.csv` (produced by `setup`).
        /// Required whenever `--summary-input` is set.
        #[arg(long)]
        patch_chrom_aliases: Option<PathBuf>,
        /// Output directory.
        #[arg(long, default_value = "data/vareffect")]
        output: PathBuf,
        /// MANE release version (e.g. "1.5").
        #[arg(long)]
        version: String,
    },
}

/// Install a compact stderr tracing subscriber.
/// `RUST_LOG=vareffect_cli=warn` (or `info`/`debug`) filters event level.
fn init_tracing() {
    use tracing_subscriber::{EnvFilter, fmt};

    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .without_time()
        .with_target(false)
        .compact()
        .init();
}

fn main() -> Result<()> {
    init_tracing();
    let cli = Cli::parse();

    match cli.command {
        Command::Init(args) => init::run(&args)?,

        Command::Check(args) => {
            if !check::run(&args)? {
                std::process::exit(1);
            }
        }

        Command::Annotate {
            input,
            output,
            fasta,
            transcripts,
            threads,
            patch_aliases,
        } => {
            annotate::run(&annotate::AnnotateConfig {
                input: &input,
                output: &output,
                fasta: &fasta,
                transcripts: &transcripts,
                threads,
                patch_aliases: patch_aliases.as_deref(),
            })?;
        }

        Command::Setup {
            config,
            output,
            fasta_only,
            models_only,
        } => {
            setup::run(
                config.as_deref(),
                fasta_only,
                models_only,
                output.as_deref(),
            )?;
        }

        Command::BuildTranscripts {
            input,
            summary_input,
            patch_chrom_aliases,
            output,
            version,
        } => {
            if summary_input.is_some() && patch_chrom_aliases.is_none() {
                bail!(
                    "--summary-input requires --patch-chrom-aliases (point it at \
                     the patch_chrom_aliases.csv produced by `vareffect-cli setup`)"
                );
            }
            let start = Instant::now();
            let (out, size) = builders::transcript_models::build(
                &input,
                summary_input.as_deref(),
                patch_chrom_aliases.as_deref(),
                &output,
                &version,
            )?;
            print_build_header();
            print_build_summary("transcript_models", &out, size, start);
        }
    }

    Ok(())
}
