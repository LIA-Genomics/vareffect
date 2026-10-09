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
mod source_file;
mod vcf;

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};

use crate::common::{print_build_header, print_build_summary};
use crate::config::SourceName;

/// Build tooling for vareffect runtime data.
///
/// Downloads and indexes a GRCh38 reference FASTA, builds transcript model
/// stores from MANE, RefSeq or GENCODE GFF3 releases, and provisions everything the `vareffect`
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
    /// Downloads the GRCh38 reference FASTA and converts it to a flat binary
    /// (`GRCh38.bin` + `.bin.idx`), downloads the NCBI assembly report and
    /// builds patch-chrom aliases, then builds the store of the transcript
    /// source selected by `[vareffect.transcripts] source`: MANE (GFF3 +
    /// summary TSV -> `transcript_models.bin`), RefSeq
    /// (`transcript_models_refseq.bin`) or GENCODE
    /// (`transcript_models_gencode.bin`). RefSeq and GENCODE downloads must
    /// match their pinned SHA-256, release and genome build.
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
        /// Number of threads for parallel annotation (0 = all logical cores).
        #[arg(long, default_value_t = 0)]
        threads: usize,
        /// Path to `patch_chrom_aliases.csv` (optional).
        #[arg(long)]
        patch_aliases: Option<PathBuf>,
    },

    /// Build a transcript model store from a MANE, RefSeq or GENCODE GFF3.
    ///
    /// Produces `transcript_models.bin` (MANE), `transcript_models_refseq.bin`
    /// or `transcript_models_gencode.bin`, with a `.excluded.tsv` listing
    /// every candidate transcript that was not admitted and why. When
    /// `--summary-input` is provided (MANE only), every transcript is
    /// cross-validated against the MANE summary TSV and the build fails on
    /// mismatches.
    ///
    /// Prefer `setup` for production MANE builds -- it downloads the source
    /// files automatically.
    BuildTranscripts {
        /// Annotation source of `--input`.
        #[arg(long, value_enum, default_value_t = SourceName::Mane)]
        source: SourceName,
        /// Path to the source GFF3: `MANE.GRCh38.vX.X.refseq_genomic.gff.gz`,
        /// `GCF_000001405.40_GRCh38.p14_genomic.gff.gz` or
        /// `gencode.vNN.annotation.gff3.gz`.
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
        /// Source release (e.g. "1.5" for MANE, "RS_2025_08" for RefSeq,
        /// "50" for GENCODE).
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
            source,
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
            let source = vareffect::TranscriptSource::from(source);
            let (out, size) =
                builders::transcript_models::build(&builders::transcript_models::BuildRequest {
                    source,
                    input: &input,
                    summary_input: summary_input.as_deref(),
                    patch_aliases_input: patch_chrom_aliases.as_deref(),
                    output_dir: &output,
                    release: &version,
                })?;
            print_build_header();
            print_build_summary(common::store_stem(source)?, &out, size, start);
        }
    }

    Ok(())
}
