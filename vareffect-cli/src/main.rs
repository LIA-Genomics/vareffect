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
use std::str::FromStr;
use std::time::Instant;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand, ValueEnum};
use vareffect::Assembly;

use crate::common::{print_build_header, print_build_summary};

/// CLI value for `--assembly`. Maps to [`vareffect::Assembly`] for the
/// runtime path; the `All` variant is `setup`-only and means "build
/// every assembly the config defines".
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum AssemblyArg {
    /// GRCh38 (NCBI MANE source).
    Grch38,
    /// GRCh37 (NCBI RefSeq Select source).
    Grch37,
    /// Build every assembly defined in the config. Setup-only — the
    /// `annotate` subcommand requires a single assembly.
    All,
}

impl AssemblyArg {
    fn to_assembly(self) -> Result<Assembly> {
        match self {
            AssemblyArg::Grch38 => Ok(Assembly::GRCh38),
            AssemblyArg::Grch37 => Ok(Assembly::GRCh37),
            AssemblyArg::All => bail!("--assembly all is only valid for `setup`"),
        }
    }
}

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
        /// Reference build to set up. Defaults to `grch38`; `all` builds
        /// every assembly defined in the config (an additional ~6 GB of
        /// raw downloads + ~3 GB of binaries on top of GRCh38).
        #[arg(long, value_enum, default_value_t = AssemblyArg::Grch38)]
        assembly: AssemblyArg,
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
        /// Reference build the input VCF was called against. Required —
        /// no default, because misrouting under the wrong chrom table
        /// silently produces off-by-millions annotations.
        ///
        /// Accepts `grch37` / `grch38`. UCSC aliases `hg19` / `hg38` are
        /// rejected with an error pointing the caller at the explicit
        /// form (UCSC `hg19` chrM differs from GRCh37 chrMT, so the
        /// alias would silently mis-annotate every chrM variant).
        #[arg(long, value_enum)]
        assembly: AssemblyArg,
        /// Path to the input VCF file (`.vcf` or `.vcf.gz`).
        #[arg(long)]
        input: PathBuf,
        /// Path to the output VCF file (`.vcf` or `.vcf.gz`).
        #[arg(long)]
        output: PathBuf,
        /// Path to the flat-binary genome (e.g. `GRCh38.bin` or
        /// `GRCh37.bin`). Must match `--assembly`.
        #[arg(long)]
        fasta: PathBuf,
        /// Path to the transcript models bin (e.g.
        /// `transcript_models_grch38.bin`).
        #[arg(long)]
        transcripts: PathBuf,
        /// Number of threads for parallel annotation.
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Path to the per-assembly `patch_chrom_aliases_*.csv` (optional).
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
        /// Reference build to drive admit-tag selection (`grch38` →
        /// MANE, `grch37` → RefSeq Select) and chrom-table translation.
        #[arg(long, value_enum)]
        assembly: AssemblyArg,
        /// Path to the GFF3 input file. For GRCh38: MANE
        /// `MANE.GRCh38.vX.X.refseq_genomic.gff.gz`. For GRCh37: NCBI
        /// `GCF_000001405.25_GRCh37.p13_genomic.gff.gz`.
        #[arg(long)]
        input: PathBuf,
        /// Optional path to a cross-validation source. For GRCh38:
        /// `MANE.GRCh38.vX.X.summary.txt.gz`. GRCh37 has no Stage A
        /// equivalent — pass `None` and the build emits a documented
        /// warning about the missing quality gate.
        #[arg(long)]
        summary_input: Option<PathBuf>,
        /// Path to the per-assembly `patch_chrom_aliases_*.csv`
        /// (produced by `setup`). Required whenever `--summary-input`
        /// is set.
        #[arg(long)]
        patch_chrom_aliases: Option<PathBuf>,
        /// Output directory.
        #[arg(long, default_value = "data/vareffect")]
        output: PathBuf,
        /// Upstream data version string (e.g. "1.5" for MANE,
        /// "GRCh37.p13" for the NCBI RefSeq GRCh37 release).
        #[arg(long)]
        version: String,
        /// File-stem (no extension) for the output bin / sha256 /
        /// manifest sibling files. Defaults to a per-assembly value.
        #[arg(long)]
        basename: Option<String>,
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
            assembly,
            input,
            output,
            fasta,
            transcripts,
            threads,
            patch_aliases,
        } => {
            let assembly = assembly.to_assembly()?;
            annotate::run(&annotate::AnnotateConfig {
                assembly,
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
            assembly,
            fasta_only,
            models_only,
        } => {
            let filter = match assembly {
                AssemblyArg::Grch38 => setup::AssemblyFilter::Single(Assembly::GRCh38),
                AssemblyArg::Grch37 => setup::AssemblyFilter::Single(Assembly::GRCh37),
                AssemblyArg::All => setup::AssemblyFilter::All,
            };
            setup::run(
                config.as_deref(),
                fasta_only,
                models_only,
                output.as_deref(),
                filter,
            )?;
        }

        Command::BuildTranscripts {
            assembly,
            input,
            summary_input,
            patch_chrom_aliases,
            output,
            version,
            basename,
        } => {
            if summary_input.is_some() && patch_chrom_aliases.is_none() {
                bail!(
                    "--summary-input requires --patch-chrom-aliases (point it at \
                     the per-assembly patch_chrom_aliases CSV produced by \
                     `vareffect-cli setup`)"
                );
            }
            let asm = assembly.to_assembly()?;
            let admit_tags = match asm {
                Assembly::GRCh38 => builders::transcript_models::ADMIT_TAGS_MANE,
                Assembly::GRCh37 => builders::transcript_models::ADMIT_TAGS_REFSEQ_SELECT,
            };
            let basename = basename.unwrap_or_else(|| match asm {
                Assembly::GRCh38 => "transcript_models_grch38".to_string(),
                Assembly::GRCh37 => "transcript_models_grch37".to_string(),
            });
            let start = Instant::now();
            let (out, size) = builders::transcript_models::build(
                &input,
                summary_input.as_deref(),
                patch_chrom_aliases.as_deref(),
                &output,
                &version,
                asm,
                admit_tags,
                &basename,
            )?;
            print_build_header();
            print_build_summary(&basename, &out, size, start);
        }
    }

    Ok(())
}

/// Bridge `AssemblyArg::Grch{37,38}` to [`vareffect::Assembly`] without
/// triggering the `--assembly all` rejection. The runtime parser is
/// unused for the non-`All` cases but kept as a shared correctness check
/// for the rejection-of-aliases path that the CLI never reaches.
#[allow(dead_code)]
fn _check_assembly_alias_rejection() {
    debug_assert!(Assembly::from_str("hg19").is_err());
    debug_assert!(Assembly::from_str("hg38").is_err());
}
