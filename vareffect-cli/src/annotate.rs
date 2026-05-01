//! `vareffect annotate` subcommand — annotate a VCF with consequence
//! predictions and write a VEP-compatible CSQ INFO field.
//!
//! Reads a VCF (plain text or gzip), annotates each variant via
//! [`vareffect::VarEffect::annotate`], and writes the annotated VCF with a
//! `CSQ` INFO field matching VEP's `--vcf` output format.
//!
//! # Threading model
//!
//! Data lines are read into chunks of [`CHUNK_SIZE`]. Each chunk is annotated
//! in parallel with `rayon::par_iter().map().collect()` (which preserves
//! input order), then written sequentially. This bounds memory to one chunk
//! while maintaining VCF output order.
//!
//! # Error handling
//!
//! Malformed lines, annotation errors, and panics inside `annotate()` are
//! logged as warnings. The original line is passed through unmodified. The
//! pipeline never drops a line or aborts for a single bad variant.

use std::io::{BufRead, Write};
use std::panic::AssertUnwindSafe;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use vareffect::{Assembly, VarEffect};

use crate::csq;
use crate::vcf;

/// Number of data lines per parallel chunk.
///
/// 10 000 lines * ~200 bytes/line = ~2 MB per chunk. Large enough to
/// amortize rayon overhead, small enough to bound memory for WGS VCFs.
const CHUNK_SIZE: usize = 10_000;

/// Configuration for the annotate pipeline.
pub struct AnnotateConfig<'a> {
    /// Reference build the input VCF was called against. Drives both
    /// the [`VarEffect`] slot routing and the chrom-table selection
    /// inside [`vareffect::FastaReader`].
    pub assembly: Assembly,
    /// Path to the input VCF (`.vcf` or `.vcf.gz`).
    pub input: &'a Path,
    /// Path to the output VCF (`.vcf` or `.vcf.gz`).
    pub output: &'a Path,
    /// Path to the flat binary genome (e.g. `GRCh38.bin` / `GRCh37.bin`).
    pub fasta: &'a Path,
    /// Path to the per-assembly transcript models bin.
    pub transcripts: &'a Path,
    /// Number of rayon worker threads.
    pub threads: usize,
    /// Optional path to the per-assembly `patch_chrom_aliases_*.csv`.
    pub patch_aliases: Option<&'a Path>,
}

/// Annotation statistics, updated atomically from rayon worker threads.
struct Counters {
    annotated: AtomicU64,
    intergenic: AtomicU64,
    skipped: AtomicU64,
    errored: AtomicU64,
}

/// Run the annotate pipeline.
///
/// Loads the vareffect engine, reads the input VCF, annotates each variant
/// in parallel chunks, and writes the output VCF with CSQ INFO fields.
///
/// # Errors
///
/// Returns an error if the input/output files cannot be opened, the
/// vareffect data fails to load, or an I/O error occurs during writing.
pub fn run(config: &AnnotateConfig<'_>) -> Result<()> {
    let start = Instant::now();

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .context("configuring rayon thread pool")?;

    tracing::info!(assembly = %config.assembly, "loading vareffect data");
    let load_start = Instant::now();
    let mut builder = VarEffect::builder();
    builder = match (config.assembly, config.patch_aliases) {
        (Assembly::GRCh38, Some(p)) => builder
            .with_grch38_and_patch_aliases(config.transcripts, config.fasta, p)
            .context("loading GRCh38 data with patch aliases")?,
        (Assembly::GRCh38, None) => builder
            .with_grch38(config.transcripts, config.fasta)
            .context("loading GRCh38 data")?,
        (Assembly::GRCh37, Some(p)) => builder
            .with_grch37_and_patch_aliases(config.transcripts, config.fasta, p)
            .context("loading GRCh37 data with patch aliases")?,
        (Assembly::GRCh37, None) => builder
            .with_grch37(config.transcripts, config.fasta)
            .context("loading GRCh37 data")?,
    };
    let ve = Arc::new(builder.build().context("assembling VarEffect")?);
    let load_elapsed = load_start.elapsed();
    tracing::info!(
        transcripts = ve.transcripts(config.assembly)?.len(),
        "vareffect loaded in {:.1}ms",
        load_elapsed.as_secs_f64() * 1000.0
    );

    let reader = vcf::open_reader(config.input)?;
    let mut writer = vcf::open_writer(config.output)?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] {msg}")
            .expect("valid progress template"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(200));
    pb.set_message("processing headers");

    let mut lines = reader.lines();
    process_headers(&mut lines, &mut writer)?;

    let annotate_start = Instant::now();
    let counters = Counters {
        annotated: AtomicU64::new(0),
        intergenic: AtomicU64::new(0),
        skipped: AtomicU64::new(0),
        errored: AtomicU64::new(0),
    };

    let mut total = 0u64;
    let mut chunk = Vec::with_capacity(CHUNK_SIZE);

    for line_result in lines {
        let line = line_result.context("reading VCF data line")?;
        chunk.push(line);

        if chunk.len() == CHUNK_SIZE {
            process_chunk(&chunk, &ve, config.assembly, &counters, &mut writer)?;
            total += chunk.len() as u64;
            pb.set_message(format!("{total} variants processed"));
            chunk.clear();
        }
    }

    // Flush remaining lines.
    if !chunk.is_empty() {
        process_chunk(&chunk, &ve, config.assembly, &counters, &mut writer)?;
        total += chunk.len() as u64;
    }

    writer.flush().context("flushing output")?;

    let annotate_elapsed = annotate_start.elapsed();
    let total_elapsed = start.elapsed();
    let annotated = counters.annotated.load(Ordering::Relaxed);
    let intergenic = counters.intergenic.load(Ordering::Relaxed);
    let skipped = counters.skipped.load(Ordering::Relaxed);
    let errored = counters.errored.load(Ordering::Relaxed);
    let annotate_secs = annotate_elapsed.as_secs_f64();
    let rate = if annotate_secs > 0.0 {
        total as f64 / annotate_secs
    } else {
        0.0
    };

    pb.finish_and_clear();
    tracing::info!(
        "{annotated} annotated, {intergenic} intergenic, {skipped} skipped, {errored} errors. \
         {total} variants in {:.1}ms ({rate:.0} variants/sec) \
         [load: {:.1}ms, total: {:.1}ms]",
        annotate_elapsed.as_secs_f64() * 1000.0,
        load_elapsed.as_secs_f64() * 1000.0,
        total_elapsed.as_secs_f64() * 1000.0,
    );

    Ok(())
}

/// Read VCF header lines, inject the CSQ header, and write them to output.
///
/// Passes all `##` meta-information lines through verbatim. Inserts
/// [`csq::CSQ_HEADER`] immediately before the `#CHROM` column header line,
/// then writes the `#CHROM` line itself.
fn process_headers(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
    writer: &mut dyn Write,
) -> Result<()> {
    for line_result in lines {
        let line = line_result.context("reading VCF header")?;
        if line.starts_with("##") {
            writeln!(writer, "{line}").context("writing meta header")?;
        } else if line.starts_with("#CHROM") {
            writeln!(writer, "{}", csq::CSQ_HEADER).context("writing CSQ header")?;
            writeln!(writer, "{line}").context("writing column header")?;
            return Ok(());
        } else {
            // Data line before #CHROM — malformed VCF but don't crash.
            // Inject CSQ header and then handle this line as data.
            tracing::warn!("data line encountered before #CHROM header");
            writeln!(writer, "{}", csq::CSQ_HEADER).context("writing CSQ header")?;
            writeln!(writer, "{line}").context("writing unexpected data line")?;
            return Ok(());
        }
    }

    // If we exhaust the iterator without finding #CHROM, this is an
    // empty or header-only VCF. Still inject the CSQ header.
    writeln!(writer, "{}", csq::CSQ_HEADER).context("writing CSQ header (empty VCF)")?;
    Ok(())
}

/// Annotate a chunk of VCF lines in parallel and write results in order.
fn process_chunk(
    chunk: &[String],
    ve: &Arc<VarEffect>,
    assembly: Assembly,
    counters: &Counters,
    writer: &mut dyn Write,
) -> Result<()> {
    let annotated: Vec<String> = chunk
        .par_iter()
        .map(|line| annotate_line(line, ve, assembly, counters))
        .collect();

    for line in &annotated {
        writeln!(writer, "{line}").context("writing annotated line")?;
    }

    Ok(())
}

/// Annotate a single VCF data line, returning the (possibly modified) line.
///
/// On any error (malformed line, annotation failure, panic), the original
/// line is returned unchanged.
fn annotate_line(line: &str, ve: &VarEffect, assembly: Assembly, counters: &Counters) -> String {
    let record = match vcf::parse_vcf_line(line) {
        Ok(r) => r,
        Err(e) => {
            tracing::warn!("skipping malformed line: {e}");
            counters.skipped.fetch_add(1, Ordering::Relaxed);
            return line.to_string();
        }
    };

    let mut all_csq_parts: Vec<String> = Vec::new();

    // VCF POS is 1-based; VarEffect::annotate() expects 0-based.
    let pos_0based = record.pos - 1;

    // Normalize chromosome name: the FASTA index expects UCSC-style "chr"
    // prefix (which it translates to NCBI RefSeq internally). ClinVar and
    // Ensembl VCFs use bare names ("1", "X", "MT").
    let chrom = normalize_chrom(record.chrom);

    for alt in &record.alt_alleles {
        // Skip empty, spanning deletions, missing alleles, and symbolic alleles.
        if alt.is_empty()
            || *alt == "*"
            || *alt == "."
            || (alt.starts_with('<') && alt.ends_with('>'))
        {
            tracing::debug!(
                chrom = %record.chrom, pos = record.pos, alt = %alt,
                "skipping non-sequence allele"
            );
            counters.skipped.fetch_add(1, Ordering::Relaxed);
            continue;
        }

        let ref_bytes = record.ref_allele.as_bytes();
        let alt_bytes = alt.as_bytes();

        // Compute the trimmed ALT for the CSQ Allele field.
        let trimmed_alt = csq::trimmed_csq_allele(ref_bytes, alt_bytes);

        // Wrap annotate() in catch_unwind to survive panics. VarEffect is
        // immutable/read-only (mmap'd FASTA + Arc'd transcript store), so
        // AssertUnwindSafe is sound — no mutable state to corrupt on unwind.
        let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
            ve.annotate(assembly, &chrom, pos_0based, ref_bytes, alt_bytes)
        }));

        match result {
            Ok(Ok(results)) => {
                let csq_str = csq::format_variant_csq(&trimmed_alt, &results.consequences);
                if !csq_str.is_empty() {
                    all_csq_parts.push(csq_str);
                    counters.annotated.fetch_add(1, Ordering::Relaxed);
                } else {
                    tracing::debug!(
                        chrom = %record.chrom, pos = record.pos, alt = %alt,
                        "intergenic, no transcript overlap"
                    );
                    counters.intergenic.fetch_add(1, Ordering::Relaxed);
                }
            }
            Ok(Err(e)) => {
                tracing::warn!(
                    chrom = %record.chrom,
                    pos = record.pos,
                    alt = %alt,
                    "annotation failed: {e}"
                );
                counters.skipped.fetch_add(1, Ordering::Relaxed);
            }
            Err(_) => {
                tracing::warn!(
                    chrom = %record.chrom,
                    pos = record.pos,
                    alt = %alt,
                    "annotation panicked, skipping"
                );
                counters.errored.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    if all_csq_parts.is_empty() {
        line.to_string()
    } else {
        let combined = all_csq_parts.join(",");
        vcf::write_annotated_line(&record, &combined)
    }
}

/// Normalize a VCF chromosome name to UCSC-style (`chr1`, `chrX`, `chrM`).
///
/// The `FastaReader` expects UCSC-style names as input (it translates them
/// to NCBI RefSeq internally). ClinVar and Ensembl VCFs use bare names
/// (`1`, `X`, `MT`). This function handles the common cases:
///
/// - Already prefixed (`chr1`) — returned as-is.
/// - Bare autosome/sex (`1`, `X`, `Y`) — prefixed with `chr`.
/// - Mitochondrial (`MT`) — mapped to `chrM`.
fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        return chrom.to_string();
    }
    if chrom == "MT" {
        return "chrM".to_string();
    }
    format!("chr{chrom}")
}
