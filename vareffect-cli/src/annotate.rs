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

use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::panic::AssertUnwindSafe;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use vareffect::{VarEffect, VarEffectError};

use crate::csq;
use crate::vcf;

/// Number of data lines per parallel chunk.
///
/// 10 000 lines * ~200 bytes/line = ~2 MB per chunk. Large enough to
/// amortize rayon overhead, small enough to bound memory for WGS VCFs.
const CHUNK_SIZE: usize = 10_000;

/// Number of chunks the reader thread may buffer ahead of annotation.
///
/// Bounds memory to ~`PIPELINE_DEPTH * CHUNK_SIZE` lines while letting input
/// reading overlap annotation and output compression.
const PIPELINE_DEPTH: usize = 3;

/// Configuration for the annotate pipeline.
pub struct AnnotateConfig<'a> {
    /// Path to the input VCF (`.vcf` or `.vcf.gz`).
    pub input: &'a Path,
    /// Path to the output VCF (`.vcf` or `.vcf.gz`).
    pub output: &'a Path,
    /// Path to the flat binary genome (`GRCh38.bin`).
    pub fasta: &'a Path,
    /// Path to `transcript_models.bin`.
    pub transcripts: &'a Path,
    /// Number of rayon worker threads.
    pub threads: usize,
    /// Optional path to `patch_chrom_aliases.csv`.
    pub patch_aliases: Option<&'a Path>,
}

/// Annotation statistics, updated atomically from rayon worker threads.
///
/// The first four count *alleles*; `malformed` counts unparseable *lines*.
/// `chrom_not_found` dedupes the missing-chromosome warning.
struct Counters {
    annotated: AtomicU64,
    intergenic: AtomicU64,
    skipped: AtomicU64,
    errored: AtomicU64,
    malformed: AtomicU64,
    chrom_not_found: Mutex<HashSet<String>>,
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

    // -- 1. Configure rayon thread pool -----------------------------------
    let threads = resolve_threads(config.threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("configuring rayon thread pool")?;

    // -- 2. Load VarEffect ------------------------------------------------
    tracing::info!("loading vareffect data");
    let load_start = Instant::now();
    let ve = Arc::new(
        VarEffect::open_with_patch_aliases(config.transcripts, config.fasta, config.patch_aliases)
            .context("loading vareffect data")?,
    );
    let load_elapsed = load_start.elapsed();
    tracing::info!(
        transcripts = ve.transcripts().len(),
        "vareffect loaded in {:.1}ms",
        load_elapsed.as_secs_f64() * 1000.0
    );

    // -- 3. Open I/O ------------------------------------------------------
    let reader = vcf::open_reader(config.input)?;
    // BGZF output compresses on its own worker pool; give it the annotation
    // budget (the stages overlap in the pipeline -- tune via --threads).
    let mut writer = vcf::open_writer(config.output, threads)?;

    // -- 4. Progress bar --------------------------------------------------
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] {msg}")
            .expect("valid progress template"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(200));
    pb.set_message("processing headers");

    // -- 5. Process headers -----------------------------------------------
    // Drop trailing '\r' so CRLF input doesn't leave it in the spliced output.
    let mut lines = reader.lines().map(|r| r.map(strip_trailing_cr));
    process_headers(&mut lines, &mut writer)?;

    // -- 6. Pipelined parallel annotation ---------------------------------
    // A reader thread prefetches chunks while the main thread annotates each
    // chunk in parallel (rayon) and writes it; BGZF compression runs on the
    // writer's own worker pool. The stages overlap. Output order is preserved
    // (FIFO channel + a single in-order annotator).
    let annotate_start = Instant::now();
    let counters = Counters {
        annotated: AtomicU64::new(0),
        intergenic: AtomicU64::new(0),
        skipped: AtomicU64::new(0),
        errored: AtomicU64::new(0),
        malformed: AtomicU64::new(0),
        chrom_not_found: Mutex::new(HashSet::new()),
    };

    let (tx, rx) = sync_channel::<Vec<String>>(PIPELINE_DEPTH);
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut chunk = Vec::with_capacity(CHUNK_SIZE);
        for line_result in lines {
            chunk.push(line_result.context("reading VCF data line")?);
            if chunk.len() == CHUNK_SIZE {
                let full = std::mem::replace(&mut chunk, Vec::with_capacity(CHUNK_SIZE));
                if tx.send(full).is_err() {
                    return Ok(()); // receiver hung up (main errored); stop quietly
                }
            }
        }
        if !chunk.is_empty() {
            let _ = tx.send(chunk);
        }
        Ok(())
    });

    let mut total = 0u64;
    for chunk in rx {
        let annotated = annotate_chunk(&chunk, &ve, &counters);
        for line in &annotated {
            writeln!(writer, "{line}").context("writing annotated line")?;
        }
        total += chunk.len() as u64;
        pb.set_message(format!("{total} variants processed"));
    }

    // Reader finished (channel closed): surface any read error before finalize.
    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    writer.finish().context("finalizing output")?;

    // -- 7. Summary -------------------------------------------------------
    let annotate_elapsed = annotate_start.elapsed();
    let total_elapsed = start.elapsed();
    let annotated = counters.annotated.load(Ordering::Relaxed);
    let intergenic = counters.intergenic.load(Ordering::Relaxed);
    let skipped = counters.skipped.load(Ordering::Relaxed);
    let errored = counters.errored.load(Ordering::Relaxed);
    let malformed = counters.malformed.load(Ordering::Relaxed);
    let annotate_secs = annotate_elapsed.as_secs_f64();
    let rate = if annotate_secs > 0.0 {
        total as f64 / annotate_secs
    } else {
        0.0
    };

    pb.finish_and_clear();
    tracing::info!(
        "{annotated} annotated, {intergenic} intergenic, {skipped} skipped, {errored} errored \
         (alleles); {malformed} malformed lines; {total} variant lines in {:.1}ms \
         ({rate:.0} lines/sec) [load: {:.1}ms, total: {:.1}ms]",
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

    // No #CHROM found (empty/truncated input): inject CSQ header and warn.
    tracing::warn!(
        "no #CHROM line found in input; wrote header only and annotated no variants \
         (is the input a complete, valid VCF?)"
    );
    writeln!(writer, "{}", csq::CSQ_HEADER).context("writing CSQ header (empty VCF)")?;
    Ok(())
}

/// Drop a single trailing `\r` (`BufRead::lines` strips `\n` but not `\r`).
fn strip_trailing_cr(mut line: String) -> String {
    if line.ends_with('\r') {
        line.pop();
    }
    line
}

/// Resolve worker-thread count: `0` means all logical cores (fallback `1`).
fn resolve_threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(std::num::NonZeroUsize::get)
            .unwrap_or(1)
    } else {
        requested
    }
}

/// Annotate a chunk of VCF lines in parallel, returning the rewritten lines in
/// input order (rayon's `collect` preserves order).
fn annotate_chunk(chunk: &[String], ve: &Arc<VarEffect>, counters: &Counters) -> Vec<String> {
    chunk
        .par_iter()
        .map(|line| annotate_line(line, ve, counters))
        .collect()
}

/// Annotate a single VCF data line, returning the (possibly modified) line.
///
/// On any error (malformed line, annotation failure, panic), the original
/// line is returned unchanged.
fn annotate_line(line: &str, ve: &VarEffect, counters: &Counters) -> String {
    let record = match vcf::parse_vcf_line(line) {
        Ok(r) => r,
        Err(e) => {
            tracing::warn!("skipping malformed line: {e}");
            counters.malformed.fetch_add(1, Ordering::Relaxed);
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
            ve.annotate(&chrom, pos_0based, ref_bytes, alt_bytes)
        }));

        match result {
            Ok(Ok(results)) => {
                let csq_str = csq::format_variant_csq(&trimmed_alt, &results);
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
                if let VarEffectError::ChromNotFound { chrom } = &e {
                    // Warn once per distinct chromosome, not once per allele.
                    let first_sight = counters
                        .chrom_not_found
                        .lock()
                        .map(|mut seen| seen.insert(chrom.clone()))
                        .unwrap_or(true);
                    if first_sight {
                        tracing::warn!(
                            chrom = %chrom,
                            "chromosome not found in genome index; variants on it cannot be \
                             annotated (is the VCF using a chromosome naming convention \
                             vareffect does not recognize?)"
                        );
                    }
                } else {
                    tracing::warn!(
                        chrom = %record.chrom,
                        pos = record.pos,
                        alt = %alt,
                        "annotation failed: {e}"
                    );
                }
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_trailing_cr_removes_one_cr() {
        assert_eq!(strip_trailing_cr("chr1\t100".to_string()), "chr1\t100");
        assert_eq!(strip_trailing_cr("chr1\t100\r".to_string()), "chr1\t100");
        assert_eq!(strip_trailing_cr(String::new()), "");
        assert_eq!(strip_trailing_cr("x\r\r".to_string()), "x\r");
    }

    #[test]
    fn resolve_threads_auto_is_at_least_one() {
        assert!(resolve_threads(0) >= 1);
        assert_eq!(resolve_threads(1), 1);
        assert_eq!(resolve_threads(8), 8);
    }

    #[test]
    fn process_headers_injects_csq_immediately_before_chrom() {
        let input: [std::io::Result<String>; 2] = [
            Ok("##fileformat=VCFv4.2".to_string()),
            Ok("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string()),
        ];
        let mut lines = input.into_iter();
        let mut out: Vec<u8> = Vec::new();
        process_headers(&mut lines, &mut out).expect("process headers");

        let text = String::from_utf8(out).expect("utf8");
        let csq_pos = text.find("ID=CSQ").expect("CSQ header written");
        let chrom_pos = text.find("#CHROM").expect("#CHROM written");
        assert!(csq_pos < chrom_pos, "CSQ header must precede #CHROM");
    }

    #[test]
    fn process_headers_tolerates_missing_chrom() {
        let input: [std::io::Result<String>; 2] = [
            Ok("##fileformat=VCFv4.2".to_string()),
            Ok("##contig=<ID=chr1,length=1000>".to_string()),
        ];
        let mut lines = input.into_iter();
        let mut out: Vec<u8> = Vec::new();
        process_headers(&mut lines, &mut out).expect("process headers");

        let text = String::from_utf8(out).expect("utf8");
        assert!(text.contains("ID=CSQ"), "CSQ header still injected");
        assert!(
            !text.contains("#CHROM"),
            "no #CHROM line in header-only input"
        );
    }
}
