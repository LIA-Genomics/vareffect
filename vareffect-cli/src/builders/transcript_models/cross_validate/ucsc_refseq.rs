//! UCSC `ncbiRefSeq.txt.gz` + `ncbiRefSeqSelect.txt.gz` parser for the
//! cross-validation pipeline.
//!
//! Used by GRCh37 builds because no NCBI tabular RefSeq Select release
//! exists. Reads the Select sidecar into a `HashSet` of accessions,
//! then streams the genePred-extended coordinates file and emits one
//! [`CrossValidationRow`] per row whose `name` is in the Select set.
//!
//! # Coordinate convention
//!
//! genePred-extended is **0-based half-open** natively, identical to
//! `TranscriptModel::tx_start` / `tx_end`. No conversion needed.
//!
//! # chrM exclusion (load-bearing)
//!
//! UCSC `hg19` chrM is **NC_001807** (the original 1981 Anderson
//! reference). NCBI's GRCh37.p13 chrMT is **NC_012920.1** (the rCRS,
//! ~10 bp shorter, multiple substitutions). Both sources normalize
//! their mitochondrial chromosome name to `chrM`, so naive
//! cross-validation would systematically false-positive on every chrM
//! coordinate comparison.
//!
//! The parser skips chrM rows on the source side and reports `chrM` in
//! `excluded_built_chroms` so the inverse "every built transcript is in
//! the source" check also skips chrM. A single `tracing::warn!`
//! documents the exclusion per build.
//!
//! # File format
//!
//! Both files are tab-separated, no header, genePred-extended (16
//! columns). `ncbiRefSeqSelect.txt.gz` is the genePred-shaped table dump
//! filtered to the RefSeq Select tier; `ncbiRefSeq.txt.gz` carries the
//! full set. The parser uses the Select dump only for its `name` column
//! (the set of accessions to admit) and pulls coordinates from the full
//! `ncbiRefSeq` file.
//!
//! genePred-extended columns used:
//!   1: `name`         — versioned RefSeq accession
//!   2: `chrom`        — UCSC contig name
//!   3: `strand`       — `+` or `-`
//!   4: `txStart`      — 0-based inclusive
//!   5: `txEnd`        — 0-based exclusive
//!  12: `name2`        — gene symbol

use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use vareffect::{Strand, TranscriptTier};

use super::{CrossValidationRow, ParsedSource};

const LABEL: &str = "UCSC ncbiRefSeq+Select";

/// Parse the UCSC ncbiRefSeq + Select pair into a [`ParsedSource`].
///
/// chrM rows are skipped during the stream with a single warning; `chrM`
/// is added to `excluded_built_chroms` so the inverse check honors the
/// same exclusion. Built transcripts on chrM are validated separately
/// (e.g. via ClinVar concordance) because the UCSC and NCBI references
/// disagree.
pub(super) fn parse(coords_path: &Path, select_path: &Path) -> Result<ParsedSource> {
    let select_set = load_select_set(select_path)?;
    let (rows, parse_errors) = parse_coords(coords_path, &select_set)?;

    let mut excluded_built_chroms = HashSet::new();
    excluded_built_chroms.insert("chrM".to_string());

    Ok(ParsedSource {
        label: LABEL,
        rows,
        excluded_built_chroms,
        parse_errors,
    })
}

/// Read `ncbiRefSeqSelect.txt.gz` into a `HashSet<String>` keyed on the
/// versioned accession in column 2. UCSC ships this file as a
/// genePred-extended dump (≥16 columns), identical in shape to
/// `ncbiRefSeq.txt.gz` but filtered to the RefSeq Select tier. We only
/// need the `name` field; coordinates come from the full `ncbiRefSeq`
/// file via [`parse_coords`].
fn load_select_set(path: &Path) -> Result<HashSet<String>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    const NAME_IDX: usize = 1;
    const MIN_FIELDS: usize = 16;

    let mut set: HashSet<String> = HashSet::new();
    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.context("reading ncbiRefSeqSelect line")?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < MIN_FIELDS {
            bail!(
                "ncbiRefSeqSelect at {} line {}: expected genePred-extended row \
                 (>= {} fields), got {} fields: {:?}",
                path.display(),
                line_no + 1,
                MIN_FIELDS,
                fields.len(),
                trimmed,
            );
        }
        let name = fields[NAME_IDX];
        if name.is_empty() {
            bail!(
                "ncbiRefSeqSelect at {} line {}: empty name field",
                path.display(),
                line_no + 1,
            );
        }
        set.insert(name.to_string());
    }

    if set.is_empty() {
        bail!(
            "ncbiRefSeqSelect at {} is empty -- expected several thousand RefSeq \
             Select accessions for hg19",
            path.display()
        );
    }
    Ok(set)
}

/// Stream `ncbiRefSeq.txt.gz` and emit a `CrossValidationRow` for each
/// row whose `name` is in the Select set, skipping chrM. Returns the
/// list of rows alongside any structural-parse errors collected
/// during the stream.
fn parse_coords(
    path: &Path,
    select_set: &HashSet<String>,
) -> Result<(Vec<CrossValidationRow>, Vec<String>)> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // genePred-extended column indices (0-based). The schema is stable
    // across UCSC's hg19/hg38 ncbiRefSeq tables — no header is provided,
    // so pinning the indices is the only option.
    const NAME_IDX: usize = 1;
    const CHROM_IDX: usize = 2;
    const STRAND_IDX: usize = 3;
    const TX_START_IDX: usize = 4;
    const TX_END_IDX: usize = 5;
    const NAME2_IDX: usize = 12;
    const MIN_FIELDS: usize = 16;

    let mut rows: Vec<CrossValidationRow> = Vec::new();
    let mut parse_errors: Vec<String> = Vec::new();
    let mut chrm_skipped: usize = 0;

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.context("reading ncbiRefSeq line")?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < MIN_FIELDS {
            parse_errors.push(format!(
                "{LABEL} line {}: short row -- {} fields, need at least {}",
                line_no + 1,
                fields.len(),
                MIN_FIELDS,
            ));
            continue;
        }

        let name = fields[NAME_IDX];
        if !select_set.contains(name) {
            continue;
        }

        let chrom = fields[CHROM_IDX];
        if chrom == "chrM" {
            chrm_skipped += 1;
            continue;
        }

        let strand = match fields[STRAND_IDX] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => {
                parse_errors.push(format!("{name}: {LABEL} bad strand {other:?}"));
                continue;
            }
        };

        let tx_start = match fields[TX_START_IDX].parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                parse_errors.push(format!(
                    "{name}: {LABEL} unparseable txStart {:?}",
                    fields[TX_START_IDX]
                ));
                continue;
            }
        };
        let tx_end = match fields[TX_END_IDX].parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                parse_errors.push(format!(
                    "{name}: {LABEL} unparseable txEnd {:?}",
                    fields[TX_END_IDX]
                ));
                continue;
            }
        };

        rows.push(CrossValidationRow {
            accession: name.to_string(),
            chrom: chrom.to_string(),
            strand,
            tx_start,
            tx_end,
            gene_symbol: fields[NAME2_IDX].to_string(),
            // UCSC's Select filter is the source of truth for the tier
            // assignment; every row that survived the filter is a
            // RefSeq Select transcript.
            expected_tier: Some(TranscriptTier::RefSeqSelect),
        });
    }

    if chrm_skipped > 0 {
        tracing::warn!(
            "chrM transcripts not cross-validated: UCSC hg19 chrM = NC_001807 \
             (Anderson 1981), GRCh37 chrMT = NC_012920.1 (rCRS) -- coordinates \
             differ by ~10 bp and indels.",
        );
        eprintln!(
            "  UCSC ncbiRefSeq: skipped {chrm_skipped} chrM row(s) due to UCSC \
             NC_001807 vs GRCh37 NC_012920.1 reference divergence",
        );
    }

    Ok((rows, parse_errors))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(suffix: &str, contents: &str) -> NamedTempFile {
        let mut file = tempfile::Builder::new()
            .suffix(suffix)
            .tempfile()
            .expect("tempfile");
        file.write_all(contents.as_bytes()).expect("write");
        file.flush().expect("flush");
        file
    }

    fn ncbi_refseq_row(
        name: &str,
        chrom: &str,
        strand: &str,
        start: u64,
        end: u64,
        name2: &str,
    ) -> String {
        format!(
            "0\t{name}\t{chrom}\t{strand}\t{start}\t{end}\t{start}\t{end}\t1\t{start},\t{end},\t0\t{name2}\tcmpl\tcmpl\t0,\n"
        )
    }

    #[test]
    fn select_set_loads_versioned_accessions() {
        let select = write_temp(
            ".select.txt",
            &(ncbi_refseq_row("NM_000546.5", "chr17", "+", 100, 1000, "TP53")
                + "# comment line\n"
                + &ncbi_refseq_row("NM_006772.2", "chr1", "+", 200, 2000, "GENE_X")
                + &ncbi_refseq_row("NM_007294.4", "chr17", "+", 300, 3000, "BRCA1")),
        );
        let set = load_select_set(select.path()).unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains("NM_000546.5"));
        assert!(set.contains("NM_006772.2"));
        assert!(set.contains("NM_007294.4"));
    }

    #[test]
    fn select_set_rejects_short_rows() {
        // UCSC ships ncbiRefSeqSelect.txt.gz as a genePred-extended dump
        // (>=16 columns). A 2-column row signals a different file and
        // should fail loudly rather than producing an empty Select set.
        let select = write_temp(".select.txt", "1\tNM_000546.5\n");
        let err = load_select_set(select.path()).expect_err("short row must be rejected");
        let msg = format!("{err:#}");
        assert!(msg.contains("genePred-extended"), "msg: {msg}");
    }

    #[test]
    fn coords_filter_to_select_set_only() {
        let select = write_temp(
            ".select.txt",
            &ncbi_refseq_row("NM_SELECT.1", "chr1", "+", 100, 1000, "GENE_A"),
        );
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &(ncbi_refseq_row("NM_SELECT.1", "chr1", "+", 100, 1000, "GENE_A")
                + &ncbi_refseq_row("NM_NOT_SELECT.1", "chr1", "+", 200, 800, "GENE_B")),
        );
        let parsed = parse(coords.path(), select.path()).unwrap();
        assert!(parsed.parse_errors.is_empty());
        assert_eq!(parsed.rows.len(), 1);
        assert_eq!(parsed.rows[0].accession, "NM_SELECT.1");
        assert_eq!(parsed.rows[0].tx_start, 100);
        assert_eq!(parsed.rows[0].tx_end, 1000);
        assert_eq!(parsed.rows[0].strand, Strand::Plus);
        assert_eq!(parsed.rows[0].gene_symbol, "GENE_A");
        assert_eq!(
            parsed.rows[0].expected_tier,
            Some(TranscriptTier::RefSeqSelect),
        );
        assert!(parsed.excluded_built_chroms.contains("chrM"));
        assert_eq!(parsed.label, LABEL);
    }

    #[test]
    fn coords_skip_chrm_with_warning() {
        let select = write_temp(
            ".select.txt",
            &ncbi_refseq_row("NM_MITO.1", "chrM", "+", 100, 200, "MT-CO1"),
        );
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &ncbi_refseq_row("NM_MITO.1", "chrM", "+", 100, 200, "MT-CO1"),
        );
        let parsed = parse(coords.path(), select.path()).unwrap();
        assert!(parsed.rows.is_empty());
        assert!(parsed.excluded_built_chroms.contains("chrM"));
    }

    #[test]
    fn coords_surface_short_rows_as_parse_errors() {
        let select = write_temp(
            ".select.txt",
            &ncbi_refseq_row("NM_SELECT.1", "chr1", "+", 100, 1000, "GENE_A"),
        );
        // 5 fields — well short of MIN_FIELDS (16).
        let coords = write_temp(".ncbiRefSeq.txt", "0\tNM_SELECT.1\tchr1\t+\t100\n");
        let parsed = parse(coords.path(), select.path()).unwrap();
        assert!(parsed.rows.is_empty());
        assert_eq!(parsed.parse_errors.len(), 1);
        assert!(parsed.parse_errors[0].contains("short row"));
    }

    #[test]
    fn coords_surface_bad_strand_as_parse_error() {
        let select = write_temp(
            ".select.txt",
            &ncbi_refseq_row("NM_SEL.1", "chr1", "+", 100, 200, "GENE"),
        );
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &ncbi_refseq_row("NM_SEL.1", "chr1", "?", 100, 200, "GENE"),
        );
        let parsed = parse(coords.path(), select.path()).unwrap();
        assert!(parsed.rows.is_empty());
        assert_eq!(parsed.parse_errors.len(), 1);
        assert!(parsed.parse_errors[0].contains("bad strand"));
    }
}
