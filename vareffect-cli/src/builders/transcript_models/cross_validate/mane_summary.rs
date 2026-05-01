//! NCBI MANE summary TSV parser for the cross-validation pipeline.
//!
//! Parses `MANE.GRCh38.vX.X.summary.txt.gz` into a [`ParsedSource`] for
//! the shared comparison loop. Column layout is resolved from the
//! header (`#GeneID...`) so future column additions don't break parsing.
//!
//! The summary's coordinates are 1-based fully-closed; `parse` converts
//! to 0-based half-open before emitting (matches `TranscriptModel`).
//!
//! GRCh38-only by construction: the source filename, the `MANE_status`
//! column, and the `GRCh38_chr` column name all hard-code GRCh38.

use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use vareffect::{Strand, TranscriptTier};

use super::{CrossValidationRow, ParsedSource};

const LABEL: &str = "MANE summary";

/// Parse a MANE summary TSV at `path` into a [`ParsedSource`].
///
/// `excluded_built_chroms` is empty: both NCBI GRCh38.p14 and MANE use
/// NC_012920.1 for chrM, so no source-vs-built reference divergence
/// applies here.
pub(super) fn parse(path: &Path) -> Result<ParsedSource> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut header_indices: Option<SummaryColumns> = None;
    let mut rows: Vec<CrossValidationRow> = Vec::new();
    let mut parse_errors: Vec<String> = Vec::new();

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.context("reading summary TSV line")?;
        if line.is_empty() {
            continue;
        }

        if header_indices.is_none() {
            // Header starts with `#` in NCBI's format.
            let header = line.trim_start_matches('#');
            header_indices = Some(SummaryColumns::resolve(header).with_context(|| {
                format!(
                    "parsing summary TSV header at {}:{}",
                    path.display(),
                    line_no + 1
                )
            })?);
            continue;
        }
        let cols = header_indices
            .as_ref()
            .expect("header line must precede data rows");

        let fields: Vec<&str> = line.split('\t').collect();
        let Some(raw) = cols.extract(&fields) else {
            parse_errors.push(format!(
                "{LABEL} line {}: short row -- {} fields, need at least {}",
                line_no + 1,
                fields.len(),
                cols.required_field_count(),
            ));
            continue;
        };

        let strand = match raw.chr_strand {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => {
                parse_errors.push(format!("{}: {LABEL} bad strand {other:?}", raw.refseq_nuc));
                continue;
            }
        };

        // 1-based fully-closed → 0-based half-open:
        //   tx_start = chr_start - 1
        //   tx_end   = chr_end       (1-based inclusive == 0-based exclusive)
        let chr_start = match raw.chr_start.parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                parse_errors.push(format!(
                    "{}: {LABEL} unparseable chr_start {:?}",
                    raw.refseq_nuc, raw.chr_start
                ));
                continue;
            }
        };
        let chr_end = match raw.chr_end.parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                parse_errors.push(format!(
                    "{}: {LABEL} unparseable chr_end {:?}",
                    raw.refseq_nuc, raw.chr_end
                ));
                continue;
            }
        };
        let tx_start = chr_start.saturating_sub(1);
        let tx_end = chr_end;

        let expected_tier = match raw.mane_status {
            "MANE Select" => Some(TranscriptTier::ManeSelect),
            "MANE Plus Clinical" => Some(TranscriptTier::ManePlusClinical),
            _ => None,
        };

        rows.push(CrossValidationRow {
            accession: raw.refseq_nuc.to_string(),
            chrom: raw.grch38_chr.to_string(),
            strand,
            tx_start,
            tx_end,
            gene_symbol: raw.symbol.to_string(),
            expected_tier,
        });
    }

    Ok(ParsedSource {
        label: LABEL,
        rows,
        excluded_built_chroms: Default::default(),
        parse_errors,
    })
}

/// Column index resolver for the MANE summary TSV header. Indices are
/// resolved from the header text so future column additions in NCBI's
/// schema don't silently break parsing.
#[derive(Debug)]
struct SummaryColumns {
    refseq_nuc: usize,
    symbol: usize,
    mane_status: usize,
    grch38_chr: usize,
    chr_start: usize,
    chr_end: usize,
    chr_strand: usize,
}

/// Borrowed-field view of a single TSV row, indexed via [`SummaryColumns`].
/// Fields are `&str` slices pointing into the caller's line buffer — no
/// allocation per row.
struct RawFields<'a> {
    refseq_nuc: &'a str,
    symbol: &'a str,
    mane_status: &'a str,
    grch38_chr: &'a str,
    chr_start: &'a str,
    chr_end: &'a str,
    chr_strand: &'a str,
}

impl SummaryColumns {
    fn resolve(header: &str) -> Result<Self> {
        let cols: Vec<&str> = header.trim().split('\t').map(str::trim).collect();
        let find = |name: &str| -> Result<usize> {
            cols.iter()
                .position(|&c| c == name)
                .with_context(|| format!("summary TSV missing `{name}` column"))
        };
        Ok(Self {
            refseq_nuc: find("RefSeq_nuc")?,
            symbol: find("symbol")?,
            mane_status: find("MANE_status")?,
            grch38_chr: find("GRCh38_chr")?,
            chr_start: find("chr_start")?,
            chr_end: find("chr_end")?,
            chr_strand: find("chr_strand")?,
        })
    }

    fn required_field_count(&self) -> usize {
        self.refseq_nuc
            .max(self.symbol)
            .max(self.mane_status)
            .max(self.grch38_chr)
            .max(self.chr_start)
            .max(self.chr_end)
            .max(self.chr_strand)
            + 1
    }

    fn extract<'a>(&self, fields: &[&'a str]) -> Option<RawFields<'a>> {
        if fields.len() < self.required_field_count() {
            return None;
        }
        Some(RawFields {
            refseq_nuc: fields[self.refseq_nuc],
            symbol: fields[self.symbol],
            mane_status: fields[self.mane_status],
            grch38_chr: fields[self.grch38_chr],
            chr_start: fields[self.chr_start],
            chr_end: fields[self.chr_end],
            chr_strand: fields[self.chr_strand],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_summary(contents: &str) -> tempfile::NamedTempFile {
        let mut file = tempfile::Builder::new()
            .suffix(".summary.tsv")
            .tempfile()
            .expect("tempfile");
        file.write_all(contents.as_bytes()).expect("write");
        file.flush().expect("flush");
        file
    }

    fn write_minimal_summary_tsv(
        refseq_nuc: &str,
        grch38_chr: &str,
        symbol: &str,
    ) -> tempfile::NamedTempFile {
        write_summary(&format!(
            "#GeneID\tEnsembl_nuc\tsymbol\tname\tRefSeq_nuc\tEnsembl_prot\tRefSeq_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n\
             123\tENSG00000099999.1\t{symbol}\tfake\t{refseq_nuc}\tENSP99999.1\tNP_999999.1\tMANE Select\t{grch38_chr}\t1000\t2000\t+\n"
        ))
    }

    #[test]
    fn parse_emits_zero_based_half_open_coordinates() {
        // Summary row carries 1-based fully-closed (chr_start=1000,
        // chr_end=2000); we expect tx_start=999, tx_end=2000.
        let summary = write_minimal_summary_tsv("NM_TEST.1", "chr1", "TESTGENE");
        let parsed = parse(summary.path()).expect("parse");
        assert!(parsed.parse_errors.is_empty(), "{:?}", parsed.parse_errors);
        assert_eq!(parsed.rows.len(), 1);
        let row = &parsed.rows[0];
        assert_eq!(row.accession, "NM_TEST.1");
        assert_eq!(row.tx_start, 999, "1-based 1000 -> 0-based 999");
        assert_eq!(
            row.tx_end, 2000,
            "1-based-inclusive 2000 == 0-based-exclusive 2000"
        );
        assert_eq!(row.strand, Strand::Plus);
        assert_eq!(row.expected_tier, Some(TranscriptTier::ManeSelect));
        assert!(parsed.excluded_built_chroms.is_empty());
        assert_eq!(parsed.label, LABEL);
    }

    #[test]
    fn parse_surfaces_bad_strand_as_parse_error() {
        // Strand "?" — must be reported, not silently dropped.
        let summary = write_summary(
            "#GeneID\tEnsembl_nuc\tsymbol\tname\tRefSeq_nuc\tEnsembl_prot\tRefSeq_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n\
             123\tENSG.1\tGENE\tfake\tNM_BAD.1\tENSP.1\tNP.1\tMANE Select\tchr1\t1000\t2000\t?\n",
        );
        let parsed = parse(summary.path()).expect("parse");
        assert!(parsed.rows.is_empty());
        assert_eq!(parsed.parse_errors.len(), 1);
        let msg = &parsed.parse_errors[0];
        assert!(msg.contains("NM_BAD.1"), "got: {msg}");
        assert!(msg.contains("bad strand"), "got: {msg}");
    }

    #[test]
    fn parse_surfaces_unparseable_coords_as_parse_error() {
        let summary = write_summary(
            "#GeneID\tEnsembl_nuc\tsymbol\tname\tRefSeq_nuc\tEnsembl_prot\tRefSeq_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n\
             123\tENSG.1\tGENE\tfake\tNM_XY.1\tENSP.1\tNP.1\tMANE Select\tchr1\tNOT_A_NUMBER\t2000\t+\n",
        );
        let parsed = parse(summary.path()).expect("parse");
        assert!(parsed.rows.is_empty());
        assert_eq!(parsed.parse_errors.len(), 1);
        assert!(parsed.parse_errors[0].contains("unparseable chr_start"));
    }
}
