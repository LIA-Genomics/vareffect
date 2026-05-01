//! Minimal VCF line parser and gzip-aware I/O helpers.
//!
//! Only the columns needed for annotation are parsed (CHROM, POS, REF, ALT,
//! INFO). The rest of the line is preserved verbatim via byte-offset
//! splicing, avoiding re-serialization of FORMAT + sample columns.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

/// Parsed VCF data line with zero-copy borrows from the raw line.
///
/// Only the fields needed for annotation are extracted. The raw line and
/// INFO byte offsets are retained so [`write_annotated_line`] can splice
/// the CSQ value into the output without re-serializing all columns.
pub struct VcfRecord<'a> {
    /// Chromosome (column 0).
    pub chrom: &'a str,
    /// VCF POS (column 1), 1-based.
    pub pos: u64,
    /// REF allele (column 3).
    pub ref_allele: &'a str,
    /// ALT alleles (column 4), split by comma.
    pub alt_alleles: Vec<&'a str>,
    /// Byte offset of the INFO field start within `raw`.
    info_start: usize,
    /// Byte offset of the INFO field end (exclusive) within `raw`.
    info_end: usize,
    /// The full original line (for pass-through and splicing).
    pub raw: &'a str,
}

/// Parse a VCF data line into a [`VcfRecord`].
///
/// Expects a tab-separated line with at least 8 columns (CHROM through
/// INFO). Lines with fewer columns or non-numeric POS are rejected.
///
/// # Errors
///
/// Returns an error if the line has fewer than 8 tab-separated columns
/// or if POS cannot be parsed as a positive integer.
pub fn parse_vcf_line(line: &str) -> Result<VcfRecord<'_>> {
    // Split into at most 9 parts: columns 0-7 and everything after (FORMAT +
    // samples preserved as part of raw).
    let mut col_starts = [0usize; 8];
    let mut col_ends = [0usize; 8];
    let mut col = 0;
    let mut prev = 0;

    for (i, b) in line.bytes().enumerate() {
        if b == b'\t' {
            if col < 8 {
                col_starts[col] = prev;
                col_ends[col] = i;
                col += 1;
            }
            prev = i + 1;
            if col >= 8 {
                break;
            }
        }
    }

    // Handle the last column if we haven't found 8 tabs yet.
    if col < 8 {
        col_starts[col] = prev;
        col_ends[col] = line.len();
        col += 1;
    }

    if col < 8 {
        bail!("VCF line has fewer than 8 columns");
    }

    let chrom = &line[col_starts[0]..col_ends[0]];
    let pos_str = &line[col_starts[1]..col_ends[1]];
    let pos: u64 = pos_str
        .parse()
        .with_context(|| format!("invalid POS: {pos_str:?}"))?;

    if pos == 0 {
        bail!("POS must be >= 1, got 0");
    }

    let ref_allele = &line[col_starts[3]..col_ends[3]];
    if ref_allele.is_empty() {
        bail!("REF allele is empty");
    }

    let alt_str = &line[col_starts[4]..col_ends[4]];
    let alt_alleles: Vec<&str> = alt_str.split(',').collect();

    let info_start = col_starts[7];
    let info_end = col_ends[7];

    Ok(VcfRecord {
        chrom,
        pos,
        ref_allele,
        alt_alleles,
        info_start,
        info_end,
        raw: line,
    })
}

/// Produce an annotated VCF line by splicing a CSQ value into the INFO field.
///
/// If the existing INFO field is `.`, it is replaced with `CSQ={csq}`.
/// Otherwise `;CSQ={csq}` is appended after the existing INFO content.
pub fn write_annotated_line(record: &VcfRecord<'_>, csq: &str) -> String {
    let info = &record.raw[record.info_start..record.info_end];

    if info == "." {
        // Replace "." with CSQ=...
        let mut out = String::with_capacity(record.raw.len() + csq.len() + 4);
        out.push_str(&record.raw[..record.info_start]);
        out.push_str("CSQ=");
        out.push_str(csq);
        out.push_str(&record.raw[record.info_end..]);
        out
    } else {
        // Append ;CSQ=... after existing INFO
        let mut out = String::with_capacity(record.raw.len() + csq.len() + 5);
        out.push_str(&record.raw[..record.info_end]);
        out.push_str(";CSQ=");
        out.push_str(csq);
        out.push_str(&record.raw[record.info_end..]);
        out
    }
}

/// Open a VCF file for reading, auto-detecting gzip by extension.
///
/// Returns a boxed `BufRead` that transparently decompresses `.gz` files.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open a VCF file for writing, auto-detecting gzip by extension.
///
/// `.gz` output uses `flate2` fast compression (standard gzip, not BGZF).
pub fn open_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufWriter::new(GzEncoder::new(
            file,
            Compression::fast(),
        ))))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_LINE: &str = "chr17\t7674221\t.\tG\tA\t.\tPASS\t.\tGT\t0/1";

    const MULTI_ALT_LINE: &str = "chr17\t7674221\trs1234\tG\tA,T\t100\tPASS\tDP=50\tGT\t0/1";

    #[test]
    fn parse_simple_snv() {
        let r = parse_vcf_line(SIMPLE_LINE).unwrap();
        assert_eq!(r.chrom, "chr17");
        assert_eq!(r.pos, 7674221);
        assert_eq!(r.ref_allele, "G");
        assert_eq!(r.alt_alleles, vec!["A"]);
        assert_eq!(&r.raw[r.info_start..r.info_end], ".");
    }

    #[test]
    fn parse_multi_allelic() {
        let r = parse_vcf_line(MULTI_ALT_LINE).unwrap();
        assert_eq!(r.alt_alleles, vec!["A", "T"]);
        assert_eq!(&r.raw[r.info_start..r.info_end], "DP=50");
    }

    #[test]
    fn parse_rejects_too_few_columns() {
        let line = "chr1\t100\t.\tA";
        assert!(parse_vcf_line(line).is_err());
    }

    #[test]
    fn parse_rejects_pos_zero() {
        let line = "chr1\t0\t.\tA\tG\t.\t.\t.";
        assert!(parse_vcf_line(line).is_err());
    }

    #[test]
    fn write_annotated_replaces_dot_info() {
        let r = parse_vcf_line(SIMPLE_LINE).unwrap();
        let out = write_annotated_line(&r, "A|missense_variant|MODERATE|TP53");
        assert!(out.contains("CSQ=A|missense_variant|MODERATE|TP53"));
        // "." should be replaced, not appended to.
        assert!(!out.contains(".;CSQ"));
        assert!(!out.contains("CSQ=."));
    }

    #[test]
    fn write_annotated_appends_to_existing_info() {
        let r = parse_vcf_line(MULTI_ALT_LINE).unwrap();
        let out = write_annotated_line(&r, "A|missense_variant|MODERATE|TP53");
        assert!(out.contains("DP=50;CSQ=A|missense_variant|MODERATE|TP53"));
    }

    #[test]
    fn write_annotated_preserves_trailing_columns() {
        let r = parse_vcf_line(SIMPLE_LINE).unwrap();
        let out = write_annotated_line(&r, "CSQ_VAL");
        // FORMAT and sample columns must be preserved.
        assert!(out.ends_with("\tGT\t0/1"));
    }
}
