//! Minimal VCF line parser and gzip-aware I/O helpers.
//!
//! Only the columns needed for annotation are parsed (CHROM, POS, REF, ALT,
//! INFO). The rest of the line is preserved verbatim via byte-offset
//! splicing, avoiding re-serialization of FORMAT + sample columns.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles_bgzf::io::multithreaded_writer::{Builder as BgzfBuilder, MultithreadedWriter};
use noodles_bgzf::io::writer::CompressionLevel;

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
/// Returns a boxed `BufRead` that transparently decompresses `.gz` files
/// (`MultiGzDecoder` handles BGZF and plain multi-member gzip). The reader is
/// `Send` so it can be moved onto the pipeline's dedicated reader thread.
pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Output sink for an annotated VCF.
///
/// `.gz` paths use a multithreaded **BGZF** writer (`noodles-bgzf`, pure-Rust
/// `zlib-rs` backend) -- a valid gzip stream that is also tabix/`bcftools`-
/// indexable. Call [`VcfWriter::finish`] when done: for BGZF it appends the
/// required EOF block and surfaces I/O errors `Drop` would otherwise discard.
pub enum VcfWriter {
    /// Uncompressed output.
    Plain(BufWriter<File>),
    /// Multithreaded BGZF (block-gzip) output.
    Bgzf(MultithreadedWriter<File>),
}

impl Write for VcfWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            VcfWriter::Plain(writer) => writer.write(buf),
            VcfWriter::Bgzf(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            VcfWriter::Plain(writer) => writer.flush(),
            VcfWriter::Bgzf(writer) => writer.flush(),
        }
    }
}

impl VcfWriter {
    /// Flush buffers and finalize the stream.
    ///
    /// For BGZF this appends the EOF block (required by tabix/htslib to detect
    /// truncation), surfacing any I/O error instead of discarding it on `Drop`.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or BGZF finalization fails (e.g. the device
    /// is full).
    pub fn finish(self) -> Result<()> {
        match self {
            VcfWriter::Plain(mut writer) => {
                writer.flush().context("flushing output")?;
            }
            VcfWriter::Bgzf(mut writer) => {
                writer.finish().context("finalizing BGZF output")?;
            }
        }
        Ok(())
    }
}

/// Open a VCF file for writing, auto-detecting gzip by extension.
///
/// `.gz` paths produce multithreaded BGZF using `bgzf_worker_count` compression
/// threads (clamped to at least 1); other paths are written uncompressed.
pub fn open_writer(path: &Path, bgzf_worker_count: usize) -> Result<VcfWriter> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;

    if path.extension().is_some_and(|e| e == "gz") {
        let workers = NonZero::new(bgzf_worker_count).unwrap_or(NonZero::<usize>::MIN);
        let writer = BgzfBuilder::default()
            .set_compression_level(CompressionLevel::FAST)
            .set_worker_count(workers)
            .build_from_writer(file);
        Ok(VcfWriter::Bgzf(writer))
    } else {
        Ok(VcfWriter::Plain(BufWriter::new(file)))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use flate2::Compression;
    use flate2::write::GzEncoder;

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

    /// Compress `content` into a single, self-contained gzip member.
    fn gzip_member(content: &str) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder
            .write_all(content.as_bytes())
            .expect("writing gzip member");
        encoder.finish().expect("finishing gzip member")
    }

    #[test]
    fn vcf_writer_bgzf_finish_round_trips() {
        let tmp = tempfile::Builder::new()
            .suffix(".vcf.gz")
            .tempfile()
            .expect("creating temp .vcf.gz");

        let mut writer = open_writer(tmp.path(), 2).expect("opening BGZF writer");
        writeln!(writer, "##fileformat=VCFv4.2").expect("write header");
        writeln!(writer, "chr1\t100\t.\tA\tG\t.\t.\tCSQ=x").expect("write record");
        writer.finish().expect("finalizing BGZF writer");

        // Reopen: succeeds only if finish() wrote a valid BGZF stream + EOF.
        let reader = open_reader(tmp.path()).expect("reopening BGZF output");
        let lines: Vec<String> = reader
            .lines()
            .collect::<std::io::Result<_>>()
            .expect("reading lines");
        assert_eq!(
            lines,
            vec!["##fileformat=VCFv4.2", "chr1\t100\t.\tA\tG\t.\t.\tCSQ=x"]
        );
    }

    #[test]
    fn write_annotated_sites_only_dot_info_replaces_cleanly() {
        // Sites-only line, INFO ".": CSQ replaces it cleanly, no stray '\r'.
        let line = "chr1\t100\t.\tA\tG\t.\tPASS\t.";
        let r = parse_vcf_line(line).unwrap();
        let out = write_annotated_line(&r, "csqval");
        assert!(out.ends_with("\tCSQ=csqval"));
        assert!(!out.contains(";CSQ"));
        assert!(!out.contains('\r'));
    }

    #[test]
    fn trailing_cr_corrupts_sites_only_info_if_not_stripped() {
        // Why annotate::run strips '\r' first: an un-stripped CR breaks INFO
        // "." detection and mis-splices CSQ after the carriage return.
        let with_cr = "chr1\t100\t.\tA\tG\t.\tPASS\t.\r";
        let r = parse_vcf_line(with_cr).unwrap();
        let out = write_annotated_line(&r, "csqval");
        assert!(
            out.contains(".\r;CSQ=csqval"),
            "unstripped CR corrupts the INFO field"
        );
    }

    #[test]
    fn open_reader_reads_all_members_of_multi_member_gzip() {
        // BGZF (the .vcf.gz format from bgzip/bcftools/tabix) is a stream of
        // concatenated gzip members. Simulate it with two independently
        // finished members so a single-member decoder would stop after the
        // first. The #CHROM line and variant data live in the SECOND member,
        // mirroring a real GRCh38 VCF whose contig-heavy header overflows the
        // first BGZF block.
        let member1 = "##fileformat=VCFv4.2\n##contig=<ID=chr16,length=90338345>\n";
        let member2 =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr16\t47435\t.\tC\tA\t.\t.\t.\n";

        let mut bytes = gzip_member(member1);
        bytes.extend_from_slice(&gzip_member(member2));

        let tmp = tempfile::Builder::new()
            .suffix(".vcf.gz")
            .tempfile()
            .expect("creating temp .vcf.gz");
        std::fs::write(tmp.path(), &bytes).expect("writing multi-member gzip");

        let reader = open_reader(tmp.path()).expect("opening multi-member gzip");
        let lines: Vec<String> = reader
            .lines()
            .collect::<std::io::Result<_>>()
            .expect("reading lines");

        // All four lines from BOTH members must be present. A single-member
        // GzDecoder would return only the two lines from member 1.
        assert_eq!(
            lines.len(),
            4,
            "expected lines from both gzip members, got {lines:?}"
        );
        assert!(
            lines.iter().any(|l| l.starts_with("#CHROM")),
            "second-member #CHROM line missing"
        );
        assert!(
            lines.iter().any(|l| l.starts_with("chr16\t47435")),
            "second-member data line missing"
        );
    }
}
