//! Cache naming and verification of a pinned RefSeq or GENCODE GFF3.
//!
//! `setup` downloads the GFF3 named by `[vareffect.refseq]` or
//! `[vareffect.gencode]` into `raw_dir` under a release-qualified name, then
//! refuses to build from it unless its SHA-256 equals the pin, its header
//! names the pinned release, and its header genome build matches
//! `fasta_build`.

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result, bail};
use vareffect::TranscriptSource;

use crate::common::sha256_file;
use crate::config::{PinnedSource, SourceName};
use crate::vcf;

/// Release and genome build read from a source GFF3's header comments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceHeader {
    /// Source release (e.g. `"RS_2025_08"`, `"50"`).
    pub release: String,
    /// Genome build (e.g. `"GRCh38.p14"`, `"GRCh38"`).
    pub genome_build: String,
}

/// Cache filename in `raw_dir` for a pinned source GFF3.
///
/// # Arguments
///
/// * `source` -- RefSeq or GENCODE.
/// * `genome_build` -- The configured `fasta_build` (e.g. `"GRCh38"`).
/// * `release` -- The pinned release.
///
/// # Returns
///
/// `refseq_<build>_<release>_genomic.gff.gz` or
/// `gencode_<build>_v<release>_annotation.gff3.gz`.
///
/// # Errors
///
/// Returns an error for MANE (its cache name comes from `[mane].input`) or
/// if the build or release is empty or holds a character other than an
/// ASCII letter, digit, `.`, `_` or `-`, which keeps the name inside
/// `raw_dir`.
pub fn cache_filename(source: SourceName, genome_build: &str, release: &str) -> Result<String> {
    for (key, value) in [("fasta_build", genome_build), ("release", release)] {
        let allowed = |c: char| c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-');
        if value.is_empty() || !value.chars().all(allowed) {
            bail!(
                "{key} {value:?} cannot name a cache file: use ASCII letters, digits, '.', '_' \
                 or '-'"
            );
        }
    }
    match source {
        SourceName::Mane => bail!("the MANE GFF3 cache name comes from [mane].input"),
        SourceName::Refseq => Ok(format!("refseq_{genome_build}_{release}_genomic.gff.gz")),
        SourceName::Gencode => Ok(format!(
            "gencode_{genome_build}_v{release}_annotation.gff3.gz"
        )),
    }
}

/// Verify a downloaded source GFF3 against its pin before building from it.
///
/// # Arguments
///
/// * `path` -- The cached GFF3.
/// * `source` -- RefSeq or GENCODE; selects the header format.
/// * `pinned` -- The configured URL, release and SHA-256.
/// * `genome_build` -- The configured `fasta_build`.
///
/// # Returns
///
/// The header read from the file.
///
/// # Errors
///
/// Returns an error if the file cannot be read, its SHA-256 differs from
/// `pinned.sha256`, its header lacks the release or genome build, the
/// release differs from `pinned.release`, or the build names another
/// assembly than `genome_build` (patch suffixes such as `.p14` ignored).
pub fn verify(
    path: &Path,
    source: SourceName,
    pinned: &PinnedSource,
    genome_build: &str,
) -> Result<SourceHeader> {
    let actual = sha256_file(path)?;
    if !actual.eq_ignore_ascii_case(pinned.sha256.trim()) {
        bail!(
            "{} has SHA-256 {actual}, but the config pins {}: delete the file to download it \
             again, or correct the pinned sha256",
            path.display(),
            pinned.sha256.trim(),
        );
    }
    let reader = vcf::open_reader(path)?;
    let header = parse_header(reader.lines(), source)
        .with_context(|| format!("reading the header of {}", path.display()))?;
    if header.release != pinned.release {
        bail!(
            "{} is release {}, but the config pins release {}",
            path.display(),
            header.release,
            pinned.release,
        );
    }
    if !build_matches(&header.genome_build, genome_build) {
        bail!(
            "{} is annotated on {}, but fasta_build is {genome_build}",
            path.display(),
            header.genome_build,
        );
    }
    Ok(header)
}

/// `true` if a header build and the configured build name one assembly,
/// ignoring a patch suffix on either side (`GRCh38.p14` matches `GRCh38`
/// both ways). Patches add scaffolds but keep the primary-assembly
/// coordinates, which are what the stores annotate.
fn build_matches(header_build: &str, genome_build: &str) -> bool {
    fn assembly(build: &str) -> &str {
        build.split_once('.').map_or(build, |(base, _)| base)
    }
    assembly(header_build) == assembly(genome_build)
}

/// Read the release and genome build from the leading `#` lines of a GFF3.
///
/// RefSeq: `#!annotation-source NCBI RefSeq GCF_000001405.40-RS_2025_08`
/// and `#!genome-build GRCh38.p14`. GENCODE: `#description: evidence-based
/// annotation of the human genome (GRCh38), version 50 (Ensembl 116)`.
/// Reading stops at the first line that is not a comment.
///
/// # Errors
///
/// Returns an error if a line cannot be read, the source is MANE, or the
/// release or genome build is not found.
fn parse_header(
    lines: impl Iterator<Item = std::io::Result<String>>,
    source: SourceName,
) -> Result<SourceHeader> {
    let mut release = None;
    let mut genome_build = None;
    for line in lines {
        let line = line.context("reading a header line")?;
        if !line.starts_with('#') {
            break;
        }
        match source {
            SourceName::Mane => bail!("MANE GFF3 headers are not verified"),
            SourceName::Refseq => {
                if let Some(rest) = line.strip_prefix("#!annotation-source ") {
                    // The last word is `<assembly accession>-<release>`.
                    release = rest
                        .split_whitespace()
                        .last()
                        .and_then(|word| word.split_once('-'))
                        .map(|(_, release)| release.to_string());
                } else if let Some(rest) = line.strip_prefix("#!genome-build ") {
                    genome_build = Some(rest.trim().to_string());
                }
            }
            SourceName::Gencode => {
                if let Some(rest) = line.strip_prefix("#description:") {
                    release = rest.split_once("version ").map(|(_, tail)| {
                        tail.chars()
                            .take_while(char::is_ascii_alphanumeric)
                            .collect::<String>()
                    });
                    genome_build = rest
                        .split_once('(')
                        .and_then(|(_, tail)| tail.split_once(')'))
                        .map(|(build, _)| build.to_string());
                }
            }
        }
    }
    match (release, genome_build) {
        (Some(release), Some(genome_build)) if !release.is_empty() && !genome_build.is_empty() => {
            Ok(SourceHeader {
                release,
                genome_build,
            })
        }
        _ => bail!(
            "no {} release and genome build in the header comments",
            TranscriptSource::from(source)
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;

    use flate2::Compression;
    use flate2::write::GzEncoder;
    use sha2::{Digest, Sha256};
    use tempfile::TempDir;

    /// Header lines as in `GCF_000001405.40_GRCh38.p14_genomic.gff.gz`.
    const REFSEQ_HEADER: &str = "##gff-version 3\n#!gff-spec-version 1.21\n\
        #!processor NCBI annotwriter\n#!genome-build GRCh38.p14\n\
        #!genome-build-accession NCBI_Assembly:GCF_000001405.40\n\
        #!annotation-date 08/01/2025\n\
        #!annotation-source NCBI RefSeq GCF_000001405.40-RS_2025_08\n\
        ##sequence-region NC_000001.11 1 248956422\n\
        NC_000001.11\tRefSeq\tregion\t1\t248956422\t.\t+\t.\tID=NC_000001.11:1..248956422\n";

    /// Header lines as in `gencode.v50.annotation.gff3.gz`.
    const GENCODE_HEADER: &str = "##gff-version 3\n\
        #description: evidence-based annotation of the human genome (GRCh38), version 50 \
        (Ensembl 116)\n#provider: GENCODE\n#contact: gencode-help@ebi.ac.uk\n#format: gff3\n\
        #date: 2026-04-08\n##sequence-region chr1 1 248956422\n\
        chr1\tHAVANA\tgene\t11121\t24894\t.\t+\t.\tID=ENSG00000290825.2\n";

    fn lines(text: &str) -> impl Iterator<Item = std::io::Result<String>> + '_ {
        text.lines().map(|line| Ok(line.to_string()))
    }

    fn header(release: &str, genome_build: &str) -> SourceHeader {
        SourceHeader {
            release: release.to_string(),
            genome_build: genome_build.to_string(),
        }
    }

    /// Gzip `text` into `dir/name` and return its path and hex SHA-256.
    fn write_gz(dir: &Path, name: &str, text: &str) -> (std::path::PathBuf, String) {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(text.as_bytes()).expect("gzip");
        let bytes = encoder.finish().expect("gzip finish");
        let path = dir.join(name);
        std::fs::write(&path, &bytes).expect("write");
        (path, hex::encode(Sha256::digest(&bytes)))
    }

    fn pinned(release: &str, sha256: &str) -> PinnedSource {
        PinnedSource {
            url: "https://example.com/source.gff.gz".to_string(),
            release: release.to_string(),
            sha256: sha256.to_string(),
        }
    }

    #[test]
    fn parse_header_reads_refseq_release_and_build() {
        let parsed = parse_header(lines(REFSEQ_HEADER), SourceName::Refseq).expect("parse");
        assert_eq!(parsed, header("RS_2025_08", "GRCh38.p14"));
    }

    #[test]
    fn parse_header_reads_gencode_release_and_build() {
        let parsed = parse_header(lines(GENCODE_HEADER), SourceName::Gencode).expect("parse");
        assert_eq!(parsed, header("50", "GRCh38"));
    }

    #[test]
    fn parse_header_fails_without_release_or_build() {
        assert!(parse_header(lines(GENCODE_HEADER), SourceName::Refseq).is_err());
        assert!(parse_header(lines(REFSEQ_HEADER), SourceName::Gencode).is_err());
        assert!(parse_header(lines(REFSEQ_HEADER), SourceName::Mane).is_err());
    }

    #[test]
    fn parse_header_stops_at_the_first_record() {
        let late = "##gff-version 3\nchr1\tHAVANA\tgene\t1\t2\t.\t+\t.\tID=g\n\
            #description: annotation of the human genome (GRCh38), version 50\n";
        assert!(parse_header(lines(late), SourceName::Gencode).is_err());
    }

    #[test]
    fn build_matches_ignores_patch_suffixes_only() {
        assert!(build_matches("GRCh38", "GRCh38"));
        assert!(build_matches("GRCh38.p14", "GRCh38"));
        assert!(build_matches("GRCh38", "GRCh38.p14"));
        assert!(!build_matches("GRCh37.p13", "GRCh38"));
        assert!(!build_matches("GRCh380", "GRCh38"));
    }

    #[test]
    fn cache_filename_names_source_build_and_release() {
        assert_eq!(
            cache_filename(SourceName::Refseq, "GRCh38", "RS_2025_08").expect("refseq"),
            "refseq_GRCh38_RS_2025_08_genomic.gff.gz"
        );
        assert_eq!(
            cache_filename(SourceName::Gencode, "GRCh38", "50").expect("gencode"),
            "gencode_GRCh38_v50_annotation.gff3.gz"
        );
        assert!(cache_filename(SourceName::Mane, "GRCh38", "1.5").is_err());
    }

    #[test]
    fn cache_filename_refuses_path_characters() {
        assert!(cache_filename(SourceName::Refseq, "GRCh38", "../x").is_err());
        assert!(cache_filename(SourceName::Gencode, "GRCh38", "").is_err());
        assert!(cache_filename(SourceName::Gencode, "GRCh 38", "50").is_err());
    }

    #[test]
    fn verify_accepts_a_file_matching_its_pin() {
        let tmp = TempDir::new().expect("tempdir");
        let (path, sha) = write_gz(tmp.path(), "gencode.gff3.gz", GENCODE_HEADER);
        let parsed = verify(
            &path,
            SourceName::Gencode,
            &pinned("50", &sha.to_uppercase()),
            "GRCh38",
        )
        .expect("verify");
        assert_eq!(parsed, header("50", "GRCh38"));
    }

    #[test]
    fn verify_refuses_a_sha256_mismatch() {
        let tmp = TempDir::new().expect("tempdir");
        let (path, _) = write_gz(tmp.path(), "refseq.gff.gz", REFSEQ_HEADER);
        let error = verify(
            &path,
            SourceName::Refseq,
            &pinned("RS_2025_08", &"0".repeat(64)),
            "GRCh38",
        )
        .expect_err("sha mismatch");
        assert!(error.to_string().contains("SHA-256"), "{error}");
    }

    #[test]
    fn verify_refuses_a_release_mismatch() {
        let tmp = TempDir::new().expect("tempdir");
        let (path, sha) = write_gz(tmp.path(), "refseq.gff.gz", REFSEQ_HEADER);
        let error = verify(
            &path,
            SourceName::Refseq,
            &pinned("RS_2024_08", &sha),
            "GRCh38",
        )
        .expect_err("release mismatch");
        assert!(error.to_string().contains("release RS_2024_08"), "{error}");
    }

    #[test]
    fn verify_refuses_a_genome_build_mismatch() {
        let tmp = TempDir::new().expect("tempdir");
        let (path, sha) = write_gz(tmp.path(), "gencode.gff3.gz", GENCODE_HEADER);
        let error = verify(&path, SourceName::Gencode, &pinned("50", &sha), "GRCh37")
            .expect_err("build mismatch");
        assert!(
            error.to_string().contains("fasta_build is GRCh37"),
            "{error}"
        );
    }
}
