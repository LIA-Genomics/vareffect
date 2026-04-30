//! Cross-validation dispatcher for built transcript models.
//!
//! Compares each built [`TranscriptModel`] against an independently-parsed
//! tabular source (MANE summary TSV for GRCh38, UCSC ncbiRefSeq +
//! ncbiRefSeqSelect for GRCh37) and fails the build on any mismatch in
//! chrom, strand, position, gene symbol, or tier.
//!
//! Backend contract:
//!
//! 1. Each backend parses its source into a [`ParsedSource`] holding
//!    rows in 0-based half-open coordinates (matching `TranscriptModel`).
//!    Backends convert from their native form (MANE summary is 1-based
//!    fully-closed; UCSC genePred is already 0-based half-open).
//! 2. Structural-parse failures (bad strand, unparseable coords, short
//!    rows) are surfaced as `parse_errors` strings — never silently
//!    dropped — so a malformed row doesn't masquerade as a built
//!    transcript that's "missing from source" further down.
//! 3. Each backend may emit `excluded_built_chroms` — chromosome names
//!    whose built transcripts skip the inverse "every built transcript
//!    is in the source" check. UCSC adds `chrM` because UCSC `hg19`
//!    chrM is NC_001807 and NCBI GRCh37 chrMT is NC_012920.1.
//! 4. The dispatcher sorts rows by accession before iterating so
//!    mismatch reports are reproducible across UCSC table refreshes.
//! 5. Comparison is exact: strict string-equality on the versioned
//!    accession; coordinate equality with no tolerance.
//!
//! On any non-empty mismatch list the dispatcher fails the build with
//! the first 20 entries plus the total count.

mod mane_summary;
mod ucsc_refseq;

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::PathBuf;

use anyhow::{Result, bail};
use vareffect::chrom::{Assembly, is_patch_sequence, refseq_to_ucsc};
use vareffect::{Strand, TranscriptModel, TranscriptTier};

/// Pointer to a tabular cross-validation source. Owns its `PathBuf`s so
/// callers can build the variant from `raw_dir.join(...)` temporaries
/// inside a match arm without lifetime threading.
#[derive(Debug, Clone)]
pub(crate) enum CrossValidationSource {
    /// NCBI MANE summary TSV (`MANE.GRCh38.vX.X.summary.txt.gz`). Single
    /// file with a `#GeneID...` header. GRCh38-only.
    ManeSummary { summary_path: PathBuf },
    /// UCSC `ncbiRefSeq.txt.gz` (genePred-extended coordinates) plus
    /// the sibling `ncbiRefSeqSelect.txt.gz` (filter to RefSeq Select).
    /// Used for GRCh37 because no NCBI tabular RefSeq Select release
    /// exists.
    UcscRefseqSelect {
        coords_path: PathBuf,
        select_path: PathBuf,
    },
}

/// Unified per-transcript validation record. Each backend produces a
/// `Vec<CrossValidationRow>` after format-specific parsing; the shared
/// comparison loop runs against this representation.
#[derive(Debug, Clone)]
pub(super) struct CrossValidationRow {
    /// Versioned RefSeq accession (e.g., `"NM_006772.2"`). Strict
    /// string-equality match against `TranscriptModel::accession`.
    pub accession: String,
    /// Source-native chromosome name, normalized to UCSC form by the
    /// comparison loop. MANE summary may use `"6"` / `"chr6"` /
    /// `"NC_000006.12"` / `"NW_*"` interchangeably; UCSC always uses
    /// `chr*` form.
    pub chrom: String,
    /// Strand from the source.
    pub strand: Strand,
    /// Transcript start, 0-based inclusive (matching
    /// `TranscriptModel::tx_start`).
    pub tx_start: u64,
    /// Transcript end, 0-based exclusive (matching
    /// `TranscriptModel::tx_end`).
    pub tx_end: u64,
    /// Gene symbol. Empty when the source doesn't include one.
    pub gene_symbol: String,
    /// Expected tier when the source explicitly states one (MANE
    /// summary's `MANE_status`, or UCSC's Select-list membership).
    pub expected_tier: Option<TranscriptTier>,
}

/// Output of a backend parse: a static label, the rows, the set of
/// built-side chromosomes the parser intentionally excluded, and any
/// structural-parse errors encountered while reading the source. The
/// dispatcher promotes `parse_errors` straight into the mismatch list
/// so a malformed row reports its actual reason rather than appearing
/// as "not present in the source" downstream.
pub(super) struct ParsedSource {
    pub label: &'static str,
    pub rows: Vec<CrossValidationRow>,
    pub excluded_built_chroms: HashSet<String>,
    pub parse_errors: Vec<String>,
}

/// Run cross-validation against the configured source.
pub(super) fn cross_validate(
    transcripts: &[TranscriptModel],
    patch_aliases: &HashMap<String, String>,
    assembly: Assembly,
    source: &CrossValidationSource,
) -> Result<()> {
    let parsed = match source {
        CrossValidationSource::ManeSummary { summary_path } => mane_summary::parse(summary_path)?,
        CrossValidationSource::UcscRefseqSelect {
            coords_path,
            select_path,
        } => ucsc_refseq::parse(coords_path, select_path)?,
    };

    compare(transcripts, patch_aliases, assembly, parsed)
}

/// Shared comparison loop.
fn compare(
    transcripts: &[TranscriptModel],
    patch_aliases: &HashMap<String, String>,
    assembly: Assembly,
    parsed: ParsedSource,
) -> Result<()> {
    let ParsedSource {
        label: source_label,
        mut rows,
        excluded_built_chroms,
        parse_errors,
    } = parsed;

    // Sort up-front so the "first 20 shown" bail! payload is reproducible
    // across rebuilds — UCSC's table refresh order is not stable.
    rows.sort_by(|a, b| a.accession.cmp(&b.accession));

    let by_acc: HashMap<&str, &TranscriptModel> = transcripts
        .iter()
        .map(|t| (t.accession.as_str(), t))
        .collect();

    let mut mismatches: Vec<String> = parse_errors;

    // Invariant: every accession the parser emitted gets recorded in
    // `source_seen` *before* any chrom/coord/tier check, so the inverse
    // "missing from source" pass below cannot false-positive on rows
    // that the comparison loop chose to flag for some other reason.
    let mut source_seen: HashSet<String> = HashSet::new();
    let mut rows_checked: usize = 0;
    let mut unknown_patch_aliases: BTreeSet<String> = BTreeSet::new();

    for row in &rows {
        source_seen.insert(row.accession.clone());

        let Some(tx) = by_acc.get(row.accession.as_str()) else {
            mismatches.push(format!(
                "{}: no matching transcript in built store \
                 (version drift between {} and the GFF3?)",
                row.accession, source_label
            ));
            continue;
        };
        rows_checked += 1;

        // --- Chromosome ---
        let raw = row.chrom.trim();
        let normalized = normalize_chrom(raw, assembly);
        let built_is_patch = is_patch_sequence(&tx.chrom);
        let source_is_patch =
            raw.starts_with("NW_") || raw.starts_with("NT_") || is_patch_sequence(raw);
        match (built_is_patch, source_is_patch) {
            (true, true) => {
                // Patch contig — patch_aliases is keyed on RefSeq form
                // (`NW_*`/`NT_*`); UCSC patch names appear as the
                // *value*, so compare directly when the source is
                // already UCSC-native.
                let resolved_match = if raw.starts_with("NW_") || raw.starts_with("NT_") {
                    match patch_aliases.get(raw) {
                        Some(expected_ucsc) => {
                            if expected_ucsc == &tx.chrom {
                                Ok(true)
                            } else {
                                Err(format!(
                                    "{}: patch chrom mismatch -- {} alias={}, built={}",
                                    row.accession, source_label, expected_ucsc, tx.chrom
                                ))
                            }
                        }
                        None => {
                            unknown_patch_aliases.insert(raw.to_string());
                            Ok(true)
                        }
                    }
                } else if raw == tx.chrom {
                    Ok(true)
                } else {
                    Err(format!(
                        "{}: patch chrom mismatch -- {}={}, built={}",
                        row.accession, source_label, raw, tx.chrom
                    ))
                };
                if let Err(msg) = resolved_match {
                    mismatches.push(msg);
                }
            }
            (false, false) => {
                if normalized != tx.chrom {
                    mismatches.push(format!(
                        "{}: chrom mismatch -- {}={}, built={}",
                        row.accession, source_label, normalized, tx.chrom
                    ));
                }
            }
            (built_patch, source_patch) => {
                mismatches.push(format!(
                    "{}: patch/primary asymmetry -- {}={} (patch={}), \
                     built={} (patch={})",
                    row.accession, source_label, raw, source_patch, tx.chrom, built_patch
                ));
            }
        }

        // --- Strand ---
        if row.strand != tx.strand {
            mismatches.push(format!(
                "{}: strand mismatch -- {}={:?}, built={:?}",
                row.accession, source_label, row.strand, tx.strand
            ));
        }

        // --- Coordinates (0-based half-open everywhere) ---
        if row.tx_start != tx.tx_start {
            mismatches.push(format!(
                "{}: tx_start mismatch -- {}={}, built={}",
                row.accession, source_label, row.tx_start, tx.tx_start
            ));
        }
        if row.tx_end != tx.tx_end {
            mismatches.push(format!(
                "{}: tx_end mismatch -- {}={}, built={}",
                row.accession, source_label, row.tx_end, tx.tx_end
            ));
        }

        // --- Gene symbol ---
        if !row.gene_symbol.is_empty() && row.gene_symbol != tx.gene_symbol {
            mismatches.push(format!(
                "{}: gene_symbol mismatch -- {}={}, built={}",
                row.accession, source_label, row.gene_symbol, tx.gene_symbol
            ));
        }

        // --- Tier ---
        if let Some(expected) = row.expected_tier
            && expected != tx.tier
        {
            mismatches.push(format!(
                "{}: tier mismatch -- {}={:?}, built={:?}",
                row.accession, source_label, expected, tx.tier
            ));
        }
    }

    // Inverse check: every built transcript must be named by the source,
    // except those on chromosomes the parser excluded (chrM under UCSC).
    for tx in transcripts {
        if excluded_built_chroms.contains(&tx.chrom) {
            continue;
        }
        if !source_seen.contains(&tx.accession) {
            mismatches.push(format!(
                "{}: built transcript is not present in the {} source",
                tx.accession, source_label
            ));
        }
    }

    if !unknown_patch_aliases.is_empty() {
        let aliases = unknown_patch_aliases
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(", ");
        tracing::warn!(
            store = "transcript_models",
            count = unknown_patch_aliases.len(),
            "{source_label} references patch accessions missing from \
             patch_chrom_aliases.csv -- the assembly report may be out of \
             date; aliases: [{aliases}]",
        );
    }

    eprintln!(
        "  {source_label} cross-validation: {} rows checked, {} mismatches, {} unknown patch aliases",
        rows_checked,
        mismatches.len(),
        unknown_patch_aliases.len(),
    );

    if !mismatches.is_empty() {
        let shown = mismatches
            .iter()
            .take(20)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join("\n");
        bail!(
            "{source_label} cross-validation failed with {} mismatches (first 20 shown):\n{shown}",
            mismatches.len(),
        );
    }
    Ok(())
}

/// Normalize a source-native chromosome value to UCSC form.
///
/// Handles every shape observed in the supported backends:
///   `"6"` -> `"chr6"`, `"chr6"` -> `"chr6"`, `"NC_000006.12"` -> `"chr6"`.
///
/// `NC_*` translation respects the assembly the validator was called
/// with — GRCh37 and GRCh38 use different patch versions for every
/// primary chromosome except chrM.
fn normalize_chrom(raw: &str, assembly: Assembly) -> String {
    let trimmed = raw.trim();
    if trimmed.starts_with("chr") {
        return trimmed.to_string();
    }
    if trimmed.starts_with("NC_") {
        return refseq_to_ucsc(assembly, trimmed).to_string();
    }
    if !trimmed.is_empty() {
        return format!("chr{trimmed}");
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    use vareffect::{Biotype, CdsSegment, Exon};

    fn write_temp(suffix: &str, contents: &str) -> NamedTempFile {
        let mut file = tempfile::Builder::new()
            .suffix(suffix)
            .tempfile()
            .expect("tempfile");
        file.write_all(contents.as_bytes()).expect("write");
        file.flush().expect("flush");
        file
    }

    fn fake_tx(accession: &str, chrom: &str, tx_start: u64, tx_end: u64) -> TranscriptModel {
        TranscriptModel {
            accession: accession.to_string(),
            protein_accession: Some("NP_999999.1".into()),
            gene_symbol: "TESTGENE".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: chrom.to_string(),
            strand: Strand::Plus,
            tx_start,
            tx_end,
            cds_genomic_start: Some(tx_start + 10),
            cds_genomic_end: Some(tx_end - 10),
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: tx_start,
                genomic_end: tx_end,
            }],
            cds_segments: vec![CdsSegment {
                exon_index: 0,
                genomic_start: tx_start + 10,
                genomic_end: tx_end - 10,
                phase: 0,
            }],
            tier: TranscriptTier::RefSeqSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 1,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    fn patch_transcript(accession: &str, chrom: &str, gene_symbol: &str) -> TranscriptModel {
        let mut tx = fake_tx(accession, chrom, 999, 2000);
        tx.gene_symbol = gene_symbol.to_string();
        tx.tier = TranscriptTier::ManeSelect;
        tx.hgnc_id = Some("HGNC:99999".into());
        tx
    }

    fn write_minimal_summary_tsv(
        refseq_nuc: &str,
        grch38_chr: &str,
        symbol: &str,
    ) -> NamedTempFile {
        write_temp(
            ".summary.tsv",
            &format!(
                "#GeneID\tEnsembl_nuc\tsymbol\tname\tRefSeq_nuc\tEnsembl_prot\tRefSeq_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n\
                 123\tENSG00000099999.1\t{symbol}\tfake\t{refseq_nuc}\tENSP99999.1\tNP_999999.1\tMANE Select\t{grch38_chr}\t1000\t2000\t+\n"
            ),
        )
    }

    fn ncbi_refseq_row(
        name: &str,
        chrom: &str,
        strand: &str,
        start: u64,
        end: u64,
        name2: &str,
    ) -> String {
        // genePred-extended: bin, name, chrom, strand, txStart, txEnd,
        // cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score,
        // name2, cdsStartStat, cdsEndStat, exonFrames.
        format!(
            "0\t{name}\t{chrom}\t{strand}\t{start}\t{end}\t{start}\t{end}\t1\t{start},\t{end},\t0\t{name2}\tcmpl\tcmpl\t0,\n"
        )
    }

    #[test]
    fn validator_resolves_patch_alias_hit() {
        let transcripts = vec![patch_transcript(
            "NM_PATCH1.1",
            "chr22_KI270879v1_alt",
            "PATCH1",
        )];
        let summary = write_minimal_summary_tsv("NM_PATCH1.1", "NT_187633.1", "PATCH1");
        let mut aliases = HashMap::new();
        aliases.insert(
            "NT_187633.1".to_string(),
            "chr22_KI270879v1_alt".to_string(),
        );
        let source = CrossValidationSource::ManeSummary {
            summary_path: summary.path().to_path_buf(),
        };
        cross_validate(&transcripts, &aliases, Assembly::GRCh38, &source)
            .expect("validator should accept a resolved patch alias");
    }

    #[test]
    fn validator_rejects_patch_alias_mismatch() {
        let transcripts = vec![patch_transcript(
            "NM_PATCH2.1",
            "chr22_KI270879v1_alt",
            "PATCH2",
        )];
        let summary = write_minimal_summary_tsv("NM_PATCH2.1", "NT_187633.1", "PATCH2");
        let mut aliases = HashMap::new();
        aliases.insert("NT_187633.1".to_string(), "chr9_KN196479v1_fix".to_string());
        let source = CrossValidationSource::ManeSummary {
            summary_path: summary.path().to_path_buf(),
        };
        let err = cross_validate(&transcripts, &aliases, Assembly::GRCh38, &source)
            .expect_err("validator should reject a patch alias mismatch");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("patch chrom mismatch"),
            "expected patch chrom mismatch, got {msg}"
        );
    }

    #[test]
    fn validator_rejects_built_patch_summary_primary() {
        let transcripts = vec![patch_transcript(
            "NM_ASYMM1.1",
            "chr22_KI270879v1_alt",
            "ASYMM1",
        )];
        let summary = write_minimal_summary_tsv("NM_ASYMM1.1", "NC_000022.11", "ASYMM1");
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::ManeSummary {
            summary_path: summary.path().to_path_buf(),
        };
        let err = cross_validate(&transcripts, &aliases, Assembly::GRCh38, &source)
            .expect_err("validator should reject patch/primary asymmetry");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("patch/primary asymmetry"),
            "expected asymmetry mismatch, got {msg}"
        );
    }

    #[test]
    fn validator_fails_on_mane_summary_version_drift() {
        // Built has NM_TEST.2, summary has NM_TEST.1 — strict equality
        // on the versioned accession; version drift is fatal.
        let mut tx = patch_transcript("NM_TEST.2", "chr1", "TEST");
        tx.chrom = "chr1".into();
        let summary = write_minimal_summary_tsv("NM_TEST.1", "chr1", "TEST");
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::ManeSummary {
            summary_path: summary.path().to_path_buf(),
        };
        let err = cross_validate(&[tx], &aliases, Assembly::GRCh38, &source)
            .expect_err("version drift must fail the build");
        let msg = format!("{err:#}");
        assert!(msg.contains("NM_TEST.1"));
        assert!(msg.contains("no matching transcript"));
    }

    #[test]
    fn ucsc_passes_for_matching_grch37_build() {
        let select = write_temp(".select.txt", "1\tNM_A.1\n1\tNM_B.1\n");
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &(ncbi_refseq_row("NM_A.1", "chr1", "+", 100, 200, "A")
                + &ncbi_refseq_row("NM_B.1", "chr2", "-", 300, 500, "B")),
        );
        let mut tx_b = fake_tx("NM_B.1", "chr2", 300, 500);
        tx_b.gene_symbol = "B".into();
        tx_b.strand = Strand::Minus;
        let mut tx_a = fake_tx("NM_A.1", "chr1", 100, 200);
        tx_a.gene_symbol = "A".into();
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::UcscRefseqSelect {
            coords_path: coords.path().to_path_buf(),
            select_path: select.path().to_path_buf(),
        };
        cross_validate(&[tx_a, tx_b], &aliases, Assembly::GRCh37, &source)
            .expect("matching build should pass cross-validation");
    }

    #[test]
    fn ucsc_skips_chrm_in_inverse_check() {
        // Built store has an MT-CO1 transcript on chrM. Source has only
        // the primary chr1 row. The inverse check would normally fail,
        // but `excluded_built_chroms` contains chrM.
        let select = write_temp(".select.txt", "1\tNM_PRIMARY.1\n1\tNM_MITO.1\n");
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &(ncbi_refseq_row("NM_PRIMARY.1", "chr1", "+", 100, 200, "GENE")
                + &ncbi_refseq_row("NM_MITO.1", "chrM", "+", 50, 150, "MT-CO1")),
        );
        let mut tx_primary = fake_tx("NM_PRIMARY.1", "chr1", 100, 200);
        tx_primary.gene_symbol = "GENE".into();
        let mut tx_mito = fake_tx("NM_MITO.1", "chrM", 50, 150);
        tx_mito.gene_symbol = "MT-CO1".into();
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::UcscRefseqSelect {
            coords_path: coords.path().to_path_buf(),
            select_path: select.path().to_path_buf(),
        };
        cross_validate(&[tx_primary, tx_mito], &aliases, Assembly::GRCh37, &source)
            .expect("chrM transcripts must be excluded from inverse check");
    }

    #[test]
    fn ucsc_fails_on_coordinate_drift() {
        let select = write_temp(".select.txt", "1\tNM_DRIFT.1\n");
        // Source says tx_end = 500; built has 600.
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &ncbi_refseq_row("NM_DRIFT.1", "chr1", "+", 100, 500, "DRIFT"),
        );
        let mut tx = fake_tx("NM_DRIFT.1", "chr1", 100, 600);
        tx.gene_symbol = "DRIFT".into();
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::UcscRefseqSelect {
            coords_path: coords.path().to_path_buf(),
            select_path: select.path().to_path_buf(),
        };
        let err = cross_validate(&[tx], &aliases, Assembly::GRCh37, &source)
            .expect_err("coordinate drift must fail");
        let msg = format!("{err:#}");
        assert!(msg.contains("tx_end mismatch"), "got: {msg}");
    }

    #[test]
    fn ucsc_fails_on_version_drift() {
        let select = write_temp(".select.txt", "1\tNM_TEST.1\n");
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &ncbi_refseq_row("NM_TEST.1", "chr1", "+", 100, 200, "TEST"),
        );
        let mut tx = fake_tx("NM_TEST.2", "chr1", 100, 200);
        tx.gene_symbol = "TEST".into();
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::UcscRefseqSelect {
            coords_path: coords.path().to_path_buf(),
            select_path: select.path().to_path_buf(),
        };
        let err = cross_validate(&[tx], &aliases, Assembly::GRCh37, &source)
            .expect_err("version drift must fail");
        let msg = format!("{err:#}");
        assert!(msg.contains("no matching transcript"), "got: {msg}");
    }

    #[test]
    fn ucsc_surfaces_bad_strand_as_parse_error_not_missing_from_source() {
        // The source row has an invalid strand. The parser must surface
        // it as a structural-parse error, not let it slip through and
        // re-appear later as "built transcript not in source".
        let select = write_temp(".select.txt", "1\tNM_BADSTRAND.1\n");
        let coords = write_temp(
            ".ncbiRefSeq.txt",
            &ncbi_refseq_row("NM_BADSTRAND.1", "chr1", "?", 100, 200, "GENE"),
        );
        let mut tx = fake_tx("NM_BADSTRAND.1", "chr1", 100, 200);
        tx.gene_symbol = "GENE".into();
        let aliases: HashMap<String, String> = HashMap::new();
        let source = CrossValidationSource::UcscRefseqSelect {
            coords_path: coords.path().to_path_buf(),
            select_path: select.path().to_path_buf(),
        };
        let err = cross_validate(&[tx], &aliases, Assembly::GRCh37, &source)
            .expect_err("bad strand must surface as a parse error");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("bad strand"),
            "expected explicit 'bad strand', got: {msg}"
        );
    }

    #[test]
    fn normalize_chrom_handles_all_formats() {
        assert_eq!(normalize_chrom("chr6", Assembly::GRCh38), "chr6");
        assert_eq!(normalize_chrom("6", Assembly::GRCh38), "chr6");
        assert_eq!(normalize_chrom("NC_000006.12", Assembly::GRCh38), "chr6");
        assert_eq!(normalize_chrom("NC_000006.11", Assembly::GRCh37), "chr6");
        assert_eq!(normalize_chrom("X", Assembly::GRCh38), "chrX");
        assert_eq!(normalize_chrom("", Assembly::GRCh38), "");
    }
}
