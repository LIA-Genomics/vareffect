//! Summary TSV cross-validation for built transcript models.
//!
//! Compares every row in the MANE summary TSV against the built
//! `Vec<TranscriptModel>`, checking chrom, strand, position, gene symbol,
//! and tier. Mismatches fail the build en masse so the implementer sees
//! every discrepancy at once.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use vareffect::chrom::{Assembly, is_patch_sequence, refseq_to_ucsc};
use vareffect::{Strand, TranscriptModel, TranscriptTier};

// Re-import types used only by tests.
#[cfg(test)]
use vareffect::{Biotype, CdsSegment, Exon};

/// Cross-validate a built transcript vec against NCBI's MANE summary TSV.
///
/// The summary TSV carries one row per transcript with columns including
/// `RefSeq_nuc`, `GRCh38_chr`, `chr_start`, `chr_end`, `chr_strand`,
/// `symbol`, and `MANE_status`. We assert that every summary row has a
/// matching built transcript with exactly-equal chrom, strand, position,
/// symbol, and tier. Mismatches are collected and the build is failed en
/// masse so the implementer sees every discrepancy.
///
/// `patch_aliases` -- RefSeq `NW_*`/`NT_*` -> UCSC contig name, used to
/// reconcile patch rows (summary uses RefSeq, GFF3 uses UCSC).
pub(super) fn cross_validate_summary(
    transcripts: &[TranscriptModel],
    summary_path: &Path,
    patch_aliases: &HashMap<String, String>,
    assembly: Assembly,
) -> Result<()> {
    let file = std::fs::File::open(summary_path)
        .with_context(|| format!("opening {}", summary_path.display()))?;
    let reader: Box<dyn BufRead> = if summary_path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Index the built transcripts by accession for O(1) lookup.
    let by_acc: HashMap<&str, &TranscriptModel> = transcripts
        .iter()
        .map(|t| (t.accession.as_str(), t))
        .collect();

    // Column indices are resolved from the header row so ordering changes
    // in future MANE releases don't silently break validation.
    let mut header_indices: Option<SummaryColumns> = None;
    let mut mismatches: Vec<String> = Vec::new();
    let mut rows_checked: usize = 0;
    // Track every accession the summary file names so we can run the
    // inverse check after the loop: any *built* transcript not in this set
    // indicates a parser bug (e.g. a stray `Note=...MANE_Select...`
    // admitting a non-MANE row).
    let mut summary_seen: HashSet<String> = HashSet::new();
    // Patch accessions in the summary that aren't in the alias map.
    // Deterministic order via `BTreeSet`; soft-warn at end of pass.
    let mut unknown_patch_aliases: std::collections::BTreeSet<String> =
        std::collections::BTreeSet::new();

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result.context("reading summary TSV line")?;
        if line.is_empty() {
            continue;
        }

        if header_indices.is_none() {
            // The header row starts with `#` in NCBI's format.
            let header = line.trim_start_matches('#');
            header_indices = Some(SummaryColumns::resolve(header).with_context(|| {
                format!(
                    "parsing summary TSV header at {}:{}",
                    summary_path.display(),
                    line_no + 1
                )
            })?);
            continue;
        }
        let cols = header_indices
            .as_ref()
            .expect("header line must precede data rows");

        let fields: Vec<&str> = line.split('\t').collect();
        let Some(row) = cols.extract(&fields) else {
            continue;
        };
        rows_checked += 1;

        summary_seen.insert(row.refseq_nuc.clone());

        let Some(tx) = by_acc.get(row.refseq_nuc.as_str()) else {
            mismatches.push(format!(
                "summary row {} ({}): no matching transcript in built store",
                line_no + 1,
                row.refseq_nuc
            ));
            continue;
        };

        // Chromosome -- summary uses `chr6` / `6` / `NC_000006.12` for primary
        // chroms and `NW_*`/`NT_*` for patches; GFF3 uses UCSC
        // `chr9_KN196479v1_fix`/`_alt`/... for the same patches. `patch_aliases`
        // resolves the patch case. A patch/primary disagreement between the
        // two sides is a real discrepancy and must fail the build.
        let summary_raw = row.grch38_chr.trim();
        let summary_chrom = normalize_summary_chrom(summary_raw, assembly);
        let built_is_patch = is_patch_sequence(&tx.chrom);
        let summary_is_patch = summary_raw.starts_with("NW_") || summary_raw.starts_with("NT_");
        match (built_is_patch, summary_is_patch) {
            (true, true) => match patch_aliases.get(summary_raw) {
                Some(expected_ucsc) if expected_ucsc == &tx.chrom => {}
                Some(expected_ucsc) => {
                    mismatches.push(format!(
                        "{}: patch chrom mismatch -- summary alias={}, built={}",
                        row.refseq_nuc, expected_ucsc, tx.chrom
                    ));
                }
                None => {
                    // Unknown patch accession -- record for an end-of-pass
                    // warning. The NCBI assembly report sometimes lags a
                    // new patch release; not a hard failure.
                    unknown_patch_aliases.insert(summary_raw.to_string());
                }
            },
            (false, false) => {
                if summary_chrom != tx.chrom {
                    mismatches.push(format!(
                        "{}: chrom mismatch -- summary={}, built={}",
                        row.refseq_nuc, summary_chrom, tx.chrom
                    ));
                }
            }
            (built_patch, summary_patch) => {
                // Exactly one side is a patch scaffold -- a genuine
                // disagreement, not something the alias map can reconcile.
                mismatches.push(format!(
                    "{}: patch/primary asymmetry -- summary={} (patch={}), built={} (patch={})",
                    row.refseq_nuc, summary_raw, summary_patch, tx.chrom, built_patch
                ));
            }
        }

        // Strand.
        let expected_strand = match row.chr_strand.as_str() {
            "+" => Some(Strand::Plus),
            "-" => Some(Strand::Minus),
            _ => None,
        };
        if let Some(strand) = expected_strand
            && strand != tx.strand
        {
            mismatches.push(format!(
                "{}: strand mismatch -- summary={:?}, built={:?}",
                row.refseq_nuc, strand, tx.strand
            ));
        }

        // Coordinate *exact* equality -- no tolerance. The conversion is
        // bijective (1-based inclusive end == 0-based exclusive end), so any
        // discrepancy is a real bug and masking it with `abs() > 1` would
        // silently ship wrong coordinates downstream to HGVS notation.
        if let (Ok(summary_start), Ok(summary_end)) =
            (row.chr_start.parse::<i64>(), row.chr_end.parse::<i64>())
        {
            let built_start = tx.tx_start as i64 + 1; // shift back to 1-based
            let built_end = tx.tx_end as i64; // 0-based exclusive == 1-based inclusive
            if summary_start != built_start {
                mismatches.push(format!(
                    "{}: tx_start mismatch -- summary={}, built={}",
                    row.refseq_nuc, summary_start, built_start
                ));
            }
            if summary_end != built_end {
                mismatches.push(format!(
                    "{}: tx_end mismatch -- summary={}, built={}",
                    row.refseq_nuc, summary_end, built_end
                ));
            }
        }

        // Gene symbol.
        if !row.symbol.is_empty() && row.symbol != tx.gene_symbol {
            mismatches.push(format!(
                "{}: gene_symbol mismatch -- summary={}, built={}",
                row.refseq_nuc, row.symbol, tx.gene_symbol
            ));
        }

        // Tier: MANE Select vs MANE Plus Clinical.
        let expected_tier = match row.mane_status.as_str() {
            "MANE Select" => Some(TranscriptTier::ManeSelect),
            "MANE Plus Clinical" => Some(TranscriptTier::ManePlusClinical),
            _ => None,
        };
        if let Some(tier) = expected_tier
            && tier != tx.tier
        {
            mismatches.push(format!(
                "{}: tier mismatch -- summary={:?}, built={:?}",
                row.refseq_nuc, tier, tx.tier
            ));
        }
    }

    // Inverse check: every built transcript must be named by the summary
    // TSV. A bug in tier detection (for example a substring false positive)
    // could admit a non-MANE transcript that would otherwise pass silently;
    // this guard fails the build loudly instead.
    for tx in transcripts {
        if !summary_seen.contains(&tx.accession) {
            mismatches.push(format!(
                "{}: built transcript is not present in the summary TSV",
                tx.accession
            ));
        }
    }

    // Unknown patch accessions -> deterministic WARN (not a hard failure).
    // Sorted via `BTreeSet` iteration order; one line per missing alias so
    // operators can `grep` the output.
    if !unknown_patch_aliases.is_empty() {
        tracing::warn!(
            store = "transcript_models",
            count = unknown_patch_aliases.len(),
            "summary TSV references patch accessions missing from patch_chrom_aliases.csv -- \
             the NCBI assembly report may be out of date; aliases: [{}]",
            unknown_patch_aliases
                .iter()
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    eprintln!(
        "  Summary TSV cross-validation: {} rows checked, {} mismatches, {} unknown patch aliases",
        rows_checked,
        mismatches.len(),
        unknown_patch_aliases.len(),
    );

    if !mismatches.is_empty() {
        // Surface the first 20 to keep error output bounded; the total
        // count is already in the header line above.
        let shown = mismatches.iter().take(20).cloned().collect::<Vec<_>>();
        bail!(
            "summary TSV cross-validation failed with {} mismatches (first 20 shown):\n{}",
            mismatches.len(),
            shown.join("\n")
        );
    }
    Ok(())
}

/// Normalize a summary TSV chromosome value to UCSC form.
///
/// Handles the three formats observed in NCBI files:
///   `"6"` -> `"chr6"`, `"chr6"` -> `"chr6"`, `"NC_000006.12"` -> `"chr6"`.
///
/// `NC_*` translation respects the assembly the validator was called with;
/// the GRCh37 and GRCh38 NC_* version tables differ for every primary
/// chromosome except chrM.
fn normalize_summary_chrom(raw: &str, assembly: Assembly) -> String {
    let trimmed = raw.trim();
    if trimmed.starts_with("chr") {
        return trimmed.to_string();
    }
    if trimmed.starts_with("NC_") {
        return refseq_to_ucsc(assembly, trimmed).to_string();
    }
    // Bare number or letter form (e.g. "6", "X").
    if !trimmed.is_empty() {
        return format!("chr{trimmed}");
    }
    String::new()
}

/// Column index resolver for the MANE summary TSV header.
///
/// The column order is documented in NCBI's release but we resolve indices
/// from the header text so future column additions don't break parsing.
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

    fn extract(&self, fields: &[&str]) -> Option<SummaryRow> {
        let max_idx = [
            self.refseq_nuc,
            self.symbol,
            self.mane_status,
            self.grch38_chr,
            self.chr_start,
            self.chr_end,
            self.chr_strand,
        ]
        .into_iter()
        .max()
        .unwrap();
        if fields.len() <= max_idx {
            return None;
        }
        Some(SummaryRow {
            refseq_nuc: fields[self.refseq_nuc].to_string(),
            symbol: fields[self.symbol].to_string(),
            mane_status: fields[self.mane_status].to_string(),
            grch38_chr: fields[self.grch38_chr].to_string(),
            chr_start: fields[self.chr_start].to_string(),
            chr_end: fields[self.chr_end].to_string(),
            chr_strand: fields[self.chr_strand].to_string(),
        })
    }
}

struct SummaryRow {
    refseq_nuc: String,
    symbol: String,
    mane_status: String,
    grch38_chr: String,
    chr_start: String,
    chr_end: String,
    chr_strand: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn normalize_summary_chrom_handles_all_formats() {
        assert_eq!(normalize_summary_chrom("chr6", Assembly::GRCh38), "chr6");
        assert_eq!(normalize_summary_chrom("6", Assembly::GRCh38), "chr6");
        assert_eq!(
            normalize_summary_chrom("NC_000006.12", Assembly::GRCh38),
            "chr6"
        );
        assert_eq!(normalize_summary_chrom("X", Assembly::GRCh38), "chrX");
        assert_eq!(normalize_summary_chrom("", Assembly::GRCh38), "");
    }

    /// Build a minimal patch-sequence `TranscriptModel` suitable for
    /// `cross_validate_summary` testing.
    fn patch_transcript(accession: &str, chrom: &str, gene_symbol: &str) -> TranscriptModel {
        TranscriptModel {
            accession: accession.to_string(),
            protein_accession: Some("NP_999999.1".into()),
            gene_symbol: gene_symbol.to_string(),
            hgnc_id: Some("HGNC:99999".into()),
            ensembl_accession: None,
            chrom: chrom.to_string(),
            strand: Strand::Plus,
            tx_start: 999,
            tx_end: 2000,
            cds_genomic_start: Some(1099),
            cds_genomic_end: Some(1900),
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: 999,
                genomic_end: 2000,
            }],
            cds_segments: vec![CdsSegment {
                exon_index: 0,
                genomic_start: 1099,
                genomic_end: 1900,
                phase: 0,
            }],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 1,
            genome_transcript_divergent: false,
            translational_exception: None,
        }
    }

    fn write_minimal_summary_tsv(
        refseq_nuc: &str,
        grch38_chr: &str,
        symbol: &str,
    ) -> tempfile::NamedTempFile {
        let mut file = tempfile::Builder::new()
            .suffix(".summary.tsv")
            .tempfile()
            .expect("create tempfile");
        writeln!(
            file,
            "#GeneID\tEnsembl_nuc\tsymbol\tname\tRefSeq_nuc\tEnsembl_prot\tRefSeq_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand"
        )
        .unwrap();
        writeln!(
            file,
            "123\tENSG00000099999.1\t{symbol}\tfake\t{refseq_nuc}\tENSP99999.1\tNP_999999.1\tMANE Select\t{grch38_chr}\t1000\t2000\t+"
        )
        .unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn cross_validate_summary_resolves_patch_alias_hit() {
        let transcripts = vec![patch_transcript(
            "NM_PATCH1.1",
            "chr22_KI270879v1_alt",
            "PATCH1",
        )];
        let summary_file = write_minimal_summary_tsv("NM_PATCH1.1", "NT_187633.1", "PATCH1");
        let mut aliases = HashMap::new();
        aliases.insert(
            "NT_187633.1".to_string(),
            "chr22_KI270879v1_alt".to_string(),
        );

        cross_validate_summary(
            &transcripts,
            summary_file.path(),
            &aliases,
            Assembly::GRCh38,
        )
        .expect("validator should accept a resolved patch alias");
    }

    #[test]
    fn cross_validate_summary_rejects_patch_alias_mismatch() {
        let transcripts = vec![patch_transcript(
            "NM_PATCH2.1",
            "chr22_KI270879v1_alt",
            "PATCH2",
        )];
        let summary_file = write_minimal_summary_tsv("NM_PATCH2.1", "NT_187633.1", "PATCH2");
        let mut aliases = HashMap::new();
        aliases.insert("NT_187633.1".to_string(), "chr9_KN196479v1_fix".to_string());

        let err = cross_validate_summary(
            &transcripts,
            summary_file.path(),
            &aliases,
            Assembly::GRCh38,
        )
        .expect_err("validator should reject a patch alias mismatch");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("patch chrom mismatch"),
            "expected patch chrom mismatch, got {msg}"
        );
        assert!(msg.contains("chr9_KN196479v1_fix"));
        assert!(msg.contains("chr22_KI270879v1_alt"));
    }

    #[test]
    fn cross_validate_summary_tolerates_unknown_patch_alias() {
        let transcripts = vec![patch_transcript(
            "NM_PATCH3.1",
            "chr22_KI270879v1_alt",
            "PATCH3",
        )];
        let summary_file = write_minimal_summary_tsv("NM_PATCH3.1", "NT_999999.1", "PATCH3");
        let aliases: HashMap<String, String> = HashMap::new();

        cross_validate_summary(
            &transcripts,
            summary_file.path(),
            &aliases,
            Assembly::GRCh38,
        )
        .expect("validator should tolerate unknown patch aliases");
    }

    #[test]
    fn cross_validate_summary_rejects_built_patch_summary_primary() {
        let transcripts = vec![patch_transcript(
            "NM_ASYMM1.1",
            "chr22_KI270879v1_alt",
            "ASYMM1",
        )];
        let summary_file = write_minimal_summary_tsv("NM_ASYMM1.1", "NC_000022.11", "ASYMM1");
        let aliases: HashMap<String, String> = HashMap::new();

        let err = cross_validate_summary(
            &transcripts,
            summary_file.path(),
            &aliases,
            Assembly::GRCh38,
        )
        .expect_err("validator should reject patch/primary asymmetry");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("patch/primary asymmetry"),
            "expected asymmetry mismatch, got {msg}"
        );
        assert!(msg.contains("NC_000022.11"));
        assert!(msg.contains("chr22_KI270879v1_alt"));
    }

    #[test]
    fn cross_validate_summary_rejects_built_primary_summary_patch() {
        let transcripts = vec![patch_transcript("NM_ASYMM2.1", "chr22", "ASYMM2")];
        let summary_file = write_minimal_summary_tsv("NM_ASYMM2.1", "NT_187633.1", "ASYMM2");
        let mut aliases = HashMap::new();
        aliases.insert(
            "NT_187633.1".to_string(),
            "chr22_KI270879v1_alt".to_string(),
        );

        let err = cross_validate_summary(
            &transcripts,
            summary_file.path(),
            &aliases,
            Assembly::GRCh38,
        )
        .expect_err("validator should reject patch/primary asymmetry");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("patch/primary asymmetry"),
            "expected asymmetry mismatch, got {msg}"
        );
        assert!(msg.contains("NT_187633.1"));
    }
}
