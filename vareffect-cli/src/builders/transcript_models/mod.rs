//! Transcript model builder for the `vareffect` crate's
//! [`TranscriptStore`](vareffect::TranscriptStore), from one of three GFF3
//! sources ([`TranscriptSource`]):
//!
//! - **MANE**: NCBI's MANE GFF3 (`MANE.GRCh38.vX.X.refseq_genomic.gff.gz`);
//!   only rows tagged `MANE Select` / `MANE Plus Clinical` are candidates.
//! - **RefSeq**: NCBI's full GRCh38 annotation
//!   (`GCF_000001405.40_GRCh38.p14_genomic.gff.gz`); every transcript row with
//!   a curated `transcript_id` (`NM_` / `NR_`) is a candidate, and only those
//!   on the primary assembly (`NC_`) are admitted.
//! - **GENCODE**: the comprehensive annotation
//!   (`gencode.vNN.annotation.gff3.gz`); every `transcript` row is a candidate.
//!
//! The models are written as a versioned store file via the shared
//! [`write_transcript_store`] helper, in 0-based half-open (BED/UCSC)
//! coordinates, with a sorted `<stem>.excluded.tsv` sidecar listing every
//! candidate that was not admitted and why.
//!
//! # GFF3 parsing strategy
//!
//! Single pass over the (optionally gzipped) GFF3 into scratch maps keyed on
//! the transcript row's `ID`: transcript rows, exon spans, CDS spans with
//! their phase and per-row notes (protein accession, RefSeq partial ranges
//! and exceptions, RefSeq and MANE `transl_except`), RefSeq `cDNA_match` gaps, and GENCODE
//! selenocysteine rows. Child rows (`start_codon`, `stop_codon`, UTRs,
//! `stop_codon_redefined_as_selenocysteine`, `cDNA_match`) are matched
//! explicitly and never folded into exons or CDS: the CDS rows already
//! contain the stop codon, so folding `stop_codon` in would add 3 bases.
//!
//! After the scan, candidates are assembled in `ID` order into
//! [`vareffect::TranscriptModel`]s: exons sorted 5'->3' on the transcript,
//! CDS segments mapped to their exon, and the per-source attributes
//! (completeness, genome mismatch, selenocysteine codons, biotype, tier)
//! applied. A candidate is either admitted or excluded with one reason;
//! anything else is a build error, and admitted + excluded must equal the
//! candidate count.
//!
//! If a MANE summary TSV path is supplied, a cross-validation pass compares
//! every row against the built vec and fails the build on mismatches.
//!
//! An accession may repeat only as one chrX and one chrY copy inside the
//! GRCh38 pseudoautosomal regions ([`check_repeated_accessions`]); any other
//! repeat fails the build.

mod cross_validate;
mod gff3_attrs;
mod source_rules;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use indicatif::{ProgressBar, ProgressStyle};
use sha2::{Digest, Sha256};
use vareffect::chrom::{canonical_order_key, is_patch_sequence, refseq_to_ucsc};
use vareffect::{
    Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptSource, TranscriptTier,
};

use cross_validate::cross_validate_summary;
use gff3_attrs::{
    extract_attr, extract_ensembl_from_dbxref, extract_hgnc_from_dbxref,
    extract_protein_id_from_cds_attrs, extract_refseq_acc_from_dbxref, extract_transcript_from_id,
};
use source_rules::{
    GapOp, gap_moves_cds_numbering, one_cds_per_exon, parse_gap_ops, parse_transl_except,
    phases_consistent, refseq_biotype,
};

use crate::common::{BuildOutput, StoreWrite, store_stem, write_transcript_store};

/// GFF3 feature types that carry the MANE tag on MANE Select / MANE Plus
/// Clinical transcripts in release v1.5. Any admitted transcript whose
/// feature type falls outside this set emits a build-log warning so future
/// MANE releases that introduce new biotypes (e.g. `vault_RNA`, `Y_RNA`)
/// are visible immediately instead of silently flowing through.
const KNOWN_TRANSCRIPT_FEATURE_TYPES: &[&str] = &[
    "mRNA",
    "lnc_RNA",
    "antisense_RNA",
    "snoRNA",
    "snRNA",
    "RNase_MRP_RNA",
    "telomerase_RNA",
];

/// What to build, from which files, and where to write it.
#[derive(Debug, Clone, Copy)]
pub struct BuildRequest<'a> {
    /// Annotation source; selects the candidate and attribute rules and the
    /// output file stem.
    pub source: TranscriptSource,
    /// The source GFF3 (`.gz` or uncompressed).
    pub input: &'a Path,
    /// Optional MANE summary TSV (`MANE.GRCh38.vX.X.summary.txt.gz`) for
    /// cross-validation. MANE only.
    pub summary_input: Option<&'a Path>,
    /// `patch_chrom_aliases.csv`, required with `summary_input`.
    pub patch_aliases_input: Option<&'a Path>,
    /// Directory for the store, its sidecars and the exclusion list.
    pub output_dir: &'a Path,
    /// Source release recorded in the store (e.g. `"1.5"`, `"RS_2025_08"`,
    /// `"50"`).
    pub release: &'a str,
}

/// Build a transcript model store from a MANE, RefSeq or GENCODE GFF3.
///
/// # Arguments
///
/// * `request` -- Source, input files, output directory and release. When
///   `summary_input` is set, every built transcript is cross-validated
///   against the MANE summary's chrom, strand, position, gene symbol and
///   tier; mismatches are fatal.
///
/// # Returns
///
/// Tuple of ([`BuildOutput`], store file size in bytes). The output
/// directory receives `<stem>.bin`, `<stem>.bin.sha256`,
/// `<stem>.manifest.json` and `<stem>.excluded.tsv` (stem from
/// [`store_stem`]).
///
/// # Errors
///
/// Returns an error if the GFF3 cannot be read or is inconsistent (a CDS
/// outside every exon, a transcript row span other than its exon span, a
/// padded CDS not marked 5'-incomplete, a malformed phase, `Gap` or
/// `transl_except`), if no transcript is admitted, if an accession repeats
/// outside a PAR pair, if a summary is given for a source other than MANE or
/// without `patch_aliases_input`, if cross-validation fails, if the built
/// store does not load, or if an output file cannot be written.
pub fn build(request: &BuildRequest<'_>) -> Result<(BuildOutput, usize)> {
    let source = request.source;
    // Also refuses a source this builder has no rules for, so the `_` arms
    // below (required by the `#[non_exhaustive]` enum) stand for MANE.
    let stem = store_stem(source)?;
    let parsed = parse_gff3(request.input, source)?;
    if parsed.transcripts.is_empty() {
        bail!(
            "no transcripts admitted from {} as a {source} GFF3 ({} candidates): check that \
             the input matches the source",
            request.input.display(),
            parsed.candidates
        );
    }

    check_repeated_accessions(&parsed.transcripts)?;

    if let Some(summary_path) = request.summary_input {
        if source != TranscriptSource::Mane {
            bail!("a MANE summary cross-validates MANE builds only, not {source}");
        }
        let aliases_path = request.patch_aliases_input.ok_or_else(|| {
            anyhow::anyhow!(
                "--summary-input requires --patch-chrom-aliases so patch-contig \
                 cross-validation can resolve `NW_*`/`NT_*` accessions"
            )
        })?;
        let aliases =
            crate::builders::patch_chrom_aliases::load(aliases_path).with_context(|| {
                format!(
                    "loading patch-chrom aliases from {}",
                    aliases_path.display()
                )
            })?;
        cross_validate_summary(&parsed.transcripts, summary_path, &aliases).with_context(|| {
            format!(
                "cross-validating GFF3 build against {}",
                summary_path.display()
            )
        })?;
    }

    let excluded_name = format!("{stem}.excluded.tsv");
    let excluded_tsv = parsed.excluded_tsv();
    let excluded_sha256 = hex::encode(Sha256::digest(excluded_tsv.as_bytes()));

    let written = write_transcript_store(&StoreWrite {
        source,
        release: request.release,
        transcripts: &parsed.transcripts,
        input: request.input,
        output_dir: request.output_dir,
        manifest_extra: parsed.manifest_report(&excluded_name, &excluded_sha256),
    })?;
    // Written only once the store is, so a failed build never leaves an
    // exclusion list that the previous store's manifest does not describe.
    let excluded_path = request.output_dir.join(&excluded_name);
    std::fs::write(&excluded_path, &excluded_tsv)
        .with_context(|| format!("writing {}", excluded_path.display()))?;
    Ok(written)
}

/// GRCh38 pseudoautosomal regions as 0-based half-open `(chrom, start, end)`
/// (GRC, https://www.ncbi.nlm.nih.gov/grc/human: PAR1 X/Y 10,001-2,781,479;
/// PAR2 X 155,701,383-156,030,895, Y 56,887,903-57,217,415, 1-based).
const GRCH38_PARS: [(&str, u64, u64); 4] = [
    ("chrX", 10_000, 2_781_479),
    ("chrX", 155_701_382, 156_030_895),
    ("chrY", 10_000, 2_781_479),
    ("chrY", 56_887_902, 57_217_415),
];

/// `true` if the transcript overlaps a GRCh38 PAR on its chromosome. Overlap,
/// not containment: RefSeq's chrX copies of SPRY3 (NM_001394353.1 and two
/// more) start 89 kb before PAR2 and end inside it.
fn in_par(tx: &TranscriptModel) -> bool {
    GRCH38_PARS
        .iter()
        .any(|&(chrom, start, end)| tx.chrom == chrom && tx.tx_start < end && start < tx.tx_end)
}

/// Refuse every repeated accession except a PAR pair: one chrX and one chrY
/// copy, each overlapping a GRCh38 pseudoautosomal region.
///
/// # Errors
///
/// Returns an error listing every other repeated accession with its
/// chromosomes.
fn check_repeated_accessions(transcripts: &[TranscriptModel]) -> Result<()> {
    let mut copies: HashMap<&str, Vec<&TranscriptModel>> = HashMap::new();
    for tx in transcripts {
        copies.entry(tx.accession.as_str()).or_default().push(tx);
    }
    let mut refused: Vec<String> = copies
        .iter()
        .filter(|(_, txs)| txs.len() > 1)
        .filter(|(_, txs)| {
            let par_pair = txs.len() == 2
                && txs.iter().all(|tx| in_par(tx))
                && txs.iter().any(|tx| tx.chrom == "chrX")
                && txs.iter().any(|tx| tx.chrom == "chrY");
            !par_pair
        })
        .map(|(accession, txs)| {
            let chroms: Vec<&str> = txs.iter().map(|tx| tx.chrom.as_str()).collect();
            format!("{accession} ({})", chroms.join(", "))
        })
        .collect();
    if refused.is_empty() {
        return Ok(());
    }
    refused.sort();
    bail!(
        "{} accession(s) repeat outside a chrX + chrY pseudoautosomal pair: {}",
        refused.len(),
        refused.join("; ")
    )
}

// ---------------------------------------------------------------------------
// Parse result
// ---------------------------------------------------------------------------

/// Why a candidate transcript was not admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum ExclusionReason {
    /// RefSeq row on an alternate locus or patch (`NT_` / `NW_`).
    NonPrimaryAssembly,
    /// MANE row whose accession could not be resolved.
    NoAccession,
    /// RefSeq row without a `gene-` parent.
    NoGeneParent,
    /// An exon or CDS row of the transcript had unreadable coordinates or
    /// strand, so its model would be missing that row.
    MalformedChildRow,
    /// No exon rows.
    NoExons,
    /// Two CDS segments in one exon (RefSeq ribosomal slippage).
    MultiCdsPerExon,
    /// An interior CDS phase disagrees with the CDS length before it.
    PhaseMismatch,
}

impl ExclusionReason {
    /// Label used in the exclusion list and the manifest.
    fn as_str(self) -> &'static str {
        match self {
            Self::NonPrimaryAssembly => "non_primary_assembly",
            Self::NoAccession => "no_accession",
            Self::NoGeneParent => "no_gene_parent",
            Self::MalformedChildRow => "malformed_child_row",
            Self::NoExons => "no_exons",
            Self::MultiCdsPerExon => "multi_cds_per_exon",
            Self::PhaseMismatch => "phase_mismatch",
        }
    }
}

/// A candidate transcript that was not admitted.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Exclusion {
    accession: String,
    chrom: String,
    /// Transcript row start, 0-based.
    start: u64,
    reason: ExclusionReason,
}

/// Transcripts admitted from one GFF3, with the accounting the manifest and
/// exclusion list report.
#[derive(Debug, Default)]
struct ParsedSource {
    /// Admitted models, sorted by accession, then canonical chromosome order,
    /// then start.
    transcripts: Vec<TranscriptModel>,
    /// Candidate transcript rows (admitted + excluded).
    candidates: u64,
    /// Excluded candidates, sorted by (accession, chrom, start, reason).
    excluded: Vec<Exclusion>,
    /// Rows skipped before they could be read, by reason (a short row,
    /// unparsable coordinates, an unknown strand, a gene without an ID).
    skipped_rows: BTreeMap<&'static str, u64>,
    /// Admitted transcripts whose source attributes contradict each other,
    /// by kind; kept, not excluded.
    anomalies: BTreeMap<&'static str, Vec<String>>,
    /// Admitted transcripts with a `transl_except` entry other than `aa:Sec`
    /// (not stored).
    other_transl_except: u64,
}

impl ParsedSource {
    /// The exclusion list: a header, then one sorted row per excluded
    /// candidate (`accession`, `chrom`, 0-based `start`, `reason`).
    fn excluded_tsv(&self) -> String {
        let mut tsv = String::from("#accession\tchrom\tstart\treason\n");
        for exclusion in &self.excluded {
            tsv.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                exclusion.accession,
                exclusion.chrom,
                exclusion.start,
                exclusion.reason.as_str()
            ));
        }
        tsv
    }

    /// Manifest fields describing this build: counts by outcome, skipped
    /// rows, anomalies and the exclusion list's name and SHA-256.
    fn manifest_report(
        &self,
        excluded_file: &str,
        excluded_sha256: &str,
    ) -> serde_json::Map<String, serde_json::Value> {
        let mut excluded_by_reason: BTreeMap<&str, u64> = BTreeMap::new();
        for exclusion in &self.excluded {
            *excluded_by_reason
                .entry(exclusion.reason.as_str())
                .or_default() += 1;
        }
        let mut report = serde_json::Map::new();
        report.insert(
            "counts".to_string(),
            serde_json::json!({
                "candidates": self.candidates,
                "admitted": self.transcripts.len(),
                "excluded": excluded_by_reason,
                "other_transl_except": self.other_transl_except,
            }),
        );
        report.insert(
            "skipped_rows".to_string(),
            serde_json::json!(self.skipped_rows),
        );
        report.insert("anomalies".to_string(), serde_json::json!(self.anomalies));
        report.insert(
            "excluded_file".to_string(),
            serde_json::Value::String(excluded_file.to_string()),
        );
        report.insert(
            "excluded_sha256".to_string(),
            serde_json::Value::String(excluded_sha256.to_string()),
        );
        report
    }
}

// ---------------------------------------------------------------------------
// GFF3 parsing
// ---------------------------------------------------------------------------

/// Scratch state filled by the single GFF3 pass.
#[derive(Debug, Default)]
struct Scan {
    gene_info: HashMap<String, GeneInfo>,
    /// Candidate transcript rows by GFF3 `ID`.
    pending: HashMap<String, PendingTranscript>,
    cds_ranges: HashMap<String, Vec<CdsRawSegment>>,
    cds_notes: HashMap<String, CdsNotes>,
    exon_ranges: HashMap<String, Vec<(u64, u64)>>,
    /// RefSeq `cDNA_match` gaps by (transcript accession, UCSC chrom):
    /// in-row `Gap` operations and, after [`add_coverage_breaks`], the
    /// transcript segments left out between blocks. An entry exists only for
    /// a gapped alignment.
    alignment_gaps: HashMap<(String, String), Vec<GapOp>>,
    /// RefSeq `cDNA_match` blocks by (transcript accession, UCSC chrom).
    match_blocks: HashMap<(String, String), Vec<MatchBlock>>,
    /// Parents of exon / CDS rows skipped for unreadable geometry.
    malformed_children: BTreeSet<String>,
    /// GENCODE `stop_codon_redefined_as_selenocysteine` starts by
    /// `transcript_id`.
    gencode_selenocysteine: HashMap<String, BTreeSet<u64>>,
    /// MANE admitted-candidate count per feature type, for the build log.
    feature_type_counts: HashMap<String, usize>,
    skipped_rows: BTreeMap<&'static str, u64>,
}

impl Scan {
    fn skip(&mut self, reason: &'static str) {
        *self.skipped_rows.entry(reason).or_default() += 1;
    }
}

/// Coordinates and strand of a GFF3 row, or the row-level skip reason.
fn row_geometry(fields: &[&str]) -> std::result::Result<(u64, u64, Strand), &'static str> {
    let (Ok(start), Ok(end)) = (fields[3].parse::<u64>(), fields[4].parse::<u64>()) else {
        return Err("unparsable_coordinates");
    };
    match fields[6] {
        "+" => Ok((start, end, Strand::Plus)),
        "-" => Ok((start, end, Strand::Minus)),
        _ => Err("unknown_strand"),
    }
}

/// Parse a GFF3 of `source` into its admitted models and accounting.
fn parse_gff3(input: &Path, source: TranscriptSource) -> Result<ParsedSource> {
    let file =
        std::fs::File::open(input).with_context(|| format!("opening {}", input.display()))?;

    // Transparently handle gzipped and uncompressed inputs.
    let reader: Box<dyn BufRead> = if input.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Spinner-with-rate rather than a bounded bar: counting the lines up
    // front would double the I/O.
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template(&format!(
            "  {{spinner}} parsing {source} GFF3: {{pos}} lines ({{per_sec}})"
        ))
        .expect("valid progress template"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(120));

    let mut scan = Scan::default();
    for line_result in reader.lines() {
        let line = line_result.context("reading GFF3 line")?;
        pb.inc(1);
        scan_line(&mut scan, &line, source)?;
    }
    pb.finish_with_message(format!("{source} GFF3 parsed"));

    assemble(scan, source)
}

/// Fold one GFF3 line into the scan.
fn scan_line(scan: &mut Scan, line: &str, source: TranscriptSource) -> Result<()> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return Ok(());
    }
    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() < 9 {
        scan.skip("short_row");
        return Ok(());
    }
    let feature_type = fields[2];
    let attrs = fields[8];

    match feature_type {
        // Gene rows: HGNC ID, symbol and gene biotype for MANE, the gene
        // biotype for RefSeq. GENCODE carries everything on its rows.
        "gene" | "pseudogene" => {
            if source == TranscriptSource::Gencode {
                return Ok(());
            }
            let Some(id) = extract_attr(attrs, "ID") else {
                scan.skip("gene_without_id");
                return Ok(());
            };
            let gene_symbol = extract_attr(attrs, "Name")
                .or_else(|| extract_attr(attrs, "gene"))
                .unwrap_or_default();
            scan.gene_info.insert(
                id,
                GeneInfo {
                    gene_symbol,
                    hgnc_id: extract_hgnc_from_dbxref(attrs),
                    gene_biotype: extract_attr(attrs, "gene_biotype"),
                },
            );
        }

        // Top-level region rows are metadata only.
        "region" => {}

        "CDS" => {
            let Some(parent) = extract_first_parent(attrs) else {
                scan.skip("child_without_parent");
                return Ok(());
            };
            let (gff_start, gff_end, _) = match row_geometry(&fields) {
                Ok(geometry) => geometry,
                Err(reason) => {
                    scan.skip(reason);
                    scan.malformed_children.insert(parent);
                    return Ok(());
                }
            };
            let refseq = source == TranscriptSource::RefSeq;
            // GFF3 column 8 is phase: `0`, `1`, `2`, or `.` (unknown).
            // Reject malformed values loudly.
            let phase = parse_phase(fields[7]).with_context(|| {
                format!(
                    "parsing CDS phase column 8 for parent {parent}: value {:?}",
                    fields[7]
                )
            })?;
            scan.cds_ranges
                .entry(parent.clone())
                .or_default()
                .push(CdsRawSegment {
                    gff_start,
                    gff_end,
                    phase,
                    start_range: refseq && extract_attr(attrs, "start_range").is_some(),
                    end_range: refseq && extract_attr(attrs, "end_range").is_some(),
                });
            let notes = scan.cds_notes.entry(parent.clone()).or_default();
            if notes.protein_id.is_none() {
                notes.protein_id = extract_protein_id_from_cds_attrs(attrs);
            }
            if refseq {
                notes
                    .read_refseq(attrs)
                    .with_context(|| format!("reading CDS attributes of {parent}"))?;
            } else if source == TranscriptSource::Mane {
                notes
                    .read_transl_except(attrs)
                    .with_context(|| format!("reading CDS attributes of {parent}"))?;
            }
        }

        "exon" => {
            let Some(parent) = extract_first_parent(attrs) else {
                scan.skip("child_without_parent");
                return Ok(());
            };
            match row_geometry(&fields) {
                Ok((start, end, _)) => scan
                    .exon_ranges
                    .entry(parent)
                    .or_default()
                    .push((start, end)),
                Err(reason) => {
                    scan.skip(reason);
                    scan.malformed_children.insert(parent);
                }
            }
        }

        // Child rows that are never exons or CDS: the CDS rows already
        // contain the stop codon (GENCODE README: "In GFF3 the start and
        // stop codon are included in the CDS").
        "start_codon" | "stop_codon" | "five_prime_UTR" | "three_prime_UTR" => {}

        "stop_codon_redefined_as_selenocysteine" => {
            if source != TranscriptSource::Gencode {
                return Ok(());
            }
            // The row's Parent is the CDS (`CDS:ENST...`); key on the
            // transcript instead.
            let Some(transcript_id) = extract_attr(attrs, "transcript_id") else {
                return Ok(());
            };
            match row_geometry(&fields) {
                Ok((start, end, _)) => {
                    // A codon split across an exon junction would need its
                    // rows grouped; GENCODE v50 has none.
                    if end + 1 != start + 3 {
                        bail!(
                            "selenocysteine row of {transcript_id} spans {start}..={end}, \
                             not one 3-base codon"
                        );
                    }
                    scan.gencode_selenocysteine
                        .entry(transcript_id)
                        .or_default()
                        .insert(start.saturating_sub(1));
                }
                Err(reason) => scan.skip(reason),
            }
        }

        "cDNA_match" => {
            if source != TranscriptSource::RefSeq {
                return Ok(());
            }
            let Some(target) = extract_attr(attrs, "Target") else {
                return Ok(());
            };
            let mut target_fields = target.split_whitespace();
            let Some(accession) = target_fields.next() else {
                return Ok(());
            };
            if !is_curated_refseq(accession) {
                return Ok(());
            }
            let (Some(Ok(target_start)), Some(Ok(target_end))) = (
                target_fields.next().map(str::parse::<u64>),
                target_fields.next().map(str::parse::<u64>),
            ) else {
                bail!("cDNA_match Target {target:?} has no numeric range");
            };
            let (start, end, strand) = match row_geometry(&fields) {
                Ok(geometry) => geometry,
                Err(reason) => bail!("cDNA_match of {accession}: {reason}"),
            };
            let key = (accession.to_string(), refseq_to_ucsc(fields[0]).to_string());
            scan.match_blocks
                .entry(key.clone())
                .or_default()
                .push(MatchBlock {
                    alignment_id: extract_attr(attrs, "ID").unwrap_or_default(),
                    target_start,
                    target_end,
                    genomic_start: start.saturating_sub(1),
                    genomic_end: end,
                    plus: matches!(strand, Strand::Plus),
                });
            if let Some(gap) = extract_attr(attrs, "Gap") {
                let ops = parse_gap_ops(&gap, start.saturating_sub(1), end, strand)
                    .with_context(|| format!("reading the cDNA_match Gap of {accession}"))?;
                if !ops.is_empty() {
                    scan.alignment_gaps.entry(key).or_default().extend(ops);
                }
            }
        }

        // Transcript-level catch-all: whether a row is a candidate is the
        // source's rule.
        _ => scan_candidate(scan, &fields, source)?,
    }
    Ok(())
}

/// One RefSeq `cDNA_match` row: a block of a transcript-to-genome alignment.
#[derive(Debug, Clone)]
struct MatchBlock {
    /// Alignment `ID`, shared by the blocks of one alignment.
    alignment_id: String,
    /// 1-based transcript range of the block (`Target`).
    target_start: u64,
    target_end: u64,
    /// 0-based half-open genomic span.
    genomic_start: u64,
    genomic_end: u64,
    plus: bool,
}

/// Add a [`GapOp::Insertion`] for every transcript segment an alignment
/// leaves out between two consecutive blocks (the next block's `Target`
/// does not start right after the previous one's): the genome-built model
/// lacks those transcript bases. The insertion sits at the junction after
/// the earlier block in transcript order.
fn add_coverage_breaks(scan: &mut Scan) {
    for (key, blocks) in &scan.match_blocks {
        let mut by_alignment: BTreeMap<&str, Vec<&MatchBlock>> = BTreeMap::new();
        for block in blocks {
            by_alignment
                .entry(block.alignment_id.as_str())
                .or_default()
                .push(block);
        }
        for blocks in by_alignment.values_mut() {
            blocks.sort_by_key(|block| block.target_start);
            for pair in blocks.windows(2) {
                let (previous, next) = (pair[0], pair[1]);
                if next.target_start != previous.target_end + 1 {
                    let at = if previous.plus {
                        previous.genomic_end
                    } else {
                        previous.genomic_start
                    };
                    scan.alignment_gaps
                        .entry(key.clone())
                        .or_default()
                        .push(GapOp::Insertion { at });
                }
            }
        }
    }
}

/// `true` for a curated RefSeq transcript accession (`NM_` / `NR_`).
fn is_curated_refseq(accession: &str) -> bool {
    accession.starts_with("NM_") || accession.starts_with("NR_")
}

/// Record a transcript row if it is a candidate of `source`.
///
/// - MANE: any row tagged `MANE Select` or `MANE Plus Clinical`, whatever
///   its feature type, so a future biotype is admitted automatically.
/// - GENCODE: every `transcript` row.
/// - RefSeq: every row with a curated `transcript_id`.
///
/// # Errors
///
/// Returns an error when a candidate `ID` repeats, which would otherwise
/// drop the earlier row from every count.
fn scan_candidate(scan: &mut Scan, fields: &[&str], source: TranscriptSource) -> Result<()> {
    let feature_type = fields[2];
    let attrs = fields[8];
    let tier = tier_from_attrs(attrs);
    let is_candidate = match source {
        TranscriptSource::Mane => matches!(
            tier,
            Some(TranscriptTier::ManeSelect | TranscriptTier::ManePlusClinical)
        ),
        TranscriptSource::Gencode => feature_type == "transcript",
        TranscriptSource::RefSeq => {
            extract_attr(attrs, "transcript_id").is_some_and(|id| is_curated_refseq(&id))
        }
        _ => false,
    };
    if !is_candidate {
        return Ok(());
    }
    let Some(id) = extract_attr(attrs, "ID") else {
        scan.skip("transcript_without_id");
        return Ok(());
    };
    let (tx_start_gff, tx_end_gff, strand) = match row_geometry(fields) {
        Ok(geometry) => geometry,
        Err(reason) => {
            scan.skip(reason);
            return Ok(());
        }
    };

    let accession = match source {
        // Every MANE transcript in v1.5 carries `Name=<accession>`; the
        // fallbacks are defensive.
        TranscriptSource::Mane => extract_attr(attrs, "Name")
            .or_else(|| extract_refseq_acc_from_dbxref(attrs))
            .or_else(|| extract_transcript_from_id(&id))
            .unwrap_or_default(),
        _ => extract_attr(attrs, "transcript_id").unwrap_or_else(|| id.clone()),
    };
    if source == TranscriptSource::Mane && !accession.is_empty() {
        *scan
            .feature_type_counts
            .entry(feature_type.to_string())
            .or_insert(0) += 1;
    }

    let tags = extract_attr(attrs, "tag").unwrap_or_default();
    let has_tag = |wanted: &str| tags.split(',').any(|tag| tag.trim() == wanted);
    let pending = PendingTranscript {
        accession,
        feature_type: feature_type.to_string(),
        parent: extract_first_parent(attrs),
        ensembl_accession: match source {
            TranscriptSource::Gencode => None,
            _ => extract_ensembl_from_dbxref(attrs),
        },
        seqid: fields[0].to_string(),
        strand,
        tx_start_gff,
        tx_end_gff,
        tier,
        gene_symbol: extract_attr(attrs, "gene_name").or_else(|| extract_attr(attrs, "gene")),
        hgnc_id: extract_attr(attrs, "hgnc_id").or_else(|| extract_hgnc_from_dbxref(attrs)),
        transcript_type: extract_attr(attrs, "transcript_type"),
        protein_id: extract_attr(attrs, "protein_id"),
        cds_start_nf: has_tag("cds_start_NF"),
        cds_end_nf: has_tag("cds_end_NF"),
        exception: extract_attr(attrs, "exception").is_some(),
    };
    if scan.pending.insert(id.clone(), pending).is_some() {
        bail!("transcript row ID {id} repeats");
    }
    Ok(())
}

/// The outcome of assembling one candidate.
enum Assembly {
    Admitted(Box<TranscriptModel>),
    Excluded(ExclusionReason),
}

/// Assemble every candidate, in `ID` order, into admitted models and
/// exclusions, and check that each candidate is accounted for once.
fn assemble(mut scan: Scan, source: TranscriptSource) -> Result<ParsedSource> {
    add_coverage_breaks(&mut scan);
    let mut parsed = ParsedSource {
        skipped_rows: std::mem::take(&mut scan.skipped_rows),
        ..ParsedSource::default()
    };
    let mut ids: Vec<&String> = scan.pending.keys().collect();
    ids.sort();
    let mut patch_seq_count: usize = 0;
    let mut hgnc_missing_count: usize = 0;
    let mut biotype_unknown_count: usize = 0;

    for id in ids {
        let pending = &scan.pending[id];
        parsed.candidates += 1;
        match assemble_one(id, pending, &scan, source, &mut parsed)? {
            Assembly::Admitted(model) => {
                if is_patch_sequence(&model.chrom) {
                    patch_seq_count += 1;
                }
                if model.hgnc_id.is_none() {
                    hgnc_missing_count += 1;
                }
                if matches!(model.biotype, Biotype::Unknown) {
                    biotype_unknown_count += 1;
                }
                parsed.transcripts.push(*model);
            }
            Assembly::Excluded(reason) => parsed.excluded.push(Exclusion {
                accession: pending.accession.clone(),
                chrom: chrom_name(source, &pending.seqid).to_string(),
                start: pending.tx_start_gff.saturating_sub(1),
                reason,
            }),
        }
    }
    let accounted = parsed.transcripts.len() as u64 + parsed.excluded.len() as u64;
    if accounted != parsed.candidates {
        bail!(
            "{accounted} transcripts admitted or excluded but {} candidates",
            parsed.candidates
        );
    }

    // Deterministic output: sort by accession, and the copies of a
    // repeated (PAR) accession canonically, so rebuilds with identical
    // inputs produce byte-identical MessagePack whatever the hash order.
    parsed.transcripts.sort_by(|a, b| {
        (
            a.accession.as_str(),
            canonical_order_key(&a.chrom),
            a.tx_start,
        )
            .cmp(&(
                b.accession.as_str(),
                canonical_order_key(&b.chrom),
                b.tx_start,
            ))
    });
    parsed.excluded.sort_by(|a, b| {
        (&a.accession, &a.chrom, a.start, a.reason).cmp(&(
            &b.accession,
            &b.chrom,
            b.start,
            b.reason,
        ))
    });
    for accessions in parsed.anomalies.values_mut() {
        accessions.sort();
    }

    if source == TranscriptSource::Mane {
        warn_on_uncommon_feature_types(&scan.feature_type_counts);
    }
    let mut excluded_by_reason: BTreeMap<&str, usize> = BTreeMap::new();
    for exclusion in &parsed.excluded {
        *excluded_by_reason
            .entry(exclusion.reason.as_str())
            .or_default() += 1;
    }
    eprintln!(
        "  Transcript models: {} admitted of {} candidates ({} patch-seq, {} without HGNC, \
         {} unknown biotype)",
        parsed.transcripts.len(),
        parsed.candidates,
        patch_seq_count,
        hgnc_missing_count,
        biotype_unknown_count,
    );
    if !excluded_by_reason.is_empty() {
        eprintln!("  Excluded: {excluded_by_reason:?}");
    }
    Ok(parsed)
}

/// UCSC chromosome name of a source seqid. MANE ships UCSC names (patch
/// contigs included); RefSeq `NC_` accessions are mapped and others pass
/// through.
fn chrom_name(source: TranscriptSource, seqid: &str) -> &str {
    match source {
        TranscriptSource::Gencode => seqid,
        _ => refseq_to_ucsc(seqid),
    }
}

/// Assemble one candidate.
///
/// # Errors
///
/// Returns an error for a GFF3 inconsistency that no exclusion reason
/// covers: more than `u16::MAX` exons, a transcript row span other than
/// the exon span, a CDS outside every exon, or a padded CDS not marked
/// 5'-incomplete.
fn assemble_one(
    id: &str,
    pending: &PendingTranscript,
    scan: &Scan,
    source: TranscriptSource,
    parsed: &mut ParsedSource,
) -> Result<Assembly> {
    let chrom = chrom_name(source, &pending.seqid);
    if source == TranscriptSource::RefSeq && chrom == pending.seqid {
        return Ok(Assembly::Excluded(ExclusionReason::NonPrimaryAssembly));
    }
    if pending.accession.is_empty() {
        return Ok(Assembly::Excluded(ExclusionReason::NoAccession));
    }
    if scan.malformed_children.contains(id) {
        return Ok(Assembly::Excluded(ExclusionReason::MalformedChildRow));
    }
    if source == TranscriptSource::RefSeq
        && !pending
            .parent
            .as_deref()
            .is_some_and(|parent| parent.starts_with("gene-"))
    {
        return Ok(Assembly::Excluded(ExclusionReason::NoGeneParent));
    }
    let Some(raw_exons) = scan.exon_ranges.get(id).filter(|exons| !exons.is_empty()) else {
        return Ok(Assembly::Excluded(ExclusionReason::NoExons));
    };
    let accession = pending.accession.as_str();

    // Sort by genomic start; minus-strand transcripts reverse so `exons[0]`
    // is the 5'-most exon on the transcript.
    let mut sorted_exons: Vec<(u64, u64)> = raw_exons.clone();
    sorted_exons.sort_by_key(|(s, _)| *s);
    if matches!(pending.strand, Strand::Minus) {
        sorted_exons.reverse();
    }
    if sorted_exons.len() > u16::MAX as usize {
        bail!(
            "transcript {accession} has {} exons (u16 overflow)",
            sorted_exons.len()
        );
    }
    let exons: Vec<Exon> = sorted_exons
        .iter()
        .enumerate()
        .map(|(i, (s, e))| Exon {
            exon_number: (i + 1) as u16,
            // GFF3 1-based fully-closed -> 0-based half-open.
            genomic_start: s.saturating_sub(1),
            genomic_end: *e,
        })
        .collect();

    // The model's span is the exon span; a transcript row that disagrees
    // is a GFF3 inconsistency, refused rather than silently resolved.
    let tx_start = pending.tx_start_gff.saturating_sub(1);
    let tx_end = pending.tx_end_gff;
    let exon_start = exons.iter().map(|e| e.genomic_start).min();
    let exon_end = exons.iter().map(|e| e.genomic_end).max();
    if (exon_start, exon_end) != (Some(tx_start), Some(tx_end)) {
        bail!(
            "transcript {accession} spans {tx_start}..{tx_end} but its exons span {exon_start:?}..{exon_end:?}"
        );
    }

    let raw_cds = scan.cds_ranges.get(id).map_or(&[][..], Vec::as_slice);
    let cds_segments =
        build_cds_segments(accession, &exons, &pending.strand, scan.cds_ranges.get(id))?;
    if !one_cds_per_exon(&cds_segments) {
        return Ok(Assembly::Excluded(ExclusionReason::MultiCdsPerExon));
    }
    if !phases_consistent(&cds_segments) {
        return Ok(Assembly::Excluded(ExclusionReason::PhaseMismatch));
    }
    let has_cds = !cds_segments.is_empty();
    let notes = scan.cds_notes.get(id);

    let (start_complete, end_complete) = match source {
        TranscriptSource::Gencode => (!pending.cds_start_nf, !pending.cds_end_nf),
        TranscriptSource::RefSeq => match notes {
            Some(notes) if has_cds => {
                let ranges = CdsRanges::of(raw_cds);
                if ranges.interior {
                    parsed
                        .anomalies
                        .entry("interior_cds_range")
                        .or_default()
                        .push(accession.to_string());
                }
                let (five, three) = if matches!(pending.strand, Strand::Plus) {
                    (ranges.low_partial, ranges.high_partial)
                } else {
                    (ranges.high_partial, ranges.low_partial)
                };
                // `partial=true` with no range on any CDS row: which end is
                // partial is unknown, so both are treated as incomplete.
                let unranged = notes.partial && !ranges.any;
                if unranged {
                    parsed
                        .anomalies
                        .entry("partial_without_range")
                        .or_default()
                        .push(accession.to_string());
                }
                (!(five || unranged), !(three || unranged))
            }
            Some(_) | None => (true, true),
        },
        _ => (true, true),
    };
    // A non-zero first phase pads codon 1: the CDS cannot start with a
    // complete start codon.
    if has_cds && cds_segments[0].phase != 0 && start_complete {
        bail!(
            "transcript {accession}: first CDS phase {} but the CDS is not marked 5'-incomplete",
            cds_segments[0].phase
        );
    }

    // An interior `start_range` / `end_range` marks a break in the
    // alignment's coverage inside the CDS (SHANK3 NM_001372044.2): the
    // transcript is not the genome there.
    let interior_range = source == TranscriptSource::RefSeq && CdsRanges::of(raw_cds).interior;
    let (genome_mismatch, cds_alignment_gap) = if source == TranscriptSource::RefSeq {
        let gaps = scan
            .alignment_gaps
            .get(&(accession.to_string(), chrom.to_string()));
        let mismatch_exception = notes.is_some_and(|notes| notes.mismatch_exception);
        let cds_gap = match (gaps, cds_segments.iter().map(|s| s.genomic_start).min()) {
            (Some(gaps), Some(cds_start)) => {
                let cds_end = cds_segments
                    .iter()
                    .map(|s| s.genomic_end)
                    .max()
                    .unwrap_or(cds_start);
                gap_moves_cds_numbering(gaps, cds_start, cds_end, &exons)
            }
            (Some(_), None) | (None, _) => false,
        };
        (
            pending.exception || mismatch_exception || gaps.is_some() || interior_range,
            cds_gap || interior_range,
        )
    } else {
        (false, false)
    };

    let selenocysteine: Vec<u64> = match source {
        TranscriptSource::Gencode => scan
            .gencode_selenocysteine
            .get(accession)
            .map(|codons| codons.iter().copied().collect())
            .unwrap_or_default(),
        TranscriptSource::RefSeq | TranscriptSource::Mane => notes
            .map(|notes| notes.selenocysteine.iter().copied().collect())
            .unwrap_or_default(),
        _ => Vec::new(),
    };
    if notes.is_some_and(|notes| notes.other_transl_except) {
        parsed.other_transl_except += 1;
    }

    let gene = pending
        .parent
        .as_deref()
        .and_then(|parent| scan.gene_info.get(parent));
    let (gene_symbol, hgnc_id, biotype) = match source {
        TranscriptSource::Gencode => (
            pending.gene_symbol.clone().unwrap_or_default(),
            pending.hgnc_id.clone(),
            Biotype::from_label(pending.transcript_type.as_deref().unwrap_or_default()),
        ),
        TranscriptSource::RefSeq => {
            if has_cds && accession.starts_with("NR_") {
                parsed
                    .anomalies
                    .entry("nr_with_cds")
                    .or_default()
                    .push(accession.to_string());
            }
            if !has_cds && accession.starts_with("NM_") {
                parsed
                    .anomalies
                    .entry("nm_without_cds")
                    .or_default()
                    .push(accession.to_string());
            }
            (
                pending
                    .gene_symbol
                    .clone()
                    .or_else(|| gene.map(|g| g.gene_symbol.clone()))
                    .unwrap_or_default(),
                pending
                    .hgnc_id
                    .clone()
                    .or_else(|| gene.and_then(|g| g.hgnc_id.clone())),
                refseq_biotype(
                    &pending.feature_type,
                    has_cds,
                    gene.and_then(|g| g.gene_biotype.as_deref()),
                ),
            )
        }
        // MANE: the gene row is authoritative for symbol, HGNC and biotype.
        _ => (
            gene.map(|g| g.gene_symbol.clone()).unwrap_or_default(),
            gene.and_then(|g| g.hgnc_id.clone()),
            biotype_for(accession, gene.and_then(|g| g.gene_biotype.as_deref())),
        ),
    };

    let protein_accession = match source {
        TranscriptSource::Gencode => pending
            .protein_id
            .clone()
            .or_else(|| notes.and_then(|notes| notes.protein_id.clone())),
        _ => notes.and_then(|notes| notes.protein_id.clone()),
    };

    let mut model = TranscriptModel::new(accession, gene_symbol, chrom, pending.strand, exons)
        .with_cds_segments(cds_segments)
        .with_protein_accession(protein_accession)
        .with_hgnc_id(hgnc_id)
        .with_ensembl_accession(pending.ensembl_accession.clone())
        .with_biotype(biotype)
        .with_cds_completeness(start_complete, end_complete)
        .with_genome_mismatch(genome_mismatch, cds_alignment_gap)
        .with_selenocysteine_codons(selenocysteine);
    if let Some(tier) = pending.tier {
        model = model.with_tier(tier);
    }
    Ok(Assembly::Admitted(Box::new(model)))
}

/// Warn on any admitted MANE feature type outside the known-common set and
/// log the per-type breakdown, sorted for stable output.
fn warn_on_uncommon_feature_types(feature_type_counts: &HashMap<String, usize>) {
    let mut breakdown: Vec<(&String, &usize)> = feature_type_counts.iter().collect();
    breakdown.sort_by(|a, b| a.0.cmp(b.0));
    for (ft, count) in &breakdown {
        if !KNOWN_TRANSCRIPT_FEATURE_TYPES.contains(&ft.as_str()) {
            tracing::warn!(
                feature_type = %ft,
                count = count,
                store = "transcript_models",
                "admitted MANE transcript(s) with uncommon feature type -- \
                 update KNOWN_TRANSCRIPT_FEATURE_TYPES if this is expected",
            );
        }
    }
    let breakdown_str = breakdown
        .iter()
        .map(|(ft, count)| format!("{ft}={count}"))
        .collect::<Vec<_>>()
        .join(", ");
    eprintln!("  By feature type: {breakdown_str}");
}

/// Intermediate gene-level info collected from `gene` rows.
#[derive(Debug)]
struct GeneInfo {
    gene_symbol: String,
    hgnc_id: Option<String>,
    gene_biotype: Option<String>,
}

/// A candidate transcript row, before exon/CDS assembly.
#[derive(Debug)]
struct PendingTranscript {
    /// Empty when a MANE row's accession cannot be resolved.
    accession: String,
    /// GFF3 column 3.
    feature_type: String,
    parent: Option<String>,
    ensembl_accession: Option<String>,
    /// GFF3 column 1, unmapped.
    seqid: String,
    strand: Strand,
    tx_start_gff: u64,
    tx_end_gff: u64,
    /// Highest tier tag on the row.
    tier: Option<TranscriptTier>,
    /// Row `gene_name` (GENCODE) or `gene` (RefSeq).
    gene_symbol: Option<String>,
    /// Row `hgnc_id` (GENCODE) or Dbxref HGNC (RefSeq).
    hgnc_id: Option<String>,
    /// GENCODE `transcript_type`.
    transcript_type: Option<String>,
    /// GENCODE row `protein_id`.
    protein_id: Option<String>,
    /// GENCODE `cds_start_NF` / `cds_end_NF` tags.
    cds_start_nf: bool,
    cds_end_nf: bool,
    /// RefSeq row `exception` attribute present.
    exception: bool,
}

/// Facts read from a transcript's CDS rows, combined over all of them.
#[derive(Debug, Default)]
struct CdsNotes {
    /// First protein accession seen.
    protein_id: Option<String>,

    /// RefSeq `partial=true` on any CDS row.
    partial: bool,
    /// A RefSeq CDS `exception` saying the protein follows transcript or
    /// proteomic data rather than the genome.
    mismatch_exception: bool,
    /// RefSeq or MANE `transl_except` selenocysteine codon starts.
    selenocysteine: BTreeSet<u64>,
    /// A RefSeq or MANE `transl_except` entry other than `aa:Sec`.
    other_transl_except: bool,
}

impl CdsNotes {
    /// Fold one RefSeq CDS row's attributes in.
    ///
    /// # Errors
    ///
    /// Returns an error for a malformed `transl_except`.
    fn read_refseq(&mut self, attrs: &str) -> Result<()> {
        self.partial |= extract_attr(attrs, "partial").as_deref() == Some("true");
        self.mismatch_exception |= extract_attr(attrs, "exception").is_some_and(|exception| {
            exception.contains("annotated by transcript or proteomic data")
        });
        self.read_transl_except(attrs)
    }

    /// Fold one CDS row's `transl_except` in: the MANE profile reads only
    /// this (the selenocysteine codons of its 25 selenoproteins).
    ///
    /// # Errors
    ///
    /// Returns an error for a malformed `transl_except`.
    fn read_transl_except(&mut self, attrs: &str) -> Result<()> {
        if let Some(raw) = extract_attr(attrs, "transl_except") {
            let parsed = parse_transl_except(&raw)?;
            self.selenocysteine.extend(parsed.selenocysteine);
            self.other_transl_except |= parsed.other > 0;
        }
        Ok(())
    }
}

/// Raw CDS row captured during the GFF3 scan, pre-coordinate-conversion.
///
/// The genomic start/end are 1-based fully-closed (GFF3 native); the phase
/// is the normalized GFF3 column 8 value (missing `.` -> 0). Converted to
/// [`CdsSegment`] in `build_cds_segments` once the containing exon index
/// can be resolved.
#[derive(Debug, Clone)]
struct CdsRawSegment {
    gff_start: u64,
    gff_end: u64,
    phase: u8,
    /// RefSeq `start_range` on this row (its low genomic end is partial).
    start_range: bool,
    /// RefSeq `end_range` on this row (its high genomic end is partial).
    end_range: bool,
}

/// Where a RefSeq CDS's partial-range markers sit.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct CdsRanges {
    /// A `start_range` on the row with the lowest start: the CDS's low
    /// genomic end is open.
    low_partial: bool,
    /// An `end_range` on the row with the highest end.
    high_partial: bool,
    /// A range on any other row: a coverage break inside the CDS, not an
    /// open end.
    interior: bool,
    /// Any range at all.
    any: bool,
}

impl CdsRanges {
    /// Classify the range markers of a transcript's CDS rows.
    fn of(rows: &[CdsRawSegment]) -> Self {
        let (Some(low), Some(high)) = (
            rows.iter().map(|r| r.gff_start).min(),
            rows.iter().map(|r| r.gff_end).max(),
        ) else {
            return Self::default();
        };
        let mut ranges = Self::default();
        for row in rows {
            if row.start_range {
                ranges.any = true;
                if row.gff_start == low {
                    ranges.low_partial = true;
                } else {
                    ranges.interior = true;
                }
            }
            if row.end_range {
                ranges.any = true;
                if row.gff_end == high {
                    ranges.high_partial = true;
                } else {
                    ranges.interior = true;
                }
            }
        }
        ranges
    }
}

/// Assemble `CdsSegment`s for a single transcript.
///
/// The returned vec is ordered 5'->3' on the transcript (reversed for
/// minus-strand genes) so `cds_segments[i]` aligns with the containing
/// exon in the caller's `exons` vec. The `exon_index` on each segment
/// points back into that `exons` vec for O(1) exon lookup downstream.
///
/// An empty input (`None` or empty vec) returns an empty `Vec<CdsSegment>`
/// -- non-coding transcripts carry no CDS rows and that's not an error.
///
/// # Errors
///
/// Returns an error if any CDS segment cannot be placed inside one of the
/// supplied exons -- this indicates either an off-by-one in the coordinate
/// conversion or a GFF3 inconsistency upstream, and we fail loudly rather
/// than silently drop the segment.
fn build_cds_segments(
    accession: &str,
    exons: &[Exon],
    strand: &Strand,
    raw: Option<&Vec<CdsRawSegment>>,
) -> Result<Vec<CdsSegment>> {
    let Some(raw) = raw else {
        return Ok(Vec::new());
    };
    if raw.is_empty() {
        return Ok(Vec::new());
    }

    // Convert to 0-based half-open and sort ascending by genomic start.
    // Minus-strand transcripts: reverse to match `exons` ordering (exon 0
    // lives at the highest genomic coordinate for minus strand). This keeps
    // `cds_segments[i].exon_index` monotonically increasing with transcript
    // 5'->3' order.
    let mut segments: Vec<CdsSegment> = raw
        .iter()
        .map(|r| CdsSegment {
            exon_index: 0, // placeholder -- resolved below
            // GFF3 is 1-based fully-closed; our types are 0-based half-open:
            //   genomic_start = gff_start - 1
            //   genomic_end   = gff_end       (already exclusive after swap)
            genomic_start: r.gff_start.saturating_sub(1),
            genomic_end: r.gff_end,
            phase: r.phase,
        })
        .collect();
    segments.sort_by_key(|s| s.genomic_start);
    if matches!(strand, Strand::Minus) {
        segments.reverse();
    }

    // Resolve the containing exon for each segment. Linear scan over exons
    // is fine -- the largest known human transcript (TTN) has ~363 exons,
    // so the per-segment cost is a few hundred comparisons. A binary search
    // would only help after the vec crosses ~1000 exons, which the human
    // transcriptome never does.
    for segment in &mut segments {
        let Some((idx, _)) = exons.iter().enumerate().find(|(_, e)| {
            e.genomic_start <= segment.genomic_start && segment.genomic_end <= e.genomic_end
        }) else {
            bail!(
                "transcript {accession}: CDS segment [{}, {}) does not fit inside any exon -- \
                 GFF3 inconsistency (check exon/CDS coordinate alignment)",
                segment.genomic_start,
                segment.genomic_end
            );
        };
        // `exon_index` is 0-based into the `exons` vec; cast is infallible
        // because the vec length already fits in u16 (verified upstream).
        segment.exon_index = idx as u16;
    }

    Ok(segments)
}

/// Parse a GFF3 column 8 phase value.
///
/// Valid values are `0`, `1`, `2`, or `.` (unknown). Anything else is
/// rejected as malformed -- a phase > 2 is logically impossible and
/// indicates upstream corruption.
fn parse_phase(raw: &str) -> Result<u8> {
    match raw.trim() {
        "." | "" => Ok(0),
        "0" => Ok(0),
        "1" => Ok(1),
        "2" => Ok(2),
        other => bail!("invalid GFF3 phase value {other:?} (expected 0, 1, 2, or .)"),
    }
}

/// Extract the first `Parent` value from a GFF3 attribute string.
///
/// GFF3 spec allows `Parent=a,b,c` for features with multiple parents.
/// Real MANE releases carry a single parent on transcript / exon / CDS
/// rows, but splitting defensively keeps a future multi-parent row from
/// silently breaking the gene-lookup chain.
fn extract_first_parent(attrs: &str) -> Option<String> {
    let raw = extract_attr(attrs, "Parent")?;
    raw.split(',')
        .next()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
}

/// Detect the highest curation tier on a transcript row: MANE Select >
/// MANE Plus Clinical > RefSeq Select.
///
/// Returns `None` for untagged rows. Parses the `tag=` attribute
/// structurally, so a `Note=...includes MANE_Select isoforms...` description
/// never matches.
///
/// GFF3 tag values are comma-separated, so a row can carry multiple tags
/// (`tag=basic,MANE Select`). Both the NCBI space form (`MANE Select`)
/// and the GENCODE underscore form (`MANE_Select`) are accepted.
fn tier_from_attrs(attrs: &str) -> Option<TranscriptTier> {
    let tag = extract_attr(attrs, "tag")?;
    let rank = |tier: TranscriptTier| match tier {
        TranscriptTier::ManeSelect => 0,
        TranscriptTier::ManePlusClinical => 1,
        TranscriptTier::RefSeqSelect => 2,
        // `#[non_exhaustive]`: a future tier ranks below the known ones.
        _ => 3,
    };
    tag.split(',')
        .filter_map(|value| match value.trim() {
            "MANE Select" | "MANE_Select" => Some(TranscriptTier::ManeSelect),
            "MANE Plus Clinical" | "MANE_Plus_Clinical" => Some(TranscriptTier::ManePlusClinical),
            "RefSeq Select" | "RefSeq_Select" => Some(TranscriptTier::RefSeqSelect),
            _ => None,
        })
        .min_by_key(|&tier| rank(tier))
}

/// Derive a biotype from the transcript accession and optional gene-level
/// `gene_biotype` attribute.
///
/// - `NM_*` -> [`Biotype::ProteinCoding`] (default)
/// - `NR_*` -> [`Biotype::NonCodingRna`] (default)
/// - Any non-empty gene-level `gene_biotype` overrides the accession-prefix
///   default and is parsed via [`Biotype::from_label`] so unknown upstream
///   labels survive as [`Biotype::Other`].
/// - If neither signal is available -> [`Biotype::Unknown`].
fn biotype_for(accession: &str, gene_biotype: Option<&str>) -> Biotype {
    if let Some(gb) = gene_biotype
        && !gb.is_empty()
    {
        return Biotype::from_label(gb);
    }
    if accession.starts_with("NM_") {
        Biotype::ProteinCoding
    } else if accession.starts_with("NR_") {
        Biotype::NonCodingRna
    } else {
        Biotype::Unknown
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// A MANE build of `input` into `output`, release 1.5, no summary.
    fn mane_request<'a>(input: &'a Path, output: &'a Path) -> BuildRequest<'a> {
        BuildRequest {
            source: TranscriptSource::Mane,
            input,
            summary_input: None,
            patch_aliases_input: None,
            output_dir: output,
            release: "1.5",
        }
    }

    /// Parse `path` as a MANE GFF3 and return the admitted models.
    fn parse_mane(path: &Path) -> Result<Vec<TranscriptModel>> {
        parse_gff3(path, TranscriptSource::Mane).map(|parsed| parsed.transcripts)
    }

    fn write_temp(contents: &str) -> tempfile::NamedTempFile {
        let mut file = tempfile::Builder::new()
            .suffix(".gff3")
            .tempfile()
            .expect("create tempfile");
        file.write_all(contents.as_bytes()).expect("write");
        file
    }

    // Global tracing subscriber for tests. A per-test thread-local
    // subscriber would be cleaner, but `tracing` caches per-callsite
    // `Interest` at the first event and short-circuits later subscribers,
    // so we install one global capture at first use via `OnceLock`.
    static TRACING_CAPTURE_BUF: std::sync::LazyLock<std::sync::Mutex<Vec<u8>>> =
        std::sync::LazyLock::new(|| std::sync::Mutex::new(Vec::new()));

    /// Install the global capture subscriber once per test binary.
    fn ensure_tracing_capture_installed() {
        use std::io::Write;
        use tracing_subscriber::fmt::MakeWriter;
        use tracing_subscriber::layer::SubscriberExt;

        static INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();

        struct CaptureWriter;
        impl Write for CaptureWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                TRACING_CAPTURE_BUF.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        impl<'a> MakeWriter<'a> for CaptureWriter {
            type Writer = Self;
            fn make_writer(&'a self) -> Self::Writer {
                CaptureWriter
            }
        }

        INIT.get_or_init(|| {
            let layer = tracing_subscriber::fmt::layer()
                .with_writer(CaptureWriter)
                .without_time()
                .with_target(false)
                .with_ansi(false);
            let subscriber = tracing_subscriber::registry().with(layer);
            // `set_global_default` can only be called once per process.
            // `try_init`-style ignore if another binary already set one.
            let _ = tracing::subscriber::set_global_default(subscriber);
        });
    }

    /// Run `f` with tracing capture enabled and return `(result, lines
    /// captured during f)`. Uses per-call sentinel markers so parallel
    /// tests isolate their output from the shared buffer.
    fn with_captured_logs<R>(f: impl FnOnce() -> R) -> (R, Vec<String>) {
        ensure_tracing_capture_installed();
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let begin = format!("CAPTURE_BEGIN_{id}");
        let end = format!("CAPTURE_END_{id}");

        tracing::warn!(sentinel = %begin, "capture begin");
        let result = f();
        tracing::warn!(sentinel = %end, "capture end");

        let full = String::from_utf8_lossy(&TRACING_CAPTURE_BUF.lock().unwrap()).into_owned();
        let mut between: Vec<String> = Vec::new();
        let mut inside = false;
        for line in full.lines() {
            if line.contains(&begin) {
                inside = true;
                continue;
            }
            if line.contains(&end) {
                inside = false;
                continue;
            }
            if inside {
                between.push(line.to_string());
            }
        }
        (result, between)
    }

    /// Minimal GFF3 fixture: two genes on opposite strands, one MANE
    /// Select (minus strand, 3 exons) and one MANE Plus Clinical (plus
    /// strand, 2 exons), one MANE Select non-coding `lnc_RNA`, one MANE
    /// Select row with an *uncommon* feature type (`vault_RNA`) to
    /// exercise the KNOWN_TRANSCRIPT_FEATURE_TYPES warning branch, and
    /// one untagged `lnc_RNA` row that must be dropped regardless of
    /// feature type. The TP53 CDS Dbxref uses `GenBank:` (capital B,
    /// matching real MANE v1.5) to exercise the case-insensitive helper.
    const MINI_GFF: &str = "\
##gff-version 3
chr6\tNCBI\tgene\t100\t2000\t.\t-\t.\tID=gene-SYNGAP1;Name=SYNGAP1;Dbxref=GeneID:8831,HGNC:HGNC:11497;gene_biotype=protein_coding
chr6\tNCBI\tmRNA\t100\t2000\t.\t-\t.\tID=rna-NM_006772.2;Parent=gene-SYNGAP1;Name=NM_006772.2;tag=MANE Select;Dbxref=GenBank:NM_006772.2,Ensembl:ENST00000418600.6
chr6\tNCBI\texon\t100\t500\t.\t-\t.\tID=e1;Parent=rna-NM_006772.2
chr6\tNCBI\texon\t600\t1000\t.\t-\t.\tID=e2;Parent=rna-NM_006772.2
chr6\tNCBI\texon\t1100\t2000\t.\t-\t.\tID=e3;Parent=rna-NM_006772.2
chr6\tNCBI\tCDS\t200\t500\t.\t-\t1\tID=c1;Parent=rna-NM_006772.2;protein_id=NP_006763.2
chr6\tNCBI\tCDS\t1100\t1800\t.\t-\t.\tID=c2;Parent=rna-NM_006772.2;protein_id=NP_006763.2
chr17\tNCBI\tgene\t5000\t8000\t.\t+\t.\tID=gene-TP53;Name=TP53;Dbxref=HGNC:HGNC:11998;gene_biotype=protein_coding
chr17\tNCBI\tmRNA\t5000\t8000\t.\t+\t.\tID=rna-NM_000546.6;Parent=gene-TP53;Name=NM_000546.6;tag=MANE Plus Clinical;Dbxref=Ensembl:ENST00000269305.9
chr17\tNCBI\texon\t5000\t6000\t.\t+\t.\tID=te1;Parent=rna-NM_000546.6
chr17\tNCBI\texon\t7000\t8000\t.\t+\t.\tID=te2;Parent=rna-NM_000546.6
chr17\tNCBI\tCDS\t5500\t6000\t.\t+\t0\tID=tc1;Parent=rna-NM_000546.6;Dbxref=GenBank:NP_000537.3
chr17\tNCBI\tCDS\t7000\t7500\t.\t+\t0\tID=tc2;Parent=rna-NM_000546.6;Dbxref=GenBank:NP_000537.3
chr11\tNCBI\tgene\t3000\t4000\t.\t+\t.\tID=gene-LNCX;Name=LNCX;Dbxref=HGNC:HGNC:99991;gene_biotype=lncRNA
chr11\tNCBI\tlnc_RNA\t3000\t4000\t.\t+\t.\tID=rna-NR_111111.1;Parent=gene-LNCX;Name=NR_111111.1;tag=MANE Select
chr11\tNCBI\texon\t3000\t3400\t.\t+\t.\tID=le1;Parent=rna-NR_111111.1
chr11\tNCBI\texon\t3700\t4000\t.\t+\t.\tID=le2;Parent=rna-NR_111111.1
chr5\tNCBI\tgene\t9000\t9500\t.\t+\t.\tID=gene-VAULT1;Name=VAULT1;Dbxref=HGNC:HGNC:99992;gene_biotype=vault_RNA
chr5\tNCBI\tvault_RNA\t9000\t9500\t.\t+\t.\tID=rna-NR_222222.1;Parent=gene-VAULT1;Name=NR_222222.1;tag=MANE Select
chr5\tNCBI\texon\t9000\t9500\t.\t+\t.\tID=ve1;Parent=rna-NR_222222.1
chr9_KN196479v1_fix\tNCBI\tgene\t500\t1500\t.\t+\t.\tID=gene-FIXED;Name=FIXED;Dbxref=HGNC:HGNC:99993;gene_biotype=protein_coding
chr9_KN196479v1_fix\tNCBI\tmRNA\t500\t1500\t.\t+\t.\tID=rna-NM_777777.1;Parent=gene-FIXED;Name=NM_777777.1;tag=MANE Select
chr9_KN196479v1_fix\tNCBI\texon\t500\t1500\t.\t+\t.\tID=fe1;Parent=rna-NM_777777.1
chr9_KN196479v1_fix\tNCBI\tCDS\t600\t1400\t.\t+\t0\tID=fc1;Parent=rna-NM_777777.1;protein_id=NP_777777.1
chr1\tNCBI\tgene\t10\t100\t.\t+\t.\tID=gene-OTHER;Name=OTHER
chr1\tNCBI\tlnc_RNA\t10\t100\t.\t+\t.\tID=rna-XR_999.1;Parent=gene-OTHER;Name=XR_999.1
chr1\tNCBI\texon\t10\t100\t.\t+\t.\tID=oe1;Parent=rna-XR_999.1
";

    #[test]
    fn build_writes_a_versioned_store_with_provenance() {
        use sha2::{Digest, Sha256};

        let input = write_temp(MINI_GFF);
        let output = tempfile::tempdir().expect("tempdir");
        let (built, size) = build(&mane_request(input.path(), output.path())).expect("build");

        let bin_path = output.path().join("transcript_models.bin");
        let bin = std::fs::read(&bin_path).expect("read bin");
        assert_eq!(size, bin.len());
        assert_eq!(
            &bin[..vareffect::transcript::TRANSCRIPT_STORE_HEADER_LEN],
            b"VEFXTXM\0\x02\x00\x00\x00"
        );
        let store = vareffect::TranscriptStore::load_from_path(&bin_path).expect("load");
        assert_eq!(store.len() as u64, built.record_count);
        assert_eq!(store.source(), Some(TranscriptSource::Mane));
        assert_eq!(store.source_release(), Some("1.5"));

        let manifest: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(output.path().join("transcript_models.manifest.json"))
                .expect("read manifest"),
        )
        .expect("parse manifest");
        let sha = |bytes: &[u8]| hex::encode(Sha256::digest(bytes));
        assert_eq!(manifest["sha256"], sha(&bin));
        assert_eq!(manifest["payload_sha256"], sha(&bin[12..]));
        assert_eq!(manifest["format_version"], 2);
        assert_eq!(
            manifest["source_file_sha256"],
            sha(&std::fs::read(input.path()).expect("read input"))
        );
    }

    /// A one-exon model of `accession` on `chrom` at `start..start + 1000`.
    fn copy_at(accession: &str, chrom: &str, start: u64) -> TranscriptModel {
        let exon = Exon {
            exon_number: 1,
            genomic_start: start,
            genomic_end: start + 1_000,
        };
        TranscriptModel::new(accession, "G", chrom, Strand::Plus, vec![exon])
    }

    #[test]
    fn repeated_accessions_are_admitted_only_as_a_par_pair() {
        let par_pair = [
            copy_at("NM_1.1", "chrX", 1_500_000),
            copy_at("NM_1.1", "chrY", 1_500_000),
            copy_at("NM_2.1", "chrX", 155_800_000),
            copy_at("NM_2.1", "chrY", 56_900_000),
            // chrX copy crossing into PAR2 (RefSeq SPRY3).
            copy_at("NM_8.1", "chrX", 155_700_800),
            copy_at("NM_8.1", "chrY", 56_950_000),
            copy_at("NM_3.1", "chr1", 1_000),
        ];
        check_repeated_accessions(&par_pair).expect("PAR pairs are admitted");

        let refused = [
            // Same chromosome twice.
            copy_at("NM_4.1", "chrX", 1_500_000),
            copy_at("NM_4.1", "chrX", 1_600_000),
            // Autosomes.
            copy_at("NM_5.1", "chr1", 1_000),
            copy_at("NM_5.1", "chr2", 1_000),
            // chrY copy outside the PAR.
            copy_at("NM_6.1", "chrX", 1_500_000),
            copy_at("NM_6.1", "chrY", 20_000_000),
            // Three copies.
            copy_at("NM_7.1", "chrX", 1_500_000),
            copy_at("NM_7.1", "chrY", 1_500_000),
            copy_at("NM_7.1", "chr1", 1_000),
        ];
        let error = check_repeated_accessions(&refused).expect_err("refused");
        let message = error.to_string();
        assert!(message.starts_with("4 accession(s)"), "{message}");
        for accession in ["NM_4.1", "NM_5.1", "NM_6.1", "NM_7.1"] {
            assert!(message.contains(accession), "{message}");
        }
    }

    #[test]
    fn build_refuses_a_gff3_that_admits_no_transcripts() {
        let input = write_temp(
            "chr1\tNCBI\tgene\t10\t100\t.\t+\t.\tID=gene-OTHER;Name=OTHER\n\
             chr1\tNCBI\tlnc_RNA\t10\t100\t.\t+\t.\tID=rna-XR_999.1;Parent=gene-OTHER;Name=XR_999.1\n",
        );
        let output = tempfile::tempdir().expect("tempdir");
        let error = build(&mane_request(input.path(), output.path())).expect_err("empty");
        assert!(
            error.to_string().contains("no transcripts admitted"),
            "{error}"
        );
        assert!(!output.path().join("transcript_models.bin").exists());
    }

    #[test]
    fn parses_mane_select_with_reversed_exons_and_cds() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");

        let syngap = transcripts
            .iter()
            .find(|t| t.accession == "NM_006772.2")
            .expect("SYNGAP1 transcript present");

        assert_eq!(syngap.chrom, "chr6");
        assert_eq!(syngap.strand, Strand::Minus);
        assert_eq!(syngap.tier, Some(TranscriptTier::ManeSelect));
        assert_eq!(syngap.gene_symbol, "SYNGAP1");
        assert_eq!(syngap.hgnc_id.as_deref(), Some("HGNC:11497"));
        assert_eq!(
            syngap.ensembl_accession.as_deref(),
            Some("ENST00000418600.6")
        );
        assert_eq!(syngap.protein_accession.as_deref(), Some("NP_006763.2"));
        assert_eq!(syngap.biotype, Biotype::ProteinCoding);
        assert_eq!(syngap.exon_count, 3);

        assert_eq!(syngap.tx_start, 99);
        assert_eq!(syngap.tx_end, 2000);

        assert_eq!(syngap.cds_genomic_start, Some(199));
        assert_eq!(syngap.cds_genomic_end, Some(1800));

        assert_eq!(syngap.cds_segments.len(), 2);
        assert_eq!(syngap.cds_segments[0].genomic_start, 1099);
        assert_eq!(syngap.cds_segments[0].genomic_end, 1800);
        assert_eq!(syngap.cds_segments[0].exon_index, 0);
        assert_eq!(syngap.cds_segments[1].genomic_start, 199);
        assert_eq!(syngap.cds_segments[1].genomic_end, 500);
        assert_eq!(syngap.cds_segments[1].exon_index, 2);

        assert_eq!(syngap.exons[0].exon_number, 1);
        assert_eq!(syngap.exons[0].genomic_start, 1099);
        assert_eq!(syngap.exons[0].genomic_end, 2000);
        assert!(
            syngap.exons[0].genomic_start > syngap.exons[1].genomic_start,
            "minus-strand exon[0] must be at highest genomic start"
        );
        assert_eq!(syngap.exons.last().unwrap().exon_number, 3);
        assert_eq!(syngap.exons.last().unwrap().genomic_start, 99);
    }

    #[test]
    fn mane_profile_reads_selenocysteine_transl_except() {
        // Shaped on MANE 1.5 NM_000581.4 (GPX1): every CDS row repeats the
        // transcript's `transl_except`, percent-encoded. An `aa:Met` entry
        // (alternative start) is counted, not stored.
        let gff = "\
##gff-version 3
chr3\tNCBI\tgene\t100\t2000\t.\t-\t.\tID=gene-GPX1;Name=GPX1;Dbxref=HGNC:HGNC:4553;gene_biotype=protein_coding
chr3\tNCBI\tmRNA\t100\t2000\t.\t-\t.\tID=rna-NM_000581.4;Parent=gene-GPX1;Name=NM_000581.4;tag=MANE Select
chr3\tNCBI\texon\t100\t500\t.\t-\t.\tID=e1;Parent=rna-NM_000581.4
chr3\tNCBI\texon\t1100\t2000\t.\t-\t.\tID=e2;Parent=rna-NM_000581.4
chr3\tNCBI\tCDS\t200\t500\t.\t-\t1\tID=c1;Parent=rna-NM_000581.4;protein_id=NP_000572.2;transl_except=(pos:complement(1200..1202)%2Caa:Sec)
chr3\tNCBI\tCDS\t1100\t1800\t.\t-\t0\tID=c1;Parent=rna-NM_000581.4;protein_id=NP_000572.2;transl_except=(pos:complement(1200..1202)%2Caa:Sec)
chr3\tNCBI\tgene\t5000\t6000\t.\t+\t.\tID=gene-ALT;Name=ALT;Dbxref=HGNC:HGNC:99994;gene_biotype=protein_coding
chr3\tNCBI\tmRNA\t5000\t6000\t.\t+\t.\tID=rna-NM_888888.1;Parent=gene-ALT;Name=NM_888888.1;tag=MANE Select
chr3\tNCBI\texon\t5000\t6000\t.\t+\t.\tID=ae1;Parent=rna-NM_888888.1
chr3\tNCBI\tCDS\t5100\t5900\t.\t+\t0\tID=ac1;Parent=rna-NM_888888.1;protein_id=NP_888888.1;transl_except=(pos:5100..5102%2Caa:Met)
";
        let file = write_temp(gff);
        let parsed = parse_gff3(file.path(), TranscriptSource::Mane).expect("parse");

        let gpx1 = parsed
            .transcripts
            .iter()
            .find(|t| t.accession == "NM_000581.4")
            .expect("GPX1 admitted");
        assert_eq!(gpx1.selenocysteine_codons, vec![1199]);
        let alt = parsed
            .transcripts
            .iter()
            .find(|t| t.accession == "NM_888888.1")
            .expect("alternative-start transcript admitted");
        assert!(alt.selenocysteine_codons.is_empty());
        assert_eq!(parsed.other_transl_except, 1);
    }

    #[test]
    fn parses_mane_plus_clinical_plus_strand() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");

        let tp53 = transcripts
            .iter()
            .find(|t| t.accession == "NM_000546.6")
            .expect("TP53 transcript present");

        assert_eq!(tp53.strand, Strand::Plus);
        assert_eq!(tp53.tier, Some(TranscriptTier::ManePlusClinical));
        assert_eq!(tp53.chrom, "chr17");

        assert_eq!(tp53.protein_accession.as_deref(), Some("NP_000537.3"));
        assert_eq!(tp53.ensembl_accession.as_deref(), Some("ENST00000269305.9"));

        assert_eq!(tp53.exons[0].genomic_start, 4999);
        assert_eq!(tp53.exons[0].genomic_end, 6000);
        assert_eq!(tp53.exons[1].genomic_start, 6999);
        assert_eq!(tp53.exons[1].genomic_end, 8000);

        assert_eq!(tp53.cds_segments.len(), 2);
        assert_eq!(tp53.cds_segments[0].genomic_start, 5499);
        assert_eq!(tp53.cds_segments[0].genomic_end, 6000);
        assert_eq!(tp53.cds_segments[0].exon_index, 0);
        assert_eq!(tp53.cds_segments[0].phase, 0);
        assert_eq!(tp53.cds_segments[1].genomic_start, 6999);
        assert_eq!(tp53.cds_segments[1].genomic_end, 7500);
        assert_eq!(tp53.cds_segments[1].exon_index, 1);
        assert_eq!(tp53.cds_segments[1].phase, 0);

        assert_eq!(tp53.cds_genomic_start, Some(5499));
        assert_eq!(tp53.cds_genomic_end, Some(7500));
    }

    #[test]
    fn non_mane_row_is_dropped_regardless_of_feature_type() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");

        assert!(
            transcripts.iter().all(|t| t.accession != "XR_999.1"),
            "untagged lnc_RNA must be dropped"
        );
        assert_eq!(transcripts.len(), 5);
    }

    /// One MANE transcript in PAR1 on chrY (listed first) and chrX.
    const PAR_GFF: &str = "\
##gff-version 3
chrY\tNCBI\tgene\t1000000\t1002000\t.\t+\t.\tID=gene-PARG-Y;Name=PARG;gene_biotype=protein_coding
chrY\tNCBI\tmRNA\t1000000\t1002000\t.\t+\t.\tID=rna-NM_123456.1-2;Parent=gene-PARG-Y;Name=NM_123456.1;tag=MANE Select
chrY\tNCBI\texon\t1000000\t1002000\t.\t+\t.\tID=ye1;Parent=rna-NM_123456.1-2
chrX\tNCBI\tgene\t1000000\t1002000\t.\t+\t.\tID=gene-PARG;Name=PARG;gene_biotype=protein_coding
chrX\tNCBI\tmRNA\t1000000\t1002000\t.\t+\t.\tID=rna-NM_123456.1;Parent=gene-PARG;Name=NM_123456.1;tag=MANE Select
chrX\tNCBI\texon\t1000000\t1002000\t.\t+\t.\tID=xe1;Parent=rna-NM_123456.1
";

    #[test]
    fn par_copies_are_built_in_canonical_order_with_stable_bytes() {
        let input = write_temp(PAR_GFF);
        let chroms: Vec<String> = parse_mane(input.path())
            .expect("parse")
            .into_iter()
            .map(|t| t.chrom)
            .collect();
        assert_eq!(chroms, vec!["chrX", "chrY"]);
        // Each build hashes `mrna_info` with a fresh random state.
        let builds: Vec<Vec<u8>> = (0..8)
            .map(|_| {
                let output = tempfile::tempdir().expect("tempdir");
                build(&mane_request(input.path(), output.path())).expect("build");
                std::fs::read(output.path().join("transcript_models.bin")).expect("read")
            })
            .collect();
        assert!(builds.windows(2).all(|pair| pair[0] == pair[1]));
    }

    #[test]
    fn output_is_sorted_by_accession_for_determinism() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");
        let accessions: Vec<&str> = transcripts.iter().map(|t| t.accession.as_str()).collect();
        assert_eq!(
            accessions,
            vec![
                "NM_000546.6",
                "NM_006772.2",
                "NM_777777.1",
                "NR_111111.1",
                "NR_222222.1",
            ]
        );
    }

    #[test]
    fn parses_non_coding_lnc_rna_transcript() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");

        let lncx = transcripts
            .iter()
            .find(|t| t.accession == "NR_111111.1")
            .expect("non-coding lnc_RNA transcript present");

        assert_eq!(lncx.chrom, "chr11");
        assert_eq!(lncx.strand, Strand::Plus);
        assert_eq!(lncx.tier, Some(TranscriptTier::ManeSelect));
        assert_eq!(lncx.gene_symbol, "LNCX");
        assert_eq!(lncx.hgnc_id.as_deref(), Some("HGNC:99991"));
        assert!(lncx.cds_genomic_start.is_none());
        assert!(lncx.cds_genomic_end.is_none());
        assert!(lncx.cds_segments.is_empty());
        assert!(lncx.protein_accession.is_none());
        assert_eq!(lncx.biotype, Biotype::LncRna);
        assert_eq!(lncx.exon_count, 2);
    }

    #[test]
    fn parses_uncommon_feature_type_vault_rna() {
        let file = write_temp(MINI_GFF);
        let (transcripts, logs) = with_captured_logs(|| parse_mane(file.path()).expect("parse"));

        let vault = transcripts
            .iter()
            .find(|t| t.accession == "NR_222222.1")
            .expect("vault_RNA transcript present");

        assert_eq!(vault.chrom, "chr5");
        assert_eq!(vault.strand, Strand::Plus);
        assert_eq!(vault.tier, Some(TranscriptTier::ManeSelect));
        assert_eq!(vault.gene_symbol, "VAULT1");
        assert!(vault.cds_genomic_start.is_none());
        assert!(vault.cds_segments.is_empty());
        assert!(vault.protein_accession.is_none());
        assert_eq!(vault.biotype, Biotype::VaultRna);
        assert_eq!(vault.exon_count, 1);

        assert!(
            logs.iter()
                .any(|l| l.contains("vault_RNA") && l.contains("uncommon feature type")),
            "expected uncommon-feature-type warning for vault_RNA, got: {logs:#?}"
        );
    }

    #[test]
    fn parses_patch_sequence_transcript() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_mane(file.path()).expect("parse");

        let fixed = transcripts
            .iter()
            .find(|t| t.accession == "NM_777777.1")
            .expect("patch-sequence transcript present");

        assert_eq!(fixed.chrom, "chr9_KN196479v1_fix");
        assert!(
            vareffect::chrom::is_patch_sequence(&fixed.chrom),
            "patch contig should be flagged by is_patch_sequence"
        );

        assert_eq!(fixed.strand, Strand::Plus);
        assert_eq!(fixed.tier, Some(TranscriptTier::ManeSelect));
        assert_eq!(fixed.gene_symbol, "FIXED");
        assert_eq!(fixed.biotype, Biotype::ProteinCoding);
        assert_eq!(fixed.exon_count, 1);

        assert_eq!(fixed.exons.len(), 1);
        assert_eq!(fixed.exons[0].genomic_start, 499);
        assert_eq!(fixed.exons[0].genomic_end, 1500);

        assert_eq!(fixed.cds_segments.len(), 1);
        assert_eq!(fixed.cds_segments[0].exon_index, 0);
        assert_eq!(fixed.cds_segments[0].phase, 0);
        assert_eq!(fixed.cds_segments[0].genomic_start, 599);
        assert_eq!(fixed.cds_segments[0].genomic_end, 1400);

        assert_eq!(fixed.cds_genomic_start, Some(599));
        assert_eq!(fixed.cds_genomic_end, Some(1400));

        assert_eq!(fixed.protein_accession.as_deref(), Some("NP_777777.1"));
    }

    #[test]
    fn biotype_derivation_from_accession_prefix() {
        assert_eq!(biotype_for("NM_006772.2", None), Biotype::ProteinCoding);
        assert_eq!(biotype_for("NR_001234.1", None), Biotype::NonCodingRna);
        assert_eq!(biotype_for("XR_999.1", None), Biotype::Unknown);
    }

    #[test]
    fn biotype_derivation_gene_biotype_overrides() {
        assert_eq!(biotype_for("NR_001234.1", Some("lncRNA")), Biotype::LncRna);
        assert_eq!(
            biotype_for("NM_006772.2", Some("protein_coding")),
            Biotype::ProteinCoding
        );
        match biotype_for("NR_001234.1", Some("misc_RNA")) {
            Biotype::Other(raw) => assert_eq!(raw, "misc_RNA"),
            other => panic!("expected Biotype::Other, got {other:?}"),
        }
    }

    #[test]
    fn tier_from_attrs_detects_both_tiers() {
        assert_eq!(
            tier_from_attrs("tag=MANE Select"),
            Some(TranscriptTier::ManeSelect)
        );
        assert_eq!(
            tier_from_attrs("tag=MANE Plus Clinical"),
            Some(TranscriptTier::ManePlusClinical)
        );
        assert_eq!(tier_from_attrs("tag=Basic"), None);
    }

    #[test]
    fn tier_from_attrs_accepts_comma_separated_tag_list() {
        assert_eq!(
            tier_from_attrs("ID=rna-x;tag=basic,MANE Select;Name=foo"),
            Some(TranscriptTier::ManeSelect)
        );
    }

    #[test]
    fn tier_from_attrs_rejects_substring_false_positive() {
        let attrs = "ID=rna-x;Note=similar to MANE_Select isoforms;tag=basic";
        assert_eq!(tier_from_attrs(attrs), None);
    }

    #[test]
    fn parse_phase_accepts_valid_values() {
        assert_eq!(parse_phase("0").unwrap(), 0);
        assert_eq!(parse_phase("1").unwrap(), 1);
        assert_eq!(parse_phase("2").unwrap(), 2);
        assert_eq!(parse_phase(".").unwrap(), 0);
        assert_eq!(parse_phase("").unwrap(), 0);
    }

    #[test]
    fn parse_phase_rejects_invalid_values() {
        assert!(parse_phase("3").is_err());
        assert!(parse_phase("-1").is_err());
        assert!(parse_phase("abc").is_err());
    }

    #[test]
    fn extract_first_parent_splits_multi_valued_parents() {
        assert_eq!(
            extract_first_parent("ID=x;Parent=gene-A").as_deref(),
            Some("gene-A")
        );
        assert_eq!(
            extract_first_parent("ID=x;Parent=gene-A,gene-B").as_deref(),
            Some("gene-A"),
            "multi-valued Parent should return the first ID, not the full string"
        );
        assert_eq!(extract_first_parent("ID=x"), None);
    }

    #[test]
    fn build_cds_segments_rejects_misaligned_cds() {
        let exons = vec![
            Exon {
                exon_number: 1,
                genomic_start: 100,
                genomic_end: 200,
            },
            Exon {
                exon_number: 2,
                genomic_start: 300,
                genomic_end: 400,
            },
        ];
        let raw = vec![CdsRawSegment {
            gff_start: 150,
            gff_end: 350,
            phase: 0,
            start_range: false,
            end_range: false,
        }];
        let err = build_cds_segments("NM_TEST.1", &exons, &Strand::Plus, Some(&raw)).unwrap_err();
        assert!(
            format!("{err:#}").contains("does not fit inside any exon"),
            "expected intron-span error, got {err:#}"
        );
    }

    /// A RefSeq CDS row at `start..=end` (1-based) with its range markers.
    fn ranged(start: u64, end: u64, start_range: bool, end_range: bool) -> CdsRawSegment {
        CdsRawSegment {
            gff_start: start,
            gff_end: end,
            phase: 0,
            start_range,
            end_range,
        }
    }

    #[test]
    fn interior_ranges_are_coverage_breaks_not_open_ends() {
        // SHANK3 NM_001372044.2 (+): `end_range` on an interior row, then a
        // `start_range` on the next one; the CDS ends are intact.
        let shank3 = [
            ranged(50_683_340, 50_683_400, false, false),
            ranged(50_694_775, 50_695_046, false, true),
            ranged(50_697_557, 50_697_561, true, false),
            ranged(50_698_690, 50_698_800, false, false),
        ];
        let ranges = CdsRanges::of(&shank3);
        assert!(ranges.interior && ranges.any);
        assert!(!ranges.low_partial && !ranges.high_partial);
        // A range on the lowest row is an open low end.
        let open = [
            ranged(100, 200, true, false),
            ranged(300, 400, false, false),
        ];
        assert_eq!(
            CdsRanges::of(&open),
            CdsRanges {
                low_partial: true,
                high_partial: false,
                interior: false,
                any: true
            }
        );
        assert_eq!(CdsRanges::of(&[]), CdsRanges::default());
    }

    /// RefSeq rows for `NM_000001.1` (+, CDS 101..=400 over two exons) with
    /// `extra` appended.
    fn refseq_gff(extra: &str) -> String {
        format!(
            "NC_000001.11\tBestRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene-G1;Name=G1;gene_biotype=protein_coding\n\
             NC_000001.11\tBestRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=rna-NM_000001.1;Parent=gene-G1;gene=G1;transcript_id=NM_000001.1\n\
             NC_000001.11\tBestRefSeq\texon\t1\t250\t.\t+\t.\tID=e1;Parent=rna-NM_000001.1\n\
             NC_000001.11\tBestRefSeq\texon\t301\t500\t.\t+\t.\tID=e2;Parent=rna-NM_000001.1\n\
             NC_000001.11\tBestRefSeq\tCDS\t101\t250\t.\t+\t0\tID=c1;Parent=rna-NM_000001.1;protein_id=NP_000001.1\n\
             NC_000001.11\tBestRefSeq\tCDS\t301\t450\t.\t+\t0\tID=c1;Parent=rna-NM_000001.1;protein_id=NP_000001.1\n\
             {extra}"
        )
    }

    fn parse_refseq(gff: &str) -> ParsedSource {
        let file = write_temp(gff);
        parse_gff3(file.path(), TranscriptSource::RefSeq).expect("parse")
    }

    #[test]
    fn a_coverage_break_between_match_blocks_flags_the_cds() {
        // Transcript bases 251..=260 are aligned nowhere: the blocks are
        // 1..=250 and 261..=460 (exons 1..=250 and 301..=500).
        let parsed = parse_refseq(&refseq_gff(
            "NC_000001.11\tRefSeq\tcDNA_match\t1\t250\t.\t+\t.\tID=aln1;Target=NM_000001.1 1 250 +\n\
             NC_000001.11\tRefSeq\tcDNA_match\t301\t500\t.\t+\t.\tID=aln1;Target=NM_000001.1 261 460 +\n",
        ));
        let model = &parsed.transcripts[0];
        assert!(model.genome_mismatch && model.cds_alignment_gap);

        let contiguous = parse_refseq(&refseq_gff(
            "NC_000001.11\tRefSeq\tcDNA_match\t1\t250\t.\t+\t.\tID=aln1;Target=NM_000001.1 1 250 +\n\
             NC_000001.11\tRefSeq\tcDNA_match\t301\t500\t.\t+\t.\tID=aln1;Target=NM_000001.1 251 450 +\n",
        ));
        let model = &contiguous.transcripts[0];
        assert!(!model.genome_mismatch && !model.cds_alignment_gap);
    }

    #[test]
    fn a_malformed_child_row_excludes_its_transcript() {
        let parsed = parse_refseq(&refseq_gff(
            "NC_000001.11\tBestRefSeq\texon\tx\t600\t.\t+\t.\tID=e3;Parent=rna-NM_000001.1\n\
             NC_000001.11\tBestRefSeq\texon\t700\t800\t.\t+\t.\tID=e4\n",
        ));
        assert!(parsed.transcripts.is_empty());
        assert_eq!(
            parsed.excluded[0].reason,
            ExclusionReason::MalformedChildRow
        );
        assert_eq!(parsed.skipped_rows.get("unparsable_coordinates"), Some(&1));
        assert_eq!(parsed.skipped_rows.get("child_without_parent"), Some(&1));
    }

    #[test]
    fn a_repeated_transcript_id_fails_the_build() {
        let file = write_temp(&refseq_gff(
            "NC_000001.11\tBestRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=rna-NM_000001.1;Parent=gene-G1;transcript_id=NM_000001.1\n",
        ));
        let error = parse_gff3(file.path(), TranscriptSource::RefSeq).expect_err("repeat");
        assert!(error.to_string().contains("repeats"), "{error}");
    }

    #[test]
    fn refseq_rules_admit_primary_assembly_rows_only() {
        let parsed = parse_refseq(&refseq_gff(
            "NT_187361.1\tBestRefSeq\tmRNA\t1\t100\t.\t+\t.\tID=rna-NM_000002.1;Parent=gene-G2;transcript_id=NM_000002.1\n",
        ));
        assert_eq!(parsed.candidates, 2);
        assert_eq!(parsed.transcripts.len(), 1);
        assert_eq!(
            parsed.excluded[0].reason,
            ExclusionReason::NonPrimaryAssembly
        );
        let model = &parsed.transcripts[0];
        assert_eq!(model.chrom, "chr1");
        assert_eq!(model.biotype, Biotype::ProteinCoding);
        assert_eq!(model.protein_accession.as_deref(), Some("NP_000001.1"));
        assert!(model.cds_start_complete && model.cds_end_complete);
        assert!(model.tier.is_none());
    }

    #[test]
    fn gencode_rules_read_tags_and_selenocysteine() {
        let gff = "\
chr1\tHAVANA\ttranscript\t1\t500\t.\t-\t.\tID=ENST1.1;Parent=ENSG1.1;transcript_id=ENST1.1;gene_name=G1;transcript_type=protein_coding;protein_id=ENSP1.1;hgnc_id=HGNC:1;tag=basic,cds_end_NF,MANE_Select
chr1\tHAVANA\texon\t1\t500\t.\t-\t.\tParent=ENST1.1
chr1\tHAVANA\tCDS\t101\t400\t.\t-\t0\tParent=ENST1.1
chr1\tHAVANA\tstop_codon\t101\t103\t.\t-\t0\tParent=ENST1.1
chr1\tHAVANA\tstop_codon_redefined_as_selenocysteine\t200\t202\t.\t-\t.\tParent=CDS:ENST1.1;transcript_id=ENST1.1
";
        let file = write_temp(gff);
        let parsed = parse_gff3(file.path(), TranscriptSource::Gencode).expect("parse");
        let model = &parsed.transcripts[0];
        assert_eq!(model.accession, "ENST1.1");
        assert_eq!(model.tier, Some(TranscriptTier::ManeSelect));
        assert!(model.cds_start_complete && !model.cds_end_complete);
        assert_eq!(model.selenocysteine_codons, vec![199]);
        assert_eq!(model.hgnc_id.as_deref(), Some("HGNC:1"));
        // The stop_codon row is not folded into the CDS.
        assert_eq!(
            (model.cds_genomic_start, model.cds_genomic_end),
            (Some(100), Some(400))
        );
        assert!(model.ensembl_accession.is_none());
    }
}
