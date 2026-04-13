//! MANE transcript model builder for the `vareffect` crate's
//! [`TranscriptStore`](vareffect::TranscriptStore).
//!
//! Parses NCBI's MANE GFF3 release
//! (`MANE.GRCh38.vX.X.refseq_genomic.gff.gz`) into
//! `Vec<vareffect::TranscriptModel>` and writes it as MessagePack via the
//! shared [`serialize_and_finalize`] helper. Every record carries full
//! exon structure, CDS bounds, protein accessions, Ensembl cross-refs,
//! biotype, and MANE tier, in 0-based half-open (BED/UCSC) coordinates.
//!
//! # GFF3 parsing strategy
//!
//! Single pass over the gzipped GFF3, building four scratch maps keyed on
//! feature `ID` attributes:
//!
//! - `gene_info` -- gene symbol, HGNC ID, gene biotype. HGNC is walked via
//!   `mrna.Parent -> gene_info` (authoritative), with the mRNA `Dbxref` as
//!   a resilience fallback.
//! - `mrna_info` -- transcript-level info; only mRNA rows carrying a
//!   `tag=MANE Select` or `tag=MANE Plus Clinical` tag are admitted.
//! - `cds_ranges` -- CDS spans per parent mRNA.
//! - `exon_ranges` -- exon spans per parent mRNA.
//!
//! After the scan, each admitted mRNA is assembled into a
//! [`vareffect::TranscriptModel`] with coordinates converted from GFF3's
//! 1-based fully-closed convention to 0-based half-open (BED/UCSC), exons
//! sorted into 5'->3' *transcript* order (reversed for minus-strand genes),
//! and biotype derived from the accession prefix and optional gene-level
//! `gene_biotype` attribute.
//!
//! If a summary TSV path is supplied, a cross-validation pass compares
//! every row against the built vec and fails the build on mismatches.

mod cross_validate;
mod gff3_attrs;

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use indicatif::{ProgressBar, ProgressStyle};
use vareffect::chrom::{is_patch_sequence, refseq_to_ucsc};
use vareffect::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};

use cross_validate::cross_validate_summary;
use gff3_attrs::{
    extract_attr, extract_ensembl_from_dbxref, extract_hgnc_from_dbxref,
    extract_protein_id_from_cds_attrs, extract_refseq_acc_from_dbxref, extract_transcript_from_id,
};

use crate::common::{BuildOutput, serialize_and_finalize};

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

/// Build the transcript model store from a MANE GFF3 file.
///
/// # Arguments
///
/// * `input` -- Path to `MANE.GRCh38.vX.X.refseq_genomic.gff.gz` (or an
///   uncompressed `.gff3`).
/// * `summary_input` -- Optional path to `MANE.GRCh38.vX.X.summary.txt.gz`.
///   When provided, every built transcript is cross-validated against the
///   summary's chrom, strand, position, gene symbol, and tier. Mismatches
///   are fatal build errors.
/// * `patch_aliases_input` -- Path to `patch_chrom_aliases.csv` produced by
///   [`crate::builders::patch_chrom_aliases::build_from_assembly_report`].
///   Required when `summary_input` is `Some`; reconciles `NW_*`/`NT_*`
///   summary accessions against UCSC patch chrom names.
/// * `output_dir` -- Directory for the output files
///   (`transcript_models.bin`, `transcript_models.bin.sha256`,
///   `transcript_models.manifest.json`).
/// * `version` -- MANE release version string (e.g. `"1.5"`).
///
/// # Returns
///
/// Tuple of ([`BuildOutput`], serialized byte count).
///
/// # Errors
///
/// Returns an error if the GFF3 cannot be read, the summary TSV
/// cross-validation fails, the patch-alias CSV cannot be read, or the
/// output files cannot be written. Also errors if `summary_input` is
/// provided without a matching `patch_aliases_input`.
pub fn build(
    input: &Path,
    summary_input: Option<&Path>,
    patch_aliases_input: Option<&Path>,
    output_dir: &Path,
    version: &str,
) -> Result<(BuildOutput, usize)> {
    let transcripts = parse_gff3(input)?;

    if let Some(summary_path) = summary_input {
        let aliases_path = patch_aliases_input.ok_or_else(|| {
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
        cross_validate_summary(&transcripts, summary_path, &aliases).with_context(|| {
            format!(
                "cross-validating GFF3 build against {}",
                summary_path.display()
            )
        })?;
    }

    let count = transcripts.len() as u64;
    serialize_and_finalize(
        "transcript_models",
        &transcripts,
        count,
        input,
        output_dir,
        version,
    )
}

// ---------------------------------------------------------------------------
// GFF3 parsing
// ---------------------------------------------------------------------------

/// Parse a MANE GFF3 file into `Vec<TranscriptModel>`, sorted by accession
/// for deterministic MessagePack output.
fn parse_gff3(input: &Path) -> Result<Vec<TranscriptModel>> {
    let file =
        std::fs::File::open(input).with_context(|| format!("opening {}", input.display()))?;

    // Transparently handle gzipped and uncompressed inputs.
    let reader: Box<dyn BufRead> = if input.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Scratch maps keyed on GFF3 feature `ID` attributes.
    let mut gene_info: HashMap<String, GeneInfo> = HashMap::new();
    let mut mrna_info: HashMap<String, PendingTranscript> = HashMap::new();
    // CDS scratch: mRNA_id -> Vec of (start, end, phase). Phase comes from
    // GFF3 column 8; missing (`.`) is normalized to 0 at parse time.
    let mut cds_ranges: HashMap<String, Vec<CdsRawSegment>> = HashMap::new();
    let mut exon_ranges: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    // mRNA_id -> first-observed protein accession (from CDS rows).
    let mut mrna_protein: HashMap<String, String> = HashMap::new();
    // Admitted-transcript count per feature type (mRNA, lnc_RNA, ...).
    // Used to surface uncommon feature types in the build log so future
    // MANE releases that introduce new biotypes are visible at a glance.
    let mut feature_type_counts: HashMap<String, usize> = HashMap::new();

    // Spinner-with-rate rather than the previous pinned-at-100 % bar: the
    // GFF3 is ~4 M lines and counting them up front (to feed a bounded
    // progress bar) would double the I/O. A spinner is honest about the
    // unbounded-length case while still giving a live read-rate indicator.
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("  {spinner} parsing MANE GFF3: {pos} lines ({per_sec})")
            .expect("valid progress template"),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(120));

    for line_result in reader.lines() {
        let line = line_result.context("reading GFF3 line")?;
        pb.inc(1);

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let raw_chrom = fields[0];
        let feature_type = fields[2];
        let Ok(gff_start) = fields[3].parse::<u64>() else {
            continue;
        };
        let Ok(gff_end) = fields[4].parse::<u64>() else {
            continue;
        };
        let strand = match fields[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => continue,
        };
        let attrs = fields[8];

        match feature_type {
            // Gene rows: authoritative source for HGNC ID and gene biotype.
            "gene" | "pseudogene" => {
                let Some(id) = extract_attr(attrs, "ID") else {
                    continue;
                };
                let gene_symbol = extract_attr(attrs, "Name")
                    .or_else(|| extract_attr(attrs, "gene"))
                    .unwrap_or_default();
                let hgnc_id = extract_hgnc_from_dbxref(attrs);
                let gene_biotype = extract_attr(attrs, "gene_biotype");
                gene_info.insert(
                    id,
                    GeneInfo {
                        gene_symbol,
                        hgnc_id,
                        gene_biotype,
                    },
                );
            }

            // Top-level region rows are metadata only; skip explicitly so
            // they never reach the transcript catch-all below.
            "region" => {}

            // CDS rows: collect spans, phase, and capture protein accession.
            "CDS" => {
                let Some(parent_mrna) = extract_first_parent(attrs) else {
                    continue;
                };
                // GFF3 column 8 is phase: `0`, `1`, `2`, or `.` (unknown).
                // Reject malformed values loudly -- a CDS row with phase `5`
                // indicates upstream corruption, not something we should
                // silently normalize.
                let phase = parse_phase(fields[7]).with_context(|| {
                    format!(
                        "parsing CDS phase column 8 for parent {parent_mrna}: value {:?}",
                        fields[7]
                    )
                })?;
                cds_ranges
                    .entry(parent_mrna.clone())
                    .or_default()
                    .push(CdsRawSegment {
                        gff_start,
                        gff_end,
                        phase,
                    });

                // Dual-source protein_id: `protein_id=NP_...` attr first,
                // then fall back to `Dbxref=Genbank:NP_...`.
                if !mrna_protein.contains_key(&parent_mrna)
                    && let Some(pid) = extract_protein_id_from_cds_attrs(attrs)
                {
                    mrna_protein.insert(parent_mrna, pid);
                }
            }

            // Exon rows: collect spans.
            "exon" => {
                if let Some(parent_mrna) = extract_first_parent(attrs) {
                    exon_ranges
                        .entry(parent_mrna)
                        .or_default()
                        .push((gff_start, gff_end));
                }
            }

            // Transcript-level catch-all. Any feature type that is not a
            // gene, CDS, exon, or region row is treated as a potential
            // transcript. The MANE tag is the authoritative "this is a MANE
            // transcript" filter -- unknown future biotypes (e.g. vault_RNA
            // in a later release) are admitted automatically without code
            // changes. Rows without a MANE tag are silently dropped.
            //
            // MANE v1.5 uses seven feature types for tagged transcripts:
            // mRNA, lnc_RNA, antisense_RNA, snoRNA, snRNA, RNase_MRP_RNA,
            // telomerase_RNA. Anything outside this set triggers a build-log
            // warning via the `feature_type_counts` loop below.
            _ => {
                let Some(tier) = tier_from_attrs(attrs) else {
                    continue;
                };
                let Some(id) = extract_attr(attrs, "ID") else {
                    continue;
                };
                // GFF3 spec allows `Parent=a,b,c`; transcript rows in real
                // MANE releases carry a single parent, but using the helper
                // keeps a future multi-parent row from silently failing the
                // gene lookup.
                let parent_gene = extract_first_parent(attrs).unwrap_or_default();

                // Chromosome: convert NC_ -> chr. Patch sequences fall through
                // unchanged and are flagged for the build summary.
                let chrom = refseq_to_ucsc(raw_chrom).to_string();

                // RefSeq transcript accession: Name attr -> GenBank Dbxref
                // fallback -> ID-prefix fallback. Every MANE transcript in
                // v1.5 carries `Name=<accession>` directly, so the
                // fallbacks are defensive but not normally exercised.
                let accession = extract_attr(attrs, "Name")
                    .or_else(|| extract_refseq_acc_from_dbxref(attrs))
                    .or_else(|| extract_transcript_from_id(&id))
                    .unwrap_or_default();

                if accession.is_empty() {
                    tracing::warn!(
                        transcript_id = %id,
                        store = "transcript_models",
                        "skipping MANE transcript: no accession could be resolved",
                    );
                    continue;
                }

                // Record the feature type for the end-of-build breakdown
                // and uncommon-type warning.
                *feature_type_counts
                    .entry(feature_type.to_string())
                    .or_insert(0) += 1;

                // Ensembl transcript accession (optional).
                let ensembl_accession = extract_ensembl_from_dbxref(attrs);

                mrna_info.insert(
                    id,
                    PendingTranscript {
                        accession,
                        parent_gene,
                        ensembl_accession,
                        chrom,
                        strand,
                        tx_start_gff: gff_start,
                        tx_end_gff: gff_end,
                        tier,
                    },
                );
            }
        }
    }

    pb.finish_with_message("MANE GFF3 parsed");

    // ---- Assembly pass ----
    let mut transcripts: Vec<TranscriptModel> = Vec::with_capacity(mrna_info.len());
    let mut patch_seq_count: usize = 0;
    let mut hgnc_fallback_count: usize = 0;
    let mut biotype_unknown_count: usize = 0;

    for (mrna_id, pending) in &mrna_info {
        // --- Exon assembly ---
        let Some(raw_exons) = exon_ranges.get(mrna_id) else {
            tracing::warn!(
                accession = %pending.accession,
                store = "transcript_models",
                "skipping MANE transcript: no exon rows found",
            );
            continue;
        };
        if raw_exons.is_empty() {
            tracing::warn!(
                accession = %pending.accession,
                store = "transcript_models",
                "skipping MANE transcript: empty exon list",
            );
            continue;
        }

        // Sort by genomic start ascending. For minus-strand transcripts we
        // reverse the sorted list so `exons[0]` is the 5'-most exon on the
        // transcript (highest genomic coordinate).
        let mut sorted_exons: Vec<(u64, u64)> = raw_exons.clone();
        sorted_exons.sort_by_key(|(s, _)| *s);
        if matches!(pending.strand, Strand::Minus) {
            sorted_exons.reverse();
        }

        // Convert to 0-based half-open and number 1..=N in transcript order.
        if sorted_exons.len() > u16::MAX as usize {
            bail!(
                "transcript {} has {} exons (u16 overflow)",
                pending.accession,
                sorted_exons.len()
            );
        }
        let exons: Vec<Exon> = sorted_exons
            .iter()
            .enumerate()
            .map(|(i, (s, e))| Exon {
                exon_number: (i + 1) as u16,
                // GFF3 uses 1-based fully-closed; BED/UCSC 0-based half-open:
                //   genomic_start = gff_start - 1
                //   genomic_end   = gff_end       (already exclusive)
                genomic_start: s.saturating_sub(1),
                genomic_end: *e,
            })
            .collect();

        let exon_count = exons.len() as u16;

        // --- CDS segments + bounds ---
        // For each CDS row, convert to 0-based half-open, sort into transcript
        // order (ascending genomic on plus strand, descending on minus to
        // match `exons`), and resolve the containing exon index. The global
        // `cds_genomic_start` / `cds_genomic_end` fields are derived from the
        // segments as `min` / `max` so they stay consistent by construction.
        let cds_segments = build_cds_segments(
            &pending.accession,
            &exons,
            &pending.strand,
            cds_ranges.get(mrna_id),
        )?;
        let (cds_genomic_start, cds_genomic_end) = if cds_segments.is_empty() {
            (None, None)
        } else {
            let min = cds_segments
                .iter()
                .map(|s| s.genomic_start)
                .min()
                .expect("cds_segments verified non-empty");
            let max = cds_segments
                .iter()
                .map(|s| s.genomic_end)
                .max()
                .expect("cds_segments verified non-empty");
            (Some(min), Some(max))
        };

        // --- tx bounds (convert 1-based -> 0-based) ---
        let tx_start = pending.tx_start_gff.saturating_sub(1);
        let tx_end = pending.tx_end_gff;

        // --- Gene-level lookup for HGNC, symbol, biotype ---
        let parent_info = gene_info.get(&pending.parent_gene);

        let gene_symbol = parent_info
            .map(|g| g.gene_symbol.clone())
            .filter(|s| !s.is_empty())
            .unwrap_or_default();

        let (hgnc_id, used_fallback) = match parent_info.and_then(|g| g.hgnc_id.clone()) {
            Some(id) => (Some(id), false),
            None => (None, true),
        };
        if used_fallback {
            hgnc_fallback_count += 1;
        }

        // Biotype: NM_/NR_ prefix is primary; gene-level `gene_biotype`
        // overrides. `Biotype::Unknown` if neither signal is present.
        let biotype = biotype_for(
            &pending.accession,
            parent_info.and_then(|g| g.gene_biotype.as_deref()),
        );
        if matches!(biotype, Biotype::Unknown) {
            biotype_unknown_count += 1;
        }

        // --- Protein accession ---
        let protein_accession = mrna_protein.get(mrna_id).cloned();

        // --- Patch-sequence flagging ---
        // We only track the count; the aggregate is reported in the final
        // build summary line. Per-transcript warnings are visual noise
        // (~64 lines) that obscure actionable build output.
        if is_patch_sequence(&pending.chrom) {
            patch_seq_count += 1;
        }

        transcripts.push(TranscriptModel {
            accession: pending.accession.clone(),
            protein_accession,
            gene_symbol,
            hgnc_id,
            ensembl_accession: pending.ensembl_accession.clone(),
            chrom: pending.chrom.clone(),
            strand: pending.strand,
            tx_start,
            tx_end,
            cds_genomic_start,
            cds_genomic_end,
            exons,
            cds_segments,
            tier: pending.tier,
            biotype,
            exon_count,
        });
    }

    // Deterministic output: sort by accession so rebuilds with identical
    // inputs produce byte-identical MessagePack (and hence stable sha256).
    transcripts.sort_by(|a, b| a.accession.cmp(&b.accession));

    // Build a deterministic per-feature-type breakdown for the build log.
    // Sorted by feature type for stable output across runs -- the warning
    // loop below reuses the same sorted view, so uncommon feature types
    // always appear in the same order across runs.
    let mut feature_type_breakdown: Vec<(&String, &usize)> = feature_type_counts.iter().collect();
    feature_type_breakdown.sort_by(|a, b| a.0.cmp(b.0));

    // Warn on any admitted feature type outside the known-common set.
    // This is the early-warning for future MANE releases introducing new
    // biotypes: the build still succeeds (the filter is inverted), but the
    // warning makes the change visually obvious.
    for (ft, count) in &feature_type_breakdown {
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
    let breakdown_str = feature_type_breakdown
        .iter()
        .map(|(ft, count)| format!("{ft}={count}"))
        .collect::<Vec<_>>()
        .join(", ");

    eprintln!(
        "  Transcript models: {} total ({} patch-seq, {} HGNC fallback, {} unknown biotype)",
        transcripts.len(),
        patch_seq_count,
        hgnc_fallback_count,
        biotype_unknown_count,
    );
    eprintln!("  By feature type: {breakdown_str}");

    Ok(transcripts)
}

/// Intermediate gene-level info collected from `gene` rows.
#[derive(Debug)]
struct GeneInfo {
    gene_symbol: String,
    hgnc_id: Option<String>,
    gene_biotype: Option<String>,
}

/// Intermediate mRNA-level info collected from `mRNA` rows before
/// exon/CDS assembly.
#[derive(Debug)]
struct PendingTranscript {
    accession: String,
    parent_gene: String,
    ensembl_accession: Option<String>,
    chrom: String,
    strand: Strand,
    tx_start_gff: u64,
    tx_end_gff: u64,
    tier: TranscriptTier,
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

/// Detect MANE tier from an mRNA row's attribute column.
///
/// Returns `None` for untagged rows. Parses the `tag=` attribute
/// structurally -- the previous implementation used `attrs.contains(...)`
/// which would match any substring anywhere in column 9 (for example a
/// `Note=...includes MANE_Select isoforms...` description), silently
/// admitting non-MANE transcripts.
///
/// GFF3 tag values are comma-separated, so a row can carry multiple tags
/// (`tag=basic,MANE Select`). Both the NCBI space form (`MANE Select`)
/// and the underscore form (`MANE_Select`) are accepted so the parser
/// stays compatible with either convention.
fn tier_from_attrs(attrs: &str) -> Option<TranscriptTier> {
    let tag = extract_attr(attrs, "tag")?;
    for value in tag.split(',').map(str::trim) {
        match value {
            "MANE Select" | "MANE_Select" => return Some(TranscriptTier::ManeSelect),
            "MANE Plus Clinical" | "MANE_Plus_Clinical" => {
                return Some(TranscriptTier::ManePlusClinical);
            }
            _ => {}
        }
    }
    None
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
chr6\tNCBI\tCDS\t200\t500\t.\t-\t.\tID=c1;Parent=rna-NM_006772.2;protein_id=NP_006763.2
chr6\tNCBI\tCDS\t1100\t1800\t.\t-\t.\tID=c2;Parent=rna-NM_006772.2;protein_id=NP_006763.2
chr17\tNCBI\tgene\t5000\t8000\t.\t+\t.\tID=gene-TP53;Name=TP53;Dbxref=HGNC:HGNC:11998;gene_biotype=protein_coding
chr17\tNCBI\tmRNA\t5000\t8000\t.\t+\t.\tID=rna-NM_000546.6;Parent=gene-TP53;Name=NM_000546.6;tag=MANE Plus Clinical;Dbxref=Ensembl:ENST00000269305.9
chr17\tNCBI\texon\t5000\t6000\t.\t+\t.\tID=te1;Parent=rna-NM_000546.6
chr17\tNCBI\texon\t7000\t8000\t.\t+\t.\tID=te2;Parent=rna-NM_000546.6
chr17\tNCBI\tCDS\t5500\t6000\t.\t+\t0\tID=tc1;Parent=rna-NM_000546.6;Dbxref=GenBank:NP_000537.3
chr17\tNCBI\tCDS\t7000\t7500\t.\t+\t2\tID=tc2;Parent=rna-NM_000546.6;Dbxref=GenBank:NP_000537.3
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
    fn parses_mane_select_with_reversed_exons_and_cds() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_gff3(file.path()).expect("parse");

        let syngap = transcripts
            .iter()
            .find(|t| t.accession == "NM_006772.2")
            .expect("SYNGAP1 transcript present");

        assert_eq!(syngap.chrom, "chr6");
        assert_eq!(syngap.strand, Strand::Minus);
        assert_eq!(syngap.tier, TranscriptTier::ManeSelect);
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
    fn parses_mane_plus_clinical_plus_strand() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_gff3(file.path()).expect("parse");

        let tp53 = transcripts
            .iter()
            .find(|t| t.accession == "NM_000546.6")
            .expect("TP53 transcript present");

        assert_eq!(tp53.strand, Strand::Plus);
        assert_eq!(tp53.tier, TranscriptTier::ManePlusClinical);
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
        assert_eq!(tp53.cds_segments[1].phase, 2);

        assert_eq!(tp53.cds_genomic_start, Some(5499));
        assert_eq!(tp53.cds_genomic_end, Some(7500));
    }

    #[test]
    fn non_mane_row_is_dropped_regardless_of_feature_type() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_gff3(file.path()).expect("parse");

        assert!(
            transcripts.iter().all(|t| t.accession != "XR_999.1"),
            "untagged lnc_RNA must be dropped"
        );
        assert_eq!(transcripts.len(), 5);
    }

    #[test]
    fn output_is_sorted_by_accession_for_determinism() {
        let file = write_temp(MINI_GFF);
        let transcripts = parse_gff3(file.path()).expect("parse");
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
        let transcripts = parse_gff3(file.path()).expect("parse");

        let lncx = transcripts
            .iter()
            .find(|t| t.accession == "NR_111111.1")
            .expect("non-coding lnc_RNA transcript present");

        assert_eq!(lncx.chrom, "chr11");
        assert_eq!(lncx.strand, Strand::Plus);
        assert_eq!(lncx.tier, TranscriptTier::ManeSelect);
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
        let (transcripts, logs) = with_captured_logs(|| parse_gff3(file.path()).expect("parse"));

        let vault = transcripts
            .iter()
            .find(|t| t.accession == "NR_222222.1")
            .expect("vault_RNA transcript present");

        assert_eq!(vault.chrom, "chr5");
        assert_eq!(vault.strand, Strand::Plus);
        assert_eq!(vault.tier, TranscriptTier::ManeSelect);
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
        let transcripts = parse_gff3(file.path()).expect("parse");

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
        assert_eq!(fixed.tier, TranscriptTier::ManeSelect);
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
        }];
        let err = build_cds_segments("NM_TEST.1", &exons, &Strand::Plus, Some(&raw)).unwrap_err();
        assert!(
            format!("{err:#}").contains("does not fit inside any exon"),
            "expected intron-span error, got {err:#}"
        );
    }
}
