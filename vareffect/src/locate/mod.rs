//! Variant locator: classify where a genomic position falls within a transcript.
//!
//! Given a 0-based half-open genomic coordinate and a [`TranscriptModel`], this
//! module determines the transcript region the variant occupies: CDS exon,
//! intron, UTR, splice site, upstream, or downstream. This is pure coordinate
//! geometry -- no FASTA access, no translation, no allele comparison.
//!
//! # VEP concordance
//!
//! This module replicates the geometric classification from:
//! - `Bio::EnsEMBL::Variation::TranscriptVariation::_calc_coords()`
//! - `Bio::EnsEMBL::Variation::Utils::VariationEffect` splice predicates
//!
//! Splice site definitions match VEP core:
//! - **Donor/Acceptor:** intronic positions +/-1-2 from exon boundary
//! - **Splice region:** intronic +3..+8 on the donor side or -3..-17 on the
//!   acceptor side (the latter covers `splice_polypyrimidine_tract_variant`
//!   since Ensembl release 105), plus exonic 1-3 bases from boundary

pub(crate) mod helpers;
mod indel;
#[cfg(test)]
mod tests;

pub use indel::{IndelLocation, IndelRegion, SpliceOverlapDetail, locate_indel};

use crate::error::VarEffectError;
use crate::types::{Strand, TranscriptModel};
use helpers::{
    ExonOrIntron, classify_intronic, compute_cds_offset, compute_utr_offset_3prime,
    compute_utr_offset_5prime, find_exon_or_intron, is_exonic_splice_region, to_u16,
};

/// VEP default `--distance` for upstream/downstream gene variants.
const UPSTREAM_DOWNSTREAM_LIMIT: u64 = 5_000;

/// Canonical splice donor/acceptor: intronic positions +1/+2 or -1/-2.
const SPLICE_CANONICAL_MAX: u64 = 2;

/// Donor-side splice region window: intronic +3..+8 from the 3' end of the
/// upstream exon. Matches VEP CORE `splice_region_variant` (SO:0001630) and
/// also covers `splice_donor_region_variant` (+3..+6) and
/// `splice_donor_5th_base_variant` (+5) after harness normalization — all
/// three are subsets of +3..+8, so vareffect emits a single collapsed
/// `splice_region_variant`.
const SPLICE_REGION_INTRON_MAX_DONOR: u64 = 8;

/// Acceptor-side splice region window: intronic -3..-17 from the 5' end of
/// the downstream exon. Matches VEP CORE `splice_polypyrimidine_tract_variant`
/// (SO:0002169), which is emitted by VEP since release 105 (December 2021,
/// commit `0dba4466`) and folds into `splice_region_variant` in the
/// concordance harness normalization.
///
/// Reference: `ensembl-variation` `BaseTranscriptVariationAllele.pm`
/// `_intron_effects`: `overlap(start, end, intron_end - 16, intron_end - 2)`.
const SPLICE_REGION_INTRON_MAX_ACCEPTOR: u64 = 17;

/// Exonic splice region: 1-3 bases from the exon/intron boundary.
const SPLICE_REGION_EXON_MAX: u64 = 3;

/// Which side of a splice junction a position falls on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceSide {
    /// Donor side: 5' end of the intron (3' end of the upstream exon in
    /// transcript order). The consensus dinucleotide is GT.
    Donor,
    /// Acceptor side: 3' end of the intron (5' end of the downstream exon
    /// in transcript order). The consensus dinucleotide is AG.
    Acceptor,
}

/// Classification of where a genomic variant falls relative to a transcript.
///
/// All offsets and distances are in transcript-order (5'->3') terms. Genomic
/// coordinates are resolved to transcript-relative positions internally.
///
/// This enum is `Copy` -- no heap allocations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantLocation {
    /// Variant is upstream of the transcript (5' of tx_start in transcript
    /// terms), within 5,000 bases.
    Upstream {
        /// Distance in bases from the transcript 5' end.
        distance: u64,
    },

    /// Variant is downstream of the transcript (3' of tx_end in transcript
    /// terms), within 5,000 bases.
    Downstream {
        /// Distance in bases from the transcript 3' end.
        distance: u64,
    },

    /// Variant is beyond 5,000 bases from the transcript.
    Distal,

    /// Variant falls in the 5' UTR (between transcript start and CDS start,
    /// on an exon).
    FivePrimeUtr {
        /// 0-based exon index into [`TranscriptModel::exons`].
        exon_index: u16,
        /// Offset from the CDS start, as a negative number. Counts only
        /// exonic bases, matching HGVS c.-N notation.
        offset_from_cds_start: i64,
        /// `true` if within 3 bases of an intron-adjacent exon boundary.
        is_splice_region: bool,
    },

    /// Variant falls in a CDS exon (coding region).
    CdsExon {
        /// 0-based exon index into [`TranscriptModel::exons`].
        exon_index: u16,
        /// 0-based CDS segment index for O(1) phase lookup.
        cds_segment_index: u16,
        /// 0-based offset within the CDS (cumulative coding bases before
        /// this position). The first coding base is offset 0.
        cds_offset: u32,
        /// 1-based codon number (`cds_offset / 3 + 1`).
        codon_number: u32,
        /// Position within the codon: 0, 1, or 2.
        codon_position: u8,
        /// `true` if within 3 bases of an intron-adjacent exon boundary.
        is_splice_region: bool,
    },

    /// Variant falls in the 3' UTR (between CDS end and transcript end,
    /// on an exon).
    ThreePrimeUtr {
        /// 0-based exon index into [`TranscriptModel::exons`].
        exon_index: u16,
        /// Positive offset from the last CDS base, counting only exonic
        /// bases. Matches HGVS c.*N notation.
        offset_from_cds_end: i64,
        /// `true` if within 3 bases of an intron-adjacent exon boundary.
        is_splice_region: bool,
    },

    /// Variant at a canonical splice donor site (intronic +1 or +2).
    SpliceDonor {
        /// 0-based intron index.
        intron_index: u16,
        /// Distance from the exon boundary: 1 or 2.
        offset: u8,
    },

    /// Variant at a canonical splice acceptor site (intronic -1 or -2).
    SpliceAcceptor {
        /// 0-based intron index.
        intron_index: u16,
        /// Distance from the exon boundary: 1 or 2.
        offset: u8,
    },

    /// Variant in the splice region but outside canonical +/-1-2: intronic
    /// **+3..+8** on the donor side or **-3..-17** on the acceptor side.
    ///
    /// The acceptor window covers VEP's `splice_polypyrimidine_tract_variant`
    /// (SO:0002169), which has been emitted by VEP core since Ensembl release
    /// 105 (2021); vareffect folds it, along with `splice_donor_region_variant`
    /// (+3..+6) and `splice_donor_5th_base_variant` (+5), into a single
    /// `splice_region_variant` output to keep the `Consequence` enum compact.
    /// Fine-grained splice impact is scored separately via SpliceAI.
    SpliceRegion {
        /// 0-based intron index.
        intron_index: u16,
        /// Which side of the splice junction.
        side: SpliceSide,
        /// Signed distance from exon boundary. Positive for donor (+3..+8),
        /// negative for acceptor (-3..-17). Maps to HGVS intronic offset.
        distance: i64,
    },

    /// Deep intron: more than 8 bases from the donor side or more than 17
    /// bases from the acceptor side of the nearest exon boundary.
    Intron {
        /// 0-based intron index.
        intron_index: u16,
        /// Signed distance to nearest exon boundary. Positive = donor side,
        /// negative = acceptor side. Donor wins on tie (VEP convention).
        distance_to_nearest_exon: i64,
    },

    /// Exon of a non-coding transcript.
    NonCodingExon {
        /// 0-based exon index.
        exon_index: u16,
        /// `true` if within 3 bases of an intron-adjacent exon boundary.
        is_splice_region: bool,
    },

    /// Deep intron of a non-coding transcript.
    NonCodingIntron {
        /// 0-based intron index.
        intron_index: u16,
        /// Signed distance to nearest exon boundary.
        distance_to_nearest_exon: i64,
    },
}

/// Format a 0-based exon index as VEP's EXON field: `"N/total"` (1-based).
///
/// # Examples
///
/// ```
/// use vareffect::locate::format_exon_number;
/// assert_eq!(format_exon_number(4, 11), "5/11");
/// ```
pub fn format_exon_number(exon_index: u16, exon_count: u16) -> String {
    format!("{}/{}", exon_index + 1, exon_count)
}

/// Format a 0-based intron index as VEP's INTRON field: `"N/total"` (1-based).
///
/// Total introns = `exon_count - 1`.
///
/// # Examples
///
/// ```
/// use vareffect::locate::format_intron_number;
/// assert_eq!(format_intron_number(3, 11), "4/10");
/// ```
pub fn format_intron_number(intron_index: u16, exon_count: u16) -> String {
    debug_assert!(exon_count > 1, "no introns in a single-exon transcript");
    format!("{}/{}", intron_index + 1, exon_count.saturating_sub(1))
}

/// Precomputed lookup data for O(1) CDS offset and O(log n) exon lookup.
///
/// Built once per transcript via [`LocateIndex::build`], then reused across
/// all variant queries against that transcript.
///
/// # Memory
///
/// For TTN (~363 exons): ~2.5 KB. For a full MANE store (~6 000 transcripts):
/// ~15 MB total.
#[derive(Debug, Clone)]
pub struct LocateIndex {
    /// `cumulative_cds[i]` = sum of CDS segment lengths for segments `0..i`.
    /// Length = `cds_segments.len() + 1` (prefix-sum convention).
    cumulative_cds: Vec<u32>,
    /// Maps `exon_index` -> CDS segment index, or `None` if the exon has no CDS.
    exon_to_cds_seg: Vec<Option<u16>>,
    /// 0-based exon index containing the first CDS base (transcript 5' end).
    cds_start_exon_idx: Option<usize>,
    /// 0-based exon index containing the last CDS base (transcript 3' end).
    cds_end_exon_idx: Option<usize>,
    /// Number of exonic bases between the transcript 5' end and the CDS start.
    /// Precomputed at build time for O(1) cDNA position computation. 0 for
    /// non-coding transcripts.
    utr5_exonic_len: u32,
}

impl LocateIndex {
    /// Build the index for a given transcript. O(n) in the number of exons.
    ///
    /// # Errors
    ///
    /// Returns [`VarEffectError::Malformed`] if a CDS boundary coordinate does
    /// not fall within any exon.
    pub fn build(transcript: &TranscriptModel) -> Result<Self, VarEffectError> {
        let malformed =
            |msg: &str| VarEffectError::Malformed(format!("{}: {}", transcript.accession, msg));

        let mut cumulative_cds = Vec::with_capacity(transcript.cds_segments.len() + 1);
        cumulative_cds.push(0u32);
        for seg in &transcript.cds_segments {
            let prev = *cumulative_cds
                .last()
                .expect("cumulative_cds initialized with 0");
            cumulative_cds.push(prev + (seg.genomic_end - seg.genomic_start) as u32);
        }

        let mut exon_to_cds_seg = vec![None; transcript.exons.len()];
        for (seg_idx, seg) in transcript.cds_segments.iter().enumerate() {
            exon_to_cds_seg[seg.exon_index as usize] = Some(to_u16(seg_idx));
        }

        let (cds_start_exon_idx, cds_end_exon_idx) = if transcript.cds_segments.is_empty() {
            (None, None)
        } else {
            // When CDS segments exist, both cds_genomic_start and
            // cds_genomic_end are guaranteed to be Some. We use
            // Option::zip to bind both values without redundant unwraps.
            let start_idx = transcript
                .cds_genomic_start
                .zip(transcript.cds_genomic_end)
                .map(|(cds_start, cds_end)| {
                    let target = match transcript.strand {
                        Strand::Plus => cds_start,
                        Strand::Minus => cds_end - 1,
                    };
                    transcript
                        .exons
                        .iter()
                        .position(|e| e.genomic_start <= target && target < e.genomic_end)
                        .ok_or_else(|| malformed("CDS start base not found in any exon"))
                })
                .transpose()?;

            let end_idx = transcript
                .cds_genomic_start
                .zip(transcript.cds_genomic_end)
                .map(|(cds_start, cds_end)| {
                    let target = match transcript.strand {
                        Strand::Plus => cds_end - 1,
                        Strand::Minus => cds_start,
                    };
                    transcript
                        .exons
                        .iter()
                        .position(|e| e.genomic_start <= target && target < e.genomic_end)
                        .ok_or_else(|| malformed("CDS end base not found in any exon"))
                })
                .transpose()?;

            (start_idx, end_idx)
        };

        // Precompute 5'UTR exonic length for O(1) cDNA position computation.
        let utr5_exonic_len = compute_utr5_exonic_length(transcript);

        Ok(Self {
            cumulative_cds,
            exon_to_cds_seg,
            cds_start_exon_idx,
            cds_end_exon_idx,
            utr5_exonic_len,
        })
    }

    /// Prefix-sum array of CDS segment lengths.
    ///
    /// Used by `consequence` module to map a CDS offset back to a genomic
    /// position via binary search.
    pub(crate) fn cumulative_cds(&self) -> &[u32] {
        &self.cumulative_cds
    }

    /// Total number of coding bases across all CDS segments.
    ///
    /// Returns 0 for non-coding transcripts.
    pub(crate) fn total_cds_length(&self) -> u32 {
        *self.cumulative_cds.last().unwrap_or(&0)
    }

    /// Maps exon index -> CDS segment index, or `None` if the exon has no CDS.
    #[allow(dead_code)]
    pub(crate) fn exon_to_cds_seg(&self) -> &[Option<u16>] {
        &self.exon_to_cds_seg
    }

    /// 0-based exon index containing the first CDS base (transcript 5' end).
    #[allow(dead_code)]
    pub(crate) fn cds_start_exon_idx(&self) -> Option<usize> {
        self.cds_start_exon_idx
    }

    /// 0-based exon index containing the last CDS base (transcript 3' end).
    #[allow(dead_code)]
    pub(crate) fn cds_end_exon_idx(&self) -> Option<usize> {
        self.cds_end_exon_idx
    }

    /// Number of exonic bases between the transcript 5' end and the CDS start.
    /// Precomputed at build time. Returns 0 for non-coding transcripts.
    pub(crate) fn utr5_exonic_len(&self) -> u32 {
        self.utr5_exonic_len
    }
}

/// Compute the 5'UTR exonic length for a transcript: number of exonic bases
/// between the transcript 5' end and the CDS start.
fn compute_utr5_exonic_length(transcript: &TranscriptModel) -> u32 {
    let cds_start = match transcript.strand {
        Strand::Plus => match transcript.cds_genomic_start {
            Some(s) => s,
            None => return 0,
        },
        Strand::Minus => match transcript.cds_genomic_end {
            Some(e) => e,
            None => return 0,
        },
    };

    let mut len = 0u32;
    for exon in &transcript.exons {
        match transcript.strand {
            Strand::Plus => {
                if exon.genomic_end <= cds_start {
                    len += (exon.genomic_end - exon.genomic_start) as u32;
                } else if exon.genomic_start < cds_start {
                    len += (cds_start - exon.genomic_start) as u32;
                }
            }
            Strand::Minus => {
                if exon.genomic_start >= cds_start {
                    len += (exon.genomic_end - exon.genomic_start) as u32;
                } else if exon.genomic_end > cds_start {
                    len += (exon.genomic_end - cds_start) as u32;
                }
            }
        }
    }
    len
}

/// Locate a variant within a transcript's structure.
///
/// Given a genomic position (0-based, half-open) and a [`TranscriptModel`],
/// returns a [`VariantLocation`] describing which transcript region the
/// variant falls in.
///
/// `end` is accepted for API forward-compatibility but only `start` is used
/// for classification.
///
/// # Errors
///
/// Returns [`VarEffectError::Malformed`] if:
/// - `chrom` does not match `transcript.chrom`
/// - `start` falls within `[tx_start, tx_end)` but does not match any
///   exon or intron (malformed transcript model)
/// - `start` is classified as CDS but has no corresponding CDS segment
pub fn locate_variant(
    chrom: &str,
    start: u64,
    _end: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<VariantLocation, VarEffectError> {
    debug_assert!(
        _end > start,
        "end ({}) must be greater than start ({})",
        _end,
        start,
    );

    if chrom != transcript.chrom {
        return Err(VarEffectError::Malformed(format!(
            "chromosome mismatch: caller passed '{}' but transcript {} is on '{}'",
            chrom, transcript.accession, transcript.chrom,
        )));
    }

    let is_coding = !transcript.cds_segments.is_empty();

    if start < transcript.tx_start {
        let distance = transcript.tx_start - start;
        if distance > UPSTREAM_DOWNSTREAM_LIMIT {
            return Ok(VariantLocation::Distal);
        }
        return Ok(match transcript.strand {
            Strand::Plus => VariantLocation::Upstream { distance },
            Strand::Minus => VariantLocation::Downstream { distance },
        });
    }
    if start >= transcript.tx_end {
        let distance = start - transcript.tx_end + 1;
        if distance > UPSTREAM_DOWNSTREAM_LIMIT {
            return Ok(VariantLocation::Distal);
        }
        return Ok(match transcript.strand {
            Strand::Plus => VariantLocation::Downstream { distance },
            Strand::Minus => VariantLocation::Upstream { distance },
        });
    }

    match find_exon_or_intron(start, transcript)? {
        ExonOrIntron::Intron {
            upstream,
            downstream,
        } => Ok(classify_intronic(
            start,
            upstream,
            downstream,
            to_u16(upstream),
            transcript,
            is_coding,
        )),
        ExonOrIntron::Exon(exon_idx) => {
            let splice_region = is_exonic_splice_region(start, exon_idx, transcript);
            let exon_index = to_u16(exon_idx);

            if !is_coding {
                return Ok(VariantLocation::NonCodingExon {
                    exon_index,
                    is_splice_region: splice_region,
                });
            }

            let cds_start = transcript.cds_genomic_start.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: coding transcript has no cds_genomic_start",
                    transcript.accession,
                ))
            })?;
            let cds_end = transcript.cds_genomic_end.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: coding transcript has no cds_genomic_end",
                    transcript.accession,
                ))
            })?;

            let is_5utr = match transcript.strand {
                Strand::Plus => start < cds_start,
                Strand::Minus => start >= cds_end,
            };
            let is_3utr = match transcript.strand {
                Strand::Plus => start >= cds_end,
                Strand::Minus => start < cds_start,
            };

            if is_5utr {
                let offset = compute_utr_offset_5prime(start, exon_idx, transcript, index)?;
                Ok(VariantLocation::FivePrimeUtr {
                    exon_index,
                    offset_from_cds_start: offset,
                    is_splice_region: splice_region,
                })
            } else if is_3utr {
                let offset = compute_utr_offset_3prime(start, exon_idx, transcript, index)?;
                Ok(VariantLocation::ThreePrimeUtr {
                    exon_index,
                    offset_from_cds_end: offset,
                    is_splice_region: splice_region,
                })
            } else {
                let cds = compute_cds_offset(start, exon_idx, transcript, index)?;
                Ok(VariantLocation::CdsExon {
                    exon_index,
                    cds_segment_index: cds.cds_segment_index,
                    cds_offset: cds.cds_offset,
                    codon_number: cds.codon_number,
                    codon_position: cds.codon_position,
                    is_splice_region: splice_region,
                })
            }
        }
    }
}
