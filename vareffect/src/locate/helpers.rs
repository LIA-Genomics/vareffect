//! Shared internal helpers for variant location classification.

use super::{
    SPLICE_CANONICAL_MAX, SPLICE_REGION_EXON_MAX, SPLICE_REGION_INTRON_MAX_ACCEPTOR,
    SPLICE_REGION_INTRON_MAX_DONOR, SpliceSide, VariantLocation,
};
use crate::error::VarEffectError;
use crate::types::{Strand, TranscriptModel};

/// Safe `usize` -> `u16` conversion for exon/intron indices. TTN (~363 exons)
/// is the largest known human transcript, so `u16` is always sufficient.
pub(super) fn to_u16(val: usize) -> u16 {
    u16::try_from(val).expect("exon/intron index exceeds u16::MAX")
}

/// Result of searching for a position within a transcript's exon/intron
/// structure.
pub(super) enum ExonOrIntron {
    /// Position falls within the exon at this 0-based index.
    Exon(usize),
    /// Position falls in the intron between `upstream` and `downstream`
    /// exons (in transcript order, both 0-based indices).
    Intron { upstream: usize, downstream: usize },
}

/// Binary search to find which exon or intron a position falls in.
/// O(log n) in the number of exons.
///
/// Exons are ordered 5'->3' on the transcript: ascending `genomic_start` for
/// plus strand, descending for minus strand.
pub(super) fn find_exon_or_intron(
    start: u64,
    transcript: &TranscriptModel,
) -> Result<ExonOrIntron, VarEffectError> {
    let exons = &transcript.exons;

    match transcript.strand {
        Strand::Plus => {
            let p = exons.partition_point(|e| e.genomic_start <= start);
            let idx = p.checked_sub(1).ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "position {} below all exon starts in {} [{}, {})",
                    start, transcript.accession, transcript.tx_start, transcript.tx_end,
                ))
            })?;
            if start < exons[idx].genomic_end {
                Ok(ExonOrIntron::Exon(idx))
            } else if idx + 1 < exons.len() {
                Ok(ExonOrIntron::Intron {
                    upstream: idx,
                    downstream: idx + 1,
                })
            } else {
                Err(VarEffectError::Malformed(format!(
                    "position {} past last exon in {} [{}, {})",
                    start, transcript.accession, transcript.tx_start, transcript.tx_end,
                )))
            }
        }
        Strand::Minus => {
            let idx = exons.partition_point(|e| e.genomic_start > start);
            if idx >= exons.len() {
                return Err(VarEffectError::Malformed(format!(
                    "position {} below all exon starts in {} [{}, {})",
                    start, transcript.accession, transcript.tx_start, transcript.tx_end,
                )));
            }
            if start < exons[idx].genomic_end {
                Ok(ExonOrIntron::Exon(idx))
            } else if idx > 0 {
                Ok(ExonOrIntron::Intron {
                    upstream: idx - 1,
                    downstream: idx,
                })
            } else {
                Err(VarEffectError::Malformed(format!(
                    "position {} past first exon in {} [{}, {})",
                    start, transcript.accession, transcript.tx_start, transcript.tx_end,
                )))
            }
        }
    }
}

/// Check if an exonic position is within [`SPLICE_REGION_EXON_MAX`] bases of
/// an exon boundary that borders an intron.
///
/// First/last exon boundaries that face the transcript edge (no adjacent
/// intron) do NOT count -- VEP does not report splice_region_variant at
/// transcript boundaries.
pub(super) fn is_exonic_splice_region(
    start: u64,
    exon_index: usize,
    transcript: &TranscriptModel,
) -> bool {
    let exon = &transcript.exons[exon_index];
    let exon_count = transcript.exon_count as usize;

    match transcript.strand {
        Strand::Plus => {
            let donor =
                exon_index < exon_count - 1 && (exon.genomic_end - start) <= SPLICE_REGION_EXON_MAX;
            let acceptor =
                exon_index > 0 && (start - exon.genomic_start + 1) <= SPLICE_REGION_EXON_MAX;
            donor || acceptor
        }
        Strand::Minus => {
            let donor = exon_index < exon_count - 1
                && (start - exon.genomic_start + 1) <= SPLICE_REGION_EXON_MAX;
            let acceptor = exon_index > 0 && (exon.genomic_end - start) <= SPLICE_REGION_EXON_MAX;
            donor || acceptor
        }
    }
}

/// CDS offset result: segment index, offset, codon number, and codon position.
pub(crate) struct CdsOffsetResult {
    pub(crate) cds_segment_index: u16,
    pub(crate) cds_offset: u32,
    pub(crate) codon_number: u32,
    pub(crate) codon_position: u8,
}

/// Compute the CDS offset for a position known to be within the CDS.
///
/// Uses `LocateIndex.cumulative_cds` for O(1) accumulated-offset lookup and
/// `exon_to_cds_seg` to jump directly from the exon index to the CDS segment.
pub(crate) fn compute_cds_offset(
    start: u64,
    exon_idx: usize,
    transcript: &TranscriptModel,
    index: &super::LocateIndex,
) -> Result<CdsOffsetResult, VarEffectError> {
    let seg_idx = index
        .exon_to_cds_seg
        .get(exon_idx)
        .copied()
        .flatten()
        .ok_or_else(|| {
            VarEffectError::Malformed(format!(
                "{}: exon {} has no CDS segment but position {} classified as CDS",
                transcript.accession, exon_idx, start,
            ))
        })? as usize;

    let segment = &transcript.cds_segments[seg_idx];
    let intra = match transcript.strand {
        Strand::Plus => (start - segment.genomic_start) as u32,
        Strand::Minus => (segment.genomic_end - 1 - start) as u32,
    };
    let cds_offset = index.cumulative_cds[seg_idx] + intra;
    Ok(CdsOffsetResult {
        cds_segment_index: to_u16(seg_idx),
        cds_offset,
        codon_number: cds_offset / 3 + 1,
        codon_position: (cds_offset % 3) as u8,
    })
}

/// Compute the 5'UTR offset (negative) counting only exonic bases from the
/// variant position to the CDS start.
///
/// The result maps to HGVS c.-N notation: c.-1 is the base immediately
/// upstream of the start codon, c.-2 is two bases upstream, etc.
pub(crate) fn compute_utr_offset_5prime(
    start: u64,
    exon_index: usize,
    transcript: &TranscriptModel,
    index: &super::LocateIndex,
) -> Result<i64, VarEffectError> {
    let exons = &transcript.exons;
    let malformed =
        |msg: &str| VarEffectError::Malformed(format!("{}: {}", transcript.accession, msg));

    match transcript.strand {
        Strand::Plus => {
            let cds_start = transcript
                .cds_genomic_start
                .ok_or_else(|| malformed("5'UTR offset requested but cds_genomic_start is None"))?;
            let cds_exon_idx = index
                .cds_start_exon_idx
                .ok_or_else(|| malformed("5'UTR offset but no cds_start_exon_idx"))?;

            if exon_index == cds_exon_idx {
                Ok(-((cds_start - start) as i64))
            } else {
                let mut acc: u64 = exons[exon_index].genomic_end - start;
                for exon in &exons[(exon_index + 1)..cds_exon_idx] {
                    acc += exon.genomic_end - exon.genomic_start;
                }
                acc += cds_start - exons[cds_exon_idx].genomic_start;
                Ok(-(acc as i64))
            }
        }
        Strand::Minus => {
            // CDS "start" in transcript terms is at cds_genomic_end.
            let cds_boundary = transcript
                .cds_genomic_end
                .ok_or_else(|| malformed("5'UTR offset requested but cds_genomic_end is None"))?;
            let cds_exon_idx = index
                .cds_start_exon_idx
                .ok_or_else(|| malformed("5'UTR offset but no cds_start_exon_idx"))?;

            if exon_index == cds_exon_idx {
                Ok(-((start - cds_boundary + 1) as i64))
            } else {
                let mut acc: u64 = start - exons[exon_index].genomic_start + 1;
                for exon in &exons[(exon_index + 1)..cds_exon_idx] {
                    acc += exon.genomic_end - exon.genomic_start;
                }
                acc += exons[cds_exon_idx].genomic_end - cds_boundary;
                Ok(-(acc as i64))
            }
        }
    }
}

/// Compute the 3'UTR offset (positive) counting only exonic bases from the
/// CDS end boundary to the variant position.
///
/// The result is 1-based: the first base after the stop codon is offset 1
/// (HGVS c.*1).
pub(crate) fn compute_utr_offset_3prime(
    start: u64,
    exon_index: usize,
    transcript: &TranscriptModel,
    index: &super::LocateIndex,
) -> Result<i64, VarEffectError> {
    let exons = &transcript.exons;
    let malformed =
        |msg: &str| VarEffectError::Malformed(format!("{}: {}", transcript.accession, msg));

    match transcript.strand {
        Strand::Plus => {
            let cds_end = transcript
                .cds_genomic_end
                .ok_or_else(|| malformed("3'UTR offset requested but cds_genomic_end is None"))?;
            let cds_end_exon_idx = index
                .cds_end_exon_idx
                .ok_or_else(|| malformed("3'UTR offset but no cds_end_exon_idx"))?;

            if exon_index == cds_end_exon_idx {
                Ok((start - cds_end + 1) as i64)
            } else {
                let mut acc: u64 = exons[cds_end_exon_idx].genomic_end - cds_end;
                for exon in &exons[(cds_end_exon_idx + 1)..exon_index] {
                    acc += exon.genomic_end - exon.genomic_start;
                }
                acc += start - exons[exon_index].genomic_start + 1;
                Ok(acc as i64)
            }
        }
        Strand::Minus => {
            // CDS end in transcript terms = cds_genomic_start.
            let cds_boundary = transcript
                .cds_genomic_start
                .ok_or_else(|| malformed("3'UTR offset requested but cds_genomic_start is None"))?;
            let cds_end_exon_idx = index
                .cds_end_exon_idx
                .ok_or_else(|| malformed("3'UTR offset but no cds_end_exon_idx"))?;

            if exon_index == cds_end_exon_idx {
                Ok((cds_boundary - start) as i64)
            } else {
                let mut acc: u64 = cds_boundary - exons[cds_end_exon_idx].genomic_start;
                for exon in &exons[(cds_end_exon_idx + 1)..exon_index] {
                    acc += exon.genomic_end - exon.genomic_start;
                }
                acc += exons[exon_index].genomic_end - start;
                Ok(acc as i64)
            }
        }
    }
}

/// Classify an intronic position as splice donor, acceptor, region, or deep
/// intron.
///
/// `upstream` and `downstream` are 0-based exon indices in transcript order.
/// `intron_index` is the 0-based intron index (= `upstream`).
pub(super) fn classify_intronic(
    start: u64,
    upstream: usize,
    downstream: usize,
    intron_index: u16,
    transcript: &TranscriptModel,
    is_coding: bool,
) -> VariantLocation {
    let up_exon = &transcript.exons[upstream];
    let down_exon = &transcript.exons[downstream];

    let (donor_dist, acceptor_dist) = match transcript.strand {
        Strand::Plus => {
            let d = start - up_exon.genomic_end + 1;
            let a = down_exon.genomic_start - start;
            (d, a)
        }
        Strand::Minus => {
            let d = up_exon.genomic_start - start;
            let a = start - down_exon.genomic_end + 1;
            (d, a)
        }
    };

    if donor_dist <= SPLICE_CANONICAL_MAX {
        return VariantLocation::SpliceDonor {
            intron_index,
            offset: donor_dist as u8,
        };
    }
    if acceptor_dist <= SPLICE_CANONICAL_MAX {
        return VariantLocation::SpliceAcceptor {
            intron_index,
            offset: acceptor_dist as u8,
        };
    }

    // Splice region window is asymmetric: donor side +3..+8 (SO:0001630),
    // acceptor side -3..-17 (SO:0002169, polypyrimidine tract).
    let donor_in_region = donor_dist <= SPLICE_REGION_INTRON_MAX_DONOR;
    let acceptor_in_region = acceptor_dist <= SPLICE_REGION_INTRON_MAX_ACCEPTOR;

    if donor_in_region && !acceptor_in_region {
        return VariantLocation::SpliceRegion {
            intron_index,
            side: SpliceSide::Donor,
            distance: donor_dist as i64,
        };
    }
    if acceptor_in_region && !donor_in_region {
        return VariantLocation::SpliceRegion {
            intron_index,
            side: SpliceSide::Acceptor,
            distance: -(acceptor_dist as i64),
        };
    }
    if donor_in_region && acceptor_in_region {
        // Both qualify -- prefer donor on tie (matches VEP).
        if donor_dist <= acceptor_dist {
            return VariantLocation::SpliceRegion {
                intron_index,
                side: SpliceSide::Donor,
                distance: donor_dist as i64,
            };
        }
        return VariantLocation::SpliceRegion {
            intron_index,
            side: SpliceSide::Acceptor,
            distance: -(acceptor_dist as i64),
        };
    }

    // Deep intron. Positive = donor side, negative = acceptor side.
    // Donor wins on tie.
    let signed_distance = if donor_dist <= acceptor_dist {
        donor_dist as i64
    } else {
        -(acceptor_dist as i64)
    };

    if is_coding {
        VariantLocation::Intron {
            intron_index,
            distance_to_nearest_exon: signed_distance,
        }
    } else {
        VariantLocation::NonCodingIntron {
            intron_index,
            distance_to_nearest_exon: signed_distance,
        }
    }
}
