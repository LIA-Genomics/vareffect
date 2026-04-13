//! Indel footprint location classification.

use super::helpers::{
    ExonOrIntron, compute_cds_offset, find_exon_or_intron, is_exonic_splice_region, to_u16,
};
use super::{
    LocateIndex, SPLICE_CANONICAL_MAX, SPLICE_REGION_INTRON_MAX_ACCEPTOR,
    SPLICE_REGION_INTRON_MAX_DONOR,
};
use crate::error::VarEffectError;
use crate::types::{Strand, TranscriptModel};

/// Classification of an indel's genomic footprint relative to a transcript.
///
/// Unlike point-based [`super::VariantLocation`] (for SNVs), this captures the
/// full range of affected bases, including whether the footprint overlaps
/// splice sites or crosses exon/intron boundaries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndelLocation {
    /// Primary region classification.
    pub region: IndelRegion,
    /// Whether any base in the footprint overlaps a canonical splice site
    /// (donor +1/+2 or acceptor -1/-2).
    pub overlaps_splice_canonical: bool,
    /// Whether any base in the footprint overlaps the splice region
    /// (+3..+8 / -3..-8 intronic, or 1-3 exonic from boundary).
    pub overlaps_splice_region: bool,
    /// Whether the indel crosses an exon/intron boundary.
    pub crosses_exon_boundary: bool,
    /// Exon index if the footprint starts in an exon (for EXON display field).
    pub exon_index: Option<u16>,
    /// Intron index if the footprint starts in an intron (for INTRON display
    /// field).
    pub intron_index: Option<u16>,
    /// Detailed donor/acceptor splice overlap, populated when
    /// `overlaps_splice_canonical` is `true`.
    pub splice_detail: Option<SpliceOverlapDetail>,
}

/// Detailed splice overlap for an indel footprint.
///
/// Identifies which specific splice sites (donor vs acceptor) are
/// overlapped and which introns are involved.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SpliceOverlapDetail {
    /// Whether any base overlaps a canonical splice donor site (+1/+2).
    pub overlaps_donor: bool,
    /// Whether any base overlaps a canonical splice acceptor site (-1/-2).
    pub overlaps_acceptor: bool,
    /// Whether any base overlaps the splice region (+3..+8 / -3..-8).
    pub overlaps_splice_region: bool,
    /// Intron indices (0-based, transcript order) whose donor sites are
    /// overlapped.
    pub donor_intron_indices: Vec<u16>,
    /// Intron indices (0-based, transcript order) whose acceptor sites are
    /// overlapped.
    pub acceptor_intron_indices: Vec<u16>,
}

/// Primary region for an indel footprint.
///
/// `Cds` uses half-open convention: `[cds_offset_start, cds_offset_end)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IndelRegion {
    /// Entirely within a CDS exon. Offsets are 0-based, half-open.
    Cds {
        /// 0-based CDS offset of the first affected coding base.
        cds_offset_start: u32,
        /// 0-based exclusive end (first unaffected offset).
        cds_offset_end: u32,
    },
    /// Entirely within the 5' UTR.
    FivePrimeUtr,
    /// Entirely within the 3' UTR.
    ThreePrimeUtr,
    /// Entirely within an intron.
    Intron,
    /// Entirely within a non-coding exon.
    NonCodingExon,
    /// Upstream of the transcript.
    Upstream,
    /// Downstream of the transcript.
    Downstream,
    /// Spans a region boundary: exon/intron, CDS/UTR, or multi-exon.
    BoundarySpanning,
}

/// Detailed splice site overlap analysis for an indel footprint.
///
/// Identifies which specific splice sites (donor vs acceptor) are overlapped
/// and records the intron indices involved.
fn check_splice_overlap_detailed(
    range_start: u64,
    range_end: u64,
    transcript: &TranscriptModel,
) -> SpliceOverlapDetail {
    let exons = &transcript.exons;
    let n = exons.len();
    let mut detail = SpliceOverlapDetail::default();
    if n < 2 {
        return detail;
    }

    // Narrow the scan using binary search for plus-strand transcripts.
    // Use the larger of the two window maxes (acceptor polypyrimidine tract)
    // so the filter cannot exclude acceptor-only candidates.
    let expanded_lo = range_start.saturating_sub(SPLICE_REGION_INTRON_MAX_ACCEPTOR);
    let expanded_hi = range_end + SPLICE_REGION_INTRON_MAX_ACCEPTOR;
    let first = match transcript.strand {
        Strand::Plus => {
            let p = exons.partition_point(|e| e.genomic_end <= expanded_lo);
            if p > 0 { p - 1 } else { 0 }
        }
        Strand::Minus => 0,
    };

    for i in first..n - 1 {
        let upstream_exon = &exons[i];
        let downstream_exon = &exons[i + 1];

        // Early exit for plus strand: if both exons are past the window.
        if matches!(transcript.strand, Strand::Plus)
            && upstream_exon.genomic_start > expanded_hi
            && downstream_exon.genomic_start > expanded_hi
        {
            break;
        }
        let intron_idx = to_u16(i);

        // Donor site: first 2 intronic bases adjacent to the 3' end of
        // the upstream exon (in transcript orientation).
        let (donor_start, donor_end) = match transcript.strand {
            Strand::Plus => (
                upstream_exon.genomic_end,
                upstream_exon.genomic_end + SPLICE_CANONICAL_MAX,
            ),
            Strand::Minus => (
                upstream_exon.genomic_start - SPLICE_CANONICAL_MAX,
                upstream_exon.genomic_start,
            ),
        };

        // Acceptor site: last 2 intronic bases adjacent to the 5' end of
        // the downstream exon (in transcript orientation).
        let (acceptor_start, acceptor_end) = match transcript.strand {
            Strand::Plus => (
                downstream_exon.genomic_start - SPLICE_CANONICAL_MAX,
                downstream_exon.genomic_start,
            ),
            Strand::Minus => (
                downstream_exon.genomic_end,
                downstream_exon.genomic_end + SPLICE_CANONICAL_MAX,
            ),
        };

        if range_start < donor_end && range_end > donor_start {
            detail.overlaps_donor = true;
            detail.donor_intron_indices.push(intron_idx);
        }

        if range_start < acceptor_end && range_end > acceptor_start {
            detail.overlaps_acceptor = true;
            detail.acceptor_intron_indices.push(intron_idx);
        }

        // Splice region (intronic only): +3..+8 on the donor side,
        // -3..-17 on the acceptor side. The acceptor window is extended to
        // cover VEP's `splice_polypyrimidine_tract_variant` (SO:0002169).
        let (donor_region_start, donor_region_end) = match transcript.strand {
            Strand::Plus => (
                upstream_exon.genomic_end + SPLICE_CANONICAL_MAX,
                upstream_exon.genomic_end + SPLICE_REGION_INTRON_MAX_DONOR,
            ),
            Strand::Minus => (
                upstream_exon
                    .genomic_start
                    .saturating_sub(SPLICE_REGION_INTRON_MAX_DONOR),
                upstream_exon
                    .genomic_start
                    .saturating_sub(SPLICE_CANONICAL_MAX),
            ),
        };
        let (acceptor_region_start, acceptor_region_end) = match transcript.strand {
            Strand::Plus => (
                downstream_exon
                    .genomic_start
                    .saturating_sub(SPLICE_REGION_INTRON_MAX_ACCEPTOR),
                downstream_exon
                    .genomic_start
                    .saturating_sub(SPLICE_CANONICAL_MAX),
            ),
            Strand::Minus => (
                downstream_exon.genomic_end + SPLICE_CANONICAL_MAX,
                downstream_exon.genomic_end + SPLICE_REGION_INTRON_MAX_ACCEPTOR,
            ),
        };

        if (range_start < donor_region_end && range_end > donor_region_start)
            || (range_start < acceptor_region_end && range_end > acceptor_region_start)
        {
            detail.overlaps_splice_region = true;
        }
    }

    detail
}

/// Classify an indel footprint against a transcript.
///
/// For **deletions**: `start` and `end` define the deleted range
/// `[start, end)` in 0-based half-open coordinates.
///
/// For **insertions**: `start == end` (zero-width). The insertion point is
/// between positions `start - 1` and `start`.
///
/// # Errors
///
/// Returns [`VarEffectError::Malformed`] on chromosome mismatch or
/// malformed transcript data.
pub fn locate_indel(
    chrom: &str,
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<IndelLocation, VarEffectError> {
    debug_assert!(end >= start, "end ({end}) must be >= start ({start})");

    if chrom != transcript.chrom {
        return Err(VarEffectError::Malformed(format!(
            "chromosome mismatch: caller passed '{}' but transcript {} is on '{}'",
            chrom, transcript.accession, transcript.chrom,
        )));
    }

    let is_coding = !transcript.cds_segments.is_empty();
    let is_insertion = start == end;

    // Upstream / downstream
    if end <= transcript.tx_start {
        let distance = transcript.tx_start - start;
        let region = if distance > super::UPSTREAM_DOWNSTREAM_LIMIT {
            return Ok(simple_indel_location(IndelRegion::Upstream));
        } else {
            match transcript.strand {
                Strand::Plus => IndelRegion::Upstream,
                Strand::Minus => IndelRegion::Downstream,
            }
        };
        return Ok(simple_indel_location(region));
    }
    if start >= transcript.tx_end {
        let distance = start - transcript.tx_end + 1;
        let region = if distance > super::UPSTREAM_DOWNSTREAM_LIMIT {
            return Ok(simple_indel_location(IndelRegion::Downstream));
        } else {
            match transcript.strand {
                Strand::Plus => IndelRegion::Downstream,
                Strand::Minus => IndelRegion::Upstream,
            }
        };
        return Ok(simple_indel_location(region));
    }

    // Splice overlap for the full footprint. A single `check_splice_overlap_detailed`
    // call produces donor/acceptor/region overlap booleans using half-open
    // interval intersection (the only correct way to match VEP's
    // `_intron_effects` — a previous `check_splice_overlap_for_range`
    // helper used a symmetric ±2 buffer around the footprint which
    // over-claimed canonical at intronic +3 and the last exonic bases).
    let (detail_start, detail_end) = if is_insertion {
        (start, start + 1)
    } else {
        (start, end)
    };
    let detail = check_splice_overlap_detailed(detail_start, detail_end, transcript);
    let splice_canonical = detail.overlaps_donor || detail.overlaps_acceptor;
    let splice_region = detail.overlaps_splice_region;
    let splice_detail = if splice_canonical { Some(detail) } else { None };

    // Insertions: classify the single insertion point
    if is_insertion {
        return locate_indel_insertion(
            start,
            transcript,
            index,
            is_coding,
            splice_canonical,
            splice_region,
            splice_detail,
        );
    }

    // Deletion/delins footprint straddles a transcript genomic edge
    // (`start < tx_start` or `end > tx_end`). `find_exon_or_intron`
    // rejects out-of-range coordinates with `Err(Malformed)`, so classify
    // the variant as `BoundarySpanning` here, anchoring the exon/intron
    // lookup at `start.max(tx_start)`. Downstream,
    // `annotate_boundary_spanning_deletion` uses min/max-clamped CDS
    // overlap math, so bases outside the transcript contribute nothing.
    if start < transcript.tx_start || end > transcript.tx_end {
        let anchor = start.max(transcript.tx_start);
        let anchor_loc = find_exon_or_intron(anchor, transcript)?;
        let (exon_index, intron_index) = match anchor_loc {
            ExonOrIntron::Exon(i) => (Some(to_u16(i)), None),
            ExonOrIntron::Intron { upstream, .. } => (None, Some(to_u16(upstream))),
        };
        return Ok(IndelLocation {
            region: IndelRegion::BoundarySpanning,
            overlaps_splice_canonical: splice_canonical,
            overlaps_splice_region: splice_region,
            crosses_exon_boundary: true,
            exon_index,
            intron_index,
            splice_detail,
        });
    }

    // Deletion: classify both endpoints
    let start_loc = find_exon_or_intron(start, transcript)?;
    let end_loc = find_exon_or_intron(end - 1, transcript)?;

    match (&start_loc, &end_loc) {
        // Both in the same intron
        (
            ExonOrIntron::Intron {
                upstream: u1,
                downstream: d1,
            },
            ExonOrIntron::Intron {
                upstream: u2,
                downstream: d2,
            },
        ) if u1 == u2 && d1 == d2 => {
            return Ok(IndelLocation {
                region: IndelRegion::Intron,
                overlaps_splice_canonical: splice_canonical,
                overlaps_splice_region: splice_region,
                crosses_exon_boundary: false,
                exon_index: None,
                intron_index: Some(to_u16(*u1)),
                splice_detail,
            });
        }

        // Both in the same exon -- fall through to region classification
        (ExonOrIntron::Exon(s), ExonOrIntron::Exon(e)) if s == e => {}

        // Different exons or exon/intron mix -- boundary crossing
        _ => {
            let (exon_idx, intron_idx) = match &start_loc {
                ExonOrIntron::Exon(i) => (Some(to_u16(*i)), None),
                ExonOrIntron::Intron { upstream, .. } => (None, Some(to_u16(*upstream))),
            };
            return Ok(IndelLocation {
                region: IndelRegion::BoundarySpanning,
                overlaps_splice_canonical: splice_canonical,
                overlaps_splice_region: splice_region,
                crosses_exon_boundary: true,
                exon_index: exon_idx,
                intron_index: intron_idx,
                splice_detail,
            });
        }
    }

    // Both endpoints in the same exon
    let start_exon = match start_loc {
        ExonOrIntron::Exon(i) => i,
        _ => unreachable!("both endpoints in same exon after match"),
    };

    let exon_index = to_u16(start_exon);
    let splice_region_exonic = is_exonic_splice_region(start, start_exon, transcript)
        || is_exonic_splice_region(end - 1, start_exon, transcript);

    if !is_coding {
        return Ok(IndelLocation {
            region: IndelRegion::NonCodingExon,
            overlaps_splice_canonical: splice_canonical,
            overlaps_splice_region: splice_region || splice_region_exonic,
            crosses_exon_boundary: false,
            exon_index: Some(exon_index),
            intron_index: None,
            splice_detail: splice_detail.clone(),
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

    let start_in_cds = start >= cds_start && start < cds_end;
    let end_in_cds = (end - 1) >= cds_start && (end - 1) < cds_end;

    if start_in_cds && end_in_cds {
        // For minus-strand, lower genomic -> higher CDS offset, so use min/max
        let cds_s = compute_cds_offset(start, start_exon, transcript, index)?;
        let cds_e = compute_cds_offset(end - 1, start_exon, transcript, index)?;
        let offset_lo = cds_s.cds_offset.min(cds_e.cds_offset);
        let offset_hi = cds_s.cds_offset.max(cds_e.cds_offset);
        return Ok(IndelLocation {
            region: IndelRegion::Cds {
                cds_offset_start: offset_lo,
                cds_offset_end: offset_hi + 1,
            },
            overlaps_splice_canonical: splice_canonical,
            overlaps_splice_region: splice_region || splice_region_exonic,
            crosses_exon_boundary: false,
            exon_index: Some(exon_index),
            intron_index: None,
            splice_detail: splice_detail.clone(),
        });
    }

    if !start_in_cds && !end_in_cds {
        let is_5utr = match transcript.strand {
            Strand::Plus => start < cds_start,
            Strand::Minus => start >= cds_end,
        };
        let region = if is_5utr {
            IndelRegion::FivePrimeUtr
        } else {
            IndelRegion::ThreePrimeUtr
        };
        return Ok(IndelLocation {
            region,
            overlaps_splice_canonical: splice_canonical,
            overlaps_splice_region: splice_region || splice_region_exonic,
            crosses_exon_boundary: false,
            exon_index: Some(exon_index),
            intron_index: None,
            splice_detail,
        });
    }

    // One endpoint in CDS, other in UTR -- CDS/UTR boundary
    Ok(IndelLocation {
        region: IndelRegion::BoundarySpanning,
        overlaps_splice_canonical: splice_canonical,
        overlaps_splice_region: splice_region || splice_region_exonic,
        crosses_exon_boundary: false,
        exon_index: Some(exon_index),
        intron_index: None,
        splice_detail,
    })
}

/// Classify an insertion point within a transcript.
fn locate_indel_insertion(
    pos: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    is_coding: bool,
    splice_canonical: bool,
    splice_region: bool,
    splice_detail: Option<SpliceOverlapDetail>,
) -> Result<IndelLocation, VarEffectError> {
    let loc = find_exon_or_intron(pos, transcript)?;

    match loc {
        ExonOrIntron::Intron { upstream, .. } => Ok(IndelLocation {
            region: IndelRegion::Intron,
            overlaps_splice_canonical: splice_canonical,
            overlaps_splice_region: splice_region,
            crosses_exon_boundary: false,
            exon_index: None,
            intron_index: Some(to_u16(upstream)),
            splice_detail,
        }),
        ExonOrIntron::Exon(exon_idx) => {
            let exon_index = to_u16(exon_idx);
            let splice_region_exonic = is_exonic_splice_region(pos, exon_idx, transcript);

            if !is_coding {
                return Ok(IndelLocation {
                    region: IndelRegion::NonCodingExon,
                    overlaps_splice_canonical: splice_canonical,
                    overlaps_splice_region: splice_region || splice_region_exonic,
                    crosses_exon_boundary: false,
                    exon_index: Some(exon_index),
                    intron_index: None,
                    splice_detail,
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

            let in_cds = pos >= cds_start && pos < cds_end;

            if in_cds {
                let cds = compute_cds_offset(pos, exon_idx, transcript, index)?;
                Ok(IndelLocation {
                    region: IndelRegion::Cds {
                        cds_offset_start: cds.cds_offset,
                        cds_offset_end: cds.cds_offset,
                    },
                    overlaps_splice_canonical: splice_canonical,
                    overlaps_splice_region: splice_region || splice_region_exonic,
                    crosses_exon_boundary: false,
                    exon_index: Some(exon_index),
                    intron_index: None,
                    splice_detail,
                })
            } else {
                let is_5utr = match transcript.strand {
                    Strand::Plus => pos < cds_start,
                    Strand::Minus => pos >= cds_end,
                };
                let region = if is_5utr {
                    IndelRegion::FivePrimeUtr
                } else {
                    IndelRegion::ThreePrimeUtr
                };
                Ok(IndelLocation {
                    region,
                    overlaps_splice_canonical: splice_canonical,
                    overlaps_splice_region: splice_region || splice_region_exonic,
                    crosses_exon_boundary: false,
                    exon_index: Some(exon_index),
                    intron_index: None,
                    splice_detail,
                })
            }
        }
    }
}

/// Build a simple `IndelLocation` with no splice overlap or boundary crossing.
fn simple_indel_location(region: IndelRegion) -> IndelLocation {
    IndelLocation {
        region,
        overlaps_splice_canonical: false,
        overlaps_splice_region: false,
        crosses_exon_boundary: false,
        exon_index: None,
        intron_index: None,
        splice_detail: None,
    }
}
