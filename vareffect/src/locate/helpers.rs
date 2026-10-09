//! Shared internal helpers for variant location classification.

use super::{
    FRAMESHIFT_INTRON_MAX_GAP, SPLICE_CANONICAL_MAX, SPLICE_REGION_EXON_MAX,
    SPLICE_REGION_INTRON_MAX_ACCEPTOR, SPLICE_REGION_INTRON_MAX_DONOR, SpliceSide, VariantLocation,
};
use crate::error::VarEffectError;
use crate::types::{Exon, Strand, TranscriptModel};

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

/// Number of genomic bases between two exons that are adjacent in
/// transcript order, on either strand. 0 for abutting exons.
pub(crate) fn intron_gap(a: &Exon, b: &Exon) -> u64 {
    if a.genomic_end <= b.genomic_start {
        b.genomic_start - a.genomic_end
    } else {
        a.genomic_start.saturating_sub(b.genomic_end)
    }
}

/// Genomic span `[start, end)` of the intron between two exons adjacent in
/// transcript order, when it is a frameshift intron (1 to
/// [`FRAMESHIFT_INTRON_MAX_GAP`] bases).
pub(crate) fn frameshift_gap(a: &Exon, b: &Exon) -> Option<(u64, u64)> {
    let gap = intron_gap(a, b);
    if gap == 0 || gap > FRAMESHIFT_INTRON_MAX_GAP {
        return None;
    }
    Some(if a.genomic_end <= b.genomic_start {
        (a.genomic_end, b.genomic_start)
    } else {
        (b.genomic_end, a.genomic_start)
    })
}

/// 0-based index of the frameshift intron that holds an end of the footprint
/// `[start, end)`, if any.
///
/// VEP treats only such a variant as `within_frameshift_intron`: its
/// `frameshift` and inframe predicates need both ends on coding bases, so a
/// deletion from an exon into the gap is `coding_sequence_variant`, while one
/// spanning the gap from exon to exon is a normal frameshift (VEP 116,
/// ENST00000400999.7 `c.228_228+1del` vs `c.226_229del`). An insertion
/// (`start == end`, between `start - 1` and `start`) counts only when both
/// flanking bases are intron bases (`BaseTranscriptVariationAllele.pm`
/// `_intron_effects`, insertion `overlap` on `r_start == r_end + 1`).
pub(crate) fn frameshift_intron_endpoint(
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
    index: &super::LocateIndex,
) -> Option<u16> {
    if !index.has_frameshift_intron() {
        return None;
    }
    let inside = |pos: u64, (gap_start, gap_end): (u64, u64)| gap_start <= pos && pos < gap_end;
    transcript
        .exons
        .windows(2)
        .position(|pair| {
            frameshift_gap(&pair[0], &pair[1]).is_some_and(|gap| {
                if start == end {
                    start > 0 && inside(start - 1, gap) && inside(start, gap)
                } else {
                    inside(start, gap) || inside(end - 1, gap)
                }
            })
        })
        .map(to_u16)
}

/// Canonical splice terms that a frameshift intron's windows give a variant
/// on the exon bases beside it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct SpliceSpill {
    /// The two-base donor window (the first intron bases in transcript
    /// order) holds the variant.
    pub(crate) donor: bool,
    /// The two-base acceptor window (the last intron bases) holds the
    /// variant.
    pub(crate) acceptor: bool,
}

/// Splice terms VEP gives a variant outside a frameshift intron whose
/// two-base windows spill past a one-base gap into the flanking exons.
///
/// VEP computes the donor and acceptor windows from the intron bounds as for
/// any intron, so on a one-base gap they cover the coding bases on either
/// side, which then get `splice_acceptor_variant` (base before the gap) or
/// `splice_donor_variant` (base after it) next to their coding terms (VEP
/// 116, ENST00000400999.7 c.228C>A, c.229G>A, c.228del, c.229del). A
/// deletion counts when it overlaps a window without touching the gap;
/// insertions use [`insertion_intron_effects`]. `None` when nothing spills,
/// which is always the case for a transcript without frameshift introns.
pub(crate) fn exonic_splice_spill(
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
    index: &super::LocateIndex,
) -> Option<SpliceSpill> {
    if !index.has_frameshift_intron() {
        return None;
    }
    let mut spill = SpliceSpill::default();
    for pair in transcript.exons.windows(2) {
        let Some((gap_start, gap_end)) = frameshift_gap(&pair[0], &pair[1]) else {
            continue;
        };
        let (donor, acceptor) = match transcript.strand {
            Strand::Plus => (
                (gap_start, gap_start + SPLICE_CANONICAL_MAX),
                (gap_end.saturating_sub(SPLICE_CANONICAL_MAX), gap_end),
            ),
            Strand::Minus => (
                (gap_end.saturating_sub(SPLICE_CANONICAL_MAX), gap_end),
                (gap_start, gap_start + SPLICE_CANONICAL_MAX),
            ),
        };
        let hits = |(window_start, window_end): (u64, u64)| {
            let touches_gap = start < gap_end && gap_start < end;
            !touches_gap && start < window_end && window_start < end
        };
        spill.donor |= hits(donor);
        spill.acceptor |= hits(acceptor);
    }
    (spill.donor || spill.acceptor).then_some(spill)
}

/// Splice and intron effects of an insertion, as VEP 116 computes them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct InsertionIntronEffects {
    /// Intron (0-based, transcript order) whose canonical donor
    /// dinucleotide holds both flanking bases.
    pub(crate) donor: Option<u16>,
    /// Intron whose canonical acceptor dinucleotide holds both flanking
    /// bases.
    pub(crate) acceptor: Option<u16>,
    /// `splice_region_variant`, VEP's donor-region and polypyrimidine terms
    /// folded in.
    pub(crate) splice_region: bool,
    /// Intron for which VEP's `intronic` test holds (`intron_variant`).
    pub(crate) intronic: Option<u16>,
}

impl InsertionIntronEffects {
    /// `true` when a canonical donor or acceptor site holds the insertion.
    pub(crate) fn canonical(&self) -> bool {
        self.donor.is_some() || self.acceptor.is_some()
    }
}

/// VEP 116's intron effects for an insertion between genomic bases
/// `pos - 1` and `pos` (0-based).
///
/// A transcription of `BaseTranscriptVariationAllele.pm::_intron_effects`
/// and `VariationEffect.pm::_intron_overlap` (ensembl-variation
/// release/116) for an insertion, where VEP's `r_end` is the left genomic
/// flank `l` and `r_start` the right one `r`. VEP's `overlap` on that
/// reversed pair holds only when *both* flanks lie in the window, except for
/// the first polypyrimidine test, which swaps the bounds and holds when
/// *either* flank does. With `lo`, `hi` an intron's first and last base:
///
/// - an intron is examined for the `intronic` and polypyrimidine tests when
///   a flank reaches `[lo - 3, hi + 3]`, and for the boundary tests when one
///   reaches `[lo - 3, lo + 7]` or `[hi - 7, hi + 3]`
///   (`BaseTranscriptVariation.pm::_create_intron_trees`); a frameshift
///   intron holding both flanks is skipped (`within_frameshift_intron`);
/// - canonical start site `[lo, lo + 1]`, end site `[hi - 1, hi]` (donor and
///   acceptor by strand);
/// - intronic: both flanks in `[lo + 2, hi - 2]`, or `r == lo + 2`, or
///   `l == hi - 2`;
/// - splice region, unless a canonical site holds the insertion: both flanks
///   in `[lo + 2, lo + 7]`, `[hi - 7, hi - 2]`, `[lo - 3, lo - 1]` or
///   `[hi + 1, hi + 3]`, or `r == lo`, `l == hi`, `r == lo + 2`,
///   `l == hi - 2`;
/// - polypyrimidine (acceptor side, `[hi - 16, hi - 2]` on the plus strand,
///   `[lo + 2, lo + 16]` on the minus strand) and donor region
///   (`[lo + 2, lo + 5]` / `[hi - 5, hi - 2]`), both folded into the splice
///   region and, as in VEP, not suppressed by a canonical site.
///
/// VEP overwrites the splice-region flag intron by intron, so on an exon
/// short enough for two introns' windows to reach the insertion its result
/// depends on the order it visits them; vareffect takes either intron's
/// region.
///
/// No effects for an insertion at or beyond either transcript end.
///
/// # Errors
///
/// Returns an error if the transcript's exons do not cover its span.
pub(crate) fn insertion_intron_effects(
    pos: u64,
    transcript: &TranscriptModel,
) -> Result<InsertionIntronEffects, VarEffectError> {
    let mut effects = InsertionIntronEffects::default();
    let exons = &transcript.exons;
    // An insertion at or beyond a transcript end has a flank outside every
    // exon and intron and is upstream or downstream of the transcript.
    if exons.len() < 2 || pos <= transcript.tx_start || pos >= transcript.tx_end {
        return Ok(effects);
    }
    let (l, r) = (pos - 1, pos);

    // The flanks can only reach the introns next to their own exon or intron.
    let mut candidates = Vec::with_capacity(4);
    for flank in [l, r] {
        let anchor = match find_exon_or_intron(flank, transcript)? {
            ExonOrIntron::Exon(i) => i,
            ExonOrIntron::Intron { upstream, .. } => upstream,
        };
        for intron in anchor.saturating_sub(1)..=(anchor + 1).min(exons.len() - 2) {
            if !candidates.contains(&intron) {
                candidates.push(intron);
            }
        }
    }

    let mut region = false;
    for intron in candidates {
        let (a, b) = (&exons[intron], &exons[intron + 1]);
        let (lo, hi_excl) = if a.genomic_end <= b.genomic_start {
            (a.genomic_end, b.genomic_start)
        } else {
            (b.genomic_end, a.genomic_start)
        };
        // Abutting exons leave no intron (VEP: `next unless start <= end`).
        if hi_excl <= lo {
            continue;
        }
        let hi = hi_excl - 1;
        let index = to_u16(intron);
        let both_in = |x: u64, y: u64| l >= x && r <= y;
        let either_in = |x: u64, y: u64| r >= x && l <= y;
        let frameshift = hi - lo < FRAMESHIFT_INTRON_MAX_GAP;
        if frameshift && both_in(lo, hi) {
            continue;
        }
        let near_intron = l <= hi + 3 && r + 3 >= lo;
        let near_boundary = (l <= lo + 7 && r + 3 >= lo) || (l <= hi + 3 && r + 7 >= hi);
        let lo3 = lo.saturating_sub(3);
        let (hi2, hi5, hi7, hi16) = (
            hi.saturating_sub(2),
            hi.saturating_sub(5),
            hi.saturating_sub(7),
            hi.saturating_sub(16),
        );

        if near_intron {
            if both_in(lo + 2, hi2) || r == lo + 2 || l == hi2 {
                effects.intronic.get_or_insert(index);
            }
            let polypyrimidine = match transcript.strand {
                Strand::Plus => either_in(hi16, hi2),
                Strand::Minus => either_in(lo + 2, lo + 16),
            };
            effects.splice_region |= polypyrimidine;
        }
        if near_boundary {
            let start_site = both_in(lo, lo + 1);
            let end_site = both_in(hi.saturating_sub(1), hi);
            let (donor, acceptor) = match transcript.strand {
                Strand::Plus => (start_site, end_site),
                Strand::Minus => (end_site, start_site),
            };
            if donor {
                effects.donor.get_or_insert(index);
            }
            if acceptor {
                effects.acceptor.get_or_insert(index);
            }
            let donor_region = match transcript.strand {
                Strand::Plus => both_in(lo + 2, lo + 5),
                Strand::Minus => both_in(hi5, hi2),
            };
            effects.splice_region |= donor_region;
            region |= both_in(lo + 2, lo + 7)
                || both_in(hi7, hi2)
                || both_in(lo3, lo.saturating_sub(1))
                || both_in(hi + 1, hi + 3)
                || r == lo
                || l == hi
                || r == lo + 2
                || l == hi2;
        }
    }
    effects.splice_region |= region && !effects.canonical();
    Ok(effects)
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

    if intron_gap(up_exon, down_exon) <= FRAMESHIFT_INTRON_MAX_GAP {
        return VariantLocation::FrameshiftIntron {
            intron_index,
            distance_to_nearest_exon: if donor_dist <= acceptor_dist {
                donor_dist as i64
            } else {
                -(acceptor_dist as i64)
            },
        };
    }

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
