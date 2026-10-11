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
pub(crate) const UPSTREAM_DOWNSTREAM_LIMIT: u64 = 5_000;

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

/// Largest exon gap, in bases, that is a frameshift intron rather than a
/// real intron. Such gaps correct the genome against the transcript; VEP
/// gives a variant inside one no splice or intron term.
///
/// Reference: `ensembl-variation` release/115 `BaseTranscriptVariation.pm`
/// `_create_intron_trees`: `abs($intron_end - $intron_start) <= 12` on
/// 1-based inclusive bounds, which is a gap of at most 13 bases (abutting
/// exons included).
pub const FRAMESHIFT_INTRON_MAX_GAP: u64 = 13;

/// Exonic splice region: 1-3 bases from the exon/intron boundary.
const SPLICE_REGION_EXON_MAX: u64 = 3;

/// Where a variant footprint lies relative to one transcript's span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Flank {
    /// The footprint overlaps `[tx_start, tx_end)`.
    Overlaps,
    /// 5' of the transcript in transcript terms, at this VEP `distance`.
    Upstream(u64),
    /// 3' of the transcript in transcript terms, at this VEP `distance`.
    Downstream(u64),
    /// Outside the transcript and past [`UPSTREAM_DOWNSTREAM_LIMIT`] (see
    /// [`classify_flank`]).
    Far,
}

/// Place a variant footprint beside a transcript, as VEP's `distance` does.
///
/// The distance counts from the nearer footprint edge to the nearer
/// transcript end in 1-based terms: a substitution or deletion next to the
/// transcript is at distance 1, an insertion between the transcript's end
/// base and the next base at distance 0 (VEP 116 ci8438 chr14:92071010
/// C>CCTG..., `downstream_gene_variant` with `distance` 0 on the minus-strand
/// ENST00000554672). This distance matched all 50,533 up/downstream rows of
/// VEP 116 on the issue-13 GENCODE corpus.
///
/// A substitution or deletion is near up to [`UPSTREAM_DOWNSTREAM_LIMIT`]
/// inclusive, an insertion only up to one less: VEP tests the limit against
/// the insertion's far flanking base (`_after_end` / `_before_start`). VEP
/// 116 on ENST00000628537 (chr2, plus strand): `2 55344525 T C` downstream at
/// 5000, `2 55344524 T TC` at 4999, `2 55344525 T TC` and `T TT` not
/// reported.
///
/// # Arguments
///
/// * `start` - 0-based footprint start
/// * `end` - 0-based exclusive footprint end; `start == end` is an insertion
///   between `start - 1` and `start`
/// * `transcript` - the transcript, on the variant's chromosome
///
/// # Returns
///
/// [`Flank::Overlaps`] when the footprint overlaps the transcript (an
/// insertion strictly inside it), otherwise the strand-aware side and
/// distance, or [`Flank::Far`] past the limit.
pub(crate) fn classify_flank(start: u64, end: u64, transcript: &TranscriptModel) -> Flank {
    let edge_base = u64::from(end > start);
    let (genomic_left, distance) = if end <= transcript.tx_start {
        (true, transcript.tx_start - end + edge_base)
    } else if start >= transcript.tx_end {
        (false, start - transcript.tx_end + edge_base)
    } else {
        return Flank::Overlaps;
    };
    if distance + u64::from(start == end) > UPSTREAM_DOWNSTREAM_LIMIT {
        return Flank::Far;
    }
    match (genomic_left, transcript.strand) {
        (true, Strand::Plus) | (false, Strand::Minus) => Flank::Upstream(distance),
        (true, Strand::Minus) | (false, Strand::Plus) => Flank::Downstream(distance),
    }
}

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
/// This enum is `Copy` -- no heap allocations. It is `#[non_exhaustive]`:
/// new transcript regions may be added in a minor release.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
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
        /// 0-based offset within the CDS in VEP's padded frame: the
        /// cumulative coding bases before this position plus
        /// [`LocateIndex::cds_start_pad`]. The first genomic coding base is
        /// offset `cds_start_pad()`, which is 0 for every complete CDS.
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

    /// Inside a frameshift intron: an exon gap of at most
    /// [`FRAMESHIFT_INTRON_MAX_GAP`] bases. Such a gap corrects the genome
    /// against the transcript sequence, so no splice or intron term applies.
    FrameshiftIntron {
        /// 0-based intron index.
        intron_index: u16,
        /// Signed distance to the nearest exon boundary, as for
        /// [`VariantLocation::Intron`]: positive on the donor side, negative
        /// on the acceptor side, donor on a tie. Maps to the HGVS intronic
        /// offset (`c.228+1`).
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
/// For TTN (~363 exons): ~2.5 KB per transcript.
#[derive(Debug, Clone)]
pub struct LocateIndex {
    /// `cumulative_cds[i]` = [`Self::cds_start_pad`] plus the sum of CDS
    /// segment lengths for segments `0..i`. Length =
    /// `cds_segments.len() + 1` (prefix-sum convention).
    cumulative_cds: Vec<u32>,
    /// Bases missing from the first codon of a 5'-incomplete CDS. See
    /// [`Self::cds_start_pad`].
    cds_start_pad: u32,
    /// `true` when some exon gap is a frameshift intron, so the per-variant
    /// frameshift-intron checks can be skipped for every other transcript
    /// (all of MANE 1.5).
    has_frameshift_intron: bool,
    /// `frameshift_junction[i]` is `true` when CDS segment `i` follows segment
    /// `i - 1` across a frameshift intron, a sequence correction where no
    /// exon-junction complex sits. Empty when `has_frameshift_intron` is
    /// false.
    frameshift_junction: Vec<bool>,
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
    /// 1-based codon numbers of the selenocysteine codons, ascending. Empty
    /// for every transcript without `selenocysteine_codons`.
    selenocysteine_codons: Vec<u32>,
}

impl LocateIndex {
    /// Build the index for a given transcript. O(n) in the number of exons.
    ///
    /// # Errors
    ///
    /// Returns [`VarEffectError::Malformed`] if a CDS phase is above 2, if a
    /// CDS boundary coordinate does
    /// not fall within any exon, if one exon holds more than one CDS
    /// segment (a ribosomal-slippage CDS, which the codon walk cannot
    /// represent), or if the first CDS phase is non-zero on a transcript
    /// marked [`TranscriptModel::cds_start_complete`] (a partial first codon
    /// cannot be a complete start), or if a selenocysteine codon base lies
    /// outside the CDS.
    pub fn build(transcript: &TranscriptModel) -> Result<Self, VarEffectError> {
        let malformed =
            |msg: &str| VarEffectError::Malformed(format!("{}: {}", transcript.accession, msg));

        // GFF3 phase on the first segment counts the bases to skip before
        // the first complete codon; Ensembl pads the 5' end with `N` so the
        // partial codon is codon 1 (`translateable_seq`, and BaseGXF.pm's
        // `_convert_phase` maps GFF3 phase p to start phase (3 - p) % 3).
        if transcript.cds_segments.iter().any(|seg| seg.phase > 2) {
            return Err(malformed("CDS phase outside 0..=2"));
        }
        let cds_start_pad = transcript
            .cds_segments
            .first()
            .map_or(0, |seg| (3 - u32::from(seg.phase)) % 3);
        if cds_start_pad != 0 && transcript.cds_start_complete {
            return Err(malformed(
                "non-zero first CDS phase on a CDS marked start-complete",
            ));
        }

        let mut cumulative_cds = Vec::with_capacity(transcript.cds_segments.len() + 1);
        cumulative_cds.push(cds_start_pad);
        for seg in &transcript.cds_segments {
            let prev = *cumulative_cds
                .last()
                .expect("cumulative_cds initialized with the pad");
            cumulative_cds.push(prev + (seg.genomic_end - seg.genomic_start) as u32);
        }

        let mut exon_to_cds_seg = vec![None; transcript.exons.len()];
        for (seg_idx, seg) in transcript.cds_segments.iter().enumerate() {
            let slot = exon_to_cds_seg
                .get_mut(seg.exon_index as usize)
                .ok_or_else(|| malformed("CDS segment refers to a missing exon"))?;
            if slot.is_some() {
                return Err(malformed("more than one CDS segment in one exon"));
            }
            *slot = Some(to_u16(seg_idx));
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

        let has_frameshift_intron = transcript
            .exons
            .windows(2)
            .any(|pair| helpers::frameshift_gap(&pair[0], &pair[1]).is_some());
        let frameshift_junction = if has_frameshift_intron {
            let segments = &transcript.cds_segments;
            (0..segments.len())
                .map(|i| {
                    i > 0 && {
                        let up = segments[i - 1].exon_index as usize;
                        let down = segments[i].exon_index as usize;
                        down == up + 1
                            && helpers::frameshift_gap(
                                &transcript.exons[up],
                                &transcript.exons[down],
                            )
                            .is_some()
                    }
                })
                .collect()
        } else {
            Vec::new()
        };

        let mut selenocysteine_codons = transcript
            .selenocysteine_codons
            .iter()
            .map(|&base| {
                selenocysteine_codon_number(transcript, &cumulative_cds, base).ok_or_else(|| {
                    malformed(&format!("selenocysteine codon at {base} outside the CDS"))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        selenocysteine_codons.sort_unstable();
        selenocysteine_codons.dedup();

        Ok(Self {
            cumulative_cds,
            cds_start_pad,
            has_frameshift_intron,
            frameshift_junction,
            exon_to_cds_seg,
            cds_start_exon_idx,
            cds_end_exon_idx,
            utr5_exonic_len,
            selenocysteine_codons,
        })
    }

    /// Prefix-sum array of CDS segment lengths.
    ///
    /// Used by `consequence` module to map a CDS offset back to a genomic
    /// position via binary search.
    pub(crate) fn cumulative_cds(&self) -> &[u32] {
        &self.cumulative_cds
    }

    /// Length of the CDS in VEP's padded frame: the coding bases plus
    /// [`Self::cds_start_pad`]. Use it for frame arithmetic (codon and stop
    /// positions); compare genomic lengths against [`Self::coding_length`].
    ///
    /// Returns 0 for non-coding transcripts.
    pub(crate) fn total_cds_length(&self) -> u32 {
        *self.cumulative_cds.last().unwrap_or(&0)
    }

    /// `true` when the transcript has a frameshift intron (an exon gap of at
    /// most [`FRAMESHIFT_INTRON_MAX_GAP`] bases).
    pub(crate) fn has_frameshift_intron(&self) -> bool {
        self.has_frameshift_intron
    }

    /// `true` when CDS segment `seg_idx` follows the previous one across a
    /// frameshift intron rather than a real exon-exon junction.
    pub(crate) fn is_frameshift_junction(&self, seg_idx: usize) -> bool {
        self.frameshift_junction
            .get(seg_idx)
            .copied()
            .unwrap_or(false)
    }

    /// Number of coding bases present on the genome, without the pad.
    ///
    /// Returns 0 for non-coding transcripts.
    pub(crate) fn coding_length(&self) -> u32 {
        self.total_cds_length() - self.cds_start_pad
    }

    /// Bases missing from the first codon of a 5'-incomplete CDS: `(3 - p) %
    /// 3` for GFF3 phase `p` on the first CDS segment, so 0, 1 or 2.
    ///
    /// CDS offsets, codon numbers and HGVS `c.` positions count these bases,
    /// as VEP does: with a pad of 2 the first genomic coding base is `c.3`
    /// and codon 1 is the partial codon. Offsets below the pad have no
    /// genomic base. 0 for every complete CDS and for non-coding transcripts.
    pub fn cds_start_pad(&self) -> u32 {
        self.cds_start_pad
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

    /// 1-based codon numbers of the transcript's selenocysteine codons,
    /// ascending (empty for almost every transcript).
    pub(crate) fn selenocysteine_codons(&self) -> &[u32] {
        &self.selenocysteine_codons
    }

    /// `true` when a selenocysteine codon precedes codon `codon_number`
    /// (1-based). VEP 116 scans the whole alternate translation for the
    /// first stop and reads such a codon as one, so it finds no new stop
    /// after the variant (`_stop_loss_extra_AA`, see [`crate::hgvs_p`]).
    pub(crate) fn has_selenocysteine_before(&self, codon_number: u32) -> bool {
        self.selenocysteine_codons
            .first()
            .is_some_and(|&first| first < codon_number)
    }
}

/// 1-based codon number of the CDS base at 0-based genomic position `base`,
/// or `None` when no CDS segment holds it. Either base of a codon gives the
/// same number, so the lowest genomic base works on both strands.
fn selenocysteine_codon_number(
    transcript: &TranscriptModel,
    cumulative_cds: &[u32],
    base: u64,
) -> Option<u32> {
    let (seg_idx, seg) = transcript
        .cds_segments
        .iter()
        .enumerate()
        .find(|(_, seg)| seg.genomic_start <= base && base < seg.genomic_end)?;
    let intra = match transcript.strand {
        Strand::Plus => base - seg.genomic_start,
        Strand::Minus => seg.genomic_end - 1 - base,
    } as u32;
    Some((cumulative_cds[seg_idx] + intra) / 3 + 1)
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
