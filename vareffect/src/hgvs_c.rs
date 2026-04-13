//! HGVS coding DNA notation (c. / n.) generation.
//!
//! Generates HGVS c. notation for coding transcripts and n. notation for
//! non-coding transcripts. All position data (CDS offsets, UTR offsets,
//! intronic distances) is already computed by the `locate` module; this
//! module is a pure formatter.
//!
//! Reference: HGVS Nomenclature v21.1 (2024), hgvs-nomenclature.org.
//!
//! # Concordance
//!
//! Output matches VEP `--hgvs --shift_hgvs` (release 115).

use std::fmt;

use crate::codon::{complement, reverse_complement};
use crate::consequence::helpers::compute_cdna_position_exonic;
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::locate::helpers::{
    compute_cds_offset, compute_utr_offset_3prime, compute_utr_offset_5prime,
};
use crate::locate::{LocateIndex, SpliceSide, VariantLocation, locate_variant};
use crate::types::{Strand, TranscriptModel};

// ---------------------------------------------------------------------------
// HgvsPosition — internal representation of a c./n. position component
// ---------------------------------------------------------------------------

/// HGVS c./n. position before string formatting.
///
/// Intronic positions are NOT represented here; they are formatted by
/// combining an anchor `HgvsPosition` with a signed offset.
#[derive(Debug, Clone, PartialEq, Eq)]
enum HgvsPosition {
    /// CDS position: 1-based (e.g., c.742).
    Cds(u32),
    /// 5'UTR position: negative offset from CDS start (e.g., c.-15).
    FivePrimeUtr(i64),
    /// 3'UTR position: positive offset from CDS end (e.g., c.*42).
    ThreePrimeUtr(i64),
    /// Non-coding transcript position: 1-based cDNA (e.g., n.76).
    NonCoding(u32),
}

impl fmt::Display for HgvsPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cds(pos) => write!(f, "{pos}"),
            Self::FivePrimeUtr(offset) => write!(f, "{offset}"),
            Self::ThreePrimeUtr(offset) => write!(f, "*{offset}"),
            Self::NonCoding(pos) => write!(f, "{pos}"),
        }
    }
}

// ---------------------------------------------------------------------------
// Intronic position formatting
// ---------------------------------------------------------------------------

/// Format an HGVS position string, including intronic offset if present.
///
/// For exonic positions: returns e.g. `"742"`, `"-15"`, `"*42"`.
/// For intronic positions: returns e.g. `"742+5"`, `"-15+2"`, `"*37-1"`.
fn format_position(anchor: &HgvsPosition, intronic_offset: Option<i64>) -> String {
    match intronic_offset {
        None => anchor.to_string(),
        Some(offset) if offset > 0 => format!("{anchor}+{offset}"),
        Some(offset) => {
            // offset is negative; Display of negative i64 includes the minus
            format!("{anchor}{offset}")
        }
    }
}

// ---------------------------------------------------------------------------
// Exon boundary anchor computation
// ---------------------------------------------------------------------------

/// Compute the HGVS anchor position for an exon boundary flanking an intron.
///
/// For donor side (`is_donor=true`): returns position of the last exonic base
/// of `exons[intron_index]` (the upstream exon in transcript order).
///
/// For acceptor side (`is_donor=false`): returns position of the first exonic
/// base of `exons[intron_index + 1]` (the downstream exon in transcript order).
///
/// # Returns
///
/// `HgvsPosition::Cds`, `FivePrimeUtr`, `ThreePrimeUtr`, or `NonCoding`
/// depending on where the boundary falls and the transcript's biotype.
fn exon_boundary_hgvs_anchor(
    intron_index: u16,
    is_donor: bool,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<HgvsPosition, VarEffectError> {
    // Non-coding transcript: use cDNA position directly.
    if transcript.cds_segments.is_empty() {
        let exon_idx = if is_donor {
            intron_index as usize
        } else {
            intron_index as usize + 1
        };
        let exon = &transcript.exons[exon_idx];
        let boundary_pos = match (is_donor, transcript.strand) {
            (true, Strand::Plus) => exon.genomic_end - 1,
            (true, Strand::Minus) => exon.genomic_start,
            (false, Strand::Plus) => exon.genomic_start,
            (false, Strand::Minus) => exon.genomic_end - 1,
        };
        let cdna = compute_cdna_position_exonic(boundary_pos, transcript).unwrap_or(1);
        return Ok(HgvsPosition::NonCoding(cdna));
    }

    // Coding transcript: determine which exon and the genomic boundary position.
    let exon_idx = if is_donor {
        intron_index as usize
    } else {
        intron_index as usize + 1
    };
    let exon = &transcript.exons[exon_idx];

    // Genomic position of the exon boundary base.
    let boundary_pos = match (is_donor, transcript.strand) {
        (true, Strand::Plus) => exon.genomic_end - 1,
        (true, Strand::Minus) => exon.genomic_start,
        (false, Strand::Plus) => exon.genomic_start,
        (false, Strand::Minus) => exon.genomic_end - 1,
    };

    // Classify boundary as CDS, 5'UTR, or 3'UTR.
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
        Strand::Plus => boundary_pos < cds_start,
        Strand::Minus => boundary_pos >= cds_end,
    };
    let is_3utr = match transcript.strand {
        Strand::Plus => boundary_pos >= cds_end,
        Strand::Minus => boundary_pos < cds_start,
    };

    if is_5utr {
        let offset = compute_utr_offset_5prime(boundary_pos, exon_idx, transcript, index)?;
        Ok(HgvsPosition::FivePrimeUtr(offset))
    } else if is_3utr {
        let offset = compute_utr_offset_3prime(boundary_pos, exon_idx, transcript, index)?;
        Ok(HgvsPosition::ThreePrimeUtr(offset))
    } else {
        let cds = compute_cds_offset(boundary_pos, exon_idx, transcript, index)?;
        Ok(HgvsPosition::Cds(cds.cds_offset + 1))
    }
}

// ---------------------------------------------------------------------------
// VariantLocation -> HGVS position string
// ---------------------------------------------------------------------------

/// Convert a `VariantLocation` to a formatted HGVS position string.
///
/// Returns `None` for `Upstream`, `Downstream`, and `Distal` variants which
/// have no HGVS c./n. notation.
fn position_for_variant_location(
    location: &VariantLocation,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    match *location {
        VariantLocation::CdsExon { cds_offset, .. } => {
            let pos = HgvsPosition::Cds(cds_offset + 1);
            Ok(Some(format_position(&pos, None)))
        }
        VariantLocation::FivePrimeUtr {
            offset_from_cds_start,
            ..
        } => {
            let pos = HgvsPosition::FivePrimeUtr(offset_from_cds_start);
            Ok(Some(format_position(&pos, None)))
        }
        VariantLocation::ThreePrimeUtr {
            offset_from_cds_end,
            ..
        } => {
            let pos = HgvsPosition::ThreePrimeUtr(offset_from_cds_end);
            Ok(Some(format_position(&pos, None)))
        }
        VariantLocation::SpliceDonor {
            intron_index,
            offset,
        } => {
            let anchor = exon_boundary_hgvs_anchor(intron_index, true, transcript, index)?;
            Ok(Some(format_position(&anchor, Some(offset as i64))))
        }
        VariantLocation::SpliceAcceptor {
            intron_index,
            offset,
        } => {
            let anchor = exon_boundary_hgvs_anchor(intron_index, false, transcript, index)?;
            Ok(Some(format_position(&anchor, Some(-(offset as i64)))))
        }
        VariantLocation::SpliceRegion {
            intron_index,
            side,
            distance,
        } => {
            let is_donor = matches!(side, SpliceSide::Donor);
            let anchor = exon_boundary_hgvs_anchor(intron_index, is_donor, transcript, index)?;
            Ok(Some(format_position(&anchor, Some(distance))))
        }
        VariantLocation::Intron {
            intron_index,
            distance_to_nearest_exon,
        } => {
            let is_donor = distance_to_nearest_exon > 0;
            let anchor = exon_boundary_hgvs_anchor(intron_index, is_donor, transcript, index)?;
            Ok(Some(format_position(
                &anchor,
                Some(distance_to_nearest_exon),
            )))
        }
        VariantLocation::NonCodingExon { exon_index, .. } => {
            let exon = &transcript.exons[exon_index as usize];
            // For NonCodingExon, we don't have the genomic position in the enum.
            // This path is only used by format_snv_hgvs which passes pos separately.
            // Return a placeholder; the caller should use position_for_genomic instead.
            // However, since this is called from format_snv_hgvs which has pos,
            // we handle it there directly.
            let _ = exon;
            Ok(None)
        }
        VariantLocation::NonCodingIntron {
            intron_index,
            distance_to_nearest_exon,
        } => {
            let is_donor = distance_to_nearest_exon > 0;
            let anchor = exon_boundary_hgvs_anchor(intron_index, is_donor, transcript, index)?;
            Ok(Some(format_position(
                &anchor,
                Some(distance_to_nearest_exon),
            )))
        }
        VariantLocation::Upstream { .. }
        | VariantLocation::Downstream { .. }
        | VariantLocation::Distal => Ok(None),
    }
}

/// Compute the HGVS position string for a genomic coordinate.
///
/// Delegates to [`locate_variant`] to classify the position, then converts
/// the result via [`position_for_variant_location`].
///
/// This is the workhorse for indel start/end positions. O(log n) in exons.
fn position_for_genomic(
    chrom: &str,
    pos: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    let loc = locate_variant(chrom, pos, pos + 1, transcript, index)?;

    // NonCodingExon needs the genomic pos for cDNA computation.
    if let VariantLocation::NonCodingExon { .. } = loc {
        let cdna = compute_cdna_position_exonic(pos, transcript).unwrap_or(1);
        return Ok(Some(HgvsPosition::NonCoding(cdna).to_string()));
    }

    position_for_variant_location(&loc, transcript, index)
}

// ---------------------------------------------------------------------------
// Duplication detection
// ---------------------------------------------------------------------------

/// Check if an insertion should be described as a duplication.
///
/// An insertion is a dup when the inserted sequence matches the reference
/// immediately 5' (upstream in transcript direction) of the insertion point.
///
/// Both FASTA bytes and VCF ALT bytes are in plus-strand orientation, so
/// no complement is needed for the comparison itself.
fn is_duplication(
    ins_pos: u64,
    inserted_bases: &[u8],
    chrom: &str,
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<bool, VarEffectError> {
    let ins_len = inserted_bases.len() as u64;
    if ins_len == 0 {
        return Ok(false);
    }

    // "5' on coding strand" direction determines where to fetch reference.
    let (fetch_start, fetch_end) = match transcript.strand {
        Strand::Plus => {
            // 5' = lower genomic coordinates
            if ins_pos < ins_len {
                return Ok(false); // out of bounds
            }
            (ins_pos - ins_len, ins_pos)
        }
        Strand::Minus => {
            // 5' = higher genomic coordinates
            let end = ins_pos + ins_len;
            if let Some(chrom_len) = fasta.chrom_length(chrom)
                && end > chrom_len
            {
                return Ok(false);
            }
            (ins_pos, end)
        }
    };

    let ref_seq = fasta.fetch_sequence(chrom, fetch_start, fetch_end)?;
    Ok(ref_seq == inserted_bases)
}

// ---------------------------------------------------------------------------
// Coding-strand allele formatting
// ---------------------------------------------------------------------------

/// Convert a single VCF (plus-strand) base to coding-strand uppercase char.
fn coding_strand_base(base: u8, strand: Strand) -> char {
    let b = match strand {
        Strand::Plus => base,
        Strand::Minus => complement(base),
    };
    b.to_ascii_uppercase() as char
}

/// Convert VCF (plus-strand) bases to coding-strand uppercase string.
fn coding_strand_seq(bases: &[u8], strand: Strand) -> String {
    match strand {
        Strand::Plus => bases
            .iter()
            .map(|&b| b.to_ascii_uppercase() as char)
            .collect(),
        Strand::Minus => reverse_complement(bases)
            .iter()
            .map(|&b| b.to_ascii_uppercase() as char)
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Public formatting functions
// ---------------------------------------------------------------------------

/// Format HGVS c./n. notation for an SNV.
///
/// Returns the full string including accession prefix (e.g.,
/// `"NM_000546.6:c.742C>T"` or `"NR_001566.3:n.76A>G"`).
///
/// Returns `Ok(None)` for upstream/downstream/distal variants.
///
/// # Arguments
///
/// * `pos` -- 0-based genomic position
/// * `ref_base` -- plus-strand VCF REF base (NOT coding-strand)
/// * `alt_base` -- plus-strand VCF ALT base (NOT coding-strand)
/// * `location` -- variant location within the transcript
/// * `transcript` -- transcript model
/// * `index` -- precomputed locate index
pub(crate) fn format_snv_hgvs(
    pos: u64,
    ref_base: u8,
    alt_base: u8,
    location: &VariantLocation,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    // Handle NonCodingExon specially (needs genomic pos for cDNA computation).
    if let VariantLocation::NonCodingExon { .. } = location {
        let cdna = compute_cdna_position_exonic(pos, transcript).unwrap_or(1);
        let hgvs_pos = HgvsPosition::NonCoding(cdna);
        let r = coding_strand_base(ref_base, transcript.strand);
        let a = coding_strand_base(alt_base, transcript.strand);
        return Ok(Some(format!(
            "{}:n.{}{}>{}",
            transcript.accession, hgvs_pos, r, a,
        )));
    }

    let pos_str = match position_for_variant_location(location, transcript, index)? {
        Some(s) => s,
        None => return Ok(None),
    };

    let prefix = if transcript.cds_segments.is_empty() {
        "n."
    } else {
        "c."
    };

    // Determine prefix from location type for intronic non-coding.
    let prefix = match location {
        VariantLocation::NonCodingIntron { .. } => "n.",
        _ => prefix,
    };

    let r = coding_strand_base(ref_base, transcript.strand);
    let a = coding_strand_base(alt_base, transcript.strand);

    Ok(Some(format!(
        "{}:{prefix}{pos_str}{r}>{a}",
        transcript.accession,
    )))
}

/// Format HGVS c./n. notation for a deletion.
///
/// Returns the full string (e.g., `"NM_000546.6:c.742del"` or
/// `"NM_000546.6:c.742_744del"`).
///
/// Returns `Ok(None)` if either endpoint is upstream/downstream/distal.
///
/// # Arguments
///
/// * `chrom` -- UCSC-style chromosome name
/// * `start` -- 0-based genomic start (inclusive, first deleted base)
/// * `end` -- 0-based genomic end (exclusive)
/// * `transcript` -- transcript model
/// * `index` -- precomputed locate index
pub(crate) fn format_deletion_hgvs(
    chrom: &str,
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    let start_pos = match position_for_genomic(chrom, start, transcript, index)? {
        Some(s) => s,
        None => return Ok(None),
    };

    let prefix = deletion_prefix(chrom, start, transcript, index)?;

    let del_len = end - start;
    if del_len == 1 {
        Ok(Some(format!(
            "{}:{prefix}{start_pos}del",
            transcript.accession,
        )))
    } else {
        // On minus strand, genomic start (lower coordinate) maps to a HIGHER
        // CDS position than genomic end-1. Swap so HGVS range is in transcript
        // order (lower c. position first), matching format_duplication().
        let (first_genomic, last_genomic) = match transcript.strand {
            Strand::Plus => (start, end - 1),
            Strand::Minus => (end - 1, start),
        };
        let first_pos = match position_for_genomic(chrom, first_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        let last_pos = match position_for_genomic(chrom, last_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        Ok(Some(format!(
            "{}:{prefix}{first_pos}_{last_pos}del",
            transcript.accession,
        )))
    }
}

/// Format HGVS c./n. notation for an insertion (including duplication
/// detection).
///
/// Returns the full string (e.g., `"NM_000546.6:c.742_743insACG"` or
/// `"NM_000546.6:c.742dup"`).
///
/// Returns `Ok(None)` if either flanking position is upstream/downstream/distal.
///
/// # Arguments
///
/// * `pos` -- 0-based genomic insertion point (insertion is between `pos-1`
///   and `pos`)
/// * `inserted_bases` -- plus-strand VCF ALT bases (NOT coding-strand)
/// * `chrom` -- UCSC-style chromosome name
/// * `transcript` -- transcript model
/// * `index` -- precomputed locate index
/// * `fasta` -- reference FASTA reader (for duplication detection)
pub(crate) fn format_insertion_hgvs(
    pos: u64,
    inserted_bases: &[u8],
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<Option<String>, VarEffectError> {
    let ins_len = inserted_bases.len() as u64;

    // Check for duplication first.
    if is_duplication(pos, inserted_bases, chrom, transcript, fasta)? {
        return format_duplication(pos, ins_len, chrom, transcript, index);
    }

    // Regular insertion: flanking positions. On minus strand, genomic pos
    // (higher coordinate) maps to a LOWER CDS position than pos-1, so swap
    // to keep HGVS range in transcript order (lower c. position first).
    if pos == 0 {
        return Ok(None);
    }
    let (left_genomic, right_genomic) = match transcript.strand {
        Strand::Plus => (pos - 1, pos),
        Strand::Minus => (pos, pos - 1),
    };

    let left_pos = match position_for_genomic(chrom, left_genomic, transcript, index)? {
        Some(s) => s,
        None => return Ok(None),
    };

    let right_pos = match position_for_genomic(chrom, right_genomic, transcript, index)? {
        Some(s) => s,
        None => return Ok(None),
    };

    let prefix = insertion_prefix(chrom, pos, transcript, index)?;
    let coding_seq = coding_strand_seq(inserted_bases, transcript.strand);

    Ok(Some(format!(
        "{}:{prefix}{left_pos}_{right_pos}ins{coding_seq}",
        transcript.accession,
    )))
}

/// Format HGVS c./n. notation for a delins (or MNV).
///
/// Returns the full string (e.g., `"NM_000546.6:c.742delinsAT"` or
/// `"NM_000546.6:c.742_743delinsACGT"`).
///
/// Returns `Ok(None)` if either endpoint is upstream/downstream/distal.
///
/// # Arguments
///
/// * `chrom` -- UCSC-style chromosome name
/// * `start` -- 0-based genomic start (inclusive)
/// * `end` -- 0-based genomic end (exclusive)
/// * `alt_bases` -- plus-strand VCF ALT bases (NOT coding-strand)
/// * `transcript` -- transcript model
/// * `index` -- precomputed locate index
pub(crate) fn format_delins_hgvs(
    chrom: &str,
    start: u64,
    end: u64,
    alt_bases: &[u8],
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    let start_pos = match position_for_genomic(chrom, start, transcript, index)? {
        Some(s) => s,
        None => return Ok(None),
    };

    let prefix = deletion_prefix(chrom, start, transcript, index)?;
    let coding_alt = coding_strand_seq(alt_bases, transcript.strand);

    let ref_len = end - start;
    if ref_len == 1 {
        Ok(Some(format!(
            "{}:{prefix}{start_pos}delins{coding_alt}",
            transcript.accession,
        )))
    } else {
        // Same strand-aware swap as format_deletion_hgvs: on minus strand,
        // genomic start maps to a higher CDS position than genomic end-1.
        let (first_genomic, last_genomic) = match transcript.strand {
            Strand::Plus => (start, end - 1),
            Strand::Minus => (end - 1, start),
        };
        let first_pos = match position_for_genomic(chrom, first_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        let last_pos = match position_for_genomic(chrom, last_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        Ok(Some(format!(
            "{}:{prefix}{first_pos}_{last_pos}delins{coding_alt}",
            transcript.accession,
        )))
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Format a duplication notation.
fn format_duplication(
    ins_pos: u64,
    ins_len: u64,
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<Option<String>, VarEffectError> {
    // Compute genomic range of the duplicated reference bases.
    let (dup_start, dup_end) = match transcript.strand {
        Strand::Plus => (ins_pos - ins_len, ins_pos),
        Strand::Minus => (ins_pos, ins_pos + ins_len),
    };

    let prefix = deletion_prefix(chrom, dup_start, transcript, index)?;

    if ins_len == 1 {
        let pos_str = match position_for_genomic(chrom, dup_start, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        Ok(Some(format!(
            "{}:{prefix}{pos_str}dup",
            transcript.accession,
        )))
    } else {
        // For multi-base dups: HGVS start < end in c. numbering.
        // For plus strand: dup_start has lower c. position.
        // For minus strand: dup_end-1 has lower c. position (higher genomic).
        let (first_genomic, last_genomic) = match transcript.strand {
            Strand::Plus => (dup_start, dup_end - 1),
            Strand::Minus => (dup_end - 1, dup_start),
        };

        let start_pos = match position_for_genomic(chrom, first_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        let end_pos = match position_for_genomic(chrom, last_genomic, transcript, index)? {
            Some(s) => s,
            None => return Ok(None),
        };
        Ok(Some(format!(
            "{}:{prefix}{start_pos}_{end_pos}dup",
            transcript.accession,
        )))
    }
}

/// Determine the HGVS prefix (c. or n.) for a deletion/delins/dup based on
/// the genomic position.
fn deletion_prefix(
    chrom: &str,
    pos: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<&'static str, VarEffectError> {
    if transcript.cds_segments.is_empty() {
        return Ok("n.");
    }
    let loc = locate_variant(chrom, pos, pos + 1, transcript, index)?;
    match loc {
        VariantLocation::NonCodingExon { .. } | VariantLocation::NonCodingIntron { .. } => Ok("n."),
        _ => Ok("c."),
    }
}

/// Determine the HGVS prefix for an insertion based on flanking positions.
fn insertion_prefix(
    chrom: &str,
    pos: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<&'static str, VarEffectError> {
    // Use the right flanking position (the insertion point itself).
    deletion_prefix(chrom, pos, transcript, index)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test HgvsPosition Display
    #[test]
    fn display_cds_position() {
        assert_eq!(HgvsPosition::Cds(742).to_string(), "742");
        assert_eq!(HgvsPosition::Cds(1).to_string(), "1");
    }

    #[test]
    fn display_five_prime_utr_position() {
        assert_eq!(HgvsPosition::FivePrimeUtr(-15).to_string(), "-15");
        assert_eq!(HgvsPosition::FivePrimeUtr(-1).to_string(), "-1");
    }

    #[test]
    fn display_three_prime_utr_position() {
        assert_eq!(HgvsPosition::ThreePrimeUtr(42).to_string(), "*42");
        assert_eq!(HgvsPosition::ThreePrimeUtr(1).to_string(), "*1");
    }

    #[test]
    fn display_non_coding_position() {
        assert_eq!(HgvsPosition::NonCoding(76).to_string(), "76");
    }

    // Test format_position with intronic offsets
    #[test]
    fn format_position_exonic() {
        let pos = HgvsPosition::Cds(742);
        assert_eq!(format_position(&pos, None), "742");
    }

    #[test]
    fn format_position_donor_intronic() {
        let pos = HgvsPosition::Cds(742);
        assert_eq!(format_position(&pos, Some(5)), "742+5");
    }

    #[test]
    fn format_position_acceptor_intronic() {
        let pos = HgvsPosition::Cds(743);
        assert_eq!(format_position(&pos, Some(-3)), "743-3");
    }

    #[test]
    fn format_position_utr_intronic() {
        let pos = HgvsPosition::FivePrimeUtr(-15);
        assert_eq!(format_position(&pos, Some(2)), "-15+2");

        let pos = HgvsPosition::ThreePrimeUtr(37);
        assert_eq!(format_position(&pos, Some(-1)), "*37-1");
    }

    #[test]
    fn format_position_noncoding_intronic() {
        let pos = HgvsPosition::NonCoding(76);
        assert_eq!(format_position(&pos, Some(1)), "76+1");
    }

    // Test coding_strand_base
    #[test]
    fn coding_strand_base_plus() {
        assert_eq!(coding_strand_base(b'A', Strand::Plus), 'A');
        assert_eq!(coding_strand_base(b'C', Strand::Plus), 'C');
    }

    #[test]
    fn coding_strand_base_minus() {
        assert_eq!(coding_strand_base(b'A', Strand::Minus), 'T');
        assert_eq!(coding_strand_base(b'G', Strand::Minus), 'C');
    }

    // Test coding_strand_seq
    #[test]
    fn coding_strand_seq_plus() {
        assert_eq!(coding_strand_seq(b"ACG", Strand::Plus), "ACG");
    }

    #[test]
    fn coding_strand_seq_minus() {
        assert_eq!(coding_strand_seq(b"ACG", Strand::Minus), "CGT");
    }
}
