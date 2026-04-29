//! HGVS c. reverse mapper: parse HGVS c. notation and resolve to genomic
//! VCF-style coordinates.
//!
//! This is the inverse of the forward mapper in [`crate::hgvs_c`]. Given a
//! string like `"NM_000546.6:c.742C>T"`, it looks up the transcript, maps the
//! c. position to a 0-based genomic coordinate, and returns a VCF-style
//! `(chrom, pos, ref, alt)` tuple suitable for [`crate::VarEffect::annotate`].
//!
//! # Supported variant types
//!
//! | Type | Example |
//! |------|---------|
//! | Substitution | `c.742C>T` |
//! | Deletion | `c.19del`, `c.19_21del` |
//! | Duplication | `c.20dup`, `c.20_23dup` |
//! | Insertion | `c.76_77insT` |
//! | Delins | `c.113delinsTACT`, `c.112_117delinsTG` |
//!
//! # Position types
//!
//! CDS (`c.742`), 5'UTR (`c.-14`), 3'UTR (`c.*32`), and intronic offsets
//! (`c.672+1`, `c.673-2`) are all handled, including combined forms like
//! `c.-15+1` (intronic position in the 5'UTR) and `c.*37-1` (intronic
//! position in the 3'UTR).

use crate::codon::{complement, reverse_complement};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::locate::LocateIndex;
use crate::transcript::TranscriptStore;
use crate::types::{Biotype, Strand, TranscriptModel};

// ---------------------------------------------------------------------------
// Public output type
// ---------------------------------------------------------------------------

/// A genomic variant in VCF-style coordinates.
///
/// All coordinates are 0-based (matching the internal convention used
/// throughout `vareffect`). Alleles are plus-strand, uppercase ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicVariant {
    /// UCSC-style chromosome name (e.g., `"chr17"`).
    pub chrom: String,
    /// 0-based genomic position of the first REF base.
    pub pos: u64,
    /// Reference allele (plus-strand, uppercase ASCII).
    pub ref_allele: Vec<u8>,
    /// Alternate allele (plus-strand, uppercase ASCII).
    pub alt_allele: Vec<u8>,
}

/// HGVS c. resolution result with transcript-version provenance.
///
/// Returned by [`crate::VarEffect::resolve_hgvs_c_with_meta`] so callers
/// can detect when the transcript version they asked for was not present
/// in the store and a different version was used instead. `resolved_accession`
/// holds the accession (with version) that was actually matched — when it
/// differs from the caller's input accession the caller is observing
/// transcript-version drift.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedHgvsC {
    /// Genomic coordinates reverse-mapped from the HGVS c. notation.
    pub variant: GenomicVariant,
    /// Accession (with version) actually used from the transcript store.
    pub resolved_accession: String,
}

// ---------------------------------------------------------------------------
// Parsed HGVS c. types
// ---------------------------------------------------------------------------

/// A single c. position component (the numeric part of an HGVS c. position).
///
/// Examples:
/// - `c.742`      → base=742, is_3utr=false, intronic_offset=None
/// - `c.-14`      → base=-14, is_3utr=false, intronic_offset=None
/// - `c.*32`      → base=32,  is_3utr=true,  intronic_offset=None
/// - `c.672+1`    → base=672, is_3utr=false, intronic_offset=Some(1)
/// - `c.673-2`    → base=673, is_3utr=false, intronic_offset=Some(-2)
/// - `c.-15+1`    → base=-15, is_3utr=false, intronic_offset=Some(1)
/// - `c.*37-1`    → base=37,  is_3utr=true,  intronic_offset=Some(-1)
#[derive(Debug, Clone, PartialEq, Eq)]
struct HgvsCPosition {
    /// Base position: positive for CDS (1-based), negative for 5'UTR.
    /// For 3'UTR, this is the offset after `*` (always positive).
    base: i64,
    /// Whether this is a 3'UTR position (c.*N).
    is_3utr: bool,
    /// Intronic offset: positive for donor side (+N), negative for acceptor
    /// side (-N). `None` for exonic positions.
    intronic_offset: Option<i64>,
}

/// The change portion of an HGVS c. variant.
#[derive(Debug, Clone, PartialEq, Eq)]
enum HgvsCChange {
    /// Single-base substitution: `C>T`.
    Substitution { ref_base: u8, alt_base: u8 },
    /// Deletion: `del` (bases determined from FASTA).
    Deletion,
    /// Insertion: `insACG`.
    Insertion { bases: Vec<u8> },
    /// Duplication: `dup` (bases determined from FASTA).
    Duplication,
    /// Deletion-insertion: `delinsTG`.
    Delins { bases: Vec<u8> },
}

/// A fully parsed HGVS c. notation string.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ParsedHgvsC {
    /// Transcript accession, possibly with version (e.g., `"NM_000546.6"`).
    accession: String,
    /// Start (or only) position.
    start: HgvsCPosition,
    /// End position for ranges. `None` for single-position variants.
    end: Option<HgvsCPosition>,
    /// The variant change type and associated data.
    change: HgvsCChange,
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse an HGVS c. notation string into structured components.
///
/// Accepts the full notation including accession prefix, e.g.
/// `"NM_000546.6:c.742C>T"`. The `c.` prefix is required.
fn parse_hgvs_c(input: &str) -> Result<ParsedHgvsC, VarEffectError> {
    let err = |msg: &str| VarEffectError::HgvsParseError(format!("{msg}: \"{input}\""));

    // Split accession from the c. description.
    let (accession, desc) = input
        .split_once(":c.")
        .ok_or_else(|| err("missing ':c.' separator"))?;

    if accession.is_empty() {
        return Err(err("empty accession"));
    }

    // Identify the change type. Order matters: "delins" before "del"/"ins".
    if let Some(idx) = desc.find("delins") {
        // Delins: positions are everything before "delins", inserted seq after.
        let pos_str = &desc[..idx];
        let inserted_str = &desc[idx + 6..]; // skip "delins"
        let inserted = parse_seq(inserted_str, input)?;
        let (start, end) = parse_position_range(pos_str, input)?;
        Ok(ParsedHgvsC {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsCChange::Delins { bases: inserted },
        })
    } else if let Some(ins_idx) = desc.find("ins")
        && desc.contains('_')
    {
        // Insertion: must have two flanking positions separated by '_'.
        // Format: {pos1}_{pos2}ins{seq}
        let pos_str = &desc[..ins_idx];
        let inserted_str = &desc[ins_idx + 3..]; // skip "ins"
        let inserted = parse_seq(inserted_str, input)?;
        let (start, end_opt) = parse_position_range(pos_str, input)?;
        let end = end_opt.ok_or_else(|| {
            VarEffectError::HgvsParseError(format!(
                "insertion requires two flanking positions: \"{input}\""
            ))
        })?;
        Ok(ParsedHgvsC {
            accession: accession.to_string(),
            start,
            end: Some(end),
            change: HgvsCChange::Insertion { bases: inserted },
        })
    } else if let Some(idx) = desc.find("del") {
        // Deletion: positions before "del", nothing meaningful after.
        let pos_str = &desc[..idx];
        let (start, end) = parse_position_range(pos_str, input)?;
        Ok(ParsedHgvsC {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsCChange::Deletion,
        })
    } else if let Some(idx) = desc.find("dup") {
        // Duplication: positions before "dup".
        let pos_str = &desc[..idx];
        let (start, end) = parse_position_range(pos_str, input)?;
        Ok(ParsedHgvsC {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsCChange::Duplication,
        })
    } else if let Some(gt_idx) = desc.find('>') {
        // Substitution: {pos}{ref}>{alt}
        // The REF base is the character immediately before '>'.
        if gt_idx == 0 {
            return Err(err("substitution missing REF base"));
        }
        let ref_base = desc.as_bytes()[gt_idx - 1];
        let pos_str = &desc[..gt_idx - 1];
        let alt_str = &desc[gt_idx + 1..];
        if alt_str.len() != 1 {
            return Err(err("substitution ALT must be a single base"));
        }
        let alt_base = alt_str.as_bytes()[0];
        validate_base(ref_base, input)?;
        validate_base(alt_base, input)?;
        let start = parse_single_position(pos_str, input)?;
        Ok(ParsedHgvsC {
            accession: accession.to_string(),
            start,
            end: None,
            change: HgvsCChange::Substitution { ref_base, alt_base },
        })
    } else {
        Err(err(
            "unrecognized change type (expected del/dup/ins/delins or >)",
        ))
    }
}

/// Parse a position range string, e.g. `"742"` or `"742_744"`.
///
/// Returns `(start, Some(end))` for ranges, `(start, None)` for single
/// positions.
fn parse_position_range(
    s: &str,
    input: &str,
) -> Result<(HgvsCPosition, Option<HgvsCPosition>), VarEffectError> {
    // Find the underscore that separates two positions. We must be careful:
    // position strings can contain '-' and '+' but never '_'. However the
    // first position may start with '*' or '-', neither of which is '_'.
    if let Some(sep_idx) = find_range_separator(s) {
        let left = parse_single_position(&s[..sep_idx], input)?;
        let right = parse_single_position(&s[sep_idx + 1..], input)?;
        Ok((left, Some(right)))
    } else {
        let pos = parse_single_position(s, input)?;
        Ok((pos, None))
    }
}

/// Find the underscore separating two positions in a range string.
///
/// The underscore must separate two positions. A position starts with an
/// optional `*` or `-`, followed by digits. After the digits there may be
/// `+` or `-` and more digits (intronic offset). The range separator `_`
/// appears between the end of the first position and the start of the second.
///
/// Strategy: walk forward past the first position, then look for `_`.
fn find_range_separator(s: &str) -> Option<usize> {
    let bytes = s.as_bytes();
    let mut i = 0;

    // Skip leading '*' or '-' of first position.
    if i < bytes.len() && (bytes[i] == b'*' || bytes[i] == b'-') {
        i += 1;
    }
    // Skip digits of base.
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    // Skip optional intronic offset (+N or -N).
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
    }
    // Now we should be at '_' if this is a range.
    if i < bytes.len() && bytes[i] == b'_' {
        Some(i)
    } else {
        None
    }
}

/// Parse a single c. position string (e.g. `"742"`, `"-14"`, `"*32"`,
/// `"672+1"`, `"-15+2"`).
fn parse_single_position(s: &str, input: &str) -> Result<HgvsCPosition, VarEffectError> {
    let err = |msg: &str| {
        VarEffectError::HgvsParseError(format!("{msg} in position \"{s}\": \"{input}\""))
    };

    if s.is_empty() {
        return Err(err("empty position"));
    }

    let bytes = s.as_bytes();
    let mut i = 0;

    // Check for 3'UTR marker '*'.
    let is_3utr = bytes[i] == b'*';
    if is_3utr {
        i += 1;
    }

    // Check for 5'UTR negative sign.
    let is_negative = !is_3utr && i < bytes.len() && bytes[i] == b'-';
    if is_negative {
        i += 1;
    }

    // Parse the base digits.
    let digit_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == digit_start {
        return Err(err("no digits in base position"));
    }
    let base_val: i64 = s[digit_start..i]
        .parse()
        .map_err(|_| err("base position integer overflow"))?;
    if base_val == 0 {
        return Err(err("c.0 is invalid (HGVS is 1-based)"));
    }
    let base = if is_negative { -base_val } else { base_val };

    // Parse optional intronic offset.
    let intronic_offset = if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        let sign_positive = bytes[i] == b'+';
        i += 1;
        let offset_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == offset_start {
            return Err(err("no digits in intronic offset"));
        }
        let offset_val: i64 = s[offset_start..i]
            .parse()
            .map_err(|_| err("intronic offset integer overflow"))?;
        Some(if sign_positive {
            offset_val
        } else {
            -offset_val
        })
    } else {
        None
    };

    if i != bytes.len() {
        return Err(err("unexpected trailing characters"));
    }

    Ok(HgvsCPosition {
        base,
        is_3utr,
        intronic_offset,
    })
}

/// Parse a nucleotide sequence string (e.g. `"ACG"`) into uppercase bytes.
fn parse_seq(s: &str, input: &str) -> Result<Vec<u8>, VarEffectError> {
    if s.is_empty() {
        return Err(VarEffectError::HgvsParseError(format!(
            "empty nucleotide sequence: \"{input}\""
        )));
    }
    let mut bases = Vec::with_capacity(s.len());
    for &b in s.as_bytes() {
        let upper = b.to_ascii_uppercase();
        validate_base(upper, input)?;
        bases.push(upper);
    }
    Ok(bases)
}

/// Validate that a byte is a standard nucleotide (A, C, G, T).
fn validate_base(b: u8, input: &str) -> Result<(), VarEffectError> {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => Ok(()),
        _ => Err(VarEffectError::HgvsParseError(format!(
            "invalid nucleotide '{}': \"{input}\"",
            b as char,
        ))),
    }
}

// ---------------------------------------------------------------------------
// Position-to-genomic mapper
// ---------------------------------------------------------------------------

/// Resolve an HGVS c. position to a 0-based genomic coordinate.
///
/// Handles CDS, 5'UTR, 3'UTR, and intronic positions. Requires a
/// protein-coding transcript (caller must validate `biotype`).
fn resolve_position(
    pos: &HgvsCPosition,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<u64, VarEffectError> {
    let accession = &transcript.accession;

    // Resolve the exonic anchor (ignoring intronic offset).
    let anchor = if pos.is_3utr {
        resolve_3utr(pos.base, transcript, index)?
    } else if pos.base < 0 {
        resolve_5utr(pos.base, transcript, index)?
    } else {
        resolve_cds(pos.base, transcript, index)?
    };

    // Apply intronic offset if present.
    match pos.intronic_offset {
        None => Ok(anchor),
        Some(offset) => {
            // Plus strand:  genomic = anchor + offset
            // Minus strand: genomic = anchor - offset
            // See plan for verification of all 4 combinations.
            let genomic = match transcript.strand {
                Strand::Plus => {
                    if offset >= 0 {
                        anchor.checked_add(offset as u64).ok_or_else(|| {
                            VarEffectError::PositionOutOfRange {
                                accession: accession.clone(),
                                detail: format!("intronic offset {offset} overflows"),
                            }
                        })?
                    } else {
                        anchor.checked_sub(offset.unsigned_abs()).ok_or_else(|| {
                            VarEffectError::PositionOutOfRange {
                                accession: accession.clone(),
                                detail: format!("intronic offset {offset} underflows"),
                            }
                        })?
                    }
                }
                Strand::Minus => {
                    if offset >= 0 {
                        anchor.checked_sub(offset as u64).ok_or_else(|| {
                            VarEffectError::PositionOutOfRange {
                                accession: accession.clone(),
                                detail: format!("intronic offset {offset} underflows"),
                            }
                        })?
                    } else {
                        anchor.checked_add(offset.unsigned_abs()).ok_or_else(|| {
                            VarEffectError::PositionOutOfRange {
                                accession: accession.clone(),
                                detail: format!("intronic offset {offset} overflows"),
                            }
                        })?
                    }
                }
            };
            Ok(genomic)
        }
    }
}

/// Resolve a CDS position (c.N, where N >= 1) to 0-based genomic.
///
/// Uses binary search on `cumulative_cds` — same pattern as
/// `consequence/helpers.rs:fetch_cds_sequence`.
fn resolve_cds(
    base: i64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<u64, VarEffectError> {
    debug_assert!(base > 0);
    let cds_offset = (base - 1) as u32;
    let cumulative = index.cumulative_cds();

    if cds_offset >= index.total_cds_length() {
        return Err(VarEffectError::PositionOutOfRange {
            accession: transcript.accession.clone(),
            detail: format!("c.{base} exceeds CDS length {}", index.total_cds_length()),
        });
    }

    // Binary search: find the CDS segment containing this offset.
    // cumulative[i] = sum of lengths of segments 0..i.
    // We want the largest i such that cumulative[i] <= cds_offset.
    let seg_idx = cumulative.partition_point(|&c| c <= cds_offset) - 1;
    let seg = &transcript.cds_segments[seg_idx];
    let local = cds_offset - cumulative[seg_idx];

    let genomic = match transcript.strand {
        Strand::Plus => seg.genomic_start + local as u64,
        Strand::Minus => seg.genomic_end - 1 - local as u64,
    };
    Ok(genomic)
}

/// Resolve a 5'UTR position (c.-N, where N >= 1) to 0-based genomic.
///
/// Inverse of `locate/helpers.rs:compute_utr_offset_5prime`. Walks exons
/// backward from the CDS start, consuming exonic bases.
fn resolve_5utr(
    base: i64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<u64, VarEffectError> {
    debug_assert!(base < 0);
    let mut remaining = base.unsigned_abs();
    let exons = &transcript.exons;
    let cds_start_exon = index.cds_start_exon_idx().ok_or_else(|| {
        VarEffectError::Malformed(format!(
            "{}: 5'UTR position but no cds_start_exon_idx",
            transcript.accession,
        ))
    })?;

    match transcript.strand {
        Strand::Plus => {
            let cds_start = transcript.cds_genomic_start.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: 5'UTR position but no cds_genomic_start",
                    transcript.accession,
                ))
            })?;

            // Bases available in the CDS-start exon's 5'UTR portion.
            let available = cds_start - exons[cds_start_exon].genomic_start;
            if remaining <= available {
                return Ok(cds_start - remaining);
            }
            remaining -= available;

            // Walk earlier exons.
            for i in (0..cds_start_exon).rev() {
                let exon_len = exons[i].genomic_end - exons[i].genomic_start;
                if remaining <= exon_len {
                    return Ok(exons[i].genomic_end - remaining);
                }
                remaining -= exon_len;
            }
        }
        Strand::Minus => {
            // On minus strand, transcript 5' = higher genomic coords.
            // CDS "start" in transcript terms = cds_genomic_end.
            let cds_boundary = transcript.cds_genomic_end.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: 5'UTR position but no cds_genomic_end",
                    transcript.accession,
                ))
            })?;

            let available = exons[cds_start_exon].genomic_end - cds_boundary;
            if remaining <= available {
                return Ok(cds_boundary + remaining - 1);
            }
            remaining -= available;

            for i in (0..cds_start_exon).rev() {
                let exon_len = exons[i].genomic_end - exons[i].genomic_start;
                if remaining <= exon_len {
                    return Ok(exons[i].genomic_start + remaining - 1);
                }
                remaining -= exon_len;
            }
        }
    }

    Err(VarEffectError::PositionOutOfRange {
        accession: transcript.accession.clone(),
        detail: format!("c.{base} exceeds 5'UTR exonic extent"),
    })
}

/// Resolve a 3'UTR position (c.*N, where N >= 1) to 0-based genomic.
///
/// Inverse of `locate/helpers.rs:compute_utr_offset_3prime`. Walks exons
/// forward from the CDS end, consuming exonic bases.
fn resolve_3utr(
    base: i64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<u64, VarEffectError> {
    debug_assert!(base > 0);
    let mut remaining = base as u64;
    let exons = &transcript.exons;
    let cds_end_exon = index.cds_end_exon_idx().ok_or_else(|| {
        VarEffectError::Malformed(format!(
            "{}: 3'UTR position but no cds_end_exon_idx",
            transcript.accession,
        ))
    })?;

    match transcript.strand {
        Strand::Plus => {
            let cds_end = transcript.cds_genomic_end.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: 3'UTR position but no cds_genomic_end",
                    transcript.accession,
                ))
            })?;

            // Bases available in the CDS-end exon's 3'UTR portion.
            let available = exons[cds_end_exon].genomic_end - cds_end;
            if remaining <= available {
                return Ok(cds_end + remaining - 1);
            }
            remaining -= available;

            // Walk later exons.
            for exon in &exons[(cds_end_exon + 1)..] {
                let exon_len = exon.genomic_end - exon.genomic_start;
                if remaining <= exon_len {
                    return Ok(exon.genomic_start + remaining - 1);
                }
                remaining -= exon_len;
            }
        }
        Strand::Minus => {
            // On minus strand, transcript 3' = lower genomic coords.
            // CDS "end" in transcript terms = cds_genomic_start.
            let cds_boundary = transcript.cds_genomic_start.ok_or_else(|| {
                VarEffectError::Malformed(format!(
                    "{}: 3'UTR position but no cds_genomic_start",
                    transcript.accession,
                ))
            })?;

            let available = cds_boundary - exons[cds_end_exon].genomic_start;
            if remaining <= available {
                return Ok(cds_boundary - remaining);
            }
            remaining -= available;

            for exon in &exons[(cds_end_exon + 1)..] {
                let exon_len = exon.genomic_end - exon.genomic_start;
                if remaining <= exon_len {
                    return Ok(exon.genomic_end - remaining);
                }
                remaining -= exon_len;
            }
        }
    }

    Err(VarEffectError::PositionOutOfRange {
        accession: transcript.accession.clone(),
        detail: format!("c.*{base} exceeds 3'UTR exonic extent"),
    })
}

// ---------------------------------------------------------------------------
// VCF coordinate construction
// ---------------------------------------------------------------------------

/// Convert a coding-strand base to plus-strand.
fn to_plus_strand_base(base: u8, strand: Strand) -> u8 {
    match strand {
        Strand::Plus => base.to_ascii_uppercase(),
        Strand::Minus => complement(base).to_ascii_uppercase(),
    }
}

/// Convert a coding-strand sequence to plus-strand.
fn to_plus_strand_seq(bases: &[u8], strand: Strand) -> Vec<u8> {
    match strand {
        Strand::Plus => bases.iter().map(|b| b.to_ascii_uppercase()).collect(),
        Strand::Minus => {
            let mut rc = reverse_complement(bases);
            rc.make_ascii_uppercase();
            rc
        }
    }
}

/// Build a VCF-style [`GenomicVariant`] for a substitution.
///
/// `hgvs` is the original input notation; it is threaded through so a
/// ref-allele disagreement can surface as [`VarEffectError::HgvsRefMismatch`]
/// with enough context to reproduce VEP's diagnostic wording.
fn build_substitution(
    parsed: &ParsedHgvsC,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    hgvs: &str,
) -> Result<GenomicVariant, VarEffectError> {
    let HgvsCChange::Substitution { ref_base, alt_base } = parsed.change else {
        unreachable!("build_substitution is only called for Substitution variants")
    };

    let genomic_pos = resolve_position(&parsed.start, transcript, index)?;
    let chrom = &transcript.chrom;

    // Convert coding-strand bases to plus-strand for VCF output.
    let plus_ref = to_plus_strand_base(ref_base, transcript.strand);
    let plus_alt = to_plus_strand_base(alt_base, transcript.strand);

    // Verify against FASTA.
    let fasta_ref = fasta.fetch_base(chrom, genomic_pos)?;
    if fasta_ref.to_ascii_uppercase() != plus_ref {
        return Err(VarEffectError::HgvsRefMismatch {
            hgvs: hgvs.to_string(),
            chrom: chrom.clone(),
            pos: genomic_pos,
            expected: String::from(fasta_ref as char),
            got: String::from(plus_ref as char),
        });
    }

    Ok(GenomicVariant {
        chrom: chrom.clone(),
        pos: genomic_pos,
        ref_allele: vec![plus_ref],
        alt_allele: vec![plus_alt],
    })
}

/// Resolve a position range to a normalized genomic interval `(gstart, gend)`
/// where `gstart <= gend` (both 0-based inclusive).
fn resolve_range(
    start: &HgvsCPosition,
    end: Option<&HgvsCPosition>,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Result<(u64, u64), VarEffectError> {
    let g_start = resolve_position(start, transcript, index)?;
    let g_end = match end {
        Some(e) => resolve_position(e, transcript, index)?,
        None => g_start,
    };
    // Normalize: on minus strand, start may have higher genomic coord than end.
    Ok((g_start.min(g_end), g_start.max(g_end)))
}

/// Build a VCF-style [`GenomicVariant`] for a deletion.
fn build_deletion(
    parsed: &ParsedHgvsC,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let (gstart, gend) = resolve_range(&parsed.start, parsed.end.as_ref(), transcript, index)?;
    let chrom = &transcript.chrom;

    // VCF convention: anchor base before the deletion.
    let anchor_pos = gstart
        .checked_sub(1)
        .ok_or_else(|| VarEffectError::PositionOutOfRange {
            accession: transcript.accession.clone(),
            detail: "deletion at position 0 has no anchor base".into(),
        })?;
    let anchor_base = fasta.fetch_base(chrom, anchor_pos)?;
    let deleted = fasta.fetch_sequence(chrom, gstart, gend + 1)?;

    let mut ref_allele = Vec::with_capacity(1 + deleted.len());
    ref_allele.push(anchor_base);
    ref_allele.extend_from_slice(&deleted);

    Ok(GenomicVariant {
        chrom: chrom.clone(),
        pos: anchor_pos,
        ref_allele,
        alt_allele: vec![anchor_base],
    })
}

/// Build a VCF-style [`GenomicVariant`] for a duplication.
///
/// A duplication is represented in VCF as an insertion anchored at the last
/// base of the duplicated range (`gend`). The anchor convention is
/// strand-agnostic: VCF always places the insertion after the anchor base
/// on the plus strand.
fn build_duplication(
    parsed: &ParsedHgvsC,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let (gstart, gend) = resolve_range(&parsed.start, parsed.end.as_ref(), transcript, index)?;
    let chrom = &transcript.chrom;

    // A duplication is an insertion immediately after the duplicated range.
    // VCF: anchor base = last duplicated base.
    let dup_seq = fasta.fetch_sequence(chrom, gstart, gend + 1)?;
    let anchor_base = fasta.fetch_base(chrom, gend)?;

    let mut alt_allele = Vec::with_capacity(1 + dup_seq.len());
    alt_allele.push(anchor_base);
    alt_allele.extend_from_slice(&dup_seq);

    Ok(GenomicVariant {
        chrom: chrom.clone(),
        pos: gend,
        ref_allele: vec![anchor_base],
        alt_allele,
    })
}

/// Build a VCF-style [`GenomicVariant`] for an insertion.
fn build_insertion(
    parsed: &ParsedHgvsC,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let HgvsCChange::Insertion { bases } = &parsed.change else {
        unreachable!("build_insertion is only called for Insertion variants")
    };

    let end_pos = parsed
        .end
        .as_ref()
        .expect("insertion must have end position");
    let g_left = resolve_position(&parsed.start, transcript, index)?;
    let g_right = resolve_position(end_pos, transcript, index)?;

    let chrom = &transcript.chrom;
    // VCF anchor = left flanking base (min genomic).
    let anchor_pos = g_left.min(g_right);
    let anchor_base = fasta.fetch_base(chrom, anchor_pos)?;

    let plus_inserted = to_plus_strand_seq(bases, transcript.strand);

    let mut alt_allele = Vec::with_capacity(1 + plus_inserted.len());
    alt_allele.push(anchor_base);
    alt_allele.extend_from_slice(&plus_inserted);

    Ok(GenomicVariant {
        chrom: chrom.clone(),
        pos: anchor_pos,
        ref_allele: vec![anchor_base],
        alt_allele,
    })
}

/// Build a VCF-style [`GenomicVariant`] for a delins.
fn build_delins(
    parsed: &ParsedHgvsC,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let HgvsCChange::Delins { bases } = &parsed.change else {
        unreachable!("build_delins is only called for Delins variants")
    };

    let (gstart, gend) = resolve_range(&parsed.start, parsed.end.as_ref(), transcript, index)?;
    let chrom = &transcript.chrom;

    let ref_allele = fasta.fetch_sequence(chrom, gstart, gend + 1)?;
    let alt_allele = to_plus_strand_seq(bases, transcript.strand);

    Ok(GenomicVariant {
        chrom: chrom.clone(),
        pos: gstart,
        ref_allele,
        alt_allele,
    })
}

// ---------------------------------------------------------------------------
// Transcript lookup
// ---------------------------------------------------------------------------

/// Look up a transcript by accession, with version-tolerant fallback.
///
/// 1. Exact match (e.g. `"NM_000546.6"`) — return immediately.
/// 2. Otherwise, derive the base accession (everything before the first
///    `.`, or the whole string if no dot) and scan the store for entries
///    whose accession starts with `"{base}."` and whose suffix parses as
///    a `u32`. Return the pair with the highest version.
/// 3. If still not found: `TranscriptNotFound`.
///
/// The fallback runs for **both** unversioned inputs (e.g. `"NM_000546"`)
/// and dotted inputs that miss the exact index (e.g. `"NM_000546.3"` when
/// the store only carries `.6`). The clinical trade-off — CDS coordinates
/// can shift between transcript versions when exon structure changes — is
/// the same one `lia-resolve-clients` already makes on the VEP REST
/// fallback path. Callers that need to surface drift to end users should
/// use [`VarEffect::resolve_hgvs_c_with_meta`], which returns the
/// resolved accession alongside the coordinates.
fn lookup_transcript<'a>(
    accession: &str,
    store: &'a TranscriptStore,
) -> Result<(&'a TranscriptModel, &'a LocateIndex), VarEffectError> {
    // Fast path: exact versioned lookup.
    if let Some(pair) = store.get_by_accession(accession) {
        return Ok(pair);
    }

    // Version-tolerant fallback: RefSeq/Ensembl accessions all use a
    // single `.` as the version separator, so splitting at the first dot
    // yields the unique base accession for the versioned input case.
    let base = accession.split_once('.').map_or(accession, |(b, _)| b);
    let prefix = format!("{base}.");
    let mut best: Option<(&TranscriptModel, u32)> = None;

    for tx in store.transcripts() {
        if let Some(suffix) = tx.accession.strip_prefix(&prefix)
            && let Ok(ver) = suffix.parse::<u32>()
            && best.is_none_or(|(_, v)| ver > v)
        {
            best = Some((tx, ver));
        }
    }

    if let Some((tx, _)) = best
        && let Some(pair) = store.get_by_accession(&tx.accession)
    {
        return Ok(pair);
    }

    Err(VarEffectError::TranscriptNotFound {
        accession: accession.to_string(),
    })
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Resolve an HGVS c. notation string to VCF-style genomic coordinates.
///
/// # Arguments
///
/// * `hgvs` -- Full HGVS c. notation (e.g., `"NM_000546.6:c.742C>T"`)
/// * `store` -- [`TranscriptStore`] for transcript lookup
/// * `fasta` -- [`FastaReader`] for reference allele verification and anchor
///   bases
///
/// # Returns
///
/// A [`GenomicVariant`] with 0-based position and plus-strand alleles,
/// ready for [`crate::VarEffect::annotate`].
///
/// # Errors
///
/// * [`VarEffectError::HgvsParseError`] -- input string cannot be parsed
/// * [`VarEffectError::TranscriptNotFound`] -- accession not in store
/// * [`VarEffectError::PositionOutOfRange`] -- c. position exceeds transcript
/// * [`VarEffectError::Malformed`] -- transcript model inconsistency
/// * [`VarEffectError::HgvsRefMismatch`] -- stated REF doesn't match FASTA
///   (includes the HGVS input string in its payload so callers can emit a
///   VEP-concordant diagnostic)
/// * [`VarEffectError::ChromNotFound`] / [`VarEffectError::CoordinateOutOfRange`]
///   -- FASTA access failures
pub(crate) fn resolve_hgvs_c_with_meta(
    hgvs: &str,
    store: &TranscriptStore,
    fasta: &FastaReader,
) -> Result<ResolvedHgvsC, VarEffectError> {
    let parsed = parse_hgvs_c(hgvs)?;
    let (transcript, index) = lookup_transcript(&parsed.accession, store)?;

    // c. notation requires a protein-coding transcript.
    if !matches!(transcript.biotype, Biotype::ProteinCoding) {
        return Err(VarEffectError::HgvsParseError(format!(
            "c. notation requires a protein-coding transcript, but {} is {:?}",
            transcript.accession, transcript.biotype,
        )));
    }

    // Snapshot the accession before the match — once we move into the
    // builders, the borrow on `transcript` is released.
    let resolved_accession = transcript.accession.clone();
    let variant = match &parsed.change {
        HgvsCChange::Substitution { .. } => {
            build_substitution(&parsed, transcript, index, fasta, hgvs)
        }
        HgvsCChange::Deletion => build_deletion(&parsed, transcript, index, fasta),
        HgvsCChange::Duplication => build_duplication(&parsed, transcript, index, fasta),
        HgvsCChange::Insertion { .. } => build_insertion(&parsed, transcript, index, fasta),
        HgvsCChange::Delins { .. } => build_delins(&parsed, transcript, index, fasta),
    }?;
    Ok(ResolvedHgvsC {
        variant,
        resolved_accession,
    })
}

/// Backward-compatible wrapper returning only the coordinates. Callers
/// that need transcript-version provenance (to detect and surface drift
/// from the requested HGVS version) should use [`resolve_hgvs_c_with_meta`].
pub(crate) fn resolve_hgvs_c(
    hgvs: &str,
    store: &TranscriptStore,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    resolve_hgvs_c_with_meta(hgvs, store, fasta).map(|r| r.variant)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::write_genome_binary;
    use crate::test_fixtures::{minus_strand_coding, plus_strand_coding};
    use tempfile::TempDir;

    // -- Parser tests (no I/O) -----------------------------------------------

    #[test]
    fn parse_cds_substitution() {
        let p = parse_hgvs_c("NM_000546.6:c.742C>T").unwrap();
        assert_eq!(p.accession, "NM_000546.6");
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: 742,
                is_3utr: false,
                intronic_offset: None
            }
        );
        assert!(p.end.is_none());
        assert_eq!(
            p.change,
            HgvsCChange::Substitution {
                ref_base: b'C',
                alt_base: b'T'
            }
        );
    }

    #[test]
    fn parse_5utr_substitution() {
        let p = parse_hgvs_c("NM_000546.6:c.-14G>A").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: -14,
                is_3utr: false,
                intronic_offset: None
            }
        );
        assert_eq!(
            p.change,
            HgvsCChange::Substitution {
                ref_base: b'G',
                alt_base: b'A'
            }
        );
    }

    #[test]
    fn parse_3utr_substitution() {
        let p = parse_hgvs_c("NM_000546.6:c.*32G>A").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: 32,
                is_3utr: true,
                intronic_offset: None
            }
        );
    }

    #[test]
    fn parse_intronic_donor() {
        let p = parse_hgvs_c("NM_000546.6:c.672+1G>A").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: 672,
                is_3utr: false,
                intronic_offset: Some(1)
            }
        );
    }

    #[test]
    fn parse_intronic_acceptor() {
        let p = parse_hgvs_c("NM_000546.6:c.673-2A>G").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: 673,
                is_3utr: false,
                intronic_offset: Some(-2)
            }
        );
    }

    #[test]
    fn parse_intronic_5utr() {
        let p = parse_hgvs_c("NM_XXX.1:c.-15+1G>A").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: -15,
                is_3utr: false,
                intronic_offset: Some(1)
            }
        );
    }

    #[test]
    fn parse_intronic_3utr() {
        let p = parse_hgvs_c("NM_XXX.1:c.*37+1G>A").unwrap();
        assert_eq!(
            p.start,
            HgvsCPosition {
                base: 37,
                is_3utr: true,
                intronic_offset: Some(1)
            }
        );
    }

    #[test]
    fn parse_single_del() {
        let p = parse_hgvs_c("NM_000546.6:c.19del").unwrap();
        assert_eq!(p.start.base, 19);
        assert!(p.end.is_none());
        assert_eq!(p.change, HgvsCChange::Deletion);
    }

    #[test]
    fn parse_range_del() {
        let p = parse_hgvs_c("NM_000546.6:c.19_21del").unwrap();
        assert_eq!(p.start.base, 19);
        assert_eq!(p.end.as_ref().unwrap().base, 21);
        assert_eq!(p.change, HgvsCChange::Deletion);
    }

    #[test]
    fn parse_single_dup() {
        let p = parse_hgvs_c("NM_000546.6:c.20dup").unwrap();
        assert_eq!(p.start.base, 20);
        assert!(p.end.is_none());
        assert_eq!(p.change, HgvsCChange::Duplication);
    }

    #[test]
    fn parse_range_dup() {
        let p = parse_hgvs_c("NM_000546.6:c.20_23dup").unwrap();
        assert_eq!(p.start.base, 20);
        assert_eq!(p.end.as_ref().unwrap().base, 23);
    }

    #[test]
    fn parse_insertion() {
        let p = parse_hgvs_c("NM_000546.6:c.76_77insT").unwrap();
        assert_eq!(p.start.base, 76);
        assert_eq!(p.end.as_ref().unwrap().base, 77);
        assert_eq!(p.change, HgvsCChange::Insertion { bases: vec![b'T'] });
    }

    #[test]
    fn parse_delins_range() {
        let p = parse_hgvs_c("NM_000546.6:c.112_117delinsTG").unwrap();
        assert_eq!(p.start.base, 112);
        assert_eq!(p.end.as_ref().unwrap().base, 117);
        assert_eq!(
            p.change,
            HgvsCChange::Delins {
                bases: vec![b'T', b'G']
            }
        );
    }

    #[test]
    fn parse_delins_single() {
        let p = parse_hgvs_c("NM_000546.6:c.113delinsTACT").unwrap();
        assert_eq!(p.start.base, 113);
        assert!(p.end.is_none());
        assert_eq!(
            p.change,
            HgvsCChange::Delins {
                bases: vec![b'T', b'A', b'C', b'T']
            }
        );
    }

    #[test]
    fn parse_invalid_missing_c_dot() {
        assert!(parse_hgvs_c("NM_000546.6:742C>T").is_err());
    }

    #[test]
    fn parse_invalid_c_zero() {
        assert!(parse_hgvs_c("NM_000546.6:c.0A>G").is_err());
    }

    #[test]
    fn parse_invalid_empty_accession() {
        assert!(parse_hgvs_c(":c.742C>T").is_err());
    }

    #[test]
    fn parse_invalid_bad_base() {
        assert!(parse_hgvs_c("NM_000546.6:c.742X>T").is_err());
    }

    #[test]
    fn parse_invalid_insertion_single_pos() {
        // Insertion requires two flanking positions.
        assert!(parse_hgvs_c("NM_000546.6:c.76insT").is_err());
    }

    // -- Position resolver tests (synthetic transcript) ----------------------

    // Helper: build a TranscriptStore from a single transcript.
    fn single_tx_store(tx: TranscriptModel) -> TranscriptStore {
        TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx])
    }

    // plus_strand_coding() geometry:
    //   exon 0: [1000,2000)  5'UTR [1000,1500) + CDS [1500,2000)  500bp CDS
    //   exon 1: [3000,3500)  CDS only                              500bp CDS
    //   exon 2: [4000,5000)  CDS [4000,4500) + 3'UTR [4500,5000)  500bp CDS
    //   cumulative_cds = [0, 500, 1000, 1500]
    //   Intron 0: [2000,3000), Intron 1: [3500,4000)

    #[test]
    fn resolve_cds_first_base_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.1 -> first CDS base = cds_genomic_start = 1500
        let pos = HgvsCPosition {
            base: 1,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 1500);
    }

    #[test]
    fn resolve_cds_first_base_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_MINUS.1").unwrap();
        // c.1 on minus -> cds_genomic_end - 1 = 19500 - 1 = 19499
        let pos = HgvsCPosition {
            base: 1,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 19499);
    }

    #[test]
    fn resolve_cds_cross_segment_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.501 -> first base of CDS segment 1 = 3000
        let pos = HgvsCPosition {
            base: 501,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 3000);
    }

    #[test]
    fn resolve_cds_last_base_of_segment_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.500 -> last base of CDS segment 0 = 1999
        let pos = HgvsCPosition {
            base: 500,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 1999);
    }

    #[test]
    fn resolve_5utr_same_exon_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.-1 -> one base before CDS start = 1500 - 1 = 1499
        let pos = HgvsCPosition {
            base: -1,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 1499);
    }

    #[test]
    fn resolve_5utr_same_exon_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_MINUS.1").unwrap();
        // Minus strand: CDS "start" in transcript = cds_genomic_end = 19500.
        // c.-1 = one base before CDS start = 19500 + 1 - 1 = 19500
        let pos = HgvsCPosition {
            base: -1,
            is_3utr: false,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 19500);
    }

    #[test]
    fn resolve_3utr_same_exon_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.*1 -> first base after CDS end = 4500
        let pos = HgvsCPosition {
            base: 1,
            is_3utr: true,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 4500);
    }

    #[test]
    fn resolve_3utr_same_exon_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_MINUS.1").unwrap();
        // Minus strand: CDS "end" in transcript = cds_genomic_start = 11000.
        // c.*1 = one base after CDS end = 11000 - 1 = 10999
        let pos = HgvsCPosition {
            base: 1,
            is_3utr: true,
            intronic_offset: None,
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 10999);
    }

    #[test]
    fn resolve_intronic_donor_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.500+5: anchor c.500 -> 1999, plus 5 -> 2004
        let pos = HgvsCPosition {
            base: 500,
            is_3utr: false,
            intronic_offset: Some(5),
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 2004);
    }

    #[test]
    fn resolve_intronic_acceptor_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // c.501-4: anchor c.501 -> 3000, minus 4 -> 2996
        let pos = HgvsCPosition {
            base: 501,
            is_3utr: false,
            intronic_offset: Some(-4),
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 2996);
    }

    #[test]
    fn resolve_intronic_donor_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_MINUS.1").unwrap();
        // c.1500+1: anchor c.1500 -> 18000, minus 1 -> 17999
        let pos = HgvsCPosition {
            base: 1500,
            is_3utr: false,
            intronic_offset: Some(1),
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 17999);
    }

    #[test]
    fn resolve_intronic_acceptor_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_MINUS.1").unwrap();
        // c.1501-2: anchor c.1501 -> 15999, minus (-2) -> 16001
        let pos = HgvsCPosition {
            base: 1501,
            is_3utr: false,
            intronic_offset: Some(-2),
        };
        assert_eq!(resolve_position(&pos, tx, idx).unwrap(), 16001);
    }

    #[test]
    fn resolve_out_of_range() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (tx, idx) = store.get_by_accession("NM_TEST_PLUS.1").unwrap();
        // CDS is 1500 bp, c.9999 is out of range.
        let pos = HgvsCPosition {
            base: 9999,
            is_3utr: false,
            intronic_offset: None,
        };
        assert!(matches!(
            resolve_position(&pos, tx, idx),
            Err(VarEffectError::PositionOutOfRange { .. })
        ));
    }

    // -- VCF construction tests (need synthetic FASTA) -----------------------

    /// Build a synthetic FASTA covering the test fixture chromosomes.
    ///
    /// chr1:  [0, 6000)  — covers plus_strand_coding (tx [1000,5000))
    /// chr17: [0, 21000) — covers minus_strand_coding (tx [10000,20000))
    fn build_test_fasta() -> (TempDir, FastaReader) {
        // Fill with a repeating ACGT pattern so we have deterministic ref bases.
        let chr1 = build_repeating_seq(6000);
        let chr17 = build_repeating_seq(21000);

        let contigs: Vec<(&str, &[u8])> = vec![("chr1", &chr1), ("chr17", &chr17)];
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("test.bin");
        let idx_path = tmp.path().join("test.bin.idx");
        write_genome_binary(&contigs, "test", &bin_path, &idx_path).unwrap();
        let reader = FastaReader::open_with_assembly(&bin_path, crate::Assembly::GRCh38).unwrap();
        (tmp, reader)
    }

    /// Build a repeating ACGTACGT... sequence of the given length.
    fn build_repeating_seq(len: usize) -> Vec<u8> {
        let pattern = b"ACGT";
        (0..len).map(|i| pattern[i % 4]).collect()
    }

    #[test]
    fn vcf_substitution_plus() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1 -> genomic 1500. FASTA at 1500 with ACGT pattern: 1500 % 4 = 0 -> 'A'
        // But the HGVS says REF is C. The FASTA has A at 1500, so this should
        // be an HgvsRefMismatch (HGVS-input path, distinct from VCF RefMismatch).
        let hgvs = "NM_TEST_PLUS.1:c.1C>T";
        let err = resolve_hgvs_c(hgvs, &store, &fasta).expect_err("expected mismatch");
        let VarEffectError::HgvsRefMismatch {
            hgvs: err_hgvs,
            chrom,
            pos,
            expected,
            got,
        } = &err
        else {
            panic!("expected HgvsRefMismatch, got: {err:?}");
        };
        assert_eq!(err_hgvs, hgvs);
        assert_eq!(chrom, "chr1");
        assert_eq!(*pos, 1500);
        assert_eq!(expected, "A");
        assert_eq!(got, "C");
        // Display matches VEP's wording so downstream callers can forward
        // it verbatim.
        assert_eq!(
            format!("{err}"),
            "ref allele mismatch at position 1500 for 'NM_TEST_PLUS.1:c.1C>T': \
             genome has 'A', HGVS states 'C'",
        );

        // Now use the correct ref base. 1500 % 4 = 0 -> A.
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.1A>T", &store, &fasta).unwrap();
        assert_eq!(result.chrom, "chr1");
        assert_eq!(result.pos, 1500);
        assert_eq!(result.ref_allele, vec![b'A']);
        assert_eq!(result.alt_allele, vec![b'T']);
    }

    #[test]
    fn hgvs_ref_mismatch_minus_strand_projects_to_plus() {
        // Minus-strand transcript: caller types a coding-strand base, but
        // the error's `got` field reports the plus-strand projection so it
        // compares apples-to-apples with `expected` (also plus-strand).
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1 on minus -> genomic 19499. FASTA: 19499 % 4 = 3 -> 'T' (plus).
        // Caller types c.1C>G (coding-strand REF=C). Plus-strand projection
        // of coding C is G, so `got = "G"`. FASTA has T, so mismatch.
        let hgvs = "NM_TEST_MINUS.1:c.1C>G";
        let err = resolve_hgvs_c(hgvs, &store, &fasta).expect_err("expected mismatch");
        let VarEffectError::HgvsRefMismatch {
            hgvs: err_hgvs,
            chrom,
            pos,
            expected,
            got,
        } = &err
        else {
            panic!("expected HgvsRefMismatch, got: {err:?}");
        };
        assert_eq!(err_hgvs, hgvs);
        assert_eq!(chrom, "chr17");
        assert_eq!(*pos, 19499);
        assert_eq!(expected, "T"); // plus-strand FASTA
        assert_eq!(got, "G"); // plus-strand projection of coding-C on minus
    }

    #[test]
    fn vcf_substitution_minus() {
        let tx = minus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1 on minus -> genomic 19499. FASTA: 19499 % 4 = 3 -> 'T'.
        // Minus strand: coding-strand REF = complement(T) = A.
        let result = resolve_hgvs_c("NM_TEST_MINUS.1:c.1A>G", &store, &fasta).unwrap();
        assert_eq!(result.chrom, "chr17");
        assert_eq!(result.pos, 19499);
        assert_eq!(result.ref_allele, vec![b'T']); // plus-strand
        assert_eq!(result.alt_allele, vec![b'C']); // complement of G
    }

    #[test]
    fn vcf_deletion_anchor_base() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.2del -> genomic 1501. Anchor at 1500.
        // FASTA: 1500=A, 1501=C
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.2del", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1500); // anchor
        assert_eq!(result.ref_allele, vec![b'A', b'C']); // anchor + deleted
        assert_eq!(result.alt_allele, vec![b'A']); // anchor only
    }

    #[test]
    fn vcf_range_deletion() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.2_4del -> genomic [1501, 1503]. Anchor at 1500.
        // FASTA: 1500=A, 1501=C, 1502=G, 1503=T
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.2_4del", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1500);
        assert_eq!(result.ref_allele, vec![b'A', b'C', b'G', b'T']);
        assert_eq!(result.alt_allele, vec![b'A']);
    }

    #[test]
    fn vcf_insertion_anchor_base() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1_2insG -> insertion between 1500 and 1501. Anchor at 1500.
        // FASTA: 1500=A
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.1_2insG", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1500);
        assert_eq!(result.ref_allele, vec![b'A']);
        assert_eq!(result.alt_allele, vec![b'A', b'G']);
    }

    #[test]
    fn vcf_duplication() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1dup -> dup base at 1500. FASTA: 1500=A.
        // VCF: pos=1500, ref=A, alt=AA
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.1dup", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1500);
        assert_eq!(result.ref_allele, vec![b'A']);
        assert_eq!(result.alt_allele, vec![b'A', b'A']);
    }

    #[test]
    fn vcf_range_duplication() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.1_3dup -> dup bases [1500,1502]. FASTA: A,C,G.
        // VCF: pos=1502 (last dup base), ref=G, alt=G+ACG
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.1_3dup", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1502);
        assert_eq!(result.ref_allele, vec![b'G']);
        assert_eq!(result.alt_allele, vec![b'G', b'A', b'C', b'G']);
    }

    #[test]
    fn vcf_delins() {
        let tx = plus_strand_coding();
        let store = single_tx_store(tx);
        let (_tmp, fasta) = build_test_fasta();

        // c.2delinsTT -> replace base at 1501 with TT.
        // FASTA: 1501=C.
        let result = resolve_hgvs_c("NM_TEST_PLUS.1:c.2delinsTT", &store, &fasta).unwrap();
        assert_eq!(result.pos, 1501);
        assert_eq!(result.ref_allele, vec![b'C']);
        assert_eq!(result.alt_allele, vec![b'T', b'T']);
    }

    // -- Transcript lookup tests ---------------------------------------------

    #[test]
    fn lookup_versioned() {
        let store = single_tx_store(plus_strand_coding());
        let (tx, _) = lookup_transcript("NM_TEST_PLUS.1", &store).unwrap();
        assert_eq!(tx.accession, "NM_TEST_PLUS.1");
    }

    #[test]
    fn lookup_not_found() {
        let store = single_tx_store(plus_strand_coding());
        assert!(matches!(
            lookup_transcript("NM_NONEXISTENT.1", &store),
            Err(VarEffectError::TranscriptNotFound { .. })
        ));
    }

    #[test]
    fn lookup_unversioned() {
        let store = single_tx_store(plus_strand_coding());
        let (tx, _) = lookup_transcript("NM_TEST_PLUS", &store).unwrap();
        assert_eq!(tx.accession, "NM_TEST_PLUS.1");
    }

    #[test]
    fn lookup_versioned_but_missing_falls_back_to_available() {
        // Store holds only `.1`; user cites `.2` (older) — fallback returns `.1`.
        // This mirrors the real-world case of NM_006772.2:c.1861_1862del
        // reported against a store carrying only NM_006772.3.
        let store = single_tx_store(plus_strand_coding());
        let (tx, _) = lookup_transcript("NM_TEST_PLUS.2", &store).unwrap();
        assert_eq!(tx.accession, "NM_TEST_PLUS.1");
    }

    #[test]
    fn lookup_versioned_missing_picks_highest_available() {
        // Two versions in the store; the fallback must pick the highest.
        let mut v1 = plus_strand_coding();
        v1.accession = "NM_TEST_PLUS.1".into();
        let mut v5 = plus_strand_coding();
        v5.accession = "NM_TEST_PLUS.5".into();
        let store = TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![v1, v5]);

        // Requested `.99` (absent) -> resolver returns `.5` (highest in store).
        let (tx, _) = lookup_transcript("NM_TEST_PLUS.99", &store).unwrap();
        assert_eq!(tx.accession, "NM_TEST_PLUS.5");
    }

    #[test]
    fn lookup_versioned_missing_unknown_base_errors() {
        // Store has only `NM_TEST_PLUS.1`. A request with a different base
        // accession must still fail — the fallback must not cross base-
        // accession boundaries.
        let store = single_tx_store(plus_strand_coding());
        assert!(matches!(
            lookup_transcript("NM_DIFFERENT.1", &store),
            Err(VarEffectError::TranscriptNotFound { .. })
        ));
    }

    #[test]
    fn resolve_hgvs_c_tolerates_missing_version() {
        // Regression for NM_006772.2:c.1861_1862del — user cites a version
        // absent from the store; resolution must succeed against the
        // available version with identical genomic coordinates.
        let store = single_tx_store(plus_strand_coding());
        let (_tmp, fasta) = build_test_fasta();

        let drifted = resolve_hgvs_c("NM_TEST_PLUS.9:c.1A>T", &store, &fasta).unwrap();
        let exact = resolve_hgvs_c("NM_TEST_PLUS.1:c.1A>T", &store, &fasta).unwrap();

        assert_eq!(drifted.chrom, exact.chrom);
        assert_eq!(drifted.pos, exact.pos);
        assert_eq!(drifted.ref_allele, exact.ref_allele);
        assert_eq!(drifted.alt_allele, exact.alt_allele);
    }

    #[test]
    fn resolve_hgvs_c_with_meta_reports_resolved_version() {
        // The meta variant must surface the accession actually used so
        // callers can detect drift. Exact-version hit reports the same
        // accession; version-drift hit reports the store's version.
        let store = single_tx_store(plus_strand_coding());
        let (_tmp, fasta) = build_test_fasta();

        let exact = resolve_hgvs_c_with_meta("NM_TEST_PLUS.1:c.1A>T", &store, &fasta).unwrap();
        assert_eq!(exact.resolved_accession, "NM_TEST_PLUS.1");

        let drifted = resolve_hgvs_c_with_meta("NM_TEST_PLUS.9:c.1A>T", &store, &fasta).unwrap();
        assert_eq!(drifted.resolved_accession, "NM_TEST_PLUS.1");
        // The variant payload must match the exact-version resolution.
        assert_eq!(drifted.variant, exact.variant);
    }

    // -- Integration tests (require real store + FASTA) ----------------------

    #[test]
    #[ignore]
    fn reverse_map_tp53_missense() {
        let store = load_store();
        let fasta = load_fasta();

        // TP53 c.742C>T (NM_000546.6) is a well-known pathogenic variant.
        // Minus-strand gene on chr17.
        let result = resolve_hgvs_c("NM_000546.6:c.742C>T", &store, &fasta).unwrap();
        assert_eq!(result.chrom, "chr17");
        // The genomic position should be known. Feed into annotate() to verify
        // consequence = missense_variant.
        assert!(!result.ref_allele.is_empty());
        assert!(!result.alt_allele.is_empty());

        // Round-trip: annotate the resolved variant and check consequence.
        let hits = store.query_overlap(&result.chrom, result.pos, result.pos + 1);
        let nm546 = hits
            .iter()
            .find(|(tx, _)| tx.accession == "NM_000546.6")
            .expect("NM_000546.6 should overlap");
        let csq = crate::annotate_snv(
            &result.chrom,
            result.pos,
            result.ref_allele[0],
            result.alt_allele[0],
            nm546.0,
            nm546.1,
            &fasta,
        )
        .unwrap();
        assert!(
            csq.consequences
                .iter()
                .any(|c| c.as_str() == "missense_variant"),
            "expected missense_variant, got {:?}",
            csq.consequences,
        );
    }

    #[test]
    #[ignore]
    fn reverse_map_tp53_splice() {
        let store = load_store();
        let fasta = load_fasta();

        // TP53 c.672+1G>A — canonical splice donor.
        let result = resolve_hgvs_c("NM_000546.6:c.672+1G>A", &store, &fasta).unwrap();
        assert_eq!(result.chrom, "chr17");

        let hits = store.query_overlap(&result.chrom, result.pos, result.pos + 1);
        let nm546 = hits
            .iter()
            .find(|(tx, _)| tx.accession == "NM_000546.6")
            .expect("NM_000546.6 should overlap");
        let csq = crate::annotate_snv(
            &result.chrom,
            result.pos,
            result.ref_allele[0],
            result.alt_allele[0],
            nm546.0,
            nm546.1,
            &fasta,
        )
        .unwrap();
        assert!(
            csq.consequences
                .iter()
                .any(|c| c.as_str() == "splice_donor_variant"),
            "expected splice_donor_variant, got {:?}",
            csq.consequences,
        );
    }

    // -- Integration test helpers --------------------------------------------

    #[cfg(test)]
    fn load_store() -> TranscriptStore {
        let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let workspace_root = manifest_dir
            .parent()
            .and_then(|p| p.parent())
            .expect("workspace root");
        let path = workspace_root.join("data/vareffect/transcript_models.bin");
        TranscriptStore::load_from_path(&path).unwrap_or_else(|e| {
            panic!(
                "failed to load transcript store from {}: {}",
                path.display(),
                e,
            )
        })
    }

    #[cfg(test)]
    fn load_fasta() -> FastaReader {
        let path = std::env::var("FASTA_PATH")
            .expect("FASTA_PATH env var must point to a GRCh38 genome binary");
        FastaReader::open_with_patch_aliases_and_assembly(
            std::path::Path::new(&path),
            Some(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .parent()
                    .and_then(|p| p.parent())
                    .expect("workspace root")
                    .join("data/vareffect/patch_chrom_aliases_grch38.csv")
                    .as_ref(),
            ),
            crate::Assembly::GRCh38,
        )
        .unwrap_or_else(|e| panic!("failed to open FASTA at {path}: {e}"))
    }
}
