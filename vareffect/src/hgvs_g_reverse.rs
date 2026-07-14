//! HGVS genomic (`g.`) reverse mapper: parse genomic-HGVS notation and resolve
//! to a left-aligned, VCF-canonical [`GenomicVariant`].
//!
//! This is the genomic sibling of [`crate::hgvs_reverse`] (which resolves `c.`
//! coding notation) and the reverse of the forward formatter
//! [`crate::hgvs_g::format_hgvs_g`]. Given a string like
//! `"NC_000017.11:g.7674220C>G"`, it maps the (already genomic, plus-strand,
//! 1-based) coordinates directly to a 0-based VCF-style `(chrom, pos, ref, alt)`.
//! Unlike the `c.` mapper there is no transcript, no strand handling, and no
//! CDS/UTR/intron math.
//!
//! # Supported variant types
//!
//! | Type | Example |
//! |------|---------|
//! | Substitution | `g.7674220C>G` |
//! | Deletion | `g.7674220del`, `g.7674220_7674225del` (optional stated `delACG`) |
//! | Duplication | `g.7674220dup`, `g.7674220_7674225dup` (optional stated `dupACG`) |
//! | Insertion | `g.7674220_7674221insATC` (flanks must be adjacent) |
//! | Deletion-insertion | `g.7674220delinsAT`, `g.7674220_7674225delinsAT` |
//! | Inversion | `g.7674220_7674230inv` (span > 1 nt) |
//!
//! # Normalization
//!
//! The output is **left-aligned to the 5'-most parsimonious form** (via
//! [`crate::left_align::left_align_indel`]), matching biocommons `hgvs_to_vcf`
//! (`shuffle_direction=5`), `bcftools norm`, gnomAD, and ClinVar's VCF row.
//! This is intentionally **different** from [`crate::VarEffect::resolve_hgvs_c`],
//! which emits the raw (un-normalized) form — the same repeat-region variant can
//! therefore resolve to different coordinates via the `c.` vs `g.` path.
//!
//! # Rejected input (fail-closed)
//!
//! Inputs that cannot be projected to exact chromosomal VCF coordinates are
//! rejected with [`VarEffectError::HgvsParseError`], never a guessed coordinate:
//! identity/no-change (`g.123=`), imprecise/templated notation
//! (`?`, `( )`, `[ ]`, `;`, `^`, `::`), and non-chromosomal reference sequences
//! (`NG_`, `LRG_`, `NW_`, `NT_` — their `g.` positions are region-relative).

use crate::chrom;
use crate::codon::reverse_complement;
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_reverse::{parse_seq, validate_base};
use crate::left_align::left_align_indel;
use crate::types::GenomicVariant;

// ---------------------------------------------------------------------------
// Parsed HGVS g. types
// ---------------------------------------------------------------------------

/// The change portion of an HGVS g. variant.
#[derive(Debug, Clone, PartialEq, Eq)]
enum HgvsGChange {
    /// Single-base substitution: `C>T`. Bases are uppercase plus-strand.
    Substitution { ref_base: u8, alt_base: u8 },
    /// Deletion. `stated` holds the explicit deleted bases when written
    /// (`delACG`), for validation against the genome; `None` for bare `del`.
    Deletion { stated: Option<Vec<u8>> },
    /// Duplication. `stated` holds the explicit duplicated bases when written.
    Duplication { stated: Option<Vec<u8>> },
    /// Insertion: `insACG` (inserted bases, uppercase).
    Insertion { bases: Vec<u8> },
    /// Deletion-insertion: `delinsTG` (inserted bases, uppercase).
    Delins { bases: Vec<u8> },
    /// Inversion: the reference span is replaced by its reverse-complement.
    Inversion,
}

/// A fully parsed HGVS g. notation string.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ParsedHgvsG {
    /// Genomic reference accession (`NC_*`) or UCSC name (`chr*`).
    accession: String,
    /// Start (or only) 1-based position.
    start: u64,
    /// End 1-based position for ranges. `None` for single-position variants.
    end: Option<u64>,
    /// The variant change type and associated data.
    change: HgvsGChange,
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse an HGVS g. (or `m.`) notation string into structured components.
///
/// Accepts the full notation including accession prefix, e.g.
/// `"NC_000017.11:g.7674220C>G"`. Mitochondrial `":m."` is accepted with the
/// same (linear, 1-based) numbering as `":g."`.
fn parse_hgvs_g(input: &str) -> Result<ParsedHgvsG, VarEffectError> {
    let err = |msg: &str| VarEffectError::HgvsParseError(format!("{msg}: \"{input}\""));

    // HGVS notation is pure ASCII. Reject non-ASCII up front so later byte-index
    // slicing (e.g. the `desc[..gt_idx - 1]` split around `>`) can never land
    // mid-multibyte-char and panic — a malformed input must be an error, not a
    // panic, per the crate's no-panic-on-bad-input policy.
    if !input.is_ascii() {
        return Err(err("HGVS notation must be ASCII"));
    }

    // Split accession from the g./m. description.
    let (accession, desc) = input
        .split_once(":g.")
        .or_else(|| input.split_once(":m."))
        .ok_or_else(|| err("missing ':g.' or ':m.' separator"))?;

    if accession.is_empty() {
        return Err(err("empty accession"));
    }

    // Non-chromosomal reference guard: NG_/LRG_ are numbered relative to the
    // gene-region/LRG sequence (not the chromosome), and NW_/NT_ are scaffolds
    // — none can be projected to chromosomal VCF coordinates here.
    if accession.starts_with("NG_")
        || accession.starts_with("LRG_")
        || accession.starts_with("NW_")
        || accession.starts_with("NT_")
    {
        return Err(err(
            "only chromosomal NC_/chr genomic references are supported; \
             gene-region/LRG/scaffold references use region-relative coordinates",
        ));
    }

    // Identity / no-change guard: `g.123=` / `g.123_456=` describe a reference
    // call, not a variant. `=` appears nowhere else in the g. grammar.
    if desc.contains('=') {
        return Err(err(
            "identity/reference-call (=) describes no sequence change, not a variant",
        ));
    }

    // Imprecision / templated / multi-allele guard: any of these tokens marks an
    // input that cannot be resolved to exact coordinates. The scan runs on
    // `desc` (post-`:g.`/`:m.` split) so the accession's single `:` cannot trip
    // it, and `::` is matched as a 2-char substring distinct from a single `:`.
    if desc.contains("::") || desc.contains(['?', '(', ')', '[', ']', ';', '^']) {
        return Err(err(
            "imprecise, templated, or multi-allele notation is not resolvable to exact coordinates",
        ));
    }

    // Identify the change type. Order matters: "delins" before "ins"/"del";
    // "inv" is checked early (it cannot collide with any keyword or sequence).
    if let Some(idx) = desc.find("delins") {
        let pos_str = &desc[..idx];
        let inserted = parse_seq(&desc[idx + 6..], input)?; // skip "delins"
        let (start, end) = parse_positions(pos_str, input)?;
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsGChange::Delins { bases: inserted },
        })
    } else if let Some(idx) = desc.find("inv") {
        let pos_str = &desc[..idx];
        if !desc[idx + 3..].is_empty() {
            return Err(err("unexpected trailing characters after 'inv'"));
        }
        let (start, end_opt) = parse_positions(pos_str, input)?;
        let end = end_opt.ok_or_else(|| err("inversion requires a range (start_end)"))?;
        if end <= start {
            return Err(err("inversion span must be > 1 nt"));
        }
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end: Some(end),
            change: HgvsGChange::Inversion,
        })
    } else if let Some(ins_idx) = desc.find("ins")
        && desc.contains('_')
    {
        let pos_str = &desc[..ins_idx];
        let inserted = parse_seq(&desc[ins_idx + 3..], input)?; // skip "ins"
        let (start, end_opt) = parse_positions(pos_str, input)?;
        let end = end_opt.ok_or_else(|| err("insertion requires two flanking positions"))?;
        if end != start + 1 {
            return Err(err("insertion flanks must be adjacent (end == start + 1)"));
        }
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end: Some(end),
            change: HgvsGChange::Insertion { bases: inserted },
        })
    } else if let Some(idx) = desc.find("del") {
        let pos_str = &desc[..idx];
        let stated = parse_optional_seq(&desc[idx + 3..], input)?; // skip "del"
        let (start, end) = parse_positions(pos_str, input)?;
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsGChange::Deletion { stated },
        })
    } else if let Some(idx) = desc.find("dup") {
        let pos_str = &desc[..idx];
        let stated = parse_optional_seq(&desc[idx + 3..], input)?; // skip "dup"
        let (start, end) = parse_positions(pos_str, input)?;
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end,
            change: HgvsGChange::Duplication { stated },
        })
    } else if let Some(gt_idx) = desc.find('>') {
        if gt_idx == 0 {
            return Err(err("substitution missing REF base"));
        }
        let ref_base = desc.as_bytes()[gt_idx - 1].to_ascii_uppercase();
        let pos_str = &desc[..gt_idx - 1];
        let alt_str = &desc[gt_idx + 1..];
        if alt_str.len() != 1 {
            return Err(err("substitution ALT must be a single base"));
        }
        let alt_base = alt_str.as_bytes()[0].to_ascii_uppercase();
        validate_base(ref_base, input)?;
        validate_base(alt_base, input)?;
        let start = parse_pos(pos_str, input)?;
        Ok(ParsedHgvsG {
            accession: accession.to_string(),
            start,
            end: None,
            change: HgvsGChange::Substitution { ref_base, alt_base },
        })
    } else {
        Err(err(
            "unrecognized change type (expected del/dup/ins/delins/inv or >)",
        ))
    }
}

/// Parse an optional stated sequence following `del`/`dup` (empty → `None`).
fn parse_optional_seq(s: &str, input: &str) -> Result<Option<Vec<u8>>, VarEffectError> {
    if s.is_empty() {
        Ok(None)
    } else {
        Ok(Some(parse_seq(s, input)?))
    }
}

/// Parse a genomic position range `"742"` or `"742_744"` into 1-based bounds.
///
/// Returns `(start, Some(end))` for ranges, `(start, None)` for single
/// positions. Rejects `end < start`.
fn parse_positions(s: &str, input: &str) -> Result<(u64, Option<u64>), VarEffectError> {
    if let Some((left, right)) = s.split_once('_') {
        let start = parse_pos(left, input)?;
        let end = parse_pos(right, input)?;
        if end < start {
            return Err(VarEffectError::HgvsParseError(format!(
                "range end {end} precedes start {start}: \"{input}\""
            )));
        }
        Ok((start, Some(end)))
    } else {
        Ok((parse_pos(s, input)?, None))
    }
}

/// Parse a single genomic position: a plain 1-based integer, no offset forms.
fn parse_pos(s: &str, input: &str) -> Result<u64, VarEffectError> {
    let err = |msg: &str| {
        VarEffectError::HgvsParseError(format!("{msg} in position \"{s}\": \"{input}\""))
    };
    if s.is_empty() {
        return Err(err("empty position"));
    }
    // Genomic positions are plain integers — reject the c.-only `+`/`-`/`*`
    // intronic/UTR offset forms explicitly (they have no genomic meaning).
    if !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(err(
            "genomic positions must be plain integers (no +/-/* offsets)",
        ));
    }
    let value: u64 = s.parse().map_err(|_| err("position integer overflow"))?;
    if value == 0 {
        return Err(err("position 0 is invalid (HGVS is 1-based)"));
    }
    Ok(value)
}

// ---------------------------------------------------------------------------
// Accession → chromosome
// ---------------------------------------------------------------------------

/// Map a genomic reference accession to a UCSC-style chromosome name.
///
/// `NC_*` accessions map via [`chrom::refseq_to_ucsc`] (an unknown `NC_*` passes
/// through unchanged and is treated as not-found); `chr*` names pass through;
/// anything else yields `None`.
fn accession_to_chrom(accession: &str) -> Option<String> {
    if accession.starts_with("NC_") {
        let ucsc = chrom::refseq_to_ucsc(accession);
        // `refseq_to_ucsc` echoes an unknown accession back unchanged.
        if ucsc == accession {
            None
        } else {
            Some(ucsc.to_string())
        }
    } else if accession.starts_with("chr") {
        Some(accession.to_string())
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// VCF coordinate construction (raw, pre-normalization forms)
// ---------------------------------------------------------------------------

/// Build the raw VCF-style substitution, verifying the stated REF.
fn build_substitution(
    pos: u64,
    ref_base: u8,
    alt_base: u8,
    chrom: &str,
    hgvs: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let gpos = pos - 1;
    let fasta_ref = fasta.fetch_base(chrom, gpos)?.to_ascii_uppercase();
    if fasta_ref != ref_base {
        return Err(VarEffectError::HgvsRefMismatch {
            hgvs: hgvs.to_string(),
            chrom: chrom.to_string(),
            pos: gpos,
            expected: String::from(fasta_ref as char),
            got: String::from(ref_base as char),
        });
    }
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: gpos,
        ref_allele: vec![ref_base],
        alt_allele: vec![alt_base],
    })
}

/// Verify a stated `del`/`dup` sequence against the genome span, if present.
fn verify_stated_span(
    stated: Option<&[u8]>,
    span: &[u8],
    chrom: &str,
    gstart0: u64,
    hgvs: &str,
) -> Result<(), VarEffectError> {
    if let Some(stated) = stated
        && stated != span
    {
        return Err(VarEffectError::HgvsRefMismatch {
            hgvs: hgvs.to_string(),
            chrom: chrom.to_string(),
            pos: gstart0,
            expected: String::from_utf8_lossy(span).into_owned(),
            got: String::from_utf8_lossy(stated).into_owned(),
        });
    }
    Ok(())
}

/// Build the raw VCF-style deletion (anchor base prepended).
fn build_deletion(
    start: u64,
    end: Option<u64>,
    stated: Option<&[u8]>,
    chrom: &str,
    hgvs: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let gstart0 = start - 1;
    let gend0 = end.unwrap_or(start) - 1;
    let deleted = fasta.fetch_sequence(chrom, gstart0, gend0 + 1)?;
    verify_stated_span(stated, &deleted, chrom, gstart0, hgvs)?;

    let anchor_pos =
        gstart0
            .checked_sub(1)
            .ok_or_else(|| VarEffectError::CoordinateOutOfRange {
                chrom: chrom.to_string(),
                start: 0,
                end: 0,
                chrom_len: fasta.chrom_length(chrom).unwrap_or(0),
            })?;
    let anchor_base = fasta.fetch_base(chrom, anchor_pos)?;

    let mut ref_allele = Vec::with_capacity(1 + deleted.len());
    ref_allele.push(anchor_base);
    ref_allele.extend_from_slice(&deleted);
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: anchor_pos,
        ref_allele,
        alt_allele: vec![anchor_base],
    })
}

/// Build the raw VCF-style duplication (an insertion after the last dup base).
fn build_duplication(
    start: u64,
    end: Option<u64>,
    stated: Option<&[u8]>,
    chrom: &str,
    hgvs: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let gstart0 = start - 1;
    let gend0 = end.unwrap_or(start) - 1;
    let dup_seq = fasta.fetch_sequence(chrom, gstart0, gend0 + 1)?;
    verify_stated_span(stated, &dup_seq, chrom, gstart0, hgvs)?;

    let anchor_base = fasta.fetch_base(chrom, gend0)?;
    let mut alt_allele = Vec::with_capacity(1 + dup_seq.len());
    alt_allele.push(anchor_base);
    alt_allele.extend_from_slice(&dup_seq);
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: gend0,
        ref_allele: vec![anchor_base],
        alt_allele,
    })
}

/// Build the raw VCF-style insertion (anchor = left flank base).
fn build_insertion(
    start: u64,
    inserted: &[u8],
    chrom: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let anchor_pos = start - 1; // left flanking base (1-based `start`)
    let anchor_base = fasta.fetch_base(chrom, anchor_pos)?;
    let mut alt_allele = Vec::with_capacity(1 + inserted.len());
    alt_allele.push(anchor_base);
    alt_allele.extend_from_slice(inserted);
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: anchor_pos,
        ref_allele: vec![anchor_base],
        alt_allele,
    })
}

/// Build the raw VCF-style deletion-insertion (no anchor; both alleles present).
fn build_delins(
    start: u64,
    end: Option<u64>,
    inserted: &[u8],
    chrom: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let gstart0 = start - 1;
    let gend0 = end.unwrap_or(start) - 1;
    let ref_allele = fasta.fetch_sequence(chrom, gstart0, gend0 + 1)?;
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: gstart0,
        ref_allele,
        alt_allele: inserted.to_vec(),
    })
}

/// Build the raw VCF-style inversion (`alt = reverse_complement(ref span)`).
fn build_inversion(
    start: u64,
    end: u64,
    chrom: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let gstart0 = start - 1;
    let gend0 = end - 1;
    let span = fasta.fetch_sequence(chrom, gstart0, gend0 + 1)?;
    let alt_allele = reverse_complement(&span);
    Ok(GenomicVariant {
        chrom: chrom.to_string(),
        pos: gstart0,
        ref_allele: span,
        alt_allele,
    })
}

/// Left-align the raw variant to the 5'-most parsimonious VCF-canonical form.
///
/// Pure indels (del/dup/ins) shift to the leftmost locus; sub/delins/inv are
/// usually unchanged, but `left_align_indel` also trims shared prefix/suffix
/// bases, so an affix-sharing delins/inv collapses to a smaller edit.
fn canonicalize(
    raw: GenomicVariant,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let shifted = {
        let ref_str =
            std::str::from_utf8(&raw.ref_allele).map_err(|_| VarEffectError::InvalidAllele)?;
        let alt_str =
            std::str::from_utf8(&raw.alt_allele).map_err(|_| VarEffectError::InvalidAllele)?;
        left_align_indel(fasta, &raw.chrom, raw.pos + 1, ref_str, alt_str)?
    };
    match shifted {
        Some((pos_1based, r, a)) => Ok(GenomicVariant {
            chrom: raw.chrom,
            pos: pos_1based - 1,
            ref_allele: r.into_bytes(),
            alt_allele: a.into_bytes(),
        }),
        None => Ok(raw),
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Resolve an HGVS g. notation string to a left-aligned VCF-canonical variant.
///
/// # Arguments
///
/// * `hgvs` — Full genomic HGVS notation (e.g. `"NC_000017.11:g.7674220C>G"`).
/// * `fasta` — [`FastaReader`] for reference-allele verification, anchor bases,
///   and left-alignment.
///
/// # Returns
///
/// A [`GenomicVariant`] with 0-based position, UCSC-style chromosome, and
/// uppercase plus-strand alleles, ready for [`crate::VarEffect::annotate`].
///
/// # Errors
///
/// * [`VarEffectError::HgvsParseError`] — unparseable, imprecise, identity, or
///   non-chromosomal-reference input.
/// * [`VarEffectError::ChromNotFound`] — accession maps to no known contig.
/// * [`VarEffectError::HgvsRefMismatch`] — stated REF / deleted / dup bases
///   disagree with the reference genome.
/// * [`VarEffectError::CoordinateOutOfRange`] — a deletion at the contig start
///   has no 5' anchor, or a span exceeds the contig.
pub(crate) fn resolve_hgvs_g(
    hgvs: &str,
    fasta: &FastaReader,
) -> Result<GenomicVariant, VarEffectError> {
    let parsed = parse_hgvs_g(hgvs)?;

    let chrom =
        accession_to_chrom(&parsed.accession).ok_or_else(|| VarEffectError::ChromNotFound {
            chrom: parsed.accession.clone(),
        })?;
    if fasta.chrom_length(&chrom).is_none() {
        return Err(VarEffectError::ChromNotFound { chrom });
    }

    let raw = match &parsed.change {
        HgvsGChange::Substitution { ref_base, alt_base } => {
            build_substitution(parsed.start, *ref_base, *alt_base, &chrom, hgvs, fasta)?
        }
        HgvsGChange::Deletion { stated } => build_deletion(
            parsed.start,
            parsed.end,
            stated.as_deref(),
            &chrom,
            hgvs,
            fasta,
        )?,
        HgvsGChange::Duplication { stated } => build_duplication(
            parsed.start,
            parsed.end,
            stated.as_deref(),
            &chrom,
            hgvs,
            fasta,
        )?,
        HgvsGChange::Insertion { bases } => build_insertion(parsed.start, bases, &chrom, fasta)?,
        HgvsGChange::Delins { bases } => {
            build_delins(parsed.start, parsed.end, bases, &chrom, fasta)?
        }
        HgvsGChange::Inversion => build_inversion(
            parsed.start,
            parsed.end.expect("inversion parse guarantees a range"),
            &chrom,
            fasta,
        )?,
    };

    // Null-variant guard: reject `ref == alt` (e.g. a palindromic inversion)
    // before left-alignment — feeding equal alleles to `left_align_indel` would
    // walk its shift loop toward position 1 (a degenerate result).
    if raw.ref_allele == raw.alt_allele {
        return Err(VarEffectError::HgvsParseError(format!(
            "no sequence change: \"{hgvs}\""
        )));
    }

    canonicalize(raw, fasta)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    //! CI-runnable exact-output tests over a tiny synthetic genome (no 3.1 GB
    //! FASTA). Coordinates are precise because we control the sequence.

    use super::*;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    /// Synthetic contig "1" (queried as "chr1" / accession `NC_000001.11`):
    ///   1-based pos: 1 2 3 4 5 6 7 8 9 10
    ///   base:        A A T G G G G T A A
    /// The `GGGG` run (pos 4-7) anchored by `T` at pos 3 exercises left-shifting;
    /// the `AT` at pos 2-3 is a reverse-complement palindrome.
    const CONTIG: &[u8] = b"AATGGGGTAA";

    fn fasta_with(contigs: &[(&str, &[u8])]) -> (TempDir, FastaReader) {
        let tmp = TempDir::new().expect("tempdir");
        let bin = tmp.path().join("g.bin");
        let idx = tmp.path().join("g.bin.idx");
        write_genome_binary(contigs, "test", &bin, &idx).expect("write synthetic genome");
        let fasta = FastaReader::open(&bin).expect("open synthetic genome");
        (tmp, fasta)
    }

    fn resolve(hgvs: &str) -> Result<GenomicVariant, VarEffectError> {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        resolve_hgvs_g(hgvs, &fasta)
    }

    fn gv(pos: u64, r: &[u8], a: &[u8]) -> GenomicVariant {
        GenomicVariant {
            chrom: "chr1".to_string(),
            pos,
            ref_allele: r.to_vec(),
            alt_allele: a.to_vec(),
        }
    }

    // -- Change types --------------------------------------------------------

    #[test]
    fn substitution() {
        assert_eq!(resolve("NC_000001.11:g.3T>A").unwrap(), gv(2, b"T", b"A"));
    }

    #[test]
    fn substitution_is_case_insensitive() {
        assert_eq!(resolve("NC_000001.11:g.3t>a").unwrap(), gv(2, b"T", b"A"));
    }

    #[test]
    fn single_base_deletion_left_aligns() {
        // Delete the 3'-most G of the GGGG run; left-aligns to the T anchor.
        assert_eq!(resolve("NC_000001.11:g.7del").unwrap(), gv(2, b"TG", b"T"));
    }

    #[test]
    fn range_deletion() {
        // Delete GG at pos 4-5; left-aligns within the run to the T anchor.
        assert_eq!(
            resolve("NC_000001.11:g.4_5del").unwrap(),
            gv(2, b"TGG", b"T")
        );
    }

    #[test]
    fn deletion_at_contig_start_errors() {
        // A deletion at 1-based position 1 has no 5' anchor base — the
        // `checked_sub(1)` guard fails closed instead of underflowing.
        assert!(matches!(
            resolve("NC_000001.11:g.1del"),
            Err(VarEffectError::CoordinateOutOfRange { .. })
        ));
    }

    #[test]
    fn stated_deletion_sequence_validates() {
        // g.7delG — stated matches the genome.
        assert_eq!(resolve("NC_000001.11:g.7delG").unwrap(), gv(2, b"TG", b"T"));
    }

    #[test]
    fn stated_deletion_sequence_mismatch_errors() {
        // g.7delA — genome has G at pos 7.
        assert!(matches!(
            resolve("NC_000001.11:g.7delA"),
            Err(VarEffectError::HgvsRefMismatch { .. })
        ));
    }

    #[test]
    fn duplication_left_aligns() {
        // Duplicate the 3'-most G; left-aligns to an insertion at the T anchor.
        assert_eq!(resolve("NC_000001.11:g.7dup").unwrap(), gv(2, b"T", b"TG"));
    }

    #[test]
    fn range_duplication_left_aligns() {
        // Duplicate the GG at pos 6-7; a 2 bp tandem copy that left-aligns to
        // the T anchor as an insertion of GG.
        assert_eq!(
            resolve("NC_000001.11:g.6_7dup").unwrap(),
            gv(2, b"T", b"TGG")
        );
    }

    #[test]
    fn stated_duplication_sequence_validates() {
        // g.7dupG — stated matches the genome (G at pos 7); left-aligns like g.7dup.
        assert_eq!(resolve("NC_000001.11:g.7dupG").unwrap(), gv(2, b"T", b"TG"));
    }

    #[test]
    fn stated_duplication_sequence_mismatch_errors() {
        // g.7dupA — genome has G at pos 7.
        assert!(matches!(
            resolve("NC_000001.11:g.7dupA"),
            Err(VarEffectError::HgvsRefMismatch { .. })
        ));
    }

    #[test]
    fn insertion() {
        // Insert C between the T (pos 8) and A (pos 9) — no repeat, no shift.
        assert_eq!(
            resolve("NC_000001.11:g.8_9insC").unwrap(),
            gv(7, b"T", b"TC")
        );
    }

    #[test]
    fn delins() {
        // TG (pos 3-4) -> CC; no shared affix, not shifted.
        assert_eq!(
            resolve("NC_000001.11:g.3_4delinsCC").unwrap(),
            gv(2, b"TG", b"CC")
        );
    }

    #[test]
    fn delins_trims_shared_prefix_to_substitution() {
        // TG -> TC shares the leading T; canonicalizes to the SNV G>C at pos 4.
        assert_eq!(
            resolve("NC_000001.11:g.3_4delinsTC").unwrap(),
            gv(3, b"G", b"C")
        );
    }

    #[test]
    fn delins_collapses_to_indel_and_left_aligns() {
        // GG (pos 4-5) -> G is really a deletion; collapses and left-aligns to
        // the T anchor.
        assert_eq!(
            resolve("NC_000001.11:g.4_5delinsG").unwrap(),
            gv(2, b"TG", b"T")
        );
    }

    #[test]
    fn inversion() {
        // Invert TGG (pos 3-5) -> reverse_complement = CCA.
        assert_eq!(
            resolve("NC_000001.11:g.3_5inv").unwrap(),
            gv(2, b"TGG", b"CCA")
        );
    }

    #[test]
    fn palindromic_inversion_is_rejected() {
        // AT (pos 2-3) is its own reverse-complement -> no sequence change.
        assert!(matches!(
            resolve("NC_000001.11:g.2_3inv"),
            Err(VarEffectError::HgvsParseError(_))
        ));
    }

    // -- Reference validation ------------------------------------------------

    #[test]
    fn substitution_ref_mismatch_errors() {
        // Genome has T at pos 3, not C.
        assert!(matches!(
            resolve("NC_000001.11:g.3C>A"),
            Err(VarEffectError::HgvsRefMismatch { .. })
        ));
    }

    // -- Rejections ----------------------------------------------------------

    fn is_parse_err(hgvs: &str) -> bool {
        matches!(resolve(hgvs), Err(VarEffectError::HgvsParseError(_)))
    }

    #[test]
    fn rejects_position_zero() {
        assert!(is_parse_err("NC_000001.11:g.0T>A"));
    }

    #[test]
    fn rejects_non_ascii_input() {
        // A multibyte char before '>' must fail closed, not panic on a
        // non-char-boundary slice.
        assert!(is_parse_err("NC_000001.11:g.3é>A"));
    }

    #[test]
    fn rejects_offset_forms() {
        assert!(is_parse_err("NC_000001.11:g.3+2C>A"));
        assert!(is_parse_err("NC_000001.11:g.3*1del"));
    }

    #[test]
    fn rejects_non_adjacent_insertion() {
        assert!(is_parse_err("NC_000001.11:g.3_5insA"));
    }

    #[test]
    fn rejects_single_position_inversion() {
        assert!(is_parse_err("NC_000001.11:g.3inv"));
    }

    #[test]
    fn rejects_reversed_range() {
        assert!(is_parse_err("NC_000001.11:g.5_3del"));
    }

    #[test]
    fn rejects_identity() {
        assert!(is_parse_err("NC_000001.11:g.3="));
        assert!(is_parse_err("NC_000001.11:g.3_4="));
    }

    #[test]
    fn rejects_imprecise_and_templated() {
        assert!(is_parse_err("NC_000001.11:g.?"));
        assert!(is_parse_err("NC_000001.11:g.(3_4)del"));
        assert!(is_parse_err("NC_000001.11:g.(3_4)_(6_7)del"));
        assert!(is_parse_err("NC_000001.11:g.3AC[4]"));
        assert!(is_parse_err("NC_000001.11:g.3_4insN[10]"));
        assert!(is_parse_err("NC_000001.11:g.[3T>A;5G>C]"));
    }

    #[test]
    fn rejects_non_chromosomal_references() {
        assert!(is_parse_err("NG_012232.1:g.3T>A"));
        assert!(is_parse_err("LRG_199:g.3T>A"));
        assert!(is_parse_err("NW_009646201.1:g.3T>A"));
        assert!(is_parse_err("NT_187633.1:g.3T>A"));
    }

    #[test]
    fn rejects_missing_separator() {
        assert!(is_parse_err("NM_000546.6:c.742C>T"));
    }

    // -- Accession mapping ---------------------------------------------------

    #[test]
    fn accepts_ucsc_chrom_accession() {
        assert_eq!(resolve("chr1:g.3T>A").unwrap(), gv(2, b"T", b"A"));
    }

    #[test]
    fn unknown_accession_is_chrom_not_found() {
        assert!(matches!(
            resolve("NC_000009.99:g.3T>A"),
            Err(VarEffectError::ChromNotFound { .. })
        ));
        assert!(matches!(
            resolve("chrZZ:g.3T>A"),
            Err(VarEffectError::ChromNotFound { .. })
        ));
    }

    #[test]
    fn accession_to_chrom_covers_edges() {
        assert_eq!(accession_to_chrom("NC_000017.11").as_deref(), Some("chr17"));
        assert_eq!(accession_to_chrom("NC_000024.10").as_deref(), Some("chrY"));
        assert_eq!(accession_to_chrom("NC_012920.1").as_deref(), Some("chrM"));
        assert_eq!(accession_to_chrom("chr1").as_deref(), Some("chr1"));
        assert_eq!(accession_to_chrom("NC_000009.99"), None);
        assert_eq!(accession_to_chrom("gibberish"), None);
    }

    #[test]
    fn mito_prefix_is_parsed() {
        // The `:m.` separator parses with the same numbering; accession maps to
        // chrM (resolution against a chrM contig is out of this synthetic genome).
        let parsed = parse_hgvs_g("NC_012920.1:m.8993T>G").unwrap();
        assert_eq!(parsed.accession, "NC_012920.1");
        assert_eq!(parsed.start, 8993);
        assert_eq!(
            accession_to_chrom(&parsed.accession).as_deref(),
            Some("chrM")
        );
    }

    // -- Round-trip against the forward formatter ----------------------------

    #[test]
    fn round_trips_through_format_hgvs_g() {
        // resolve (left-align) then format (3'-shift) reconciles for canonical
        // 3'-most inputs.
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        for input in [
            "NC_000001.11:g.3T>A",
            "NC_000001.11:g.7del",
            "NC_000001.11:g.7dup",
            "NC_000001.11:g.8_9insC",
        ] {
            let gvar = resolve_hgvs_g(input, &fasta).unwrap();
            let formatted = crate::hgvs_g::format_hgvs_g(
                &gvar.chrom,
                gvar.pos,
                &gvar.ref_allele,
                &gvar.alt_allele,
                &fasta,
            )
            .unwrap();
            assert_eq!(formatted.as_deref(), Some(input), "round-trip for {input}");
        }
    }
}
