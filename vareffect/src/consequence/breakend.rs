//! VCF 4.x breakend (`<BND>`) ALT-string parser.
//!
//! A breakend ALT encodes one end of a structural rearrangement (translocation,
//! large inversion, complex event) as a local reference base joined to a mate
//! locus, using the bracket grammar from the VCF specification:
//!
//! ```text
//!   t[p[   t]p]   ]p]t   [p[t          (t = local ref base(s), p = chrom:pos)
//!   t.     .t                          (single breakend: no mate locus)
//! ```
//!
//! Parsing yields the mate locus so a caller can annotate the mate breakpoint
//! via [`crate::VarEffect::annotate_breakend`]. The bracket orientation is
//! preserved for callers that need join direction (fusion strand); it does not
//! affect which transcript a breakpoint truncates.

use crate::error::VarEffectError;

/// Orientation of a parsed breakend, named by which side the local base joins
/// and which side of the mate locus the partner sequence extends from.
///
/// Only the mate locus is used for consequence annotation; this enum records
/// the join geometry for callers reconstructing fusion strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BreakendOrientation {
    /// `t[p[` — local base precedes the join; mate sequence extends to the
    /// right of `p` (forward strand).
    JoinAfterMateRight,
    /// `t]p]` — local base precedes the join; mate sequence extends to the
    /// left of `p` (reverse-complemented).
    JoinAfterMateLeft,
    /// `]p]t` — local base follows the join; mate sequence extends to the
    /// left of `p` (forward strand).
    JoinBeforeMateLeft,
    /// `[p[t` — local base follows the join; mate sequence extends to the
    /// right of `p` (reverse-complemented).
    JoinBeforeMateRight,
}

/// The mate locus of a breakend.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BreakendMate {
    /// Mate chromosome, verbatim from the ALT (e.g. `"12"`, `"chr12"`).
    pub chrom: String,
    /// Mate breakpoint position, **1-based**, verbatim from the ALT (VCF
    /// convention). Convert to 0-based (`pos - 1`) before passing to
    /// [`crate::VarEffect::annotate_breakend`].
    pub pos: u64,
    /// Join orientation.
    pub orientation: BreakendOrientation,
}

/// A parsed VCF breakend ALT.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Breakend {
    /// The local reference base(s) (`t` in the grammar), uppercase ASCII.
    pub ref_bases: String,
    /// The mate locus, or `None` for a single breakend (`t.` / `.t`).
    pub mate: Option<BreakendMate>,
}

impl Breakend {
    /// Parse a VCF breakend ALT string.
    ///
    /// Handles the four mated forms (`t[p[`, `t]p]`, `]p]t`, `[p[t`) and the
    /// two single-breakend forms (`t.`, `.t`).
    ///
    /// # Errors
    ///
    /// Returns [`VarEffectError::InvalidBreakend`] if `alt` does not match the
    /// breakend grammar, the mate position is not a valid integer, or the local
    /// base string is empty or non-`ACGTN`.
    pub fn parse(alt: &str) -> Result<Self, VarEffectError> {
        let invalid = || VarEffectError::InvalidBreakend(alt.to_string());

        // Single breakends carry no mate locus.
        if let Some(t) = alt.strip_suffix('.') {
            validate_ref_bases(t, alt)?;
            return Ok(Self {
                ref_bases: t.to_string(),
                mate: None,
            });
        }
        if let Some(t) = alt.strip_prefix('.') {
            validate_ref_bases(t, alt)?;
            return Ok(Self {
                ref_bases: t.to_string(),
                mate: None,
            });
        }

        let bytes = alt.as_bytes();
        let (first, last) = match (bytes.first(), bytes.last()) {
            (Some(&f), Some(&l)) => (f, l),
            _ => return Err(invalid()),
        };

        // Either the bracket pair leads (`]p]t` / `[p[t`, local base after) or
        // trails (`t[p[` / `t]p]`, local base before).
        let (bracket, local_after) = if first == b'[' || first == b']' {
            (first, true)
        } else if last == b'[' || last == b']' {
            (last, false)
        } else {
            return Err(invalid());
        };

        // Splitting on the bracket yields exactly three parts for a well-formed
        // breakend (the two brackets bound the mate locus).
        let parts: Vec<&str> = alt.split(bracket as char).collect();
        if parts.len() != 3 {
            return Err(invalid());
        }
        let (ref_bases, locus) = if local_after {
            (parts[2], parts[1])
        } else {
            (parts[0], parts[1])
        };
        validate_ref_bases(ref_bases, alt)?;

        // Mate locus is `chrom:pos`; split on the last ':' so contig names that
        // themselves contain a colon survive.
        let (chrom, pos_str) = locus.rsplit_once(':').ok_or_else(invalid)?;
        let pos: u64 = pos_str.parse().map_err(|_| invalid())?;
        if chrom.is_empty() {
            return Err(invalid());
        }

        let orientation = match (bracket, local_after) {
            (b'[', false) => BreakendOrientation::JoinAfterMateRight,
            (b']', false) => BreakendOrientation::JoinAfterMateLeft,
            (b']', true) => BreakendOrientation::JoinBeforeMateLeft,
            (b'[', true) => BreakendOrientation::JoinBeforeMateRight,
            _ => return Err(invalid()),
        };

        Ok(Self {
            ref_bases: ref_bases.to_string(),
            mate: Some(BreakendMate {
                chrom: chrom.to_string(),
                pos,
                orientation,
            }),
        })
    }
}

/// Validate that the local base string is non-empty and `ACGTN` (any case).
fn validate_ref_bases(t: &str, alt: &str) -> Result<(), VarEffectError> {
    if t.is_empty()
        || !t.bytes().all(|b| {
            matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
            )
        })
    {
        return Err(VarEffectError::InvalidBreakend(alt.to_string()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_t_bracket_right() {
        // t[p[ : local base precedes, mate extends right.
        let bnd = Breakend::parse("G[17:198982[").unwrap();
        assert_eq!(bnd.ref_bases, "G");
        let mate = bnd.mate.unwrap();
        assert_eq!(mate.chrom, "17");
        assert_eq!(mate.pos, 198_982);
        assert_eq!(mate.orientation, BreakendOrientation::JoinAfterMateRight);
    }

    #[test]
    fn parses_t_bracket_left() {
        let bnd = Breakend::parse("T]chr12:58877476]").unwrap();
        assert_eq!(bnd.ref_bases, "T");
        let mate = bnd.mate.unwrap();
        assert_eq!(mate.chrom, "chr12");
        assert_eq!(mate.pos, 58_877_476);
        assert_eq!(mate.orientation, BreakendOrientation::JoinAfterMateLeft);
    }

    #[test]
    fn parses_bracket_left_t() {
        // ]p]t : local base follows.
        let bnd = Breakend::parse("]13:123456]A").unwrap();
        assert_eq!(bnd.ref_bases, "A");
        let mate = bnd.mate.unwrap();
        assert_eq!(mate.chrom, "13");
        assert_eq!(mate.pos, 123_456);
        assert_eq!(mate.orientation, BreakendOrientation::JoinBeforeMateLeft);
    }

    #[test]
    fn parses_bracket_right_t() {
        let bnd = Breakend::parse("[2:321682[C").unwrap();
        assert_eq!(bnd.ref_bases, "C");
        let mate = bnd.mate.unwrap();
        assert_eq!(mate.chrom, "2");
        assert_eq!(mate.pos, 321_682);
        assert_eq!(mate.orientation, BreakendOrientation::JoinBeforeMateRight);
    }

    #[test]
    fn parses_multibase_local_sequence() {
        // Inserted sequence between the breakpoints is allowed in `t`.
        let bnd = Breakend::parse("GACGT[17:198982[").unwrap();
        assert_eq!(bnd.ref_bases, "GACGT");
        assert_eq!(bnd.mate.unwrap().pos, 198_982);
    }

    #[test]
    fn parses_single_breakend_suffix() {
        let bnd = Breakend::parse("G.").unwrap();
        assert_eq!(bnd.ref_bases, "G");
        assert!(bnd.mate.is_none());
    }

    #[test]
    fn parses_single_breakend_prefix() {
        let bnd = Breakend::parse(".G").unwrap();
        assert_eq!(bnd.ref_bases, "G");
        assert!(bnd.mate.is_none());
    }

    #[test]
    fn rejects_plain_allele() {
        assert!(matches!(
            Breakend::parse("ACGT"),
            Err(VarEffectError::InvalidBreakend(_))
        ));
    }

    #[test]
    fn rejects_missing_position() {
        assert!(matches!(
            Breakend::parse("G[17:[",),
            Err(VarEffectError::InvalidBreakend(_))
        ));
    }

    #[test]
    fn rejects_empty_local_base() {
        assert!(matches!(
            Breakend::parse("[2:321682["),
            Err(VarEffectError::InvalidBreakend(_))
        ));
    }

    #[test]
    fn rejects_non_acgtn_local_base() {
        assert!(matches!(
            Breakend::parse("X[2:321682["),
            Err(VarEffectError::InvalidBreakend(_))
        ));
    }
}
