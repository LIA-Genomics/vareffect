//! HGVS protein notation (`p.`) for variant consequences.
//!
//! Generates VEP-compatible HGVS protein notation for SNVs, indels, and
//! complex (delins / MNV) variants. Covers:
//!
//! SNV:
//! - Missense:       `p.Arg248Trp`
//! - Synonymous:     `p.Arg248=`
//! - Nonsense:       `p.Arg196Ter`
//! - Stop lost:      `p.Ter394CysextTer9` (computed via 3'UTR stop-scan)
//! - Start lost:     `p.Met1?`
//! - Start retained: `p.Met1=`
//! - Stop retained:  `p.Ter394=`
//!
//! Indel / complex:
//! - Frameshift:        `p.Glu23ValfsTer17`
//! - Inframe deletion:  `p.Phe508del`
//! - Inframe insertion: `p.Lys2_Met3insGlnSerLys`
//! - Duplication:       `p.Gly4_Gln6dup`
//! - Delins:            `p.Cys28_Lys29delinsTrp`
//! - Extension:         `p.Ter394CysextTer9`
//!
//! VEP does NOT use parentheses on predicted protein consequences (known
//! VEP issue #498). This module matches VEP's default behaviour.
//!
//! VEP uses the three-letter code `Ter` for stop codons in p. notation,
//! which matches `codon::aa_three_letter(b'*')`.

use crate::codon::{aa_three_letter, reverse_complement, translate_codon_for_transcript};
use crate::consequence::Consequence;
use crate::consequence::helpers::fetch_cds_sequence;
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::locate::LocateIndex;
use crate::types::{Strand, TranscriptModel};

/// Generate HGVS protein notation for an SNV consequence.
///
/// Returns the full `p.` string (e.g., `"p.Arg248Trp"`) or `None` if the
/// variant does not affect the protein sequence (UTR, intronic, splice, etc.).
///
/// # Arguments
///
/// * `consequences` — SO consequence terms from `ConsequenceResult`
/// * `amino_acids`  — VEP-style amino acid string: `"R/W"` for non-synonymous,
///   `"R"` for synonymous (produced by `codon::format_amino_acids`)
/// * `protein_start` — 1-based protein position
///
/// # Format details
///
/// No parentheses (matches VEP). Stop codons rendered as `Ter` (matches VEP).
/// Stop-lost renders as `extTer?` here — callers that have access to a
/// FASTA reader can compute the exact extension distance via the stop-lost
/// helper and substitute the concrete distance.
pub(crate) fn format_hgvs_p_snv(
    consequences: &[Consequence],
    amino_acids: Option<&str>,
    protein_start: Option<u32>,
) -> Option<String> {
    let aa_str = amino_acids?;
    let pos = protein_start?;

    // Parse "R/W" (non-synonymous) or "R" (synonymous) format.
    // `format_amino_acids` in codon.rs guarantees at least one character.
    let (ref_aa, alt_aa) = match aa_str.split_once('/') {
        Some((r, a)) if !r.is_empty() && !a.is_empty() => (r.as_bytes()[0], Some(a.as_bytes()[0])),
        None if !aa_str.is_empty() => (aa_str.as_bytes()[0], None),
        _ => return None,
    };

    // A partial codon at either end of an incomplete CDS translates to X:
    // the residue is unknown, so the protein effect is too. VEP prints
    // `p.TerN=` here, which names a stop that is not there (recorded in
    // VEP_DIVERGENCES.md).
    if ref_aa == b'X' {
        return Some("p.?".to_string());
    }

    // Consequence-priority order (defensive — consequences are mutually
    // exclusive in annotate_cds_snv, but the ordering is safe regardless).
    if consequences.contains(&Consequence::StartLost) {
        // Start codon disrupted — downstream consequence unknown.
        // VEP: p.Met1?  (no parentheses, position always 1)
        return Some("p.Met1?".to_string());
    }

    if consequences.contains(&Consequence::StopLost) {
        // Stop codon replaced by coding amino acid — extension with an
        // unknown new stop position (the caller's stop-lost helper
        // substitutes the computed distance).
        let alt = aa_three_letter(alt_aa?);
        return Some(format!("p.Ter{pos}{alt}extTer?"));
    }

    if consequences.contains(&Consequence::StopGained) {
        // Nonsense — amino acid replaced by stop codon.
        let r = aa_three_letter(ref_aa);
        return Some(format!("p.{r}{pos}Ter"));
    }

    if consequences.contains(&Consequence::MissenseVariant) {
        let r = aa_three_letter(ref_aa);
        let a = aa_three_letter(alt_aa?);
        return Some(format!("p.{r}{pos}{a}"));
    }

    if consequences.contains(&Consequence::SynonymousVariant)
        || consequences.contains(&Consequence::StartRetainedVariant)
    {
        let r = aa_three_letter(ref_aa);
        return Some(format!("p.{r}{pos}="));
    }

    if consequences.contains(&Consequence::StopRetainedVariant) {
        // ref_aa is '*', aa_three_letter returns "Ter"
        return Some(format!("p.Ter{pos}="));
    }

    // Non-CDS consequence with amino_acids/protein_start populated
    // (shouldn't happen for SNVs, but defensive).
    None
}

// ---------------------------------------------------------------------------
// Indel / complex-variant primitives
// ---------------------------------------------------------------------------

/// Concatenate three-letter amino acid codes for a slice of one-letter codes.
///
/// # Examples
///
/// `[b'G', b'S', b'K'] -> "GlySerLys"`, `[b'*'] -> "Ter"`.
fn format_aa_sequence(aas: &[u8]) -> String {
    let mut out = String::with_capacity(aas.len() * 3);
    for &aa in aas {
        out.push_str(aa_three_letter(aa));
    }
    out
}

/// Return the 0-based index of the first position where `ref_aas` and
/// `alt_aas` differ, or `None` if they are identical up to the shorter
/// length.
fn find_first_changed_aa(ref_aas: &[u8], alt_aas: &[u8]) -> Option<usize> {
    ref_aas.iter().zip(alt_aas.iter()).position(|(r, a)| r != a)
}

/// The first residue a frameshift changes, or `None` when none is observable.
///
/// Where the two translations agree through the shorter one, the change sits
/// at its end (the alternate frame ran out): `None` only when that end is the
/// end of the reference too.
fn first_changed_residue(ref_protein: &[u8], alt_protein: &[u8]) -> Option<usize> {
    find_first_changed_aa(ref_protein, alt_protein).or_else(|| {
        let idx = ref_protein.len().min(alt_protein.len());
        (idx < ref_protein.len()).then_some(idx)
    })
}

/// Append the transcript's 3'UTR to a frameshift's reference and alternate
/// coding sequences.
///
/// # Returns
///
/// `true` when a 3'UTR was appended, `false` for a transcript without one.
///
/// # Errors
///
/// Propagates FASTA read errors.
fn append_3prime_utr(
    ref_seq: &mut Vec<u8>,
    alt_seq: &mut Vec<u8>,
    chrom: &str,
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<bool, VarEffectError> {
    let utr_seq = fetch_3prime_utr_coding_seq(chrom, transcript, fasta)?;
    if utr_seq.is_empty() {
        return Ok(false);
    }
    ref_seq.extend_from_slice(&utr_seq);
    alt_seq.extend_from_slice(&utr_seq);
    Ok(true)
}

/// Translate a frameshift's reference and alternate sequences, each
/// truncated to whole codons.
///
/// # Arguments
///
/// * `ref_seq`, `alt_seq` - coding-strand sequences from the same codon start
/// * `first_codon` - 1-based CDS codon number of that start, for
///   selenocysteine codons in the reference
/// * `index` - the transcript's locate index
/// * `is_mito` - translate with the vertebrate mitochondrial code
///
/// # Returns
///
/// The reference and alternate one-letter translations; the alternate is
/// empty when it has no whole codon.
///
/// # Errors
///
/// Propagates translation errors.
fn translate_frameshift_pair(
    ref_seq: &[u8],
    alt_seq: &[u8],
    first_codon: u32,
    index: &LocateIndex,
    is_mito: bool,
) -> Result<(Vec<u8>, Vec<u8>), VarEffectError> {
    let ref_complete_len = ref_seq.len() - (ref_seq.len() % 3);
    let ref_protein = crate::codon::translate_reference_sequence(
        &ref_seq[..ref_complete_len],
        first_codon,
        index.selenocysteine_codons(),
        is_mito,
    )?;
    let alt_complete_len = alt_seq.len() - (alt_seq.len() % 3);
    let alt_protein = if alt_complete_len == 0 {
        Vec::new()
    } else {
        crate::codon::translate_sequence(&alt_seq[..alt_complete_len], is_mito)?
    };
    Ok((ref_protein, alt_protein))
}

/// Shift a deletion or duplication to the most C-terminal position within
/// `ref_protein` where the same amino acid(s) could be removed/inserted.
///
/// This implements the HGVS protein-level 3' rule: in amino acid repeats,
/// the most C-terminal residue is assigned as deleted/duplicated. Operates
/// on the codon-aligned window from the annotator, matching VEP's
/// `_protein_3prime_shift` behavior.
///
/// # Arguments
///
/// * `ref_protein` — full `ref_aas` window from the annotator
/// * `del_len` — length of the deleted/inserted amino acid stretch
/// * `initial_pos` — 0-based position within `ref_protein`
///
/// # Returns
///
/// Adjusted 0-based position (>= `initial_pos`).
fn apply_protein_3prime_rule(ref_protein: &[u8], del_len: usize, initial_pos: usize) -> usize {
    if del_len == 0 {
        return initial_pos;
    }
    let mut pos = initial_pos;
    // Shift right while the AA at `pos` matches the AA that would rotate
    // in from beyond the deletion range.
    while pos + del_len < ref_protein.len() && ref_protein[pos] == ref_protein[pos + del_len] {
        pos += 1;
    }
    pos
}

/// Translate a coding-strand sequence codon-by-codon and return the
/// 1-based position of the first stop codon, or `None` if no stop is
/// found before the sequence ends.
///
/// Incomplete trailing bases (`seq.len() % 3` remainder) are discarded.
/// Position 1 = the first codon of the input. The stop codon itself is
/// counted in the return value (i.e., if the 8th codon is a stop, returns
/// `Some(8)`).
fn scan_for_stop_codon(seq: &[u8], is_mito: bool) -> Option<u32> {
    for (i, codon) in seq.as_chunks::<3>().0.iter().enumerate() {
        let aa = translate_codon_for_transcript(codon, is_mito);
        if aa == b'*' {
            return Some(i as u32 + 1);
        }
    }
    None
}

/// Fetch the coding-strand 3'UTR exonic sequence starting from the first
/// base after the stop codon through all remaining exonic bases.
///
/// Returns bases in transcript 5'→3' order (coding strand). For
/// plus-strand transcripts, these are the plus-strand bases after
/// `cds_genomic_end`. For minus-strand transcripts, the reverse-
/// complemented bases before `cds_genomic_start`.
///
/// Returns an empty `Vec` for non-coding transcripts or transcripts with
/// no annotated 3'UTR exonic sequence.
fn fetch_3prime_utr_coding_seq(
    chrom: &str,
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<Vec<u8>, VarEffectError> {
    let cds_end = match transcript.strand {
        // Plus: CDS ends at cds_genomic_end (highest CDS coord).
        Strand::Plus => match transcript.cds_genomic_end {
            Some(e) => e,
            None => return Ok(Vec::new()),
        },
        // Minus: CDS ends at cds_genomic_start (lowest CDS coord =
        // transcript 3' boundary).
        Strand::Minus => match transcript.cds_genomic_start {
            Some(s) => s,
            None => return Ok(Vec::new()),
        },
    };

    let mut utr_seq = Vec::new();

    match transcript.strand {
        Strand::Plus => {
            // Walk exons in transcript order (ascending genomic).
            // 3'UTR bases are those after cds_genomic_end.
            for exon in &transcript.exons {
                if exon.genomic_end <= cds_end {
                    continue; // entirely before or at CDS end
                }
                let start = exon.genomic_start.max(cds_end);
                let end = exon.genomic_end;
                if start < end {
                    let bases = fasta.fetch_sequence(chrom, start, end)?;
                    utr_seq.extend_from_slice(&bases);
                }
            }
        }
        Strand::Minus => {
            // Exons are stored 5'→3' on transcript (descending genomic).
            // 3'UTR exons have the lowest genomic coords and appear at
            // the END of the exons array. For each exon whose genomic
            // range overlaps with [tx_start, cds_genomic_start), fetch
            // plus-strand bases and reverse-complement per chunk.
            for exon in &transcript.exons {
                if exon.genomic_start >= cds_end {
                    continue; // entirely within or before CDS
                }
                let start = exon.genomic_start;
                let end = exon.genomic_end.min(cds_end);
                if start < end {
                    let plus_bases = fasta.fetch_sequence(chrom, start, end)?;
                    // Plus-strand ascending genomic = 3'→5' on minus
                    // transcript. Reverse-complement to get coding strand.
                    let coding = reverse_complement(&plus_bases);
                    utr_seq.extend_from_slice(&coding);
                }
            }
        }
    }

    Ok(utr_seq)
}

// ---------------------------------------------------------------------------
// Indel / complex-variant formatters
// ---------------------------------------------------------------------------

/// The reference-CDS edit a frameshift describes, after HGVS 3' normalization.
///
/// Anchors the alternate-to-reference coordinate mapping: exon-exon junctions
/// downstream of `start` shift by [`CdsEdit::delta`] in the alternate
/// transcript, junctions at or before it do not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CdsEdit {
    /// 0-based CDS offset where the edit begins (3'-normalized).
    pub(crate) start: u32,
    /// Reference CDS bases replaced. Zero for a pure insertion.
    pub(crate) del_len: u32,
    /// Bases inserted at `start`. Zero for a pure deletion.
    pub(crate) ins_len: u32,
}

impl CdsEdit {
    /// Net CDS length change. Positive lengthens the alternate transcript.
    pub(crate) fn delta(&self) -> i64 {
        i64::from(self.ins_len) - i64::from(self.del_len)
    }

    /// 1-based CDS position of the last reference base the edit replaces.
    /// Equals `start` for a pure insertion, which sits between `start` and
    /// `start + 1`.
    pub(crate) fn last_replaced_pos(&self) -> u32 {
        self.start + self.del_len
    }
}

/// Where a frameshift's new termination codon sits, in **alternate**-protein
/// residue numbering.
///
/// Each variant names one outcome of the stop scan so a caller cannot confuse
/// "no stop exists" with "the stop could not be located" — they carry opposite
/// clinical meanings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PtcSite {
    /// The termination codon could not be located: the alternate allele yields
    /// no complete codon, produces no observable change, or the alternate
    /// protein is shorter than the first changed residue. Missing data, never
    /// an assertion that no stop exists.
    Indeterminate,
    /// The alternate reading frame runs to the end of the 3'UTR without
    /// reaching a stop codon (HGVS `fsTer?`).
    NoStopCodon,
    /// 1-based residue of the termination codon in the alternate protein.
    Residue(u32),
    /// The alternate frame keeps the reference stop codon in place: the
    /// protein is unchanged and no termination codon is created.
    ReferenceStop,
}

/// Protein-level result of a frameshift: the HGVS notation, the located
/// termination codon, and the edit that anchors its coordinate space.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FrameshiftProtein {
    /// HGVS p. notation, or `None` when no protein change is describable.
    pub(crate) notation: Option<String>,
    /// The new termination codon, in alternate-protein residue numbering.
    pub(crate) ptc: PtcSite,
    /// The 3'-normalized CDS edit these positions are relative to.
    pub(crate) edit: CdsEdit,
}

/// Generate HGVS p. notation for a frameshift variant.
///
/// Fetches the remaining CDS from the variant site to the end, builds the
/// alt reading frame, and scans for a new stop codon. If no stop is found
/// in the CDS, continues into the 3'UTR.
///
/// # HGVS rules
///
/// 1. The first amino acid *changed* (not first codon affected at DNA level)
///    is the starting position.
/// 2. If the first changed AA is itself a stop → nonsense (`p.Tyr4Ter`),
///    NOT `fsTer1`.
/// 3. Stop position counted from the first changed AA (position 1).
/// 4. `fsTer?` if no stop found.
///
/// # Returns
///
/// A [`FrameshiftProtein`] carrying the notation, the located termination
/// codon in alternate-protein numbering, and the 3'-normalized [`CdsEdit`]
/// those positions are relative to. The caller translates the termination
/// codon into an NMD verdict via `consequence::nmd::resolve_ptc`; this
/// function deliberately makes no NMD prediction of its own.
///
/// # Arguments
///
/// * `cds_offset_start` — 0-based CDS offset where the variant begins
/// * `cds_offset_end` — 0-based exclusive end (equal to start for insertions)
/// * `inserted_coding_bases` — coding-strand bases inserted at the variant
///   site (empty slice for pure deletions)
/// * `chrom`, `transcript`, `index`, `fasta` — for downstream CDS/UTR fetch
/// * `dna_shift` — HGVS 3' normalization shift in CDS bases (0 for
///   deletions / delins; `> 0` for insertions in repeat regions)
#[allow(clippy::too_many_arguments)]
pub(crate) fn format_hgvs_p_frameshift(
    cds_offset_start: u32,
    cds_offset_end: u32,
    inserted_coding_bases: &[u8],
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    dna_shift: u32,
) -> Result<FrameshiftProtein, VarEffectError> {
    let total_cds = index.total_cds_length();
    let is_mito = transcript.chrom == "chrM";

    // Apply HGVS 3' normalization shift when computing the protein effect.
    // VEP shifts the variant to its most 3' equivalent position before
    // determining the first changed amino acid. The shift is constrained
    // to exon boundaries by the caller, so genomic shift == CDS shift.
    let is_insertion = cds_offset_end == cds_offset_start;
    let (eff_start, eff_end) = if dna_shift > 0 {
        let shifted = cds_offset_start.saturating_add(dna_shift);
        if shifted < total_cds {
            // For insertions (zero-width), shift both start and end
            // to maintain the zero-width invariant.
            if is_insertion {
                (shifted, shifted)
            } else {
                (shifted, cds_offset_end.saturating_add(dna_shift))
            }
        } else {
            (cds_offset_start, cds_offset_end)
        }
    } else {
        (cds_offset_start, cds_offset_end)
    };

    // The 3'-normalized edit anchors every position this function returns.
    // `eff_end == eff_start` for insertions, so `del_len` is 0 there.
    let edit = CdsEdit {
        start: eff_start,
        del_len: eff_end - eff_start,
        ins_len: inserted_coding_bases.len() as u32,
    };

    let codon_start = (eff_start / 3) * 3;

    // Fetch full ref CDS from codon boundary to end.
    let mut ref_seq = fetch_cds_sequence(codon_start, total_cds, chrom, transcript, index, fasta)?;

    // Build alt CDS from pieces: prefix + inserted bases + suffix.
    let local_start = (eff_start - codon_start) as usize;
    let local_end = (eff_end - codon_start) as usize;
    let mut alt_seq =
        Vec::with_capacity(ref_seq.len() - (local_end - local_start) + inserted_coding_bases.len());
    alt_seq.extend_from_slice(&ref_seq[..local_start]);
    alt_seq.extend_from_slice(inserted_coding_bases);
    alt_seq.extend_from_slice(&ref_seq[local_end..]);

    let first_codon = codon_start / 3 + 1;
    let (mut ref_protein, mut alt_protein) =
        translate_frameshift_pair(&ref_seq, &alt_seq, first_codon, index, is_mito)?;
    let mut change = first_changed_residue(&ref_protein, &alt_protein);

    // Read on into the 3'UTR, as VEP does, when the coding sequence alone
    // shows no change or no stop after it. A deletion just before the stop
    // can leave the alternate frame ending before its first changed residue,
    // also when that residue is the stop (VEP 116 ci1766 ENST00001110657.1
    // c.417_420del p.Ter142ValextTer16; ci7204 NM_000190.4 c.1084del,
    // 3'-shifted from c.1083, p.Ter362AsnextTer?).
    let needs_utr =
        change.is_none_or(|idx| idx >= alt_protein.len() || !alt_protein[idx..].contains(&b'*'));
    if needs_utr && append_3prime_utr(&mut ref_seq, &mut alt_seq, chrom, transcript, fasta)? {
        (ref_protein, alt_protein) =
            translate_frameshift_pair(&ref_seq, &alt_seq, first_codon, index, is_mito)?;
        change = first_changed_residue(&ref_protein, &alt_protein);
    }

    if alt_protein.is_empty() {
        // Entire alt is incomplete — no protein produced.
        return Ok(FrameshiftProtein {
            notation: None,
            ptc: PtcSite::Indeterminate,
            edit,
        });
    }
    let Some(change_idx) = change else {
        // No observable change.
        return Ok(FrameshiftProtein {
            notation: None,
            ptc: PtcSite::Indeterminate,
            edit,
        });
    };

    // Read on into the 3'UTR, the first change can lie past the reference
    // stop: the frameshift keeps the stop codon and the protein (VEP 116
    // NM_000969.5 c.890_893dup and NM_000179.3 c.4069_4082dup: `p.Ter298=`,
    // `p.Ter1361=`).
    if let Some(stop_idx) = ref_protein[..change_idx].iter().position(|&aa| aa == b'*') {
        let stop_pos = (codon_start / 3) + 1 + stop_idx as u32;
        return Ok(FrameshiftProtein {
            notation: Some(format!("p.Ter{stop_pos}=")),
            ptc: PtcSite::ReferenceStop,
            edit,
        });
    }

    let protein_pos = (codon_start / 3) + 1 + change_idx as u32;
    let ref_aa = ref_protein[change_idx];
    let Some(&alt_aa) = alt_protein.get(change_idx) else {
        // The alternate frame ends before the change, with no 3'UTR to read
        // on into. Report unknown.
        return Ok(FrameshiftProtein {
            notation: Some(format!("p.{}{protein_pos}?fsTer?", aa_three_letter(ref_aa))),
            ptc: PtcSite::Indeterminate,
            edit,
        });
    };

    // Check if first changed AA in alt is a stop → nonsense, not fs.
    if alt_aa == b'*' {
        // The first changed residue is itself a stop, so it *is* the
        // termination codon — there is no `stop_n` offset to add.
        return Ok(FrameshiftProtein {
            notation: Some(format!("p.{}{}Ter", aa_three_letter(ref_aa), protein_pos)),
            ptc: PtcSite::Residue(protein_pos),
            edit,
        });
    }

    // VEP 116 counts to the first stop of the whole alternate translation,
    // which reads an upstream selenocysteine TGA as one, so it writes
    // `fsTer?` (NM_000581.4 c.149del: p.Gly50AlafsTer?). The PTC keeps the
    // real termination codon.
    let ter = |stop_n: u32| {
        if index.has_selenocysteine_before(protein_pos) {
            "?".to_string()
        } else {
            stop_n.to_string()
        }
    };
    let alt3 = aa_three_letter(alt_aa);
    let Some(star_pos) = alt_protein[change_idx..].iter().position(|&a| a == b'*') else {
        // No stop before the end of the 3'UTR.
        let notation = if ref_aa == b'*' {
            format!("p.Ter{protein_pos}{alt3}extTer?")
        } else {
            format!("p.{}{protein_pos}{alt3}fsTer?", aa_three_letter(ref_aa))
        };
        return Ok(FrameshiftProtein {
            notation: Some(notation),
            ptc: PtcSite::NoStopCodon,
            edit,
        });
    };
    // `star_pos` is 0-based from the first changed residue; HGVS `fsTer`
    // counts that residue as 1. A frameshift whose first changed residue is
    // the reference stop is an extension, counted from the first codon after
    // the old stop (VEP 116 ci8341 NM_001321821.2 c.1275_1276dup
    // p.Ter426CysextTer4).
    let stop_n = star_pos as u32 + 1;
    let notation = if ref_aa == b'*' {
        format!("p.Ter{protein_pos}{alt3}extTer{}", ter(stop_n - 1))
    } else {
        format!(
            "p.{}{protein_pos}{alt3}fsTer{}",
            aa_three_letter(ref_aa),
            ter(stop_n)
        )
    };
    Ok(FrameshiftProtein {
        notation: Some(notation),
        ptc: PtcSite::Residue(protein_pos + stop_n - 1),
        edit,
    })
}

/// Generate HGVS p. notation for an inframe deletion.
///
/// Compares ref and alt amino acid sequences to determine which amino
/// acids are deleted, applies the protein-level 3' rule for repeat
/// regions, and formats the notation.
///
/// # Arguments
///
/// * `ref_aas` — reference amino acids from codon-expanded region
/// * `alt_aas` — alternate amino acids (shorter than ref)
/// * `protein_start` — 1-based position of the first codon in the window
/// * `consequences` — SO consequence terms (checked for StartLost/StopLost)
pub(crate) fn format_hgvs_p_inframe_del(
    ref_aas: &[u8],
    alt_aas: &[u8],
    protein_start: u32,
    consequences: &[Consequence],
) -> Option<String> {
    // StartLost takes precedence (compound StartLost + InframeDeletion).
    if consequences.contains(&Consequence::StartLost) {
        return Some("p.Met1?".to_string());
    }
    // StopLost requires FASTA for extension scan — caller overrides with
    // `format_hgvs_p_del_extension` after this returns None.
    if consequences.contains(&Consequence::StopLost) {
        return None;
    }

    // Strip common prefix.
    let prefix_len = ref_aas
        .iter()
        .zip(alt_aas.iter())
        .take_while(|(r, a)| r == a)
        .count();

    // Strip common suffix from remaining portions.
    let ref_tail = &ref_aas[prefix_len..];
    let alt_tail = &alt_aas[prefix_len..];
    let suffix_len = ref_tail
        .iter()
        .rev()
        .zip(alt_tail.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    let remaining_ref = &ref_aas[prefix_len..ref_aas.len() - suffix_len];
    let remaining_alt = &alt_aas[prefix_len..alt_aas.len() - suffix_len];

    if remaining_ref.is_empty() {
        // No ref AAs removed — shouldn't happen for inframe_deletion,
        // but defensive.
        return None;
    }

    if remaining_alt.is_empty() {
        // Pure deletion.
        let del_len = remaining_ref.len();
        let adjusted_pos = apply_protein_3prime_rule(ref_aas, del_len, prefix_len);
        let start_pos = protein_start + adjusted_pos as u32;

        if del_len == 1 {
            let aa = ref_aas[adjusted_pos];
            Some(format!("p.{}{}del", aa_three_letter(aa), start_pos))
        } else {
            let first_aa = ref_aas[adjusted_pos];
            let last_aa = ref_aas[adjusted_pos + del_len - 1];
            let end_pos = start_pos + del_len as u32 - 1;
            Some(format!(
                "p.{}{}_{}{}del",
                aa_three_letter(first_aa),
                start_pos,
                aa_three_letter(last_aa),
                end_pos,
            ))
        }
    } else {
        // Deletion changed flanking codons → delins at protein level.
        let start_pos = protein_start + prefix_len as u32;
        let end_pos = start_pos + remaining_ref.len() as u32 - 1;

        // Truncate alt at stop codon if present.
        let alt_display = truncate_at_stop(remaining_alt);

        if remaining_ref.len() == 1 {
            Some(format!(
                "p.{}{}delins{}",
                aa_three_letter(remaining_ref[0]),
                start_pos,
                format_aa_sequence(alt_display),
            ))
        } else {
            Some(format!(
                "p.{}{}_{}{}delins{}",
                aa_three_letter(remaining_ref[0]),
                start_pos,
                aa_three_letter(
                    *remaining_ref
                        .last()
                        .expect("remaining_ref verified non-empty")
                ),
                end_pos,
                format_aa_sequence(alt_display),
            ))
        }
    }
}

/// Reference protein residues read on demand, for placing an in-frame
/// insertion's HGVS p. against the whole protein rather than a codon window.
pub(crate) trait ReferenceProtein {
    /// Number of residues before the terminal stop codon (VEP's `_peptide`
    /// length).
    fn residue_count(&self) -> u32;

    /// Residues `[start, end)` by 0-based index, `end <= residue_count()`.
    ///
    /// # Errors
    ///
    /// Returns an error if the coding sequence cannot be read.
    fn residues(&self, start: u32, end: u32) -> Result<Vec<u8>, VarEffectError>;
}

/// The peptides of an in-frame insertion in VEP's codon window
/// (`TranscriptVariationAllele::codon`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct InsertionPeptides {
    /// Reference peptide; empty for an insertion between two codons.
    pub ref_aas: Vec<u8>,
    /// Alternate peptide.
    pub alt_aas: Vec<u8>,
    /// 1-based residue VEP's `translation_start` gives: the codon holding the
    /// insertion, or the codon after it at a codon boundary.
    pub tl_start: u32,
    /// 1-based residue VEP's `translation_end` gives: `tl_start` inside a
    /// codon, `tl_start - 1` at a codon boundary.
    pub tl_end: u32,
}

/// Residues read per step while shifting an inserted peptide 3'.
const SHIFT_CHUNK: u32 = 32;

/// Generate HGVS p. notation for an in-frame insertion, as VEP 116's
/// `hgvs_protein` does.
///
/// Port of `TranscriptVariationAllele.pm` (release/116): `_clip_alleles`
/// (common prefix and suffix, left unclipped when both peptides start with a
/// stop), `_get_hgvs_protein_type`, `_check_peptides_post_var` and
/// `_shift_3prime` (the inserted peptide shifts C-terminal across the whole
/// protein), `_check_for_peptide_duplication`, `_get_surrounding_peptides`
/// and `_get_hgvs_protein_format`. VEP 116 writes `p.Glu199_Tyr200insTer`
/// (NM_001625.4 c.597_599dup), `p.Ter394delinsTer` and
/// `p.Ter394delinsTrpGly` (NM_000546.6 c.1180_1181insAAG / insGGG), and no
/// HGVS p. for a non-duplicating insertion just before the stop codon
/// (NM_000546.6 c.1179_1180insCTGTAA).
///
/// # Arguments
///
/// * `peptides` -- reference and alternate peptides at the HGVS (3'-shifted)
///   insertion point
/// * `protein` -- reference protein of the transcript
///
/// # Returns
///
/// The notation without accession, or `None` where VEP writes none.
///
/// # Errors
///
/// Returns an error if reference residues cannot be read.
pub(crate) fn format_hgvs_p_inframe_ins(
    peptides: &InsertionPeptides,
    protein: &impl ReferenceProtein,
) -> Result<Option<String>, VarEffectError> {
    let original_ref = peptides.ref_aas.as_slice();
    let mut ref_aas = original_ref;
    let mut alt_aas = peptides.alt_aas.as_slice();
    let mut start = i64::from(peptides.tl_start);
    let mut end = i64::from(peptides.tl_end);

    // `_clip_alleles` returns before clipping anything when the next residues
    // are both stops.
    let both_stop = |r: &[u8], a: &[u8]| r.first() == Some(&b'*') && a.first() == Some(&b'*');
    let mut unclipped = false;
    while let (Some(&r), Some(&a)) = (ref_aas.first(), alt_aas.first()) {
        if both_stop(ref_aas, alt_aas) {
            unclipped = true;
            break;
        }
        if r != a {
            break;
        }
        ref_aas = &ref_aas[1..];
        alt_aas = &alt_aas[1..];
        start += 1;
    }
    if unclipped {
        ref_aas = original_ref;
        alt_aas = &peptides.alt_aas;
        start = i64::from(peptides.tl_start);
    } else {
        while let (Some(&r), Some(&a)) = (ref_aas.last(), alt_aas.last()) {
            if r != a {
                break;
            }
            ref_aas = &ref_aas[..ref_aas.len() - 1];
            alt_aas = &alt_aas[..alt_aas.len() - 1];
            end -= 1;
        }
    }

    if !ref_aas.is_empty() {
        // `delins`: one residue after clipping, or the unclipped stop codon.
        let alt = format_aa_sequence(truncate_at_stop(alt_aas));
        let first = aa_three_letter(ref_aas[0]);
        if start == end {
            return Ok(Some(format!("p.{first}{start}delins{alt}")));
        }
        let (lo, hi) = (start.min(end), start.max(end));
        let last = aa_three_letter(*ref_aas.last().expect("ref_aas checked non-empty"));
        return Ok(Some(format!("p.{first}{lo}_{last}{hi}delins{alt}")));
    }

    // `ins`. VEP turns the first stop into `X` before shifting, so the
    // inserted peptide never rotates past it and never matches as a dup; a
    // marker that matches no residue does the same without being confused
    // with an unknown amino acid.
    const STOP_MARKER: u8 = 0;
    let mut inserted = alt_aas.to_vec();
    if let Some(stop) = inserted.iter().position(|&aa| aa == b'*') {
        inserted[stop] = STOP_MARKER;
    }
    let count = i64::from(protein.residue_count());
    // `_get_surrounding_peptides` reads the protein with the original
    // reference peptide appended when that starts with a stop.
    let stop_tail: &[u8] = if original_ref.first() == Some(&b'*') {
        original_ref
    } else {
        &[]
    };
    let surrounding_len = count + stop_tail.len() as i64;
    let surrounding = |from: i64, to: i64| -> Result<Vec<u8>, VarEffectError> {
        let mut out = Vec::new();
        if from < count {
            out = protein.residues(from as u32, to.min(count) as u32)?;
        }
        for i in from.max(count)..to.min(surrounding_len) {
            out.push(stop_tail[(i - count) as usize]);
        }
        Ok(out)
    };

    // `_check_peptides_post_var` reads from residue `end + 1` (0-based index
    // `end`) and shifts nothing unless the protein extends past it.
    let post_start = end;
    if surrounding_len > end + 1 {
        let post_len = surrounding_len - post_start;
        let max_steps = post_len - inserted.len() as i64;
        let mut step = 0;
        'shift: while step <= max_steps {
            let chunk_end = (post_start + step + i64::from(SHIFT_CHUNK)).min(surrounding_len);
            let chunk = surrounding(post_start + step, chunk_end)?;
            for &next in &chunk {
                if step > max_steps || inserted[0] != next {
                    break 'shift;
                }
                inserted.rotate_left(1);
                start += 1;
                end += 1;
                step += 1;
            }
        }
    }

    // `_check_for_peptide_duplication`.
    let len = inserted.len() as i64;
    let dup_from = start - len - 1;
    if dup_from >= 0 && !inserted.contains(&b'*') {
        let upstream = surrounding(dup_from, dup_from + len)?;
        if upstream == inserted {
            let dup_start = start - len;
            let dup_end = start - 1;
            let first = aa_three_letter(inserted[0]);
            if dup_start < dup_end {
                let last = aa_three_letter(inserted[inserted.len() - 1]);
                return Ok(Some(format!("p.{first}{dup_start}_{last}{dup_end}dup")));
            }
            return Ok(Some(format!("p.{first}{dup_start}dup")));
        }
    }

    // `p.{R1}{m}_{R2}{m+1}ins...` with `m = min(start, end)`; none when the
    // protein has no residue after `m`.
    let m = start.min(end);
    if m < 1 || surrounding_len <= m {
        return Ok(None);
    }
    let flank = surrounding(m - 1, m + 1)?;
    for aa in &mut inserted {
        if *aa == STOP_MARKER {
            *aa = b'*';
        }
    }
    Ok(Some(format!(
        "p.{}{}_{}{}ins{}",
        aa_three_letter(flank[0]),
        m,
        aa_three_letter(flank[1]),
        m + 1,
        format_aa_sequence(truncate_at_stop(&inserted)),
    )))
}

/// Generate HGVS p. notation for an inframe delins or MNV.
///
/// Compares ref and alt amino acid sequences after stripping common
/// prefix/suffix. Produces missense notation for 1:1 substitution,
/// or delins notation for more complex changes.
///
/// # Arguments
///
/// * `ref_aas` — reference amino acids from codon-expanded region
/// * `alt_aas` — alternate amino acids
/// * `protein_start` — 1-based position of the first codon in the window
/// * `protein_end` — 1-based position of the last codon in the window
/// * `consequences` — SO consequence terms
pub(crate) fn format_hgvs_p_delins(
    ref_aas: &[u8],
    alt_aas: &[u8],
    protein_start: u32,
    _protein_end: u32,
    consequences: &[Consequence],
) -> Option<String> {
    if consequences.contains(&Consequence::StartLost) {
        return Some("p.Met1?".to_string());
    }
    if consequences.contains(&Consequence::StopLost) {
        return None;
    }

    // Strip common prefix.
    let prefix_len = ref_aas
        .iter()
        .zip(alt_aas.iter())
        .take_while(|(r, a)| r == a)
        .count();

    // Strip common suffix from remaining portions.
    let ref_tail = &ref_aas[prefix_len..];
    let alt_tail = &alt_aas[prefix_len..];
    let suffix_len = ref_tail
        .iter()
        .rev()
        .zip(alt_tail.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    let remaining_ref = &ref_aas[prefix_len..ref_aas.len() - suffix_len];
    let remaining_alt = &alt_aas[prefix_len..alt_aas.len() - suffix_len];

    if remaining_ref.is_empty() && remaining_alt.is_empty() {
        // Entirely synonymous across all codons.
        if ref_aas.first() == Some(&b'M') && protein_start == 1 {
            return Some("p.Met1=".to_string());
        }
        if ref_aas.last() == Some(&b'*') {
            let pos = protein_start + ref_aas.len() as u32 - 1;
            return Some(format!("p.Ter{}=", pos));
        }
        // Synonymous — report the first AA at the window start.
        let aa = ref_aas[0];
        return Some(format!("p.{}{}=", aa_three_letter(aa), protein_start));
    }

    // 1:1 amino acid substitution → missense/nonsense notation.
    if remaining_ref.len() == 1 && remaining_alt.len() == 1 {
        let r = remaining_ref[0];
        let a = remaining_alt[0];
        let pos = protein_start + prefix_len as u32;

        if r == a {
            if r == b'*' {
                return Some(format!("p.Ter{}=", pos));
            }
            return Some(format!("p.{}{}=", aa_three_letter(r), pos));
        }
        if a == b'*' {
            return Some(format!("p.{}{}Ter", aa_three_letter(r), pos));
        }
        if r == b'*' {
            // Stop lost → needs extension scan. Return None here.
            return None;
        }
        return Some(format!(
            "p.{}{}{}",
            aa_three_letter(r),
            pos,
            aa_three_letter(a),
        ));
    }

    // Multi-AA delins.
    if remaining_ref.is_empty() {
        // Pure insertion at protein level (rare for delins call path).
        return None;
    }

    let start_pos = protein_start + prefix_len as u32;
    let end_pos = start_pos + remaining_ref.len() as u32 - 1;
    let alt_display = truncate_at_stop(remaining_alt);

    if remaining_ref.len() == 1 {
        Some(format!(
            "p.{}{}delins{}",
            aa_three_letter(remaining_ref[0]),
            start_pos,
            format_aa_sequence(alt_display),
        ))
    } else {
        Some(format!(
            "p.{}{}_{}{}delins{}",
            aa_three_letter(remaining_ref[0]),
            start_pos,
            aa_three_letter(
                *remaining_ref
                    .last()
                    .expect("remaining_ref verified non-empty")
            ),
            end_pos,
            format_aa_sequence(alt_display),
        ))
    }
}

/// Compute the actual extension distance for a stop_lost variant by
/// scanning into the 3'UTR for a new in-frame stop codon.
///
/// # Arguments
///
/// * `alt_aa` — the amino acid replacing the original stop codon
/// * `protein_pos` — 1-based position of the original stop codon
/// * `chrom`, `transcript`, `index`, `fasta` — for 3'UTR fetch
///
/// # Returns
///
/// The full `p.Ter{pos}{AltAA}extTer{N}` string, or `extTer?` if no
/// in-frame stop is found in the 3'UTR. A transcript with a
/// selenocysteine codon always gets `extTer?`: VEP 116 counts to the first
/// stop of the whole alternate translation, which reads the selenocysteine
/// `TGA` as one (NM_000581.4 c.610T>C: `p.Ter204GlnextTer?`).
pub(crate) fn format_hgvs_p_extension(
    alt_aa: u8,
    protein_pos: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<String, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";
    let alt_name = aa_three_letter(alt_aa);
    if index.has_selenocysteine_before(protein_pos) {
        return Ok(format!("p.Ter{protein_pos}{alt_name}extTer?"));
    }

    let utr_seq = fetch_3prime_utr_coding_seq(chrom, transcript, fasta)?;

    if let Some(utr_stop) = scan_for_stop_codon(&utr_seq, is_mito) {
        // Extension distance = position of the new stop codon in the
        // 3'UTR (1-based). The replacement amino acid at the original
        // stop position is already captured by `{alt_name}` in the
        // notation and is NOT counted in the distance.
        Ok(format!("p.Ter{protein_pos}{alt_name}extTer{utr_stop}"))
    } else {
        Ok(format!("p.Ter{protein_pos}{alt_name}extTer?"))
    }
}

/// Compute the extension distance for a stop codon deletion.
///
/// When the stop codon is deleted (inframe deletion of the entire stop
/// codon), the reading frame extends into the 3'UTR until a new in-frame
/// stop is found. Unlike SNV stop-lost (where the stop is replaced by a
/// coding amino acid), there is no replacement residue — the extension
/// distance counts only 3'UTR codons.
///
/// VEP format: `p.Ter{pos}delextTer{N}` (e.g., `p.Ter394delextTer8`).
///
/// # Arguments
///
/// * `protein_pos` — 1-based position of the original stop codon
/// * `chrom`, `transcript`, `fasta` — for 3'UTR sequence fetch
/// * `selenocysteine_upstream` — a selenocysteine codon precedes the first
///   deleted residue ([`LocateIndex::has_selenocysteine_before`])
///
/// With a selenocysteine codon left upstream the extension is `delextTer?`,
/// as in [`format_hgvs_p_extension`] (NM_000581.4 c.610_612del:
/// `p.Ter204delextTer?`); VEP 116 counts when the deletion removes it
/// (NM_001093771.3 c.1942_1950del: `p.Sec648_Ter650delextTer53`).
pub(crate) fn format_hgvs_p_del_extension(
    protein_pos: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    selenocysteine_upstream: bool,
    fasta: &FastaReader,
) -> Result<String, VarEffectError> {
    if selenocysteine_upstream {
        return Ok(format!("p.Ter{protein_pos}delextTer?"));
    }
    let is_mito = transcript.chrom == "chrM";
    let utr_seq = fetch_3prime_utr_coding_seq(chrom, transcript, fasta)?;

    if let Some(utr_stop) = scan_for_stop_codon(&utr_seq, is_mito) {
        // For deletion, the stop codon is removed entirely (no
        // replacement AA occupies the position). The first 3'UTR codon
        // directly fills position `protein_pos`, so the extension
        // distance is one less than the UTR stop position.
        let ext_distance = utr_stop - 1;
        Ok(format!("p.Ter{protein_pos}delextTer{ext_distance}"))
    } else {
        Ok(format!("p.Ter{protein_pos}delextTer?"))
    }
}

/// Truncate an amino acid slice at the first stop codon (`b'*'`),
/// including the stop in the returned slice. Used for display in delins
/// notation where AAs after a stop are not listed.
fn truncate_at_stop(aas: &[u8]) -> &[u8] {
    match aas.iter().position(|&a| a == b'*') {
        Some(pos) => &aas[..=pos],
        None => aas,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a consequence slice from a single variant.
    fn csq(c: Consequence) -> Vec<Consequence> {
        vec![c]
    }

    #[test]
    fn partial_codon_has_unknown_protein_effect() {
        // VEP prints `p.Ter1=` / `p.Ter309=` for amino acid X; vareffect
        // reports the unknown effect instead (VEP_DIVERGENCES.md).
        let five_prime =
            format_hgvs_p_snv(&csq(Consequence::CodingSequenceVariant), Some("X"), Some(1));
        assert_eq!(five_prime.as_deref(), Some("p.?"));
        let three_prime = format_hgvs_p_snv(
            &[
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant,
            ],
            Some("X"),
            Some(309),
        );
        assert_eq!(three_prime.as_deref(), Some("p.?"));
    }

    #[test]
    fn missense_basic() {
        let result = format_hgvs_p_snv(&csq(Consequence::MissenseVariant), Some("R/W"), Some(248));
        assert_eq!(result.as_deref(), Some("p.Arg248Trp"));
    }

    #[test]
    fn missense_braf() {
        let result = format_hgvs_p_snv(&csq(Consequence::MissenseVariant), Some("V/E"), Some(600));
        assert_eq!(result.as_deref(), Some("p.Val600Glu"));
    }

    #[test]
    fn synonymous() {
        // amino_acids is "L" (single letter, no slash) for synonymous
        let result = format_hgvs_p_snv(&csq(Consequence::SynonymousVariant), Some("L"), Some(54));
        assert_eq!(result.as_deref(), Some("p.Leu54="));
    }

    #[test]
    fn nonsense() {
        let result = format_hgvs_p_snv(&csq(Consequence::StopGained), Some("W/*"), Some(26));
        assert_eq!(result.as_deref(), Some("p.Trp26Ter"));
    }

    #[test]
    fn stop_lost() {
        let result = format_hgvs_p_snv(&csq(Consequence::StopLost), Some("*/R"), Some(328));
        assert_eq!(result.as_deref(), Some("p.Ter328ArgextTer?"));
    }

    #[test]
    fn start_lost() {
        let result = format_hgvs_p_snv(&csq(Consequence::StartLost), Some("M/V"), Some(1));
        assert_eq!(result.as_deref(), Some("p.Met1?"));
    }

    #[test]
    fn start_retained() {
        // amino_acids is "M" (single letter) for start-retained synonymous
        let result = format_hgvs_p_snv(&csq(Consequence::StartRetainedVariant), Some("M"), Some(1));
        assert_eq!(result.as_deref(), Some("p.Met1="));
    }

    #[test]
    fn stop_retained() {
        // amino_acids is "*" (single letter) for stop-retained synonymous
        let result =
            format_hgvs_p_snv(&csq(Consequence::StopRetainedVariant), Some("*"), Some(394));
        assert_eq!(result.as_deref(), Some("p.Ter394="));
    }

    #[test]
    fn no_amino_acids_returns_none() {
        let result = format_hgvs_p_snv(&csq(Consequence::MissenseVariant), None, Some(248));
        assert_eq!(result, None);
    }

    #[test]
    fn no_protein_start_returns_none() {
        let result = format_hgvs_p_snv(&csq(Consequence::MissenseVariant), Some("R/W"), None);
        assert_eq!(result, None);
    }

    #[test]
    fn all_20_amino_acids() {
        // Verify all 20 standard amino acids produce correct three-letter
        // codes when used in a missense context.
        let cases: &[(u8, &str)] = &[
            (b'A', "Ala"),
            (b'C', "Cys"),
            (b'D', "Asp"),
            (b'E', "Glu"),
            (b'F', "Phe"),
            (b'G', "Gly"),
            (b'H', "His"),
            (b'I', "Ile"),
            (b'K', "Lys"),
            (b'L', "Leu"),
            (b'M', "Met"),
            (b'N', "Asn"),
            (b'P', "Pro"),
            (b'Q', "Gln"),
            (b'R', "Arg"),
            (b'S', "Ser"),
            (b'T', "Thr"),
            (b'V', "Val"),
            (b'W', "Trp"),
            (b'Y', "Tyr"),
        ];
        for &(one, three) in cases {
            let aa_str = format!("{}/A", one as char);
            let result =
                format_hgvs_p_snv(&csq(Consequence::MissenseVariant), Some(&aa_str), Some(100));
            let expected = format!("p.{three}100Ala");
            assert_eq!(
                result.as_deref(),
                Some(expected.as_str()),
                "failed for amino acid '{}'",
                one as char,
            );
        }
    }

    #[test]
    fn non_cds_consequence_returns_none() {
        // IntronVariant with no amino_acids/protein_start → None
        let result = format_hgvs_p_snv(&csq(Consequence::IntronVariant), None, None);
        assert_eq!(result, None);
    }

    // -----------------------------------------------------------------------
    // Primitives
    // -----------------------------------------------------------------------

    #[test]
    fn format_aa_sequence_basic() {
        assert_eq!(format_aa_sequence(b"GSK"), "GlySerLys");
        assert_eq!(format_aa_sequence(b"*"), "Ter");
        assert_eq!(format_aa_sequence(&[]), "");
        assert_eq!(format_aa_sequence(b"M"), "Met");
    }

    #[test]
    fn find_first_changed_basic() {
        assert_eq!(find_first_changed_aa(b"MRK", b"MWK"), Some(1));
        assert_eq!(find_first_changed_aa(b"MRK", b"MRK"), None);
        assert_eq!(find_first_changed_aa(b"MRK", b"WRK"), Some(0));
        assert_eq!(find_first_changed_aa(b"MR", b"MRK"), None); // identical up to shorter
    }

    #[test]
    fn scan_for_stop_basic() {
        // TAA = stop codon (standard code)
        assert_eq!(scan_for_stop_codon(b"ATGTAA", false), Some(2));
        // First codon is stop
        assert_eq!(scan_for_stop_codon(b"TAA", false), Some(1));
        // No stop
        assert_eq!(scan_for_stop_codon(b"ATGATG", false), None);
        // Incomplete trailing ignored
        assert_eq!(scan_for_stop_codon(b"ATGTAAG", false), Some(2));
        // TAG
        assert_eq!(scan_for_stop_codon(b"ATGTAGATG", false), Some(2));
        // Empty
        assert_eq!(scan_for_stop_codon(b"", false), None);
        // Less than 3 bases
        assert_eq!(scan_for_stop_codon(b"AT", false), None);
    }

    #[test]
    fn protein_3prime_rule_shift() {
        // ref: M W S S S H D — delete one S at position 2 (0-based).
        // Shifts to position 4 (most C-terminal S).
        let ref_prot = b"MWSSSHD";
        assert_eq!(apply_protein_3prime_rule(ref_prot, 1, 2), 4);

        // No shift when no repeat
        let ref_prot2 = b"MWFHD";
        assert_eq!(apply_protein_3prime_rule(ref_prot2, 1, 2), 2);

        // Two-AA deletion in repeat: SS at position 2.
        // ref: MWSSSSHD — delete SS. Shifts from 2 to 4.
        let ref_prot3 = b"MWSSSSHD";
        assert_eq!(apply_protein_3prime_rule(ref_prot3, 2, 2), 4);
    }

    // -----------------------------------------------------------------------
    // Inframe deletion
    // -----------------------------------------------------------------------

    #[test]
    fn inframe_del_single() {
        // CFTR deltaF508: ref=[I,F], alt=[I], protein_start=507
        let result = format_hgvs_p_inframe_del(b"IF", b"I", 507, &[Consequence::InframeDeletion]);
        assert_eq!(result.as_deref(), Some("p.Phe508del"));
    }

    #[test]
    fn inframe_del_range() {
        // EGFR E746_A750del: ref=[K,E,L,R,E,A], alt=[K], protein_start=745
        let result =
            format_hgvs_p_inframe_del(b"KELREA", b"K", 745, &[Consequence::InframeDeletion]);
        assert_eq!(result.as_deref(), Some("p.Glu746_Ala750del"));
    }

    #[test]
    fn inframe_del_3prime_rule() {
        // ref: M W S S S H D, alt: M W S S H D → p.Ser5del (3' rule)
        let result =
            format_hgvs_p_inframe_del(b"MWSSSHD", b"MWSSHD", 1, &[Consequence::InframeDeletion]);
        assert_eq!(result.as_deref(), Some("p.Ser5del"));
    }

    #[test]
    fn inframe_del_becomes_delins() {
        // ref=[R,K], alt=[W] → p.Arg10_Lys11delinsTrp
        let result = format_hgvs_p_inframe_del(b"RK", b"W", 10, &[Consequence::InframeDeletion]);
        assert_eq!(result.as_deref(), Some("p.Arg10_Lys11delinsTrp"));
    }

    #[test]
    fn start_lost_returns_met1() {
        let result = format_hgvs_p_inframe_del(
            b"M",
            &[],
            1,
            &[Consequence::StartLost, Consequence::InframeDeletion],
        );
        assert_eq!(result.as_deref(), Some("p.Met1?"));
    }

    #[test]
    fn start_lost_plus_inframe_del() {
        // Compound consequence: StartLost takes priority.
        let result = format_hgvs_p_inframe_del(
            b"MR",
            b"R",
            1,
            &[Consequence::StartLost, Consequence::InframeDeletion],
        );
        assert_eq!(result.as_deref(), Some("p.Met1?"));
    }

    // -----------------------------------------------------------------------
    // Inframe insertion
    // -----------------------------------------------------------------------

    /// An in-memory reference protein without its terminal stop.
    struct TestProtein(&'static [u8]);

    impl ReferenceProtein for TestProtein {
        fn residue_count(&self) -> u32 {
            self.0.len() as u32
        }

        fn residues(&self, start: u32, end: u32) -> Result<Vec<u8>, VarEffectError> {
            Ok(self.0[start as usize..end as usize].to_vec())
        }
    }

    fn ins_p(
        protein: &'static [u8],
        ref_aas: &[u8],
        alt_aas: &[u8],
        tl: (u32, u32),
    ) -> Option<String> {
        let peptides = InsertionPeptides {
            ref_aas: ref_aas.to_vec(),
            alt_aas: alt_aas.to_vec(),
            tl_start: tl.0,
            tl_end: tl.1,
        };
        format_hgvs_p_inframe_ins(&peptides, &TestProtein(protein)).expect("in-memory protein")
    }

    #[test]
    fn inframe_ins_between_codons() {
        // VEP 116 NM_000059.4 CDS 42/43 shape: `-/L` between Ile and Phe.
        assert_eq!(
            ins_p(b"MPIFST", b"", b"L", (4, 3)).as_deref(),
            Some("p.Ile3_Phe4insLeu")
        );
    }

    #[test]
    fn inframe_ins_inside_codon_clips_the_reference() {
        // NM_000546.6 c.525_526insAGC shape: `R/RS` -> p.Arg175_Cys176insSer.
        assert_eq!(
            ins_p(b"MARCD", b"R", b"RS", (3, 3)).as_deref(),
            Some("p.Arg3_Cys4insSer")
        );
    }

    #[test]
    fn inframe_ins_inside_codon_without_flank_is_delins() {
        // NM_000546.6 c.523_524insTTT: `R/LC` -> p.Arg175delinsLeuCys.
        assert_eq!(
            ins_p(b"MARCD", b"R", b"LC", (3, 3)).as_deref(),
            Some("p.Arg3delinsLeuCys")
        );
    }

    #[test]
    fn inframe_ins_dup_shifts_across_the_whole_repeat() {
        assert_eq!(
            ins_p(b"MAKKKKD", b"", b"K", (3, 2)).as_deref(),
            Some("p.Lys6dup")
        );
        assert_eq!(
            ins_p(b"MAGSGSD", b"", b"GS", (3, 2)).as_deref(),
            Some("p.Gly5_Ser6dup")
        );
    }

    #[test]
    fn inframe_ins_of_a_stop_is_an_insertion() {
        // VEP 116 NM_206933.4 c.987_988insTAA (`-/*`) and NM_001625.4
        // c.597_599dup (`E/E*`): p.Asp329_Asn330insTer, p.Glu199_Tyr200insTer.
        assert_eq!(
            ins_p(b"MEYQ", b"", b"*", (3, 2)).as_deref(),
            Some("p.Glu2_Tyr3insTer")
        );
        assert_eq!(
            ins_p(b"MEYQ", b"E", b"E*", (2, 2)).as_deref(),
            Some("p.Glu2_Tyr3insTer")
        );
        // NM_000374.5 c.385_386insAGT (`S/*S`): p.Ala128_Ser129insTer.
        assert_eq!(
            ins_p(b"MASQ", b"S", b"*S", (3, 3)).as_deref(),
            Some("p.Ala2_Ser3insTer")
        );
    }

    #[test]
    fn inframe_ins_keeps_unknown_residues_apart_from_the_stop() {
        assert_eq!(
            ins_p(b"MEYQ", b"", b"GX*", (3, 2)).as_deref(),
            Some("p.Glu2_Tyr3insGlyXaaTer")
        );
    }

    #[test]
    fn inframe_ins_before_the_stop_codon_has_no_notation() {
        // VEP 116 NM_000546.6 c.1179_1180insCTGTAA (`-/L*`): no HGVS p.
        assert_eq!(ins_p(b"MKD", b"", b"L*", (4, 3)), None);
    }

    #[test]
    fn inframe_ins_does_not_shift_onto_the_last_residue_alone() {
        // `_get_surrounding_peptides` gives no post sequence for the last
        // residue, so the insertion stays where it is.
        assert_eq!(
            ins_p(b"MKAD", b"", b"D", (4, 3)).as_deref(),
            Some("p.Ala3_Asp4insAsp")
        );
    }

    #[test]
    fn inframe_ins_inside_the_stop_codon() {
        // VEP 116 NM_000546.6 c.1180_1181insAAG (`*/*G`) and insGGG (`*/WG`).
        assert_eq!(
            ins_p(b"MKD", b"*", b"*G", (4, 4)).as_deref(),
            Some("p.Ter4delinsTer")
        );
        assert_eq!(
            ins_p(b"MKD", b"*", b"WG", (4, 4)).as_deref(),
            Some("p.Ter4delinsTrpGly")
        );
        // The flank reaches the stop codon through the reference peptide.
        assert_eq!(
            ins_p(b"MKD", b"*", b"L*", (4, 4)).as_deref(),
            Some("p.Asp3_Ter4insLeu")
        );
    }

    // -----------------------------------------------------------------------
    // Delins / MNV
    // -----------------------------------------------------------------------

    #[test]
    fn delins_shrink() {
        // ref=[C,K], alt=[W] → p.Cys28_Lys29delinsTrp
        let result =
            format_hgvs_p_delins(b"CK", b"W", 28, 29, &[Consequence::ProteinAlteringVariant]);
        assert_eq!(result.as_deref(), Some("p.Cys28_Lys29delinsTrp"));
    }

    #[test]
    fn delins_grow() {
        // ref=[C], alt=[W,V] → p.Cys28delinsTrpVal
        let result =
            format_hgvs_p_delins(b"C", b"WV", 28, 28, &[Consequence::ProteinAlteringVariant]);
        assert_eq!(result.as_deref(), Some("p.Cys28delinsTrpVal"));
    }

    #[test]
    fn delins_single_missense() {
        // ref=[R], alt=[W] → p.Arg248Trp (missense, not delins)
        let result = format_hgvs_p_delins(b"R", b"W", 248, 248, &[Consequence::MissenseVariant]);
        assert_eq!(result.as_deref(), Some("p.Arg248Trp"));
    }

    #[test]
    fn delins_synonymous() {
        // ref=[R], alt=[R] → p.Arg248=
        let result = format_hgvs_p_delins(b"R", b"R", 248, 248, &[Consequence::SynonymousVariant]);
        assert_eq!(result.as_deref(), Some("p.Arg248="));
    }

    #[test]
    fn delins_stop_gained() {
        // ref=[P,K], alt=[L,*] → p.Pro578_Lys579delinsLeuTer
        let result = format_hgvs_p_delins(
            b"PK",
            b"L*",
            578,
            579,
            &[Consequence::StopGained, Consequence::ProteinAlteringVariant],
        );
        assert_eq!(result.as_deref(), Some("p.Pro578_Lys579delinsLeuTer"));
    }

    #[test]
    fn truncate_at_stop_helper() {
        assert_eq!(truncate_at_stop(b"L*K"), b"L*");
        assert_eq!(truncate_at_stop(b"LK"), b"LK");
        assert_eq!(truncate_at_stop(b"*"), b"*");
    }
}
