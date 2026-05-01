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
    for (i, codon_bytes) in seq.chunks_exact(3).enumerate() {
        // SAFETY: chunks_exact(3) guarantees exactly 3 bytes.
        let codon: [u8; 3] = [codon_bytes[0], codon_bytes[1], codon_bytes[2]];
        let aa = translate_codon_for_transcript(&codon, is_mito);
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
) -> Result<Option<String>, VarEffectError> {
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

    // When the variant overlaps the stop codon, extend both ref and alt
    // into the 3'UTR before translation. Without this, the alt protein
    // is too short to identify the first changed amino acid and the new
    // stop position, producing "p.TerN?fsTer?" instead of the correct
    // "p.TerNXxxfsTerM" notation that VEP emits.
    let stop_codon_start = total_cds.saturating_sub(3);
    let extended_into_utr = if eff_start >= stop_codon_start {
        let utr_seq = fetch_3prime_utr_coding_seq(chrom, transcript, fasta)?;
        if !utr_seq.is_empty() {
            ref_seq.extend_from_slice(&utr_seq);
            alt_seq.extend_from_slice(&utr_seq);
            true
        } else {
            false
        }
    } else {
        false
    };

    // Translate ref: truncate to complete codons (always exact when
    // UTR extension is not applied; may have trailing bases after UTR).
    let ref_complete_len = ref_seq.len() - (ref_seq.len() % 3);
    let ref_protein = crate::codon::translate_sequence(&ref_seq[..ref_complete_len], is_mito)?;

    // Translate alt: truncate to complete codons, keep trailing bases
    // for CDS/UTR junction carryover.
    let alt_complete_len = alt_seq.len() - (alt_seq.len() % 3);
    let trailing_bases = if extended_into_utr {
        // UTR already included; no separate continuation needed.
        Vec::new()
    } else {
        alt_seq[alt_complete_len..].to_vec()
    };
    if alt_complete_len == 0 {
        // Entire alt is incomplete — no protein produced.
        return Ok(None);
    }
    let alt_protein = crate::codon::translate_sequence(&alt_seq[..alt_complete_len], is_mito)?;

    // Find first changed amino acid.
    let change_idx = match find_first_changed_aa(&ref_protein, &alt_protein) {
        Some(idx) => idx,
        None => {
            // ref and alt identical through min length. The change is at
            // the truncation point (alt is shorter due to deletion).
            let idx = ref_protein.len().min(alt_protein.len());
            if idx >= ref_protein.len() {
                return Ok(None); // no observable change
            }
            idx
        }
    };

    let protein_pos = (codon_start / 3) + 1 + change_idx as u32;
    let ref_aa = ref_protein[change_idx];

    // Check if first changed AA in alt is a stop → nonsense, not fs.
    if let Some(&alt_aa) = alt_protein.get(change_idx) {
        if alt_aa == b'*' {
            return Ok(Some(format!(
                "p.{}{}Ter",
                aa_three_letter(ref_aa),
                protein_pos,
            )));
        }

        // Scan for stop in the alt protein from the change point onward.
        // First check within the already-translated alt protein.
        if let Some(star_pos) = alt_protein[change_idx..].iter().position(|&a| a == b'*') {
            // star_pos is 0-based from change_idx. HGVS position is
            // 1-based.
            let stop_n = star_pos as u32 + 1;
            return Ok(Some(format!(
                "p.{}{}{}fsTer{}",
                aa_three_letter(ref_aa),
                protein_pos,
                aa_three_letter(alt_aa),
                stop_n,
            )));
        }

        // No stop in CDS portion. Continue into 3'UTR with trailing
        // byte carryover — unless UTR was already spliced into the
        // sequences above (avoid double-counting).
        if !extended_into_utr {
            let codons_past_change = alt_protein.len() - change_idx;

            let utr_seq = fetch_3prime_utr_coding_seq(chrom, transcript, fasta)?;
            let mut continuation = trailing_bases;
            continuation.extend_from_slice(&utr_seq);

            if let Some(utr_stop) = scan_for_stop_codon(&continuation, is_mito) {
                let stop_n = codons_past_change as u32 + utr_stop;
                return Ok(Some(format!(
                    "p.{}{}{}fsTer{}",
                    aa_three_letter(ref_aa),
                    protein_pos,
                    aa_three_letter(alt_aa),
                    stop_n,
                )));
            }
        }

        // No stop found anywhere.
        return Ok(Some(format!(
            "p.{}{}{}fsTer?",
            aa_three_letter(ref_aa),
            protein_pos,
            aa_three_letter(alt_aa),
        )));
    }

    // alt_protein shorter than change_idx — the change is beyond the alt
    // protein. This can happen for very large deletions. Report unknown.
    Ok(Some(format!(
        "p.{}{}?fsTer?",
        aa_three_letter(ref_aa),
        protein_pos,
    )))
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

/// Generate HGVS p. notation for an inframe insertion.
///
/// Checks for duplication first (inserted AAs match the immediately
/// preceding reference), then formats as insertion or duplication.
///
/// # Arguments
///
/// * `ref_aas` — reference amino acid(s) at the affected codon
/// * `alt_aas` — alternate amino acids (longer than ref)
/// * `protein_start` — 1-based position of the affected codon
/// * `consequences` — SO consequence terms
/// * `right_flanking_aa` — amino acid at `protein_start + 1` (next codon,
///   needed for HGVS flanking position notation)
pub(crate) fn format_hgvs_p_inframe_ins(
    ref_aas: &[u8],
    alt_aas: &[u8],
    protein_start: u32,
    consequences: &[Consequence],
    right_flanking_aa: u8,
) -> Option<String> {
    if consequences.contains(&Consequence::StartLost) {
        return Some("p.Met1?".to_string());
    }
    if consequences.contains(&Consequence::StopLost) {
        return None;
    }

    // Strip common prefix and suffix to isolate inserted AAs.
    let prefix_len = ref_aas
        .iter()
        .zip(alt_aas.iter())
        .take_while(|(r, a)| r == a)
        .count();

    let ref_tail = &ref_aas[prefix_len..];
    let alt_tail = &alt_aas[prefix_len..];
    let suffix_len = ref_tail
        .iter()
        .rev()
        .zip(alt_tail.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    let remaining_ref = &ref_aas[prefix_len..ref_aas.len() - suffix_len];
    let inserted_aas = &alt_aas[prefix_len..alt_aas.len() - suffix_len];

    if inserted_aas.is_empty() {
        return None; // no actual insertion
    }

    // If ref portion is non-empty, this is really a delins.
    if !remaining_ref.is_empty() {
        let start_pos = protein_start + prefix_len as u32;
        let end_pos = start_pos + remaining_ref.len() as u32 - 1;
        let alt_display = truncate_at_stop(inserted_aas);
        if remaining_ref.len() == 1 {
            return Some(format!(
                "p.{}{}delins{}",
                aa_three_letter(remaining_ref[0]),
                start_pos,
                format_aa_sequence(alt_display),
            ));
        }
        return Some(format!(
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
        ));
    }

    // Pure insertion. Check if stop gained (alone).
    if consequences.contains(&Consequence::StopGained) && inserted_aas == [b'*'] {
        // The insertion introduced a premature stop at the ref codon position.
        let pos = protein_start + prefix_len as u32;
        if prefix_len > 0 {
            let ref_aa = ref_aas[prefix_len - 1];
            return Some(format!("p.{}{}Ter", aa_three_letter(ref_aa), pos));
        }
    }

    // Duplication check: do the inserted AAs match the immediately
    // preceding (N-terminal) reference AAs?
    let ins_len = inserted_aas.len();
    if prefix_len >= ins_len {
        let preceding = &ref_aas[prefix_len - ins_len..prefix_len];
        if preceding == inserted_aas {
            // Duplication. Apply protein 3' rule.
            let dup_start_0based = prefix_len - ins_len;
            let adjusted = apply_protein_3prime_rule(ref_aas, ins_len, dup_start_0based);
            let dup_start = protein_start + adjusted as u32;

            if ins_len == 1 {
                let aa = ref_aas[adjusted];
                return Some(format!("p.{}{}dup", aa_three_letter(aa), dup_start));
            }
            let first_aa = ref_aas[adjusted];
            let last_aa = ref_aas[adjusted + ins_len - 1];
            let dup_end = dup_start + ins_len as u32 - 1;
            return Some(format!(
                "p.{}{}_{}{}dup",
                aa_three_letter(first_aa),
                dup_start,
                aa_three_letter(last_aa),
                dup_end,
            ));
        }
    }

    // Insertion notation: p.{Left}_{Right}ins{Inserted}
    let left_aa;
    let left_pos;
    let right_aa;
    let right_pos;

    if prefix_len > 0 && prefix_len <= ref_aas.len() {
        left_aa = ref_aas[prefix_len - 1];
        left_pos = protein_start + prefix_len as u32 - 1;
        // Right flanking: if the insertion is at the end of the ref
        // window, use the caller-provided next codon AA.
        if prefix_len < ref_aas.len() {
            right_aa = ref_aas[prefix_len];
            right_pos = left_pos + 1;
        } else {
            right_aa = right_flanking_aa;
            right_pos = left_pos + 1;
        }
    } else {
        // Insertion at position 0 within the codon window — edge case.
        left_aa = ref_aas[0];
        left_pos = protein_start;
        right_aa = right_flanking_aa;
        right_pos = left_pos + 1;
    }

    // Truncate inserted AAs at stop if present.
    let ins_display = truncate_at_stop(inserted_aas);

    Some(format!(
        "p.{}{}_{}{}ins{}",
        aa_three_letter(left_aa),
        left_pos,
        aa_three_letter(right_aa),
        right_pos,
        format_aa_sequence(ins_display),
    ))
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
/// in-frame stop is found in the 3'UTR.
pub(crate) fn format_hgvs_p_extension(
    alt_aa: u8,
    protein_pos: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    _index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<String, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";
    let alt_name = aa_three_letter(alt_aa);

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
pub(crate) fn format_hgvs_p_del_extension(
    protein_pos: u32,
    chrom: &str,
    transcript: &TranscriptModel,
    fasta: &FastaReader,
) -> Result<String, VarEffectError> {
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

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a consequence slice from a single variant.
    fn csq(c: Consequence) -> Vec<Consequence> {
        vec![c]
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

    #[test]
    fn inframe_ins_simple() {
        // ref=[K], alt=[K,Q,S,K], protein_start=2, right_flanking=M
        let result =
            format_hgvs_p_inframe_ins(b"K", b"KQSK", 2, &[Consequence::InframeInsertion], b'M');
        assert_eq!(result.as_deref(), Some("p.Lys2_Met3insGlnSerLys"));
    }

    #[test]
    fn inframe_ins_dup_single() {
        // ref=[K], alt=[K,K], protein_start=5, right_flanking=M
        // Inserted K matches preceding K → p.Lys5dup
        let result =
            format_hgvs_p_inframe_ins(b"K", b"KK", 5, &[Consequence::InframeInsertion], b'M');
        assert_eq!(result.as_deref(), Some("p.Lys5dup"));
    }

    #[test]
    fn inframe_ins_with_stop() {
        // Inserted sequence contains stop → include Ter
        let result = format_hgvs_p_inframe_ins(
            b"P",
            b"PG*",
            2,
            &[Consequence::InframeInsertion, Consequence::StopGained],
            b'I',
        );
        assert_eq!(result.as_deref(), Some("p.Pro2_Ile3insGlyTer"));
    }

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
