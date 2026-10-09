//! Indel (insertion and deletion) consequence annotation.

use super::helpers::{
    add_splice_spill, build_base_result, compute_cdna_position_exonic,
    compute_cdna_position_for_cds, fetch_cds_sequence, fetch_ref_codon, finalize_consequences,
    has_complete_start_codon, is_incomplete_terminal_codon,
};
use super::{Consequence, ConsequenceResult};
use crate::codon::{
    format_amino_acids_indel, format_codons_indel, reverse_complement,
    translate_codon_for_transcript, translate_sequence,
};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_c;
use crate::locate::helpers::{
    exonic_splice_spill, frameshift_intron_endpoint, insertion_intron_effects,
};
use crate::locate::{
    IndelLocation, IndelRegion, LocateIndex, format_exon_number, format_intron_number, locate_indel,
};
use crate::types::{Strand, TranscriptModel, TranscriptTier};

/// Annotate a deletion against one transcript.
///
/// `start` and `end` define the deleted range `[start, end)` in 0-based
/// half-open coordinates (post-trimming). `deleted_bases` is the plus-strand
/// sequence of deleted bases (trimmed REF).
///
/// # Consequence assignment priority
///
/// 1. Splice canonical overlap -> splice_donor_variant / splice_acceptor_variant
/// 2. Exon/intron boundary crossing -> boundary-spanning analysis
/// 3. CDS: frameshift (`len%3 != 0`) or inframe_deletion (`len%3 == 0`)
/// 4. Non-CDS: same mapping as SNVs
pub fn annotate_deletion(
    chrom: &str,
    start: u64,
    end: u64,
    deleted_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let mut result = annotate_deletion_core(
        chrom,
        start,
        end,
        deleted_bases,
        transcript,
        locate_index,
        fasta,
    )?;
    add_splice_spill(
        &mut result,
        exonic_splice_spill(start, end, transcript, locate_index),
    );
    Ok(result)
}

/// [`annotate_deletion`] without the splice terms a one-base frameshift intron spills
/// onto exon bases.
fn annotate_deletion_core(
    chrom: &str,
    start: u64,
    end: u64,
    _deleted_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let location = locate_indel(chrom, start, end, transcript, locate_index)?;

    // HGVS 3' normalization: shift deletion to the most 3' equivalent
    // position on the coding strand for HGVS notation. Skip boundary-spanning
    // deletions (multi-exon shifts are not meaningful for HGVS c.).
    let shift = if !location.crosses_exon_boundary {
        crate::normalize::compute_3prime_shift_deletion(chrom, start, end, transcript, fasta)?
    } else {
        0
    };
    let (hgvs_start, hgvs_end) = match transcript.strand {
        Strand::Plus => (start + shift, end + shift),
        Strand::Minus => (start - shift, end - shift),
    };
    let hgvs = hgvs_c::format_deletion_hgvs(chrom, hgvs_start, hgvs_end, transcript, locate_index)?;

    if let Some(intron) = frameshift_intron_endpoint(start, end, transcript, locate_index) {
        let mut result = super::complex::frameshift_intron_indel_result(
            start,
            end,
            intron,
            &location,
            transcript,
            locate_index,
        );
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    if location.overlaps_splice_canonical
        && exonic_splice_spill(start, end, transcript, locate_index).is_none()
    {
        let mut result = build_splice_indel_result(transcript, &location, start, end)?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    if location.crosses_exon_boundary {
        let mut result = super::complex::annotate_boundary_spanning_deletion(
            chrom,
            start,
            end,
            transcript,
            locate_index,
            fasta,
            &location,
        )?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    // CDS/UTR boundary within the same exon
    if matches!(location.region, IndelRegion::BoundarySpanning) {
        let mut result = super::complex::annotate_boundary_spanning_deletion(
            chrom,
            start,
            end,
            transcript,
            locate_index,
            fasta,
            &location,
        )?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    let mut result = match &location.region {
        IndelRegion::Cds {
            cds_offset_start,
            cds_offset_end,
        } => {
            let del_cds_len = cds_offset_end - cds_offset_start;
            let exon_index = location.exon_index.unwrap_or(0);

            if is_incomplete_terminal_codon(*cds_offset_start, locate_index) {
                partial_terminal_codon_indel_result(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_end,
                    &[],
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                )
            } else if is_incomplete_terminal_codon(cds_offset_end - 1, locate_index) {
                deletion_into_partial_codon_result(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_end,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                )
            } else if !del_cds_len.is_multiple_of(3) {
                annotate_cds_frameshift(
                    chrom,
                    *cds_offset_start,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                    *cds_offset_end,
                    None,
                    None,
                    0,
                )
            } else {
                annotate_cds_inframe_deletion(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_end,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                )
            }
        }
        _ => build_noncds_indel_result(transcript, &location, start, end),
    }?;

    result.hgvs_c = hgvs;
    Ok(result)
}

/// Translate a coding sequence whose last codon may be incomplete: complete
/// codons as usual, a trailing partial codon as `X` unless the complete
/// codons are a lone stop (`TranscriptVariationAllele.pm::peptide`,
/// ensembl-variation release/116: VEP 116 c.658_659insGA `X/*`, edge
/// insGGGTGAG `-/G*X`).
fn translate_with_partial_codon(seq: &[u8], is_mito: bool) -> Result<Vec<u8>, VarEffectError> {
    let whole = seq.len() - seq.len() % 3;
    let mut aas = translate_sequence(&seq[..whole], is_mito)?;
    if whole < seq.len() && aas != b"*" {
        aas.push(b'X');
    }
    Ok(aas)
}

/// Result for a deletion that starts in a complete codon and ends in the
/// incomplete terminal codon of a CDS whose length is not a whole number of
/// codons.
///
/// Ending inside the partial codon, it shifts or removes complete codons:
/// `frameshift_variant` (or `inframe_deletion`), as VEP 116 gives for
/// NM_001424184.1 `c.657_658del` (codons `tcATa/tca`, amino acids `SX/S`).
/// Removing the CDS end itself leaves no reading frame to shift and no stop
/// codon to lose: `incomplete_terminal_codon_variant` +
/// `coding_sequence_variant`. VEP's `frameshift_variant` + `stop_lost` there
/// (`c.656_659del`) names a stop the transcript does not have
/// (VEP_DIVERGENCES.md), and has no HGVS p. Otherwise HGVS p. is `p.?` and
/// the termination codon is undetermined: the new frame runs past the
/// annotated sequence.
#[allow(clippy::too_many_arguments)]
fn deletion_into_partial_codon_result(
    chrom: &str,
    cds_offset_start: u32,
    cds_offset_end: u32,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
    is_splice_region: bool,
) -> Result<ConsequenceResult, VarEffectError> {
    let total = index.total_cds_length();
    let codon_start = (cds_offset_start / 3) * 3;
    let ref_seq = fetch_cds_sequence(codon_start, total, chrom, transcript, index, fasta)?;
    let local_start = (cds_offset_start - codon_start) as usize;
    let local_end = (cds_offset_end - codon_start) as usize;
    let mut alt_seq = ref_seq[..local_start].to_vec();
    alt_seq.extend_from_slice(&ref_seq[local_end..]);

    let mut consequences = if cds_offset_end >= total {
        vec![
            Consequence::IncompleteTerminalCodonVariant,
            Consequence::CodingSequenceVariant,
        ]
    } else if (cds_offset_end - cds_offset_start).is_multiple_of(3) {
        vec![Consequence::InframeDeletion]
    } else {
        vec![Consequence::FrameshiftVariant]
    };
    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }
    let impact = finalize_consequences(&mut consequences);

    let is_mito = transcript.chrom == "chrM";
    let ref_aas = translate_with_partial_codon(&ref_seq, is_mito)?;
    let alt_aas = translate_with_partial_codon(&alt_seq, is_mito)?;
    let protein_start = codon_start / 3 + 1;

    let mut result = build_base_result(transcript, consequences);
    result.impact = impact;
    result.protein_start = Some(protein_start);
    result.protein_end = Some(total.div_ceil(3));
    result.codons = Some(format_codons_indel(
        &ref_seq,
        &alt_seq,
        local_start,
        local_end,
    ));
    result.amino_acids = Some(format_amino_acids_indel(&ref_aas, &alt_aas));
    result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
    result.cds_position = Some(cds_offset_start + 1);
    result.cds_position_end = Some(cds_offset_end);
    result.cdna_position = Some(compute_cdna_position_for_cds(cds_offset_start, index));
    result.cdna_position_end = Some(compute_cdna_position_for_cds(cds_offset_end - 1, index));
    // A deletion removing the CDS end has no protein change to name (VEP 116
    // NM_001424184.1 c.656_659del: no HGVS p.); one ending inside the partial
    // codon has an unknown one.
    result.hgvs_p = (cds_offset_end < total).then(|| "p.?".to_string());
    if result
        .consequences
        .contains(&Consequence::FrameshiftVariant)
    {
        result.ptc = super::PtcStatus::Indeterminate;
    }
    Ok(result)
}

/// Result for a deletion or insertion that starts in the incomplete terminal
/// codon of a CDS whose length is not a whole number of codons.
///
/// VEP's `frameshift` and inframe predicates do not apply to a partial codon
/// (`VariationEffect.pm::partial_codon`). VEP 116: removing the whole partial
/// codon is `protein_altering_variant` + `incomplete_terminal_codon_variant`
/// (ENST00000379099.3 `c.923del`, NM_001424184.1 `c.658_659del`: codons
/// `TA/-`, amino acids `X/-`); removing part of it, or inserting without
/// completing a codon, is `incomplete_terminal_codon_variant` +
/// `coding_sequence_variant` (`c.658del`: `Ta/a`, `X`); an insertion that
/// completes a stop codon adds `stop_gained` (`c.658_659insG`: `ta/tGa`,
/// `X/*`). HGVS p. is `p.?`: the residue is unknown.
#[allow(clippy::too_many_arguments)]
fn partial_terminal_codon_indel_result(
    chrom: &str,
    cds_offset_start: u32,
    cds_offset_end: u32,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
) -> Result<ConsequenceResult, VarEffectError> {
    let total = index.total_cds_length();
    let codon_start = total - total % 3;
    let ref_bases = fetch_cds_sequence(codon_start, total, chrom, transcript, index, fasta)?;
    let is_insertion = cds_offset_end == cds_offset_start;
    // An insertion goes before `cds_offset_start` on both strands.
    let edit_start = cds_offset_start;
    let edit_end = if is_insertion {
        edit_start
    } else {
        cds_offset_end.min(total)
    };
    let local_start = (edit_start - codon_start) as usize;
    let local_end = (edit_end - codon_start) as usize;
    let coding_inserted = match transcript.strand {
        Strand::Plus => inserted_bases.to_vec(),
        Strand::Minus => reverse_complement(inserted_bases),
    };
    let mut alt_bases = Vec::with_capacity(ref_bases.len() + coding_inserted.len());
    alt_bases.extend_from_slice(&ref_bases[..local_start]);
    alt_bases.extend_from_slice(&coding_inserted);
    alt_bases.extend_from_slice(&ref_bases[local_end..]);

    let is_mito = transcript.chrom == "chrM";
    // Every completed codon is translated, as VEP does: `X/*` with
    // `stop_gained` (NM_001424184.1 c.658_659insG), `X/W*X` for
    // c.658_659insGGTGA, `X/R` otherwise (ENST00000518351.2 c.544_545insG).
    let translation = if alt_bases.len() >= 3 {
        Some(translate_with_partial_codon(&alt_bases, is_mito)?)
    } else {
        None
    };
    let stop_index = translation
        .as_ref()
        .and_then(|aas| aas.iter().position(|&aa| aa == b'*'));
    let (mut consequences, amino_acids) = match (alt_bases.is_empty(), &translation) {
        (true, _) => (
            vec![
                Consequence::ProteinAlteringVariant,
                Consequence::IncompleteTerminalCodonVariant,
            ],
            "X/-".to_string(),
        ),
        (false, Some(aas)) => {
            let mut terms = vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant,
            ];
            if stop_index.is_some() {
                terms.push(Consequence::StopGained);
            }
            let peptide: String = aas.iter().map(|&b| char::from(b)).collect();
            (terms, format!("X/{peptide}"))
        }
        (false, None) => (
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant,
            ],
            "X".to_string(),
        ),
    };
    let impact = finalize_consequences(&mut consequences);
    let mut result = build_base_result(transcript, consequences);
    result.impact = impact;
    let protein = codon_start / 3 + 1;
    result.protein_start = Some(protein);
    result.protein_end = Some(protein);
    result.codons = Some(format_codons_indel(
        &ref_bases,
        &alt_bases,
        local_start,
        local_end,
    ));
    result.amino_acids = Some(amino_acids);
    result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
    // An insertion reports the base before and the base after the insertion
    // point, as the other insertion paths do (VEP 116 NM_001424184.1
    // c.658_659insG: CDS 658); a deletion its first and last removed base.
    let (first, last) = if is_insertion {
        (edit_start.saturating_sub(1), edit_start.min(total - 1))
    } else {
        (edit_start, edit_end - 1)
    };
    result.cds_position = Some(first + 1);
    result.cds_position_end = Some(last + 1);
    result.cdna_position = Some(compute_cdna_position_for_cds(first, index));
    result.cdna_position_end = Some(compute_cdna_position_for_cds(last, index));
    result.hgvs_p = Some("p.?".to_string());
    if let Some(stop) = stop_index {
        // The new stop sits at or after the original terminus, in the last
        // exon: resolve it like any other termination codon.
        result.ptc = super::nmd::resolve_ptc(
            crate::hgvs_p::PtcSite::Residue(protein + stop as u32),
            &crate::hgvs_p::CdsEdit {
                start: edit_start,
                del_len: edit_end - edit_start,
                ins_len: coding_inserted.len() as u32,
            },
            index,
        );
    }
    Ok(result)
}

/// Where a CDS insertion sits relative to the incomplete terminal codon of a
/// CDS whose length is not a whole number of codons.
#[derive(Debug, PartialEq)]
pub(super) enum PartialCodonInsertion {
    /// Between two bases of the partial codon.
    Inside,
    /// Between the last complete codon and the partial codon, at this
    /// padded CDS offset (the insertion goes before it).
    FivePrimeEdge(u32),
}

/// Classify a CDS insertion against the incomplete terminal codon.
///
/// `insert_before` is the CDS offset the insertion goes before
/// ([`IndelRegion::Cds`] for an insertion). `None` when the CDS has no
/// partial codon or the insertion is elsewhere.
pub(super) fn partial_codon_insertion_site(
    insert_before: u32,
    index: &LocateIndex,
) -> Option<PartialCodonInsertion> {
    let total = index.total_cds_length();
    if total.is_multiple_of(3) {
        return None;
    }
    let partial_start = total - total % 3;
    if insert_before == partial_start && partial_start > 0 {
        // With no complete codon before it (a CDS shorter than three bases)
        // there is no edge; the insertion is inside the only codon.
        Some(PartialCodonInsertion::FivePrimeEdge(insert_before))
    } else if insert_before >= partial_start {
        Some(PartialCodonInsertion::Inside)
    } else {
        None
    }
}

/// Result for an insertion between the last complete codon and the
/// incomplete terminal codon.
///
/// VEP 116 translates the inserted bases on their own
/// ([`translate_with_partial_codon`]). Every such insertion is `inframe_insertion` +
/// `incomplete_terminal_codon_variant`; `stop_gained` is added when the
/// translation has a stop and `coding_sequence_variant` when it has an `X`
/// (NM_001424184.1 c.657_658insG `-/X`, insGGG `-/G`, insTGA and insTGAG
/// `-/*` and insGGGTGAG `-/G*X` with `stop_gained`; ENST00000518351.2 the
/// same on the minus strand).
/// Codons are `-/<inserted>` and the positions are the two flanking bases.
/// A stop is a termination codon at the original terminus, in the last
/// exon, and gets HGVS p. `p.?` as on the rest of the partial codon (VEP
/// writes `p.Ter220=` or nothing); without a stop there is no HGVS p.
fn partial_codon_edge_insertion_result(
    insert_before: u32,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    index: &LocateIndex,
    exon_index: u16,
) -> Result<ConsequenceResult, VarEffectError> {
    let coding_inserted = match transcript.strand {
        Strand::Plus => inserted_bases.to_vec(),
        Strand::Minus => reverse_complement(inserted_bases),
    };
    let is_mito = transcript.chrom == "chrM";
    let translation = translate_with_partial_codon(&coding_inserted, is_mito)?;
    let stop_codon_index = translation
        .iter()
        .position(|&aa| aa == b'*')
        .map(|i| i as u32);

    let mut consequences = vec![
        Consequence::InframeInsertion,
        Consequence::IncompleteTerminalCodonVariant,
    ];
    if stop_codon_index.is_some() {
        consequences.push(Consequence::StopGained);
    }
    // VEP's `coding_unknown`: the peptide has an `X`.
    if translation.contains(&b'X') {
        consequences.push(Consequence::CodingSequenceVariant);
    }
    let impact = finalize_consequences(&mut consequences);
    let mut result = build_base_result(transcript, consequences);
    result.impact = impact;
    let inserted: String = coding_inserted.iter().map(|&b| char::from(b)).collect();
    let peptide: String = translation.iter().map(|&b| char::from(b)).collect();
    result.codons = Some(format!("-/{inserted}"));
    result.amino_acids = Some(format!("-/{peptide}"));
    result.protein_start = Some(insert_before / 3);
    result.protein_end = Some(insert_before / 3 + 1);
    result.cds_position = Some(insert_before);
    result.cds_position_end = Some(insert_before + 1);
    // `insert_before` is at least 3: the 5' edge exists only after a
    // complete codon (`partial_codon_insertion_site`).
    result.cdna_position = Some(compute_cdna_position_for_cds(insert_before - 1, index));
    result.cdna_position_end = Some(compute_cdna_position_for_cds(insert_before, index));
    result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
    if let Some(stop) = stop_codon_index {
        let residue = insert_before / 3 + 1 + stop;
        result.hgvs_p = Some("p.?".to_string());
        result.ptc = super::nmd::resolve_ptc(
            crate::hgvs_p::PtcSite::Residue(residue),
            &crate::hgvs_p::CdsEdit {
                start: insert_before,
                del_len: 0,
                ins_len: coding_inserted.len() as u32,
            },
            index,
        );
    }
    Ok(result)
}

/// Annotate an insertion against one transcript.
///
/// `pos` is the 0-based insertion point (post-trimming). The insertion
/// occurs between positions `pos - 1` and `pos`. `inserted_bases` is the
/// plus-strand sequence of inserted bases (trimmed ALT).
///
/// The region comes from [`locate_indel`] (classified at the insertion
/// point, both strands), the splice and intron terms from VEP's intron
/// predicates for an insertion: a canonical-site insertion is the canonical
/// term alone (VEP 116 NM_000546.6 `c.782+1_782+2insC`:
/// `splice_donor_variant`), and an insertion at an exon/intron junction is
/// the coding change plus `splice_region_variant` (`c.782_782+1insT`).
pub fn annotate_insertion(
    chrom: &str,
    pos: u64,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let location = locate_indel(chrom, pos, pos, transcript, locate_index)?;
    let effects = insertion_intron_effects(pos, transcript)?;

    // HGVS 3' normalization: shift the insertion to the most 3' equivalent
    // position on the coding strand, across exon/intron junctions as VEP's
    // genomic shift does, with the inserted bases rotated to match. At the
    // shifted position, is_duplication() in format_insertion_hgvs re-checks
    // the 5' flank, so an insertion shifted to a repeat boundary is a dup.
    let shift = crate::normalize::compute_3prime_shift_insertion(
        chrom,
        pos,
        inserted_bases,
        transcript.strand,
        fasta,
    )?;
    let shifted_pos = match transcript.strand {
        Strand::Plus => pos + shift,
        Strand::Minus => pos - shift,
    };
    let shifted_bases =
        crate::normalize::rotate_inserted_bases(inserted_bases, shift, transcript.strand);
    let hgvs = hgvs_c::format_insertion_hgvs(
        shifted_pos,
        &shifted_bases,
        chrom,
        transcript,
        locate_index,
        fasta,
    )?;
    // VEP writes HGVS p. at the shifted insertion point (`hgvs_protein`
    // translates the shifted allele), so it needs that point's region.
    let shifted_region = if shift == 0 {
        location.region.clone()
    } else {
        locate_indel(chrom, shifted_pos, shifted_pos, transcript, locate_index)?.region
    };

    if let Some(intron) = frameshift_intron_endpoint(pos, pos, transcript, locate_index) {
        let mut result = super::complex::frameshift_intron_indel_result(
            pos,
            pos,
            intron,
            &location,
            transcript,
            locate_index,
        );
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    let mut result = match &location.region {
        // The intron terms are added below from VEP's predicates: a
        // canonical-site insertion is not `intronic`.
        IndelRegion::Intron => {
            let mut result = build_base_result(transcript, Vec::new());
            result.intron = location
                .intron_index
                .map(|idx| format_intron_number(idx, transcript.exon_count));
            Ok(result)
        }
        IndelRegion::Cds {
            cds_offset_start, ..
        } => {
            let ins_len = inserted_bases.len() as u32;
            // A junction insertion has no exon number; the placeholder is
            // cleared below.
            let exon_index = location.exon_index.unwrap_or(0);

            match partial_codon_insertion_site(*cds_offset_start, locate_index) {
                Some(PartialCodonInsertion::Inside) => partial_terminal_codon_indel_result(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_start,
                    inserted_bases,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                ),
                Some(PartialCodonInsertion::FivePrimeEdge(edit)) => {
                    partial_codon_edge_insertion_result(
                        edit,
                        inserted_bases,
                        transcript,
                        locate_index,
                        exon_index,
                    )
                }
                None if !ins_len.is_multiple_of(3) => annotate_cds_frameshift(
                    chrom,
                    *cds_offset_start,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    false,
                    *cds_offset_start,
                    Some(inserted_bases),
                    None,
                    0, // Pattern G deferred: applying the genomic shift
                       // as a CDS offset shift causes regressions (~60 cases
                       // where the shifted CDS position produces a different
                       // reading frame than VEP expects). Root cause: VEP's
                       // internal protein-level shift differs from the
                       // DNA-level shift. Needs per-variant investigation.
                ),
                None => annotate_cds_inframe_insertion(
                    chrom,
                    *cds_offset_start,
                    inserted_bases,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    false,
                    match (shift, &shifted_region) {
                        (
                            1..,
                            IndelRegion::Cds {
                                cds_offset_start, ..
                            },
                        ) => Some(ShiftedCdsInsertion {
                            cds_offset: *cds_offset_start,
                            inserted_bases: &shifted_bases,
                        }),
                        _ => None,
                    },
                ),
            }
        }
        _ => build_noncds_indel_result(transcript, &location, pos, pos),
    }?;
    // A coding insertion shifted out of the CDS has no HGVS p. in VEP 116,
    // whose `hgvs_protein` needs translation coordinates at the shifted
    // point: into the intron (NM_004341.5 c.4728+1_4728+2dup) or past the
    // stop codon (NM_003036.4 c.2152_*3dup). The consequences stay on the
    // unshifted insertion.
    if matches!(location.region, IndelRegion::Cds { .. })
        && !matches!(shifted_region, IndelRegion::Cds { .. })
    {
        result.hgvs_p = None;
    }

    for (hit, term) in [
        (effects.donor.is_some(), Consequence::SpliceDonorVariant),
        (
            effects.acceptor.is_some(),
            Consequence::SpliceAcceptorVariant,
        ),
        (effects.splice_region, Consequence::SpliceRegionVariant),
        (effects.intronic.is_some(), Consequence::IntronVariant),
    ] {
        if hit && !result.consequences.contains(&term) {
            result.consequences.push(term);
        }
    }
    result.impact = finalize_consequences(&mut result.consequences);
    if location.exon_index.is_none() {
        result.exon = None;
    }
    result.hgvs_c = hgvs;
    Ok(result)
}

/// Build the result for a CDS frameshift (shared by deletion, insertion,
/// and complex delins).
///
/// Reports `protein_start = protein_end = first_affected_codon` (1-based).
/// If the frameshift fully removes the start codon (codon 1), reports
/// `start_lost` instead. If it partially overlaps codon 1, reports
/// `frameshift_variant` only (VEP behavior).
#[allow(clippy::too_many_arguments)]
pub(super) fn annotate_cds_frameshift(
    chrom: &str,
    cds_offset_start: u32,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
    is_splice_region: bool,
    cds_offset_end: u32,
    inserted_bases: Option<&[u8]>,
    delins_alt: Option<&[u8]>,
    dna_shift: u32,
) -> Result<ConsequenceResult, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";
    let is_insertion = cds_offset_end == cds_offset_start;

    // For insertions at a codon boundary, the affected protein position is
    // the preceding codon (the one being disrupted), matching VEP's
    // `genomic2pep()` convention. Non-boundary insertions and deletions are
    // unaffected because integer division floors to the same codon.
    let codon_number = if is_insertion && cds_offset_start > 0 && cds_offset_start.is_multiple_of(3)
    {
        cds_offset_start / 3
    } else {
        cds_offset_start / 3 + 1
    };

    // VEP's start_lost() checks _overlaps_start_codon() AND
    // _ins_del_start_altered() (methionine actually changed).
    let overlaps_start_codon = codon_number == 1 && has_complete_start_codon(index);
    let start_lost = if overlaps_start_codon {
        if cds_offset_start == 0 && cds_offset_end >= 3 {
            true
        } else if is_insertion && cds_offset_start < 3 {
            // Insertion within codon 1 shifts reading frame, destroying start
            true
        } else if !is_insertion && cds_offset_start < 3 {
            // Partial deletion within codon 1 -- check if Met is destroyed
            let ref_codon = fetch_ref_codon(0, chrom, transcript, index, fasta)?;
            let del_start = cds_offset_start as usize;
            let del_end = (cds_offset_end.min(3)) as usize;
            let mut remaining: Vec<u8> = Vec::with_capacity(3);
            remaining.extend_from_slice(&ref_codon[..del_start]);
            remaining.extend_from_slice(&ref_codon[del_end..]);
            remaining.len() < 3
                || translate_codon_for_transcript(
                    remaining[..3]
                        .try_into()
                        .expect("length >= 3 guaranteed by short-circuit above"),
                    is_mito,
                ) != b'M'
        } else {
            false
        }
    } else {
        false
    };

    let predicts_nmd = !start_lost && super::nmd::predicts_nmd(cds_offset_start + 1, index);
    // VEP's `frameshift` and `start_lost` predicates fire independently —
    // there is no suppression of frameshift when the start codon is also
    // destroyed. Emit both. `finalize_consequences` (consequence/helpers.rs)
    // sorts by severity_rank so the insertion order here doesn't matter.
    let mut consequences = vec![Consequence::FrameshiftVariant];
    if start_lost {
        consequences.push(Consequence::StartLost);
    }

    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }

    let codon_start = (cds_offset_start / 3) * 3;

    // Compute coding-strand inserted bases for hgvs_p frameshift formatter.
    let coding_inserted_for_hgvs: Vec<u8> = if let Some(bases) = inserted_bases {
        match transcript.strand {
            Strand::Plus => bases.to_vec(),
            Strand::Minus => reverse_complement(bases),
        }
    } else if let Some(alt_bases) = delins_alt {
        match transcript.strand {
            Strand::Plus => alt_bases.to_vec(),
            Strand::Minus => reverse_complement(alt_bases),
        }
    } else {
        Vec::new() // pure deletion
    };

    // Insertions: single-codon context with "-/BASE" format.
    // Deletions: full codon-aligned range matching VEP's multi-codon display.
    let (codons, amino_acids, protein_end) = if is_insertion {
        let ref_codon = fetch_ref_codon(codon_start, chrom, transcript, index, fasta)?;
        let ref_aa = translate_codon_for_transcript(&ref_codon, is_mito);
        let codons = inserted_bases.map(|bases| {
            let coding_bases = match transcript.strand {
                Strand::Plus => bases.to_vec(),
                Strand::Minus => reverse_complement(bases),
            };
            let ins_str: String = coding_bases.iter().map(|&b| b as char).collect();
            format!("-/{ins_str}")
        });
        let aa_char = if ref_aa == b'X' { 'X' } else { ref_aa as char };
        (codons, format!("{aa_char}/X"), codon_number)
    } else {
        // Expand to codon boundaries
        let codon_end = ((cds_offset_end - 1) / 3 + 1) * 3;
        let ref_seq = fetch_cds_sequence(codon_start, codon_end, chrom, transcript, index, fasta)?;

        // Build alt by removing deleted bases (and inserting replacement
        // bases for delins variants)
        let local_del_start = (cds_offset_start - codon_start) as usize;
        let local_del_end = (cds_offset_end - codon_start) as usize;
        let replacement_len = delins_alt.map_or(0, |b| b.len());
        let mut alt_seq =
            Vec::with_capacity(ref_seq.len() - (local_del_end - local_del_start) + replacement_len);
        alt_seq.extend_from_slice(&ref_seq[..local_del_start]);
        if let Some(alt_bases) = delins_alt {
            let coding_alt = match transcript.strand {
                Strand::Plus => alt_bases.to_vec(),
                Strand::Minus => reverse_complement(alt_bases),
            };
            alt_seq.extend_from_slice(&coding_alt);
        }
        alt_seq.extend_from_slice(&ref_seq[local_del_end..]);

        // Translate both. Truncate alt to complete codons -- the trailing
        // 1-2 bases are the start of the shifted reading frame.
        let ref_aas = translate_sequence(&ref_seq, is_mito)?;
        let translate_len = alt_seq.len() - (alt_seq.len() % 3);
        let mut alt_aas = translate_sequence(&alt_seq[..translate_len], is_mito)?;
        alt_aas.push(b'X'); // frameshift indicator

        let codons_str = format_codons_indel(&ref_seq, &alt_seq, local_del_start, local_del_end);
        let aa_str = format_amino_acids_indel(&ref_aas, &alt_aas);
        (Some(codons_str), aa_str, codon_end / 3)
    };

    let impact = finalize_consequences(&mut consequences);

    // Deletions: 1-based inclusive [first_deleted, last_deleted]
    // Insertions: flanking positions [before, after] (VEP convention)
    let (cds_pos, cds_pos_end, cdna_position, cdna_end) = if is_insertion {
        let cdna_after = compute_cdna_position_for_cds(cds_offset_start, index);
        let cdna_before = if cds_offset_start > 0 {
            compute_cdna_position_for_cds(cds_offset_start - 1, index)
        } else {
            cdna_after
        };
        (
            cds_offset_start,
            cds_offset_start + 1,
            cdna_before.min(cdna_after),
            cdna_before.max(cdna_after),
        )
    } else {
        let cdna_start = compute_cdna_position_for_cds(cds_offset_start, index);
        let last_offset = cds_offset_end - 1;
        let cdna_end = compute_cdna_position_for_cds(last_offset, index);
        (cds_offset_start + 1, last_offset + 1, cdna_start, cdna_end)
    };

    // Generate HGVS protein notation for the frameshift, and with it the
    // termination codon the formatter locates on the way.
    //
    // When the start codon is destroyed there is no defined reading frame, so
    // the formatter is skipped and the termination codon is genuinely
    // undetermined -- consistent with `predicts_nmd` being forced false above.
    let (hgvs_p, ptc) = if start_lost {
        (Some("p.Met1?".to_string()), super::PtcStatus::Indeterminate)
    } else {
        let frameshift = crate::hgvs_p::format_hgvs_p_frameshift(
            cds_offset_start,
            cds_offset_end,
            &coding_inserted_for_hgvs,
            chrom,
            transcript,
            index,
            fasta,
            dna_shift,
        )?;
        let ptc = super::nmd::resolve_ptc(frameshift.ptc, &frameshift.edit, index);
        (frameshift.notation, ptc)
    };

    Ok(ConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        protein_accession: transcript.protein_accession.clone(),
        consequences,
        impact,
        protein_start: Some(codon_number),
        protein_end: Some(protein_end),
        codons,
        amino_acids: Some(amino_acids),
        exon: Some(format_exon_number(exon_index, transcript.exon_count)),
        intron: None,
        cds_position: Some(cds_pos),
        cds_position_end: Some(cds_pos_end),
        cdna_position: Some(cdna_position),
        cdna_position_end: Some(cdna_end),
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
        is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
        is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
        hgvs_c: None,
        hgvs_p,
        predicts_nmd,
        ptc,
    })
}

/// Build the result for an inframe CDS deletion.
///
/// Expands to codon boundaries, fetches ref codons, removes deleted bases,
/// retranslates, and compares amino acids.
#[allow(clippy::too_many_arguments)]
fn annotate_cds_inframe_deletion(
    chrom: &str,
    cds_offset_start: u32,
    cds_offset_end: u32,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
    is_splice_region: bool,
) -> Result<ConsequenceResult, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";

    let total_cds = index.total_cds_length();
    let codon_start = (cds_offset_start / 3) * 3;
    let codon_end = ((cds_offset_end - 1) / 3 + 1) * 3;

    let ref_seq = fetch_cds_sequence(codon_start, codon_end, chrom, transcript, index, fasta)?;

    let local_del_start = (cds_offset_start - codon_start) as usize;
    let local_del_end = (cds_offset_end - codon_start) as usize;
    let mut alt_seq = Vec::with_capacity(ref_seq.len() - (local_del_end - local_del_start));
    alt_seq.extend_from_slice(&ref_seq[..local_del_start]);
    alt_seq.extend_from_slice(&ref_seq[local_del_end..]);

    let ref_aas = translate_sequence(&ref_seq, is_mito)?;
    let alt_aas = translate_sequence(&alt_seq, is_mito)?;

    let first_codon = codon_start / 3;
    let mut consequences = Vec::new();

    if first_codon == 0
        && has_complete_start_codon(index)
        && ((cds_offset_start == 0 && cds_offset_end >= 3) || alt_aas.first() != Some(&b'M'))
    {
        consequences.push(Consequence::StartLost);
    }

    // The `!ref_aas.contains(b'*')` half is load-bearing, not symmetry: an
    // inframe deletion overlapping the reference terminator can fuse the
    // flanking bases into a fresh stop at essentially the same place, which is
    // not a premature termination. VEP guards it the same way -- its
    // `stop_gained` predicate is `($alt_pep =~ /\*/) and ($ref_pep !~ /\*/)`
    // (`ensembl-variation`, `VariationEffect.pm::stop_gained`).
    if alt_aas.contains(&b'*') && !ref_aas.contains(&b'*') {
        consequences.push(Consequence::StopGained);
    }

    if ref_aas.contains(&b'*') && !alt_aas.contains(&b'*') {
        consequences.push(Consequence::StopLost);
    }

    consequences.push(Consequence::InframeDeletion);

    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }

    let impact = finalize_consequences(&mut consequences);
    let protein_start = first_codon + 1;
    let protein_end = codon_end / 3;
    let cdna_start = compute_cdna_position_for_cds(cds_offset_start, index);
    let cdna_end = compute_cdna_position_for_cds(cds_offset_end - 1, index);

    // For the protein-level 3' rule, extend the AA window with
    // downstream context so deletions in amino acid repeats are shifted
    // to the most C-terminal position.
    let ext_codon_end = (codon_end + 90).min(total_cds);
    let (ref_aas_ext, alt_aas_ext) = if ext_codon_end > codon_end {
        let ext_seq =
            fetch_cds_sequence(codon_end, ext_codon_end, chrom, transcript, index, fasta)?;
        let ext_aas = translate_sequence(&ext_seq, is_mito)?;
        let mut r = ref_aas.clone();
        r.extend_from_slice(&ext_aas);
        let mut a = alt_aas.clone();
        a.extend_from_slice(&ext_aas);
        (r, a)
    } else {
        (ref_aas.clone(), alt_aas.clone())
    };

    let mut hgvs_p = crate::hgvs_p::format_hgvs_p_inframe_del(
        &ref_aas_ext,
        &alt_aas_ext,
        protein_start,
        &consequences,
    );

    // Stop_lost: compute extension distance via 3'UTR stop-scan.
    // format_hgvs_p_inframe_del returns None for StopLost because it lacks
    // FASTA access; override here with the actual extension notation.
    if hgvs_p.is_none() && consequences.contains(&Consequence::StopLost) {
        let stop_protein_pos = ref_aas
            .iter()
            .position(|&a| a == b'*')
            .map(|i| protein_start + i as u32)
            .unwrap_or(protein_end);
        hgvs_p = Some(crate::hgvs_p::format_hgvs_p_del_extension(
            stop_protein_pos,
            chrom,
            transcript,
            fasta,
        )?);
    }

    let predicts_nmd = consequences.contains(&Consequence::StopGained)
        && super::nmd::predicts_nmd(cds_offset_start + 1, index);

    // The new stop can sit anywhere in the retained window, so locate it
    // rather than assuming it coincides with the variant.
    let ptc = super::nmd::ptc_for_window(
        &consequences,
        &alt_aas,
        &ref_aas,
        first_codon,
        &crate::hgvs_p::CdsEdit {
            start: cds_offset_start,
            del_len: cds_offset_end - cds_offset_start,
            ins_len: 0,
        },
        index,
    );

    Ok(ConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        protein_accession: transcript.protein_accession.clone(),
        consequences,
        impact,
        protein_start: Some(protein_start),
        protein_end: Some(protein_end),
        codons: Some(format_codons_indel(
            &ref_seq,
            &alt_seq,
            local_del_start,
            local_del_end,
        )),
        amino_acids: Some(format_amino_acids_indel(&ref_aas, &alt_aas)),
        exon: Some(format_exon_number(exon_index, transcript.exon_count)),
        intron: None,
        cds_position: Some(cds_offset_start + 1),
        cds_position_end: Some(cds_offset_end),
        cdna_position: Some(cdna_start),
        cdna_position_end: Some(cdna_end),
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
        is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
        is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
        hgvs_c: None,
        hgvs_p,
        predicts_nmd,
        ptc,
    })
}

/// An insertion's CDS point and plus-strand bases after HGVS 3'
/// normalization.
struct ShiftedCdsInsertion<'a> {
    /// 0-based CDS offset of the base the shifted insertion precedes.
    cds_offset: u32,
    /// Plus-strand inserted bases, rotated by the shift.
    inserted_bases: &'a [u8],
}

/// Build the result for an inframe CDS insertion.
///
/// Expands the insertion-point codon, splices in inserted bases
/// (reverse-complemented for minus-strand), retranslates, and compares.
///
/// `shifted` is the insertion after HGVS 3' normalization when that moved it
/// to another CDS point: the HGVS p. notation is computed there, with the
/// rotated bases, so protein-level dup detection works at the shifted
/// position (VEP `hgvs_protein`). All other result fields remain at the
/// original `cds_offset`.
#[allow(clippy::too_many_arguments)]
fn annotate_cds_inframe_insertion(
    chrom: &str,
    cds_offset: u32,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
    is_splice_region: bool,
    shifted: Option<ShiftedCdsInsertion<'_>>,
) -> Result<ConsequenceResult, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";

    // Expand the codon window to include the preceding codon when the
    // insertion is at a codon boundary (cds_offset % 3 == 0). VEP uses
    // both flanking codons for the HGVS p. insertion notation; without
    // this expansion, the left flanking AA would be off by one position.
    let at_codon_boundary = cds_offset.is_multiple_of(3) && cds_offset >= 3;
    let codon_start = if at_codon_boundary {
        (cds_offset / 3 - 1) * 3
    } else {
        (cds_offset / 3) * 3
    };
    let codon_end = codon_start + if at_codon_boundary { 6 } else { 3 };
    let codon_end = codon_end.min(index.total_cds_length());

    let ref_seq = fetch_cds_sequence(codon_start, codon_end, chrom, transcript, index, fasta)?;

    let local_ins_pos = (cds_offset - codon_start) as usize;
    let coding_inserted = match transcript.strand {
        Strand::Plus => inserted_bases.to_vec(),
        Strand::Minus => reverse_complement(inserted_bases),
    };
    let mut alt_seq = Vec::with_capacity(ref_seq.len() + coding_inserted.len());
    alt_seq.extend_from_slice(&ref_seq[..local_ins_pos]);
    alt_seq.extend_from_slice(&coding_inserted);
    alt_seq.extend_from_slice(&ref_seq[local_ins_pos..]);

    let ref_aas = translate_sequence(&ref_seq, is_mito)?;
    let alt_aas = translate_sequence(&alt_seq, is_mito)?;

    let first_codon = codon_start / 3;
    let mut consequences = Vec::new();

    if first_codon == 0 && has_complete_start_codon(index) && alt_aas.first() != Some(&b'M') {
        consequences.push(Consequence::StartLost);
    }

    // Premature stop in alt (only if ref had no stop)
    if alt_aas.contains(&b'*') && !ref_aas.contains(&b'*') {
        consequences.push(Consequence::StopGained);
    }

    if ref_aas.contains(&b'*') && !alt_aas.contains(&b'*') {
        consequences.push(Consequence::StopLost);
    }

    if !consequences
        .iter()
        .any(|c| matches!(c, Consequence::StartLost))
    {
        consequences.push(Consequence::InframeInsertion);
    }

    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }

    let impact = finalize_consequences(&mut consequences);
    let protein_pos = first_codon + 1;
    // VEP reports flanking positions for insertions
    let cdna_after = compute_cdna_position_for_cds(cds_offset, index);
    let cdna_before = if cds_offset > 0 {
        compute_cdna_position_for_cds(cds_offset - 1, index)
    } else {
        cdna_after
    };

    // HGVS 3' normalization: when the insertion shifted to another CDS
    // point, recompute the codon window there with the rotated bases so
    // protein-level dup detection sees the shifted flanking amino acids. All
    // non-HGVS fields (consequences, protein_start, codons, amino_acids,
    // etc.) stay at the original position.
    //
    // The codon window must extend far enough upstream to include the
    // preceding amino acids that the dup detection compares against. For
    // an insertion of N amino acids, we need at least N preceding codons
    // (= 3*N CDS bases) before the shifted insertion point.
    //
    // A shift that leaves the CDS (into the 3' UTR, e.g. SKI
    // NM_003036.4:c.2152_*3dup, or an intron) has no HGVS p. in VEP; the
    // caller clears it, and the unshifted formatter only fills the field
    // meanwhile.
    let ins_aa_count = (coding_inserted.len() / 3) as u32;
    // HGVS p. at an insertion point, from a window reaching `ins_aa_count + 1`
    // codons upstream, so the residue before the insertion and any duplicated
    // residues are inside it. A window starting at a split codon would leave
    // the left flank out (VEP 116 NM_006231.4 c.2465_2467dup:
    // p.Lys822_Gly823insGlu).
    let format_at = |offset: u32, inserted: &[u8]| -> Result<Option<String>, VarEffectError> {
        let total = index.total_cds_length();
        let w_start = (offset / 3).saturating_sub(ins_aa_count + 1) * 3;
        let w_end = ((offset / 3 + 1) * 3).min(total);
        let w_ref = fetch_cds_sequence(w_start, w_end, chrom, transcript, index, fasta)?;
        let w_local = (offset - w_start) as usize;
        let mut w_alt = Vec::with_capacity(w_ref.len() + inserted.len());
        w_alt.extend_from_slice(&w_ref[..w_local]);
        w_alt.extend_from_slice(inserted);
        w_alt.extend_from_slice(&w_ref[w_local..]);
        let w_right_flank = if w_end < total {
            let next = fetch_ref_codon(w_end, chrom, transcript, index, fasta)?;
            translate_codon_for_transcript(&next, is_mito)
        } else {
            b'*'
        };
        Ok(crate::hgvs_p::format_hgvs_p_inframe_ins(
            &translate_sequence(&w_ref, is_mito)?,
            &translate_sequence(&w_alt, is_mito)?,
            w_start / 3 + 1,
            &consequences,
            w_right_flank,
        ))
    };
    let hgvs_p = if let Some(shifted) = shifted {
        let shifted_coding = match transcript.strand {
            Strand::Plus => shifted.inserted_bases.to_vec(),
            Strand::Minus => reverse_complement(shifted.inserted_bases),
        };
        format_at(shifted.cds_offset, &shifted_coding)?
    } else {
        // Without a DNA-level shift the insertion may still be a
        // protein-level duplication (codon degeneracy can produce identical
        // amino acids from different codons): probe the next codon boundary
        // and keep that result only when it is a duplication; otherwise the
        // insertion point keeps the mid-codon amino acid identity.
        let aligned_cds = (cds_offset.div_ceil(3) * 3).min(index.total_cds_length());
        let aligned_dup = if aligned_cds != cds_offset {
            format_at(aligned_cds, &coding_inserted)?.filter(|p| p.contains("dup"))
        } else {
            None
        };
        match aligned_dup {
            Some(dup) => Some(dup),
            None => format_at(cds_offset, &coding_inserted)?,
        }
    };

    let predicts_nmd = consequences.contains(&Consequence::StopGained)
        && super::nmd::predicts_nmd(cds_offset + 1, index);

    // An in-frame insertion can carry its own stop codon, so the termination
    // codon may lie inside the inserted sequence rather than at the variant.
    let ptc = super::nmd::ptc_for_window(
        &consequences,
        &alt_aas,
        &ref_aas,
        first_codon,
        &crate::hgvs_p::CdsEdit {
            start: cds_offset,
            del_len: 0,
            ins_len: inserted_bases.len() as u32,
        },
        index,
    );

    Ok(ConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        protein_accession: transcript.protein_accession.clone(),
        consequences,
        impact,
        protein_start: Some(protein_pos),
        protein_end: Some(protein_pos),
        codons: Some(format_codons_indel(
            &ref_seq,
            &alt_seq,
            local_ins_pos,
            local_ins_pos,
        )),
        amino_acids: Some(format_amino_acids_indel(&ref_aas, &alt_aas)),
        exon: Some(format_exon_number(exon_index, transcript.exon_count)),
        intron: None,
        cds_position: Some(cds_offset),
        cds_position_end: Some(cds_offset + 1),
        cdna_position: Some(cdna_before.min(cdna_after)),
        cdna_position_end: Some(cdna_before.max(cdna_after)),
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
        is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
        is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
        hgvs_c: None,
        hgvs_p,
        predicts_nmd,
        ptc,
    })
}

/// Build a splice consequence result for an indel that overlaps a canonical
/// splice site (+/-1-2 of an exon boundary).
///
/// Uses [`SpliceOverlapDetail`](crate::locate::SpliceOverlapDetail) to
/// distinguish donor vs acceptor. If a deletion overlaps both, both
/// consequences are reported. Secondary SO terms (`IntronVariant`,
/// `CodingSequenceVariant`, `SpliceRegionVariant`) are added when the
/// variant's genomic footprint also overlaps intron, CDS, or splice-region
/// bases — matching VEP's independent predicate evaluation.
///
/// # Arguments
///
/// * `var_start` — 0-based inclusive start of the variant's genomic footprint
/// * `var_end` — 0-based exclusive end of the variant's genomic footprint
///   (for insertions, use `(pos, pos)` since the footprint is zero-width)
pub(super) fn build_splice_indel_result(
    transcript: &TranscriptModel,
    location: &IndelLocation,
    var_start: u64,
    var_end: u64,
) -> Result<ConsequenceResult, VarEffectError> {
    // Invariant: `build_splice_indel_result` is only called when
    // `overlaps_splice_canonical` is true, which implies
    // `splice_detail.overlaps_donor || splice_detail.overlaps_acceptor`.
    // If that invariant is ever broken, surface it as an error rather than
    // silently emitting `SpliceDonorVariant` (which previously masked the
    // canonical-overlap calculation bug).
    let detail = location.splice_detail.as_ref().ok_or_else(|| {
        VarEffectError::Malformed(format!(
            "{}: build_splice_indel_result called without splice_detail",
            transcript.accession,
        ))
    })?;

    let mut consequences = Vec::with_capacity(4);
    if detail.overlaps_donor {
        consequences.push(Consequence::SpliceDonorVariant);
    }
    if detail.overlaps_acceptor {
        consequences.push(Consequence::SpliceAcceptorVariant);
    }

    if consequences.is_empty() {
        return Err(VarEffectError::Malformed(format!(
            "{}: splice canonical flag set but detail reports neither donor nor acceptor",
            transcript.accession,
        )));
    }

    // Secondary terms: VEP independently evaluates intron, CDS, and
    // splice-region overlap even when a canonical splice site is hit.
    let has_intron =
        matches!(location.region, IndelRegion::Intron) || location.crosses_exon_boundary;
    if has_intron {
        consequences.push(Consequence::IntronVariant);
    }

    // CDS overlap: check the variant's genomic footprint against the
    // transcript's CDS range. The `if let` naturally skips non-coding
    // transcripts (where cds_genomic_start/end are None).
    if let (Some(cds_start), Some(cds_end)) =
        (transcript.cds_genomic_start, transcript.cds_genomic_end)
        && var_start < cds_end
        && var_end > cds_start
    {
        consequences.push(Consequence::CodingSequenceVariant);
    }

    if location.overlaps_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }

    let impact = finalize_consequences(&mut consequences);
    let mut result = build_base_result(transcript, consequences);
    result.impact = impact;

    if let Some(&idx) = detail
        .donor_intron_indices
        .first()
        .or(detail.acceptor_intron_indices.first())
    {
        result.intron = Some(format_intron_number(idx, transcript.exon_count));
    } else if let Some(intron_idx) = location.intron_index {
        result.intron = Some(format_intron_number(intron_idx, transcript.exon_count));
    }
    if let Some(exon_idx) = location.exon_index {
        result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
    }
    Ok(result)
}

/// cDNA positions of the two bases either side of an exonic insertion
/// between `pos - 1` and `pos`, lower first.
///
/// At an exon/intron junction one flank is intronic; VEP reports the mRNA
/// point next to the exonic flank (VEP 116 NM_000546.6 `c.-28-1_-28insT`:
/// cDNA 114-115).
fn insertion_cdna_positions(pos: u64, transcript: &TranscriptModel) -> (Option<u32>, Option<u32>) {
    let right = compute_cdna_position_exonic(pos, transcript);
    let left = pos
        .checked_sub(1)
        .and_then(|p| compute_cdna_position_exonic(p, transcript));
    match (left, right) {
        (Some(a), Some(b)) => (Some(a.min(b)), Some(a.max(b))),
        (Some(c), None) | (None, Some(c)) => {
            // The exonic flank is the 5' one when it is the left base on the
            // plus strand or the right base on the minus strand.
            let exonic_is_five_prime = match transcript.strand {
                Strand::Plus => left.is_some(),
                Strand::Minus => right.is_some(),
            };
            if exonic_is_five_prime {
                (Some(c), Some(c + 1))
            } else {
                (Some(c.saturating_sub(1)), Some(c))
            }
        }
        (None, None) => (None, None),
    }
}

/// Build a non-CDS consequence result for an indel.
pub(super) fn build_noncds_indel_result(
    transcript: &TranscriptModel,
    location: &IndelLocation,
    start: u64,
    end: u64,
) -> Result<ConsequenceResult, VarEffectError> {
    let is_insertion = end == start;
    let mut consequences = match &location.region {
        IndelRegion::FivePrimeUtr => vec![Consequence::FivePrimeUtrVariant],
        IndelRegion::ThreePrimeUtr => vec![Consequence::ThreePrimeUtrVariant],
        IndelRegion::Intron => vec![Consequence::IntronVariant],
        // An insertion at an exon junction of a non-coding transcript is in
        // no exon (VEP 116 ENST00000424518.6 n.197_197+1insT:
        // `non_coding_transcript_variant`).
        IndelRegion::NonCodingExon if is_insertion && location.exon_index.is_none() => {
            vec![Consequence::NonCodingTranscriptVariant]
        }
        IndelRegion::NonCodingExon => vec![Consequence::NonCodingTranscriptExonVariant],
        IndelRegion::Upstream => vec![Consequence::UpstreamGeneVariant],
        IndelRegion::Downstream => vec![Consequence::DownstreamGeneVariant],
        _ => vec![Consequence::CodingSequenceVariant],
    };
    if location.overlaps_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }
    let mut result = build_base_result(transcript, consequences);
    if let Some(exon_idx) = location.exon_index {
        result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        if !is_insertion {
            result.cdna_position = compute_cdna_position_exonic(start, transcript);
            result.cdna_position_end = compute_cdna_position_exonic(end - 1, transcript);
        }
    }
    let exonic = matches!(
        location.region,
        IndelRegion::FivePrimeUtr | IndelRegion::ThreePrimeUtr | IndelRegion::NonCodingExon
    );
    if is_insertion && exonic {
        (result.cdna_position, result.cdna_position_end) =
            insertion_cdna_positions(start, transcript);
    }
    if let Some(intron_idx) = location.intron_index {
        result.intron = Some(format_intron_number(intron_idx, transcript.exon_count));
    }
    Ok(result)
}
