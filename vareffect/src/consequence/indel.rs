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
    exonic_splice_spill, frameshift_intron_endpoint, frameshift_junction_insertion_point,
    is_frameshift_junction_insertion,
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
    // An insertion is classified by the CDS offset of genomic base `pos`. On
    // the minus strand that base precedes the insertion point in coding
    // order, so the bases go in after it (VEP 116 ENST00000518351.2
    // c.544_545insG).
    let edit_start = if is_insertion && matches!(transcript.strand, Strand::Minus) {
        cds_offset_start + 1
    } else {
        cds_offset_start
    };
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
    /// After the last CDS base, which only the minus strand classifies as a
    /// CDS insertion.
    AfterCds,
}

/// Classify a CDS insertion against the incomplete terminal codon.
///
/// An insertion is classified by the CDS offset of genomic base `pos`; on
/// the minus strand that base precedes the insertion point in coding order,
/// so the coding insertion point is one later. `None` when the CDS has no
/// partial codon or the insertion is elsewhere.
pub(super) fn partial_codon_insertion_site(
    classified_offset: u32,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Option<PartialCodonInsertion> {
    let total = index.total_cds_length();
    if total.is_multiple_of(3) {
        return None;
    }
    let partial_start = total - total % 3;
    let insert_before = match transcript.strand {
        Strand::Plus => classified_offset,
        Strand::Minus => classified_offset + 1,
    };
    if insert_before >= total {
        Some(PartialCodonInsertion::AfterCds)
    } else if insert_before == partial_start && partial_start > 0 {
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
pub fn annotate_insertion(
    chrom: &str,
    pos: u64,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let mut result =
        annotate_insertion_core(chrom, pos, inserted_bases, transcript, locate_index, fasta)?;
    add_splice_spill(
        &mut result,
        exonic_splice_spill(pos, pos, transcript, locate_index),
    );
    // An insertion at the edge of a frameshift intron whose canonical windows
    // do not spill (gaps of 2 to 13 bases) also sits next to the gap's
    // intronic bases: VEP 116 ENST00000436261.7 c.1475_1475+1insT is
    // frameshift_variant + splice_region_variant + intron_variant.
    if is_frameshift_junction_insertion(pos, transcript, locate_index)
        && exonic_splice_spill(pos, pos, transcript, locate_index).is_none()
    {
        for term in [Consequence::SpliceRegionVariant, Consequence::IntronVariant] {
            if !result.consequences.contains(&term) {
                result.consequences.push(term);
            }
        }
        result.impact = finalize_consequences(&mut result.consequences);
    }
    Ok(result)
}

/// [`annotate_insertion`] without the splice terms a one-base frameshift intron spills
/// onto exon bases.
fn annotate_insertion_core(
    chrom: &str,
    pos: u64,
    inserted_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let location_pos = frameshift_junction_insertion_point(pos, transcript, locate_index);
    let location = locate_indel(chrom, location_pos, location_pos, transcript, locate_index)?;

    // HGVS 3' normalization: shift insertion to the most 3' equivalent
    // position on the coding strand for HGVS notation. At the shifted
    // position, is_duplication() in format_insertion_hgvs naturally re-checks
    // the 5' flanking — if the shift moved the insertion point to a repeat
    // boundary, dup is detected.
    let shift = if !location.crosses_exon_boundary {
        crate::normalize::compute_3prime_shift_insertion(
            chrom,
            pos,
            inserted_bases,
            transcript,
            fasta,
        )?
    } else {
        0
    };
    let shifted_pos = match transcript.strand {
        Strand::Plus => pos + shift,
        Strand::Minus => pos - shift,
    };
    let hgvs = hgvs_c::format_insertion_hgvs(
        shifted_pos,
        inserted_bases,
        chrom,
        transcript,
        locate_index,
        fasta,
    )?;

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

    if location.overlaps_splice_canonical
        && exonic_splice_spill(pos, pos, transcript, locate_index).is_none()
    {
        let mut result = build_splice_indel_result(transcript, &location, pos, pos)?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    // Boundary crossing: insertions at exact exon/intron junctions are rare.
    if location.crosses_exon_boundary {
        let mut consequences = vec![Consequence::IntronVariant];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        if let Some(intron_idx) = location.intron_index {
            result.intron = Some(format_intron_number(intron_idx, transcript.exon_count));
        }
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    let mut result = match &location.region {
        IndelRegion::Cds {
            cds_offset_start, ..
        } => {
            let ins_len = inserted_bases.len() as u32;
            let exon_index = location.exon_index.unwrap_or(0);

            match partial_codon_insertion_site(*cds_offset_start, transcript, locate_index) {
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
                Some(PartialCodonInsertion::AfterCds) => {
                    // Only the minus strand gets here, where the 3'UTR lies
                    // below the CDS. Without one the insertion is past the
                    // transcript end, as on the plus strand (VEP 116
                    // ENST00000518351.2: `downstream_gene_variant`; with a
                    // 3'UTR, ENST00000654733.1 c.*1dup:
                    // `3_prime_UTR_variant`).
                    let has_utr3 = transcript
                        .cds_genomic_start
                        .is_some_and(|cds_start| transcript.tx_start < cds_start);
                    if !has_utr3 {
                        return Ok(build_base_result(
                            transcript,
                            vec![Consequence::DownstreamGeneVariant],
                        ));
                    }
                    let mut utr_location = location.clone();
                    utr_location.region = IndelRegion::ThreePrimeUtr;
                    build_noncds_indel_result(transcript, &utr_location, pos, pos)
                }
                None if !ins_len.is_multiple_of(3) => annotate_cds_frameshift(
                    chrom,
                    *cds_offset_start,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
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
                    location.overlaps_splice_region,
                    shift as u32,
                ),
            }
        }
        _ => build_noncds_indel_result(transcript, &location, pos, pos),
    }?;

    // An insertion between an exon and a frameshift intron has no exon
    // number in VEP 116 (ENST00000400999.7 c.228_228+1insA).
    if is_frameshift_junction_insertion(pos, transcript, locate_index) {
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

/// Build the result for an inframe CDS insertion.
///
/// Expands the insertion-point codon, splices in inserted bases
/// (reverse-complemented for minus-strand), retranslates, and compares.
///
/// `dna_shift` is the HGVS 3' normalization shift in genomic bases. When
/// `> 0`, the HGVS p. notation is recomputed at the shifted CDS offset so
/// protein-level dup detection works at the shifted position. All other
/// result fields remain at the original `cds_offset`.
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
    dna_shift: u32,
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

    // Fetch the right flanking AA (next codon) for HGVS insertion notation.
    let right_flanking_aa = if codon_end < index.total_cds_length() {
        let next_codon = fetch_ref_codon(codon_end, chrom, transcript, index, fasta)?;
        translate_codon_for_transcript(&next_codon, is_mito)
    } else {
        b'*' // insertion before stop codon
    };

    // HGVS 3' normalization: when the DNA-level shift is non-zero and the
    // shifted CDS offset stays strictly within the CDS, recompute the codon
    // window at the shifted CDS offset so protein-level dup detection sees
    // the shifted flanking amino acids. All non-HGVS fields (consequences,
    // protein_start, codons, amino_acids, etc.) stay at the original
    // position.
    //
    // The codon window must extend far enough upstream to include the
    // preceding amino acids that the dup detection compares against. For
    // an insertion of N amino acids, we need at least N preceding codons
    // (= 3*N CDS bases) before the shifted insertion point.
    //
    // If the 3' shift reaches or crosses the stop codon (insertion
    // normalizes at/into the 3' UTR, e.g. HGVS c.*Ndup for last-exon
    // duplications like SKI NM_003036.4:c.2152_*3dup), the shifted window
    // is not meaningful at the protein level: the inserted codons would be
    // placed after the translated stop, producing a bogus
    // stop-straddling amino acid sequence. Fall back to the unshifted
    // formatter in that case (HGVS c. still carries the UTR-shifted
    // notation via the outer `shifted_pos`). `saturating_add` guards
    // against a pathological u32 overflow before the bounds check runs.
    let shifted_cds = cds_offset.saturating_add(dna_shift);
    let hgvs_p = if dna_shift > 0 && shifted_cds < index.total_cds_length() {
        // Extend the codon window upstream by the inserted AA count so
        // format_hgvs_p_inframe_ins can detect duplications at the
        // shifted position.
        let ins_aa_count = (coding_inserted.len() / 3) as u32;
        let upstream_codons = ins_aa_count + 1; // +1 for the boundary codon
        let s_codon_start = (shifted_cds / 3).saturating_sub(upstream_codons) * 3;
        let s_codon_end = ((shifted_cds / 3 + 1) * 3).min(index.total_cds_length());
        let s_ref =
            fetch_cds_sequence(s_codon_start, s_codon_end, chrom, transcript, index, fasta)?;
        let s_local = (shifted_cds - s_codon_start) as usize;
        let mut s_alt = Vec::with_capacity(s_ref.len() + coding_inserted.len());
        s_alt.extend_from_slice(&s_ref[..s_local]);
        s_alt.extend_from_slice(&coding_inserted);
        s_alt.extend_from_slice(&s_ref[s_local..]);
        let s_ref_aas = translate_sequence(&s_ref, is_mito)?;
        let s_alt_aas = translate_sequence(&s_alt, is_mito)?;
        let s_protein_pos = s_codon_start / 3 + 1;
        let s_right_flank = if s_codon_end < index.total_cds_length() {
            let nc = fetch_ref_codon(s_codon_end, chrom, transcript, index, fasta)?;
            translate_codon_for_transcript(&nc, is_mito)
        } else {
            b'*'
        };
        crate::hgvs_p::format_hgvs_p_inframe_ins(
            &s_ref_aas,
            &s_alt_aas,
            s_protein_pos,
            &consequences,
            s_right_flank,
        )
    } else {
        // When DNA-level shift is 0, the insertion may still be a
        // protein-level duplication (codon degeneracy can produce
        // identical AAs from different DNA codons). Try dup detection
        // at the codon-aligned position first; fall back to the
        // original window for standard insertion notation.
        let ins_aa_count = (coding_inserted.len() / 3) as u32;
        let upstream_codons = ins_aa_count + 1;

        // Dup probe: align the insertion to the next codon boundary
        // so the inserted bases translate cleanly. Only used when the
        // aligned result is actually a dup; non-dup results must use
        // the original (unaligned) position to preserve the correct
        // mid-codon amino acid identity.
        let aligned_cds = cds_offset.div_ceil(3) * 3;
        let aligned_cds = aligned_cds.min(index.total_cds_length());
        let aligned_dup = if aligned_cds != cds_offset {
            let a_start = (aligned_cds / 3).saturating_sub(upstream_codons) * 3;
            let a_end = ((aligned_cds / 3 + 1) * 3).min(index.total_cds_length());
            let a_ref = fetch_cds_sequence(a_start, a_end, chrom, transcript, index, fasta)?;
            let a_local = (aligned_cds - a_start) as usize;
            let mut a_alt = Vec::with_capacity(a_ref.len() + coding_inserted.len());
            a_alt.extend_from_slice(&a_ref[..a_local]);
            a_alt.extend_from_slice(&coding_inserted);
            a_alt.extend_from_slice(&a_ref[a_local..]);
            let a_ref_aas = translate_sequence(&a_ref, is_mito)?;
            let a_alt_aas = translate_sequence(&a_alt, is_mito)?;
            let a_pos = a_start / 3 + 1;
            let a_flank = if a_end < index.total_cds_length() {
                let nc = fetch_ref_codon(a_end, chrom, transcript, index, fasta)?;
                translate_codon_for_transcript(&nc, is_mito)
            } else {
                b'*'
            };
            let r = crate::hgvs_p::format_hgvs_p_inframe_ins(
                &a_ref_aas,
                &a_alt_aas,
                a_pos,
                &consequences,
                a_flank,
            );
            // Only adopt the aligned result when it detected a dup.
            match &r {
                Some(s) if s.contains("dup") => Some(r),
                _ => None,
            }
        } else {
            None
        };

        if let Some(dup_result) = aligned_dup {
            dup_result
        } else {
            // No dup detected — use the original window.
            crate::hgvs_p::format_hgvs_p_inframe_ins(
                &ref_aas,
                &alt_aas,
                protein_pos,
                &consequences,
                right_flanking_aa,
            )
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

/// Build a non-CDS consequence result for an indel.
pub(super) fn build_noncds_indel_result(
    transcript: &TranscriptModel,
    location: &IndelLocation,
    start: u64,
    end: u64,
) -> Result<ConsequenceResult, VarEffectError> {
    let mut consequences = match &location.region {
        IndelRegion::FivePrimeUtr => vec![Consequence::FivePrimeUtrVariant],
        IndelRegion::ThreePrimeUtr => vec![Consequence::ThreePrimeUtrVariant],
        IndelRegion::Intron => vec![Consequence::IntronVariant],
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
        let is_insertion = end == start;
        if is_insertion {
            // VEP reports flanking cDNA positions for insertions,
            // always with start <= end regardless of strand.
            let cdna_at = compute_cdna_position_exonic(start, transcript);
            let cdna_adj = if start > 0 {
                compute_cdna_position_exonic(start - 1, transcript)
            } else {
                cdna_at
            };
            // Both values are Option<u32>; use flatten for min/max.
            match (cdna_at, cdna_adj) {
                (Some(a), Some(b)) => {
                    result.cdna_position = Some(a.min(b));
                    result.cdna_position_end = Some(a.max(b));
                }
                (v @ Some(_), None) | (None, v @ Some(_)) => {
                    result.cdna_position = v;
                    result.cdna_position_end = v;
                }
                (None, None) => {}
            }
        } else {
            result.cdna_position = compute_cdna_position_exonic(start, transcript);
            result.cdna_position_end = compute_cdna_position_exonic(end - 1, transcript);
        }
    }
    if let Some(intron_idx) = location.intron_index {
        result.intron = Some(format_intron_number(intron_idx, transcript.exon_count));
    }
    Ok(result)
}
