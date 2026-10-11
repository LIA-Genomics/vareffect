//! Indel (insertion and deletion) consequence annotation.

use super::helpers::{
    AltSite, add_splice_spill, build_base_result, complete_codons_end,
    compute_cdna_position_exonic, compute_cdna_position_for_cds, ends_in_stop_codon,
    fetch_cds_sequence, fetch_ref_codon, finalize_consequences, has_complete_start_codon,
    is_incomplete_terminal_codon, without_gapped_numbering, without_unknown_alt_effects,
};
use super::{Consequence, ConsequenceResult};
use crate::codon::{
    format_amino_acids_indel, format_codons_indel, reverse_complement,
    translate_reference_sequence, translate_sequence,
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
    Ok(without_gapped_numbering(transcript, result))
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

    // HGVS 3' normalization: shift the deletion to the most 3' equivalent
    // position on the coding strand, on the genome as VEP 116 does, across
    // exon/intron junctions. The consequences stay on the unshifted deletion.
    let shift = crate::normalize::compute_3prime_shift_deletion(
        chrom,
        start,
        end,
        transcript.strand,
        fasta,
    )?;
    let (hgvs_start, hgvs_end) = match transcript.strand {
        Strand::Plus => (start + shift, end + shift),
        Strand::Minus => (start - shift, end - shift),
    };
    // A deletion with no HGVS c. at its own position keeps none: its shifted
    // range is not reproduced the way VEP writes it (deletions from the CDS
    // across a 5'UTR intron, NM_025150.5; VEP_DIVERGENCES.md 2.3).
    let hgvs = if shift > 0
        && hgvs_c::format_deletion_hgvs(chrom, start, end, transcript, locate_index)?.is_none()
    {
        None
    } else {
        hgvs_c::format_deletion_hgvs(chrom, hgvs_start, hgvs_end, transcript, locate_index)?
    };
    let shifted_region = if shift > 0 {
        Some(locate_indel(chrom, hgvs_start, hgvs_end, transcript, locate_index)?.region)
    } else {
        None
    };

    let mut result = 'annotate: {
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
            break 'annotate result;
        }

        if location.overlaps_splice_canonical
            && exonic_splice_spill(start, end, transcript, locate_index).is_none()
        {
            let mut result = build_splice_indel_result(transcript, &location, start, end)?;
            result.hgvs_c = hgvs;
            break 'annotate result;
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
            break 'annotate result;
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
            break 'annotate result;
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
                    indel_into_partial_codon_result(
                        chrom,
                        *cds_offset_start,
                        *cds_offset_end,
                        &[],
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
        result
    };

    // VEP writes HGVS p. from the shifted deletion (`hgvs_protein` translates
    // the shifted allele): none for a coding deletion shifted out of the CDS
    // (VEP 116 NM_000632.4 c.1497+1del), and the shifted deletion's for one
    // that reaches the CDS from an intron and shifts wholly into it
    // (NM_000540.3 c.958_960del `p.Glu320del`). A purely intronic deletion
    // shifted into an exon has none (`c.410-1del` written `c.410del`), nor
    // has one from the 5'UTR (VEP_DIVERGENCES.md 2.3).
    if let Some(shifted_region) = &shifted_region {
        let shifted_in_cds = matches!(shifted_region, IndelRegion::Cds { .. });
        if matches!(location.region, IndelRegion::Cds { .. }) {
            if !shifted_in_cds {
                result.hgvs_p = None;
            }
        } else if shifted_in_cds
            && location.crosses_exon_boundary
            && transcript
                .cds_segments
                .iter()
                .any(|segment| start < segment.genomic_end && segment.genomic_start < end)
        {
            result.hgvs_p = annotate_deletion_core(
                chrom,
                hgvs_start,
                hgvs_end,
                &[],
                transcript,
                locate_index,
                fasta,
            )?
            .hgvs_p;
        }
    }
    Ok(result)
}

/// Translate a coding sequence whose last codon may be incomplete: complete
/// codons as usual, a trailing partial codon as `X` unless the complete
/// codons are a lone stop (`TranscriptVariationAllele.pm::peptide`,
/// ensembl-variation release/116).
pub(super) fn translate_with_partial_codon(
    seq: &[u8],
    is_mito: bool,
) -> Result<Vec<u8>, VarEffectError> {
    translate_reference_with_partial_codon(seq, 1, &[], is_mito)
}

/// [`translate_with_partial_codon`] for a reference sequence starting at
/// 1-based codon `first_codon_number`, reading `TGA` at a selenocysteine
/// codon as `U` ([`translate_reference_sequence`]).
///
/// # Errors
///
/// Never in practice: the complete codons are a whole number of bases.
pub(super) fn translate_reference_with_partial_codon(
    seq: &[u8],
    first_codon_number: u32,
    selenocysteine_codons: &[u32],
    is_mito: bool,
) -> Result<Vec<u8>, VarEffectError> {
    let whole = seq.len() - seq.len() % 3;
    let mut aas = translate_reference_sequence(
        &seq[..whole],
        first_codon_number,
        selenocysteine_codons,
        is_mito,
    )?;
    if whole < seq.len() && aas != b"*" {
        aas.push(b'X');
    }
    Ok(aas)
}

/// Result for a deletion, deletion-insertion or MNV that starts in a complete
/// codon and ends in the incomplete terminal codon of a CDS, or at its end.
///
/// Codons and amino acids run from the first changed codon to the CDS end,
/// the partial codon read as `X`. HGVS p. is `p.?` (none for a deletion that
/// removes the CDS end). A stop created in a complete codon gets its PTC as
/// an SNV stop there does; a frameshift without one is
/// [`super::PtcStatus::Indeterminate`].
///
/// Terms follow VEP 116 (`VEP_DIVERGENCES.md` 3.1): a pure deletion is
/// `frameshift_variant` or `inframe_deletion`, or, removing the CDS end,
/// `incomplete_terminal_codon_variant` + `coding_sequence_variant`. A
/// deletion-insertion or MNV gets VEP's frameshift, `inframe_insertion`
/// (the alternate peptide starts or ends with the reference one),
/// `protein_altering_variant`, `stop_gained` and `missense_variant`
/// predicates, with `coding_sequence_variant` beside every result but a
/// frameshift or `protein_altering_variant` (VEP's `coding_unknown`).
///
/// # Arguments
///
/// * `cds_offset_start`, `cds_offset_end` - Padded CDS offsets of the
///   replaced bases, half-open
/// * `replacement` - The inserted bases in plus-strand genomic order (empty
///   for a pure deletion)
///
/// # Errors
///
/// Returns an error if the CDS bases cannot be read from the reference.
#[allow(clippy::too_many_arguments)]
pub(super) fn indel_into_partial_codon_result(
    chrom: &str,
    cds_offset_start: u32,
    cds_offset_end: u32,
    replacement: &[u8],
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
    let coding_replacement = match transcript.strand {
        Strand::Plus => replacement.to_vec(),
        Strand::Minus => reverse_complement(replacement),
    };
    let mut alt_seq = ref_seq[..local_start].to_vec();
    alt_seq.extend_from_slice(&coding_replacement);
    alt_seq.extend_from_slice(&ref_seq[local_end..]);

    let is_mito = transcript.chrom == "chrM";
    let ref_aas = translate_reference_with_partial_codon(
        &ref_seq,
        codon_start / 3 + 1,
        index.selenocysteine_codons(),
        is_mito,
    )?;
    let alt_aas = translate_with_partial_codon(&alt_seq, is_mito)?;

    let deleted = cds_offset_end - cds_offset_start;
    let mut consequences = if replacement.is_empty() {
        if cds_offset_end >= total {
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant,
            ]
        } else if deleted.is_multiple_of(3) {
            vec![Consequence::InframeDeletion]
        } else {
            vec![Consequence::FrameshiftVariant]
        }
    } else {
        let inserted = coding_replacement.len() as u32;
        let mut terms = Vec::new();
        if !inserted.abs_diff(deleted).is_multiple_of(3) {
            terms.push(Consequence::FrameshiftVariant);
        } else if inserted > deleted && alt_extends_ref(&ref_aas, &alt_aas) {
            terms.push(Consequence::InframeInsertion);
        } else if inserted != deleted
            && alt_aas.first() != Some(&b'*')
            && !alt_aas.starts_with(&ref_aas)
            && !alt_aas.ends_with(&ref_aas)
        {
            // VEP's `protein_altering_variant` matches the whole alternate
            // peptide, stop included: `c.140_142delinsATTAATTGA` (`LX/Y*LX`)
            // is `stop_gained` + `coding_sequence_variant`.
            terms.push(Consequence::ProteinAlteringVariant);
        }
        let stop_gained = alt_aas.contains(&b'*') && !ref_aas.contains(&b'*');
        if stop_gained {
            terms.push(Consequence::StopGained);
        } else if inserted == deleted && residue_changed(&ref_aas, &alt_aas) {
            terms.push(Consequence::MissenseVariant);
        }
        // VEP's `coding_unknown`: the reference peptide ends in `X`, so
        // `coding_sequence_variant` joins every term but frameshift and
        // protein_altering_variant.
        if !terms.iter().any(|t| {
            matches!(
                t,
                Consequence::FrameshiftVariant | Consequence::ProteinAlteringVariant
            )
        }) {
            terms.push(Consequence::CodingSequenceVariant);
        }
        terms
    };
    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }
    let impact = finalize_consequences(&mut consequences);
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
    // NM_001424184.1 c.656_659del: no HGVS p.); every other change reaching
    // the partial codon has an unknown one.
    result.hgvs_p = (!replacement.is_empty() || cds_offset_end < total).then(|| "p.?".to_string());
    // A new stop in the complete codons is located as for an SNV at that
    // codon; a frameshift without one runs past the annotated sequence.
    result.ptc = if result.consequences.contains(&Consequence::StopGained) {
        super::nmd::ptc_for_window(
            &result.consequences,
            &alt_aas,
            &ref_aas,
            codon_start / 3,
            &crate::hgvs_p::CdsEdit {
                start: cds_offset_start,
                del_len: deleted,
                ins_len: coding_replacement.len() as u32,
            },
            index,
        )
    } else if result
        .consequences
        .contains(&Consequence::FrameshiftVariant)
    {
        super::PtcStatus::Indeterminate
    } else {
        super::PtcStatus::NotApplicable
    };
    Ok(result)
}

/// `true` when the alternate peptide, cut after its first stop, starts or
/// ends with the reference peptide: VEP's `inframe_insertion` test for a
/// gain of whole codons (`VariationEffect.pm::inframe_insertion`, release
/// 116).
///
/// # Arguments
///
/// * `ref_aas`, `alt_aas` - Reference and alternate peptides
pub(super) fn alt_extends_ref(ref_aas: &[u8], alt_aas: &[u8]) -> bool {
    let alt = match alt_aas.iter().position(|&a| a == b'*') {
        Some(stop) => &alt_aas[..=stop],
        None => alt_aas,
    };
    alt.starts_with(ref_aas) || alt.ends_with(ref_aas)
}

/// `true` when two equal-length peptides differ at a residue both of them
/// know: the partial-codon `X` on either side is not a residue.
///
/// # Arguments
///
/// * `ref_aas`, `alt_aas` - Reference and alternate peptides
pub(super) fn residue_changed(ref_aas: &[u8], alt_aas: &[u8]) -> bool {
    ref_aas
        .iter()
        .zip(alt_aas)
        .any(|(r, a)| r != a && *r != b'X' && *a != b'X')
}

/// Result for a deletion or insertion that starts in the incomplete terminal
/// codon of a CDS whose length is not a whole number of codons.
///
/// VEP's `frameshift` and inframe predicates do not apply there
/// (`VariationEffect.pm::partial_codon`). Removing the whole partial codon is
/// `protein_altering_variant` + `incomplete_terminal_codon_variant`; anything
/// else is `incomplete_terminal_codon_variant` + `coding_sequence_variant`,
/// plus `stop_gained` when an insertion completes a stop. HGVS p. is `p.?`.
/// Examples: `VEP_DIVERGENCES.md` 3.1.
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
/// incomplete terminal codon: `inframe_insertion` +
/// `incomplete_terminal_codon_variant`, the inserted bases translated on
/// their own ([`translate_with_partial_codon`]), with `coding_sequence_variant`
/// for an `X` and `stop_gained` for a stop. Codons are `-/<inserted>`, the
/// positions the two flanking bases. A stop gets HGVS p. `p.?`; without one
/// there is no HGVS p. Examples: `VEP_DIVERGENCES.md` 3.1.
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
        return without_unknown_alt_effects(
            &AltSite {
                alt: inserted_bases,
                is_indel: true,
                chrom,
                transcript,
                index: locate_index,
                fasta,
            },
            without_gapped_numbering(transcript, result),
        );
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
    without_unknown_alt_effects(
        &AltSite {
            alt: inserted_bases,
            is_indel: true,
            chrom,
            transcript,
            index: locate_index,
            fasta,
        },
        without_gapped_numbering(transcript, result),
    )
}

/// Build the result for a CDS frameshift (shared by deletion, insertion,
/// and complex delins).
///
/// Reports `protein_start = protein_end = first_affected_codon` (1-based).
/// A frameshift touching a complete start codon follows VEP 116's
/// `start_lost` / `start_retained_variant` (see [`frameshift_start_terms`]).
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

    let mut consequences = vec![Consequence::FrameshiftVariant];
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

    // The codons and peptides of VEP's codon window
    // (`TranscriptVariationAllele::codon` / `peptide`, release/116): the
    // reference codons the variant touches and the same window with the
    // allele applied. An insertion between two codons touches none, so its
    // reference side is empty (`-/T`, `-/X`); one inside a codon shows that
    // codon (`gat/gaTt`, `D/DX`).
    let (codons, ref_aas, alt_aas, protein_end) = if is_insertion {
        let (ref_seq, alt_seq, local_start) = if cds_offset_start.is_multiple_of(3) {
            (Vec::new(), coding_inserted_for_hgvs.clone(), 0)
        } else {
            let ref_codon = fetch_ref_codon(codon_start, chrom, transcript, index, fasta)?;
            let local_start = (cds_offset_start - codon_start) as usize;
            let mut alt_seq = Vec::with_capacity(3 + coding_inserted_for_hgvs.len());
            alt_seq.extend_from_slice(&ref_codon[..local_start]);
            alt_seq.extend_from_slice(&coding_inserted_for_hgvs);
            alt_seq.extend_from_slice(&ref_codon[local_start..]);
            (ref_codon.to_vec(), alt_seq, local_start)
        };
        let codons = format_codons_indel(&ref_seq, &alt_seq, local_start, local_start);
        let ref_aas = translate_reference_sequence(
            &ref_seq,
            codon_start / 3 + 1,
            index.selenocysteine_codons(),
            is_mito,
        )?;
        let alt_aas = translate_with_partial_codon(&alt_seq, is_mito)?;
        (Some(codons), ref_aas, alt_aas, codon_number)
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

        let ref_aas = translate_reference_sequence(
            &ref_seq,
            codon_start / 3 + 1,
            index.selenocysteine_codons(),
            is_mito,
        )?;
        let alt_aas = translate_with_partial_codon(&alt_seq, is_mito)?;
        let codons_str = format_codons_indel(&ref_seq, &alt_seq, local_del_start, local_del_end);
        (Some(codons_str), ref_aas, alt_aas, codon_end / 3)
    };
    let amino_acids = format_amino_acids_indel(&ref_aas, &alt_aas);

    // VEP's `frameshift` and the start predicates fire independently: a
    // frameshift touching the start codon also gets `start_lost` and / or
    // `start_retained_variant`.
    // `_overlaps_start_codon` compares unshifted cDNA bounds; an insertion's
    // bounds are reversed (start = end + 1), so it overlaps only from inside
    // the codon (offsets 1 and 2), a deletion when it starts in codon 1.
    let overlaps_start_codon = has_complete_start_codon(transcript, index)
        && if is_insertion {
            (1..3).contains(&cds_offset_start)
        } else {
            cds_offset_start < 3
        };
    let (start_lost, start_retained) = if overlaps_start_codon {
        let start = check_start_codon(
            chrom,
            transcript,
            index,
            fasta,
            (cds_offset_start, cds_offset_end),
            &coding_inserted_for_hgvs,
        )?;
        frameshift_start_terms(start, &ref_aas, &alt_aas)
    } else {
        (false, false)
    };
    if start_lost {
        consequences.push(Consequence::StartLost);
    }
    if start_retained {
        consequences.push(Consequence::StartRetainedVariant);
    }
    let predicts_nmd = !start_lost && super::nmd::predicts_nmd(cds_offset_start + 1, index);
    // VEP's `stop_gained` holds for a frameshift too: the window's alternate
    // peptide reads a stop the reference peptide lacks (VEP 116
    // NM_013322.3 c.ins `frameshift_variant,stop_gained`). Not when an
    // insertion just before the stop codon opens with a stop, which then
    // sits where the reference stop was (VEP 116 NM_000546.6
    // c.1179_1180insTAAG: `inframe_insertion,stop_retained_variant`,
    // VEP_DIVERGENCES.md).
    let stop_in_place = is_insertion
        && cds_offset_start + 3 == index.total_cds_length()
        && alt_aas.first() == Some(&b'*');
    if alt_aas.contains(&b'*') && !ref_aas.contains(&b'*') && !stop_in_place {
        consequences.push(Consequence::StopGained);
    }

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
        is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
        is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
        is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
        hgvs_c: None,
        hgvs_p,
        predicts_nmd,
        ptc,
    })
}

/// What an indel does to the start codon (`VariationEffect.pm`,
/// release/116).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StartCheck {
    /// `_ins_del_start_altered` is false: the transcript has a 5'UTR and
    /// `ATG` still opens the CDS, or the edited 5'UTR + CDS still ends with
    /// the whole CDS (the indel only shifts bases across the CDS start,
    /// NM_000546.6 c.1dup).
    keeps_start: bool,
    /// `_inv_start_altered` is true: the transcript has a 5'UTR and the
    /// edited CDS no longer opens with `ATG`.
    start_moved: bool,
}

/// Check an indel at CDS offsets `start..end` (`start == end` for an
/// insertion), with `alt` (coding strand) in place of the removed bases,
/// against the start codon as VEP edits the 5'UTR plus CDS.
///
/// # Errors
///
/// Returns an error if the reference CDS or 5'UTR cannot be read.
fn check_start_codon(
    chrom: &str,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    (start, end): (u32, u32),
    alt: &[u8],
) -> Result<StartCheck, VarEffectError> {
    let (start, end) = (start as usize, end as usize);
    let window_end = (end as u32 + 3).min(index.total_cds_length());
    let cds = fetch_cds_sequence(0, window_end, chrom, transcript, index, fasta)?;
    let utr5_len = index.utr5_exonic_len() as usize;
    let mut edited = Vec::with_capacity(cds.len() + alt.len());
    edited.extend_from_slice(&cds[..start]);
    edited.extend_from_slice(alt);
    edited.extend_from_slice(&cds[end..]);
    let starts_atg = edited.len() >= 3 && edited[..3].eq_ignore_ascii_case(b"ATG");
    let start_moved = utr5_len > 0 && !starts_atg;
    if utr5_len > 0 && starts_atg {
        return Ok(StartCheck {
            keeps_start: true,
            start_moved,
        });
    }
    // The edited 5'UTR + CDS is shorter than the CDS alone.
    if utr5_len + alt.len() < end - start {
        return Ok(StartCheck {
            keeps_start: false,
            start_moved,
        });
    }
    // It ends with the CDS when the bases before `cds[end..]` (5'UTR tail,
    // kept CDS prefix, inserted bases) end with `cds[..end]`.
    let mut head = utr5_tail(end, transcript, index, fasta)?;
    head.extend_from_slice(&cds[..start]);
    head.extend_from_slice(alt);
    Ok(StartCheck {
        keeps_start: head.ends_with(&cds[..end]),
        start_moved,
    })
}

/// The last `n` exonic 5'UTR bases on the coding strand (fewer when the
/// 5'UTR is shorter).
///
/// # Errors
///
/// Returns an error if the reference cannot be read.
pub(super) fn utr5_tail(
    n: usize,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<Vec<u8>, VarEffectError> {
    let Some(first_cds) = transcript.cds_segments.first() else {
        return Ok(Vec::new());
    };
    if n == 0 || index.utr5_exonic_len() == 0 {
        return Ok(Vec::new());
    }
    let plus = matches!(transcript.strand, Strand::Plus);
    // 5'UTR exon pieces in transcript order.
    let pieces: Vec<(u64, u64)> = transcript
        .exons
        .iter()
        .filter_map(|exon| {
            let (piece_start, piece_end) = if plus {
                (
                    exon.genomic_start,
                    exon.genomic_end.min(first_cds.genomic_start),
                )
            } else {
                (
                    exon.genomic_start.max(first_cds.genomic_end),
                    exon.genomic_end,
                )
            };
            (piece_start < piece_end).then_some((piece_start, piece_end))
        })
        .collect();
    let mut tail: Vec<u8> = Vec::with_capacity(n);
    for &(piece_start, piece_end) in pieces.iter().rev() {
        let missing = (n - tail.len()) as u64;
        let (fetch_start, fetch_end) = if plus {
            (
                piece_end.saturating_sub(missing).max(piece_start),
                piece_end,
            )
        } else {
            (piece_start, (piece_start + missing).min(piece_end))
        };
        let chunk = fasta.fetch_sequence(&transcript.chrom, fetch_start, fetch_end)?;
        let chunk = if plus {
            chunk
        } else {
            reverse_complement(&chunk)
        };
        let mut joined = chunk;
        joined.extend_from_slice(&tail);
        tail = joined;
        if tail.len() >= n {
            break;
        }
    }
    Ok(tail)
}

/// VEP 116's start terms for a frameshift touching a complete start codon,
/// as `(start_lost, start_retained_variant)` (`VariationEffect.pm`,
/// release/116):
///
/// - `start_retained_variant` when the indel keeps the start codon
///   ([`StartCheck::keeps_start`]; NM_021922.3 c.3_12del rebuilds `ATG`);
/// - `start_lost` when it does not, when `ATG` no longer opens the CDS of a
///   transcript with a 5'UTR ([`StartCheck::start_moved`]), or when the
///   alternate peptide of the window (which begins at codon 1) neither
///   starts nor ends with the reference peptide and is not just `X`. Both
///   can hold: NM_000546.6 c.1dup and KRAS NM_004985.5 c.1del, which only
///   shift the start codon, are `start_lost,start_retained_variant`.
///
/// # Arguments
///
/// * `start` - The indel's effect on the start codon
/// * `ref_aas`, `alt_aas` - The window's peptides
fn frameshift_start_terms(start: StartCheck, ref_aas: &[u8], alt_aas: &[u8]) -> (bool, bool) {
    let peptide_lost = !ref_aas.is_empty()
        && !alt_aas.is_empty()
        && alt_aas != b"X"
        && !alt_aas.starts_with(ref_aas)
        && !alt_aas.ends_with(ref_aas);
    (
        !start.keeps_start || start.start_moved || peptide_lost,
        start.keeps_start,
    )
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

    let codon_start = (cds_offset_start / 3) * 3;
    let codon_end = ((cds_offset_end - 1) / 3 + 1) * 3;

    let ref_seq = fetch_cds_sequence(codon_start, codon_end, chrom, transcript, index, fasta)?;

    let local_del_start = (cds_offset_start - codon_start) as usize;
    let local_del_end = (cds_offset_end - codon_start) as usize;
    let mut alt_seq = Vec::with_capacity(ref_seq.len() - (local_del_end - local_del_start));
    alt_seq.extend_from_slice(&ref_seq[..local_del_start]);
    alt_seq.extend_from_slice(&ref_seq[local_del_end..]);

    let ref_aas = translate_reference_sequence(
        &ref_seq,
        codon_start / 3 + 1,
        index.selenocysteine_codons(),
        is_mito,
    )?;
    let alt_aas = translate_sequence(&alt_seq, is_mito)?;

    let first_codon = codon_start / 3;
    let mut consequences = Vec::new();

    if first_codon == 0
        && has_complete_start_codon(transcript, index)
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
    let ext_codon_end = (codon_end + 90).min(complete_codons_end(index));
    let (ref_aas_ext, alt_aas_ext) = if ext_codon_end > codon_end {
        let ext_seq =
            fetch_cds_sequence(codon_end, ext_codon_end, chrom, transcript, index, fasta)?;
        let ext_aas = translate_reference_sequence(
            &ext_seq,
            codon_end / 3 + 1,
            index.selenocysteine_codons(),
            is_mito,
        )?;
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
            index.has_selenocysteine_before(protein_start),
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
        is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
        is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
        is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
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

/// VEP's codon window of an in-frame insertion
/// (`TranscriptVariationAllele::codon`, release/116).
///
/// Between two codons VEP's `translation_end` is the codon before the
/// insertion and `translation_start` the codon after it, so the reference
/// codon is empty and the alternate is the inserted bases alone (VEP 116
/// NM_000059.4 CDS 42/43: `-/TTA`, `-/L`). Inside a codon both are that codon
/// (NM_000546.6 c.523_524insTTT: `cgc/cTTTgc`, `R/LC`).
#[derive(Debug)]
pub(super) struct InsertionCodons {
    /// Reference codon on the coding strand; empty between two codons.
    ref_codon: Vec<u8>,
    /// `ref_codon` with the inserted bases spliced in.
    alt_codon: Vec<u8>,
    /// Index in `ref_codon` the bases are inserted at.
    local: usize,
    /// 0-based index of the codon holding the insertion, or of the codon
    /// after it at a codon boundary.
    codon_index: u32,
    /// Translations of both codons and VEP's translation coordinates.
    pub(super) peptides: crate::hgvs_p::InsertionPeptides,
}

impl InsertionCodons {
    /// Build the window for `coding_inserted` (coding strand) inserted before
    /// CDS offset `cds_offset`.
    ///
    /// # Errors
    ///
    /// Returns an error if the reference codon cannot be read.
    fn new(
        cds_offset: u32,
        coding_inserted: &[u8],
        chrom: &str,
        transcript: &TranscriptModel,
        index: &LocateIndex,
        fasta: &FastaReader,
    ) -> Result<Self, VarEffectError> {
        let codon_index = cds_offset / 3;
        let local = (cds_offset % 3) as usize;
        let ref_codon = if local == 0 {
            Vec::new()
        } else {
            fetch_cds_sequence(
                codon_index * 3,
                codon_index * 3 + 3,
                chrom,
                transcript,
                index,
                fasta,
            )?
        };
        Self::from_codon(
            ref_codon,
            coding_inserted,
            local,
            codon_index,
            index.selenocysteine_codons(),
            transcript.chrom == "chrM",
        )
    }

    /// Build the window from the reference codon (empty between two codons).
    ///
    /// # Arguments
    ///
    /// * `ref_codon` -- coding-strand reference codon, empty when `local` is 0
    /// * `coding_inserted` -- inserted bases on the coding strand
    /// * `local` -- index in `ref_codon` the bases go before
    /// * `codon_index` -- 0-based codon holding the insertion, or after it
    /// * `selenocysteine_codons` -- 1-based selenocysteine codon numbers,
    ///   read as `U` in the reference codon only
    /// * `is_mito` -- translate with NCBI table 2
    ///
    /// # Errors
    ///
    /// Returns an error if a codon is not a whole number of bases.
    pub(super) fn from_codon(
        ref_codon: Vec<u8>,
        coding_inserted: &[u8],
        local: usize,
        codon_index: u32,
        selenocysteine_codons: &[u32],
        is_mito: bool,
    ) -> Result<Self, VarEffectError> {
        let mut alt_codon = Vec::with_capacity(ref_codon.len() + coding_inserted.len());
        alt_codon.extend_from_slice(&ref_codon[..local]);
        alt_codon.extend_from_slice(coding_inserted);
        alt_codon.extend_from_slice(&ref_codon[local..]);
        let peptides = crate::hgvs_p::InsertionPeptides {
            ref_aas: translate_reference_sequence(
                &ref_codon,
                codon_index + 1,
                selenocysteine_codons,
                is_mito,
            )?,
            alt_aas: translate_sequence(&alt_codon, is_mito)?,
            tl_start: codon_index + 1,
            tl_end: if local == 0 {
                codon_index
            } else {
                codon_index + 1
            },
        };
        Ok(Self {
            ref_codon,
            alt_codon,
            local,
            codon_index,
            peptides,
        })
    }
}

/// The reference protein of a transcript, translated from the genome on
/// demand.
pub(super) struct CdsProtein<'a> {
    chrom: &'a str,
    transcript: &'a TranscriptModel,
    index: &'a LocateIndex,
    fasta: &'a FastaReader,
    residue_count: u32,
}

impl<'a> CdsProtein<'a> {
    /// Reader for `transcript`'s protein; the residue count excludes a
    /// terminal stop codon, as VEP's `_peptide` does.
    ///
    /// # Errors
    ///
    /// Returns an error if the last codon cannot be read.
    pub(super) fn new(
        chrom: &'a str,
        transcript: &'a TranscriptModel,
        index: &'a LocateIndex,
        fasta: &'a FastaReader,
    ) -> Result<Self, VarEffectError> {
        let codons = index.total_cds_length() / 3;
        let residue_count = if ends_in_stop_codon(chrom, transcript, index, fasta)? {
            codons - 1
        } else {
            codons
        };
        Ok(Self {
            chrom,
            transcript,
            index,
            fasta,
            residue_count,
        })
    }
}

impl crate::hgvs_p::ReferenceProtein for CdsProtein<'_> {
    fn residue_count(&self) -> u32 {
        self.residue_count
    }

    fn residues(&self, start: u32, end: u32) -> Result<Vec<u8>, VarEffectError> {
        // The genome reader refuses an empty range.
        if start >= end {
            return Ok(Vec::new());
        }
        let seq = fetch_cds_sequence(
            start * 3,
            end * 3,
            self.chrom,
            self.transcript,
            self.index,
            self.fasta,
        )?;
        translate_reference_sequence(
            &seq,
            start + 1,
            self.index.selenocysteine_codons(),
            self.transcript.chrom == "chrM",
        )
    }
}

/// VEP 116 consequence terms of an in-frame CDS insertion
/// (`VariationEffect.pm`, release/116), without `splice_region_variant`.
///
/// - `inframe_insertion`: the alternate peptide, cut after its first stop,
///   starts or ends with the reference peptide (always between codons).
/// - `protein_altering_variant`: the uncut alternate peptide neither starts
///   nor ends with the reference peptide, and neither starts with a stop
///   (NM_000546.6 c.523_524insTTT, `R/LC`).
/// - `stop_retained_variant` (`ref_eq_alt_sequence`): the protein with the
///   alternate peptide spliced in is the reference followed by a stop
///   (NM_003036.4 c.2152_*3dup), or the stop stays at the same index of
///   both peptides (NM_000546.6 c.1180_1181insAAG, `*/*G`).
/// - `stop_gained` / `stop_lost`: a stop only in the alternate / reference
///   peptide, unless the stop is retained.
/// - Codon 1 (`start_lost`, `start_retained_variant`): VEP edits the 5'UTR
///   plus CDS; `_inv_start_altered` calls `start_lost` whenever a transcript
///   with a 5'UTR no longer has `ATG` at the start codon position
///   (NM_000546.6 c.1_2insCGA: `start_lost,start_retained_variant`;
///   c.2_3insGAT: `start_retained_variant` alone).
///
/// An inserted `N` (amino acid `X`) takes VEP's sequence-based stop checks,
/// which are not ported; no stop term is given then.
///
/// # Arguments
///
/// * `codons` -- VEP's codon window at the (unshifted) insertion point
/// * `inserted_len` -- number of inserted bases
/// * `complete_start` -- codon 1 is a complete start codon (no `cds_start_NF`)
/// * `has_utr5` -- the transcript has 5'UTR exonic bases
/// * `protein` -- reference protein of the transcript
///
/// # Errors
///
/// Returns an error if reference residues cannot be read.
pub(super) fn inframe_insertion_terms(
    codons: &InsertionCodons,
    inserted_len: usize,
    complete_start: bool,
    has_utr5: bool,
    protein: &impl crate::hgvs_p::ReferenceProtein,
) -> Result<Vec<Consequence>, VarEffectError> {
    let ref_aas = codons.peptides.ref_aas.as_slice();
    let alt_aas = codons.peptides.alt_aas.as_slice();

    // Only an insertion inside codon 1 overlaps the start codon.
    let overlaps_start = codons.codon_index == 0 && codons.local != 0 && complete_start;
    let start_is_atg = codons.alt_codon[..3].eq_ignore_ascii_case(b"ATG");
    let start_altered = overlaps_start
        && !(has_utr5 && start_is_atg)
        && codons.alt_codon[inserted_len..] != codons.ref_codon[..];
    let start_retained = overlaps_start && !start_altered;

    let alt_to_stop = match alt_aas.iter().position(|&aa| aa == b'*') {
        Some(stop) => &alt_aas[..=stop],
        None => alt_aas,
    };
    let inframe = !(start_retained && alt_aas.ends_with(ref_aas))
        && (alt_to_stop.starts_with(ref_aas) || alt_to_stop.ends_with(ref_aas));
    let start_lost = overlaps_start
        && ((start_altered && !inframe)
            || (has_utr5 && !start_is_atg)
            || (!ref_aas.is_empty()
                && !alt_aas.is_empty()
                && alt_aas != b"X"
                && !alt_aas.starts_with(ref_aas)
                && !alt_aas.ends_with(ref_aas)));

    let alt_unknown = alt_aas.contains(&b'X');
    let ref_stop = ref_aas.iter().position(|&aa| aa == b'*');
    let alt_stop = alt_aas.iter().position(|&aa| aa == b'*');
    let stop_lost = !alt_unknown && ref_stop.is_some() && alt_stop.is_none();
    let stop_retained = !stop_lost && !alt_aas.is_empty() && !alt_unknown && {
        // The spliced protein matches the reference for its first
        // `residue_count` residues and has a stop next: the stop must come
        // from the alternate peptide at index `k`.
        let first = codons.peptides.tl_start - 1;
        let k = protein.residue_count().saturating_sub(first) as usize;
        let elongated = alt_aas.get(k) == Some(&b'*')
            && alt_aas[..k] == protein.residues(first, first + k as u32)?[..];
        elongated || (ref_stop.is_some() && ref_stop == alt_stop)
    };
    let stop_gained = !stop_retained && !stop_lost && alt_stop.is_some() && ref_stop.is_none();
    let protein_altering = !start_lost
        && ref_aas.len() != alt_aas.len()
        && ref_aas.first() != Some(&b'*')
        && alt_aas.first() != Some(&b'*')
        && !alt_aas.starts_with(ref_aas)
        && !alt_aas.ends_with(ref_aas);

    Ok([
        (start_lost, Consequence::StartLost),
        (start_retained, Consequence::StartRetainedVariant),
        (stop_gained, Consequence::StopGained),
        (stop_lost, Consequence::StopLost),
        (stop_retained, Consequence::StopRetainedVariant),
        (inframe && !start_lost, Consequence::InframeInsertion),
        (protein_altering, Consequence::ProteinAlteringVariant),
    ]
    .into_iter()
    .filter_map(|(hit, term)| hit.then_some(term))
    .collect())
}

/// Build the result for an inframe CDS insertion.
///
/// Consequences, codons, amino acids and positions come from VEP's codon
/// window ([`InsertionCodons`]) at the original insertion point; HGVS p.
/// from the window at `shifted`, the insertion after HGVS 3' normalization,
/// when that moved it to another CDS point (VEP `hgvs_protein` translates the
/// shifted allele), placed against the whole protein.
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
    let coding = |bases: &[u8]| match transcript.strand {
        Strand::Plus => bases.to_vec(),
        Strand::Minus => reverse_complement(bases),
    };
    let coding_inserted = coding(inserted_bases);
    let codons = InsertionCodons::new(
        cds_offset,
        &coding_inserted,
        chrom,
        transcript,
        index,
        fasta,
    )?;
    let protein = CdsProtein::new(chrom, transcript, index, fasta)?;

    // VEP 116 gives an allele with a non-ACGT base no peptide: inside codon 1
    // it is `start_retained_variant` (its start checks skip ambiguous
    // sequence), elsewhere `coding_unknown` (`coding_sequence_variant`), with
    // no amino acids and no HGVS p. (every N-containing in-frame insertion in
    // the ClinVar sample checked against VEP REST). VEP's strand-dependent
    // stop-codon overlap for such alleles is not ported (VEP_DIVERGENCES.md).
    let unambiguous = coding_inserted
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'));
    let mut consequences = if unambiguous {
        inframe_insertion_terms(
            &codons,
            coding_inserted.len(),
            has_complete_start_codon(transcript, index),
            index.utr5_exonic_len() > 0,
            &protein,
        )?
    } else if codons.codon_index == 0
        && codons.local != 0
        && has_complete_start_codon(transcript, index)
    {
        vec![Consequence::StartRetainedVariant]
    } else {
        vec![Consequence::CodingSequenceVariant]
    };
    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }
    let impact = finalize_consequences(&mut consequences);

    // VEP reports flanking positions for insertions
    let cdna_after = compute_cdna_position_for_cds(cds_offset, index);
    let cdna_before = if cds_offset > 0 {
        compute_cdna_position_for_cds(cds_offset - 1, index)
    } else {
        cdna_after
    };

    let hgvs_p = if !unambiguous {
        None
    } else if consequences.contains(&Consequence::StartLost) {
        Some("p.Met1?".to_string())
    } else {
        let at_hgvs_point = match shifted {
            // A shift into the incomplete terminal codon leaves no whole codon
            // to translate; VEP's partial-codon peptide is not ported, so no
            // HGVS p. is written there.
            Some(shifted)
                if partial_codon_insertion_site(shifted.cds_offset, index)
                    == Some(PartialCodonInsertion::Inside) =>
            {
                None
            }
            Some(shifted) => Some(
                InsertionCodons::new(
                    shifted.cds_offset,
                    &coding(shifted.inserted_bases),
                    chrom,
                    transcript,
                    index,
                    fasta,
                )?
                .peptides,
            ),
            None => Some(codons.peptides.clone()),
        };
        match at_hgvs_point {
            Some(peptides) => crate::hgvs_p::format_hgvs_p_inframe_ins(&peptides, &protein)?,
            None => None,
        }
    };

    let predicts_nmd = consequences.contains(&Consequence::StopGained)
        && super::nmd::predicts_nmd(cds_offset + 1, index);

    // An in-frame insertion can carry its own stop codon, so the termination
    // codon may lie inside the inserted sequence rather than at the variant.
    let peptides = &codons.peptides;
    let ptc = super::nmd::ptc_for_window(
        &consequences,
        &peptides.alt_aas,
        &peptides.ref_aas,
        codons.codon_index,
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
        // Between two codons VEP reports both flanking residues.
        protein_start: Some(peptides.tl_end.min(peptides.tl_start)),
        protein_end: Some(peptides.tl_start),
        codons: Some(format_codons_indel(
            &codons.ref_codon,
            &codons.alt_codon,
            codons.local,
            codons.local,
        )),
        amino_acids: unambiguous
            .then(|| format_amino_acids_indel(&peptides.ref_aas, &peptides.alt_aas)),
        exon: Some(format_exon_number(exon_index, transcript.exon_count)),
        intron: None,
        cds_position: Some(cds_offset),
        cds_position_end: Some(cds_offset + 1),
        cdna_position: Some(cdna_before.min(cdna_after)),
        cdna_position_end: Some(cdna_before.max(cdna_after)),
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
        is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
        is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
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
