//! Boundary-spanning deletions, complex delins, and MNV annotation.

use super::helpers::{
    build_base_result, compute_cdna_position_for_cds, fetch_cds_sequence, fetch_ref_codon,
    finalize_consequences, is_incomplete_terminal_codon,
};
use super::indel::{annotate_cds_frameshift, build_noncds_indel_result, build_splice_indel_result};
use super::{Consequence, ConsequenceResult, Impact};
use crate::codon::{
    format_amino_acids_indel, format_codons_indel, reverse_complement,
    translate_codon_for_transcript, translate_sequence,
};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_c;
use crate::locate::{
    IndelLocation, IndelRegion, LocateIndex, format_exon_number, format_intron_number, locate_indel,
};
use crate::types::{Strand, TranscriptModel, TranscriptTier};

/// Compute the number of CDS bases deleted by a genomic range `[del_start, del_end)`.
///
/// Walks all CDS segments and sums the overlap with the deletion footprint.
fn compute_cds_projected_length(del_start: u64, del_end: u64, transcript: &TranscriptModel) -> u32 {
    let mut cds_bases = 0u32;
    for seg in &transcript.cds_segments {
        let overlap_start = del_start.max(seg.genomic_start);
        let overlap_end = del_end.min(seg.genomic_end);
        if overlap_start < overlap_end {
            cds_bases += (overlap_end - overlap_start) as u32;
        }
    }
    cds_bases
}

/// For each CDS segment overlapping `[del_start, del_end)`, compute the
/// CDS offset range of the overlap. Returns `(cds_offset_start, cds_offset_end)`
/// pairs in half-open convention, ordered by CDS position.
fn compute_cds_overlap_offsets(
    del_start: u64,
    del_end: u64,
    transcript: &TranscriptModel,
    index: &LocateIndex,
) -> Vec<(u32, u32)> {
    let cumulative = index.cumulative_cds();
    let mut offsets = Vec::new();
    for (seg_idx, seg) in transcript.cds_segments.iter().enumerate() {
        let overlap_start = del_start.max(seg.genomic_start);
        let overlap_end = del_end.min(seg.genomic_end);
        if overlap_start >= overlap_end {
            continue;
        }
        let (intra_start, intra_end) = match transcript.strand {
            Strand::Plus => (
                (overlap_start - seg.genomic_start) as u32,
                (overlap_end - seg.genomic_start) as u32,
            ),
            Strand::Minus => (
                (seg.genomic_end - overlap_end) as u32,
                (seg.genomic_end - overlap_start) as u32,
            ),
        };
        let seg_cds_start = cumulative[seg_idx];
        offsets.push((seg_cds_start + intra_start, seg_cds_start + intra_end));
    }
    offsets
}

/// Annotate a deletion that crosses an exon/intron boundary or spans
/// multiple exons.
///
/// **Priority**:
/// 1. Splice canonical overlap -> splice donor/acceptor consequence.
/// 2. CDS-projected deletion -> frame analysis (frameshift or inframe).
/// 3. No CDS overlap -> intronic/UTR consequence.
#[allow(clippy::too_many_arguments)]
pub(super) fn annotate_boundary_spanning_deletion(
    chrom: &str,
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
    location: &IndelLocation,
) -> Result<ConsequenceResult, VarEffectError> {
    let is_mito = transcript.chrom == "chrM";

    if location.overlaps_splice_canonical {
        return build_splice_indel_result(transcript, location, start, end);
    }

    let cds_del_len = compute_cds_projected_length(start, end, transcript);
    let total_cds = locate_index.total_cds_length();

    // No CDS bases deleted -> purely intronic/UTR deletion
    if cds_del_len == 0 {
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
        return Ok(result);
    }

    // Entire CDS deleted -> transcript ablation
    if cds_del_len >= total_cds {
        let consequences = vec![Consequence::TranscriptAblation];
        let mut result = build_base_result(transcript, consequences);
        result.impact = Impact::High;
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }

    let cds_offsets = compute_cds_overlap_offsets(start, end, transcript, locate_index);
    if cds_offsets.is_empty() {
        return Ok(build_base_result(
            transcript,
            vec![Consequence::CodingSequenceVariant],
        ));
    }

    let global_cds_start = cds_offsets
        .first()
        .expect("cds_offsets verified non-empty")
        .0;
    let global_cds_end = cds_offsets
        .last()
        .expect("cds_offsets verified non-empty")
        .1;

    // Start codon removal
    if global_cds_start == 0 && cds_del_len >= 3 {
        let mut consequences = vec![Consequence::StartLost];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        result.protein_start = Some(1);
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }

    // Stop-codon-spanning deletion regression. The deletion
    // overlaps the last CDS codon (the stop) AND extends genomically into
    // the 3' UTR. VEP emits `{stop_lost, 3_prime_UTR_variant}` without
    // `frameshift_variant` because its `frameshift` predicate short-circuits
    // when `cds_end` is undefined (see `ensembl-variation`
    // `VariationEffect.pm::frameshift`: `return 0 unless defined $bvfo->cds_start
    // && defined $bvfo->cds_end;`). vareffect replicates that behavior here
    // by detecting the geometry and emitting the stop_lost/UTR pair
    // directly, skipping both the frameshift and inframe branches below.
    //
    // `touches_stop`: deletion overlaps the last CDS codon `[total_cds-3, total_cds)`.
    // `extends_into_3utr`: genomic footprint crosses past `cds_genomic_end`
    // (plus strand) or before `cds_genomic_start` (minus strand).
    let cds_genomic_start = transcript.cds_genomic_start.ok_or_else(|| {
        VarEffectError::Malformed(format!(
            "{}: coding transcript has no cds_genomic_start",
            transcript.accession,
        ))
    })?;
    let cds_genomic_end = transcript.cds_genomic_end.ok_or_else(|| {
        VarEffectError::Malformed(format!(
            "{}: coding transcript has no cds_genomic_end",
            transcript.accession,
        ))
    })?;
    let touches_stop = total_cds >= 3 && global_cds_end + 3 > total_cds;
    let extends_into_3utr = match transcript.strand {
        Strand::Plus => end > cds_genomic_end,
        Strand::Minus => start < cds_genomic_start,
    };
    if touches_stop && extends_into_3utr {
        let mut consequences = vec![Consequence::StopLost, Consequence::ThreePrimeUtrVariant];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        // Stop codon is the last CDS codon; protein position = total_cds / 3.
        let stop_codon_number = (total_cds / 3).max(1);
        result.protein_start = Some(stop_codon_number);
        result.protein_end = Some(stop_codon_number);
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }

    if !cds_del_len.is_multiple_of(3) {
        // Frameshift
        let first_codon = global_cds_start / 3;
        let mut consequences = vec![Consequence::FrameshiftVariant];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);

        let codon_start_offset = first_codon * 3;
        let ref_aa = if codon_start_offset + 3 <= total_cds {
            let ref_codon =
                fetch_ref_codon(codon_start_offset, chrom, transcript, locate_index, fasta)?;
            translate_codon_for_transcript(&ref_codon, is_mito)
        } else {
            b'X'
        };

        let hgvs_p = crate::hgvs_p::format_hgvs_p_frameshift(
            global_cds_start,
            global_cds_end,
            &[], // pure deletion, no inserted bases
            chrom,
            transcript,
            locate_index,
            fasta,
            0, // no 3' shift for boundary-spanning deletions
        )?;

        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        result.protein_start = Some(first_codon + 1);
        result.protein_end = Some(first_codon + 1);
        result.amino_acids = Some(format!("{}/X", char::from(ref_aa)));
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        result.hgvs_p = hgvs_p;
        result.predicts_nmd = super::nmd::predicts_nmd(global_cds_start + 1, locate_index);
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }

    // Inframe deletion: fetch, remove deleted CDS bases, translate, compare
    let codon_start = (global_cds_start / 3) * 3;
    let codon_end = ((global_cds_end - 1) / 3 + 1) * 3;
    let codon_end = codon_end.min(total_cds);

    let ref_seq = fetch_cds_sequence(
        codon_start,
        codon_end,
        chrom,
        transcript,
        locate_index,
        fasta,
    )?;

    // The deletion may span non-contiguous CDS regions (across introns),
    // but in the CDS sequence they are contiguous.
    let local_del_start = (global_cds_start - codon_start) as usize;
    let local_del_end = (global_cds_end - codon_start) as usize;
    let mut alt_seq = Vec::with_capacity(ref_seq.len() - (local_del_end - local_del_start));
    alt_seq.extend_from_slice(&ref_seq[..local_del_start]);
    alt_seq.extend_from_slice(&ref_seq[local_del_end..]);

    let ref_aas = translate_sequence(&ref_seq, is_mito)?;
    let alt_aas = if alt_seq.is_empty() || !alt_seq.len().is_multiple_of(3) {
        Vec::new()
    } else {
        translate_sequence(&alt_seq, is_mito)?
    };

    let first_codon = codon_start / 3;
    let mut consequences = Vec::new();

    if first_codon == 0 && alt_aas.first() != Some(&b'M') {
        consequences.push(Consequence::StartLost);
    }
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
        consequences.push(Consequence::InframeDeletion);
    }
    if location.overlaps_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }

    let impact = finalize_consequences(&mut consequences);
    let protein_start = first_codon + 1;
    let protein_end = codon_end / 3;
    let cdna_start = compute_cdna_position_for_cds(global_cds_start, locate_index);
    let cdna_end = compute_cdna_position_for_cds(global_cds_end.saturating_sub(1), locate_index);

    // Extend AA window with downstream context for the protein 3' rule.
    let ext_codon_end = (codon_end + 90).min(total_cds);
    let (ref_aas_ext, alt_aas_ext) = if ext_codon_end > codon_end {
        let ext_seq = fetch_cds_sequence(
            codon_end,
            ext_codon_end,
            chrom,
            transcript,
            locate_index,
            fasta,
        )?;
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
        && super::nmd::predicts_nmd(global_cds_start + 1, locate_index);

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
        exon: location
            .exon_index
            .map(|idx| format_exon_number(idx, transcript.exon_count)),
        intron: location
            .intron_index
            .map(|idx| format_intron_number(idx, transcript.exon_count)),
        cds_position: Some(global_cds_start + 1),
        cds_position_end: Some(global_cds_end),
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
    })
}

/// Annotate a complex delins (deletion + insertion where both trimmed_ref
/// and trimmed_alt are non-empty with different lengths).
///
/// For CDS variants: expand to codon boundaries, replace deleted CDS bases
/// with inserted alt bases, then check frame. If the resulting CDS alt
/// length is not a multiple of 3 -> `FrameshiftVariant`. Otherwise ->
/// `ProteinAlteringVariant` (VEP convention for delins that change both
/// sequence and length).
#[allow(clippy::too_many_arguments)]
pub(super) fn annotate_complex_delins(
    chrom: &str,
    start: u64,
    trimmed_ref: &[u8],
    trimmed_alt: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let ref_len = trimmed_ref.len() as u64;
    let del_end = start + ref_len;
    let is_mito = transcript.chrom == "chrM";

    let location = locate_indel(chrom, start, del_end, transcript, locate_index)?;
    let hgvs =
        hgvs_c::format_delins_hgvs(chrom, start, del_end, trimmed_alt, transcript, locate_index)?;

    if location.overlaps_splice_canonical {
        let mut result = build_splice_indel_result(transcript, &location, start, del_end)?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    if location.crosses_exon_boundary {
        let mut result = annotate_boundary_spanning_deletion(
            chrom,
            start,
            del_end,
            transcript,
            locate_index,
            fasta,
            &location,
        )?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    match &location.region {
        IndelRegion::Cds {
            cds_offset_start,
            cds_offset_end,
        } => {
            let cds_ref_len = cds_offset_end - cds_offset_start;
            let cds_alt_len = trimmed_alt.len() as u32;
            let net_change = cds_alt_len as i64 - cds_ref_len as i64;
            let exon_index = location.exon_index.unwrap_or(0);

            if is_incomplete_terminal_codon(*cds_offset_start, locate_index) {
                let consequences = vec![Consequence::IncompleteTerminalCodonVariant];
                let mut result = build_base_result(transcript, consequences);
                result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
                result.cds_position = Some(cds_offset_start + 1);
                result.hgvs_c = hgvs;
                return Ok(result);
            }

            if net_change % 3 != 0 {
                let mut result = annotate_cds_frameshift(
                    chrom,
                    *cds_offset_start,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                    *cds_offset_end,
                    None,
                    Some(trimmed_alt),
                    0,
                )?;
                result.hgvs_c = hgvs;
                Ok(result)
            } else {
                // Inframe delins
                let codon_start = (cds_offset_start / 3) * 3;
                let codon_end = ((cds_offset_end - 1) / 3 + 1) * 3;
                let total_cds = locate_index.total_cds_length();
                let codon_end = codon_end.min(total_cds);

                let ref_seq = fetch_cds_sequence(
                    codon_start,
                    codon_end,
                    chrom,
                    transcript,
                    locate_index,
                    fasta,
                )?;

                let local_del_start = (cds_offset_start - codon_start) as usize;
                let local_del_end = (cds_offset_end - codon_start) as usize;
                let coding_alt = match transcript.strand {
                    Strand::Plus => trimmed_alt.to_vec(),
                    Strand::Minus => reverse_complement(trimmed_alt),
                };
                let mut alt_seq = Vec::with_capacity(
                    ref_seq.len() - (local_del_end - local_del_start) + coding_alt.len(),
                );
                alt_seq.extend_from_slice(&ref_seq[..local_del_start]);
                alt_seq.extend_from_slice(&coding_alt);
                alt_seq.extend_from_slice(&ref_seq[local_del_end..]);

                let ref_aas = translate_sequence(&ref_seq, is_mito)?;
                let alt_aas = if alt_seq.len().is_multiple_of(3) {
                    translate_sequence(&alt_seq, is_mito)?
                } else {
                    Vec::new()
                };

                let first_codon = codon_start / 3;
                let mut consequences = Vec::new();

                if first_codon == 0 && alt_aas.first() != Some(&b'M') {
                    consequences.push(Consequence::StartLost);
                }
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
                    consequences.push(Consequence::ProteinAlteringVariant);
                }
                if location.overlaps_splice_region {
                    consequences.push(Consequence::SpliceRegionVariant);
                }

                let impact = finalize_consequences(&mut consequences);
                let protein_start = first_codon + 1;
                let protein_end = codon_end / 3;
                let cdna_start = compute_cdna_position_for_cds(*cds_offset_start, locate_index);
                let cdna_end =
                    compute_cdna_position_for_cds(cds_offset_end.saturating_sub(1), locate_index);

                let hgvs_p = crate::hgvs_p::format_hgvs_p_delins(
                    &ref_aas,
                    &alt_aas,
                    protein_start,
                    protein_end,
                    &consequences,
                );

                let nmd = consequences.contains(&Consequence::StopGained)
                    && super::nmd::predicts_nmd(cds_offset_start + 1, locate_index);
                let mut result = ConsequenceResult {
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
                    cds_position_end: Some(*cds_offset_end),
                    cdna_position: Some(cdna_start),
                    cdna_position_end: Some(cdna_end),
                    strand: transcript.strand,
                    biotype: transcript.biotype.clone(),
                    is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
                    is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
                    is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
                    hgvs_c: None,
                    hgvs_p,
                    predicts_nmd: nmd,
                };
                result.hgvs_c = hgvs;
                Ok(result)
            }
        }
        _ => {
            let del_end = start + trimmed_ref.len() as u64;
            let mut result = build_noncds_indel_result(transcript, &location, start, del_end)?;
            result.hgvs_c = hgvs;
            Ok(result)
        }
    }
}

/// Annotate a multi-nucleotide variant (same-length substitution > 1bp).
///
/// MNVs are treated as a single mutational event: the entire ref region
/// is replaced by the alt. For CDS variants, all affected codons are
/// expanded, translated, and compared. VEP reports the most severe
/// consequence across affected codons.
#[allow(clippy::too_many_arguments)]
pub(super) fn annotate_mnv(
    chrom: &str,
    start: u64,
    ref_bases: &[u8],
    alt_bases: &[u8],
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let n = ref_bases.len() as u64;
    let end = start + n;
    let is_mito = transcript.chrom == "chrM";

    let location = locate_indel(chrom, start, end, transcript, locate_index)?;
    let hgvs = hgvs_c::format_delins_hgvs(chrom, start, end, alt_bases, transcript, locate_index)?;

    if location.overlaps_splice_canonical {
        let mut result = build_splice_indel_result(transcript, &location, start, end)?;
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    // Boundary crossing: MNV spanning exon/intron is very rare.
    // Treat as coding_sequence_variant since the intronic bases are
    // substituted (no splice site destroyed by same-length change unless
    // canonical overlap, which is caught above).
    if location.crosses_exon_boundary {
        let mut consequences = vec![Consequence::CodingSequenceVariant];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let mut result = build_base_result(transcript, consequences);
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        if let Some(intron_idx) = location.intron_index {
            result.intron = Some(format_intron_number(intron_idx, transcript.exon_count));
        }
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    match &location.region {
        IndelRegion::Cds {
            cds_offset_start,
            cds_offset_end,
        } => {
            let exon_index = location.exon_index.unwrap_or(0);

            if is_incomplete_terminal_codon(*cds_offset_start, locate_index) {
                let consequences = vec![Consequence::IncompleteTerminalCodonVariant];
                let mut result = build_base_result(transcript, consequences);
                result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
                result.cds_position = Some(cds_offset_start + 1);
                result.hgvs_c = hgvs;
                return Ok(result);
            }

            let codon_start = (cds_offset_start / 3) * 3;
            let codon_end = ((*cds_offset_end - 1) / 3 + 1) * 3;
            let total_cds = locate_index.total_cds_length();
            let codon_end = codon_end.min(total_cds);

            let ref_seq = fetch_cds_sequence(
                codon_start,
                codon_end,
                chrom,
                transcript,
                locate_index,
                fasta,
            )?;

            // VCF alt_bases are in ascending genomic order (plus strand).
            // For minus-strand transcripts, the lowest genomic position maps
            // to the HIGHEST CDS offset. Reverse-complement the entire alt
            // block to get coding-strand bases in 5'-to-3' transcript order.
            let coding_alt = match transcript.strand {
                Strand::Plus => alt_bases.to_vec(),
                Strand::Minus => reverse_complement(alt_bases),
            };
            let mut alt_seq = ref_seq.clone();
            for i in 0..n as u32 {
                let cds_offset = cds_offset_start + i;
                if cds_offset >= codon_end {
                    break;
                }
                let local_pos = (cds_offset - codon_start) as usize;
                alt_seq[local_pos] = coding_alt[i as usize];
            }

            let ref_aas = translate_sequence(&ref_seq, is_mito)?;
            let alt_aas = translate_sequence(&alt_seq, is_mito)?;

            let first_codon = codon_start / 3;
            let mut consequences = Vec::new();

            let all_same = ref_aas.iter().zip(alt_aas.iter()).all(|(r, a)| r == a);

            if all_same {
                if first_codon == 0 && ref_aas.first() == Some(&b'M') {
                    consequences.push(Consequence::StartRetainedVariant);
                } else if ref_aas.last() == Some(&b'*') {
                    consequences.push(Consequence::StopRetainedVariant);
                } else {
                    consequences.push(Consequence::SynonymousVariant);
                }
            } else if first_codon == 0
                && ref_aas.first() == Some(&b'M')
                && alt_aas.first() != Some(&b'M')
            {
                consequences.push(Consequence::StartLost);
            } else if alt_aas.contains(&b'*') && !ref_aas.contains(&b'*') {
                consequences.push(Consequence::StopGained);
            } else if ref_aas.contains(&b'*') && !alt_aas.contains(&b'*') {
                consequences.push(Consequence::StopLost);
            } else {
                consequences.push(Consequence::MissenseVariant);
            }

            if location.overlaps_splice_region {
                consequences.push(Consequence::SpliceRegionVariant);
            }

            let impact = finalize_consequences(&mut consequences);
            let protein_start = first_codon + 1;
            let protein_end = codon_end / 3;
            let cdna_start = compute_cdna_position_for_cds(*cds_offset_start, locate_index);
            let cdna_end =
                compute_cdna_position_for_cds(cds_offset_end.saturating_sub(1), locate_index);

            let hgvs_p = crate::hgvs_p::format_hgvs_p_delins(
                &ref_aas,
                &alt_aas,
                protein_start,
                protein_end,
                &consequences,
            );

            let predicts_nmd = consequences.contains(&Consequence::StopGained)
                && super::nmd::predicts_nmd(cds_offset_start + 1, locate_index);
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
                    (cds_offset_start - codon_start) as usize,
                    (cds_offset_end - codon_start) as usize,
                )),
                amino_acids: Some(format_amino_acids_indel(&ref_aas, &alt_aas)),
                exon: Some(format_exon_number(exon_index, transcript.exon_count)),
                intron: None,
                cds_position: Some(cds_offset_start + 1),
                cds_position_end: Some(*cds_offset_end),
                cdna_position: Some(cdna_start),
                cdna_position_end: Some(cdna_end),
                strand: transcript.strand,
                biotype: transcript.biotype.clone(),
                is_mane_select: transcript.tier == TranscriptTier::ManeSelect,
                is_mane_plus_clinical: transcript.tier == TranscriptTier::ManePlusClinical,
                is_refseq_select: transcript.tier == TranscriptTier::RefSeqSelect,
                hgvs_c: hgvs,
                hgvs_p,
                predicts_nmd,
            })
        }
        _ => {
            let mnv_end = start + ref_bases.len() as u64;
            let mut result = build_noncds_indel_result(transcript, &location, start, mnv_end)?;
            result.hgvs_c = hgvs;
            Ok(result)
        }
    }
}
