//! Indel (insertion and deletion) consequence annotation.

use super::helpers::{
    build_base_result, compute_cdna_position_exonic, compute_cdna_position_for_cds,
    fetch_cds_sequence, fetch_ref_codon, finalize_consequences,
};
use super::{Consequence, ConsequenceResult};
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

    if location.overlaps_splice_canonical {
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

            if !del_cds_len.is_multiple_of(3) {
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
    let location = locate_indel(chrom, pos, pos, transcript, locate_index)?;

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

    if location.overlaps_splice_canonical {
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

            if !ins_len.is_multiple_of(3) {
                annotate_cds_frameshift(
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
                )
            } else {
                annotate_cds_inframe_insertion(
                    chrom,
                    *cds_offset_start,
                    inserted_bases,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                    shift as u32,
                )
            }
        }
        _ => build_noncds_indel_result(transcript, &location, pos, pos),
    }?;

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
    let overlaps_start_codon = codon_number == 1;
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

    // Generate HGVS protein notation for the frameshift.
    let hgvs_p = if start_lost {
        Some("p.Met1?".to_string())
    } else {
        crate::hgvs_p::format_hgvs_p_frameshift(
            cds_offset_start,
            cds_offset_end,
            &coding_inserted_for_hgvs,
            chrom,
            transcript,
            index,
            fasta,
            dna_shift,
        )?
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
        && ((cds_offset_start == 0 && cds_offset_end >= 3) || alt_aas.first() != Some(&b'M'))
    {
        consequences.push(Consequence::StartLost);
    }

    if alt_aas.contains(&b'*') {
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

    if first_codon == 0 && alt_aas.first() != Some(&b'M') {
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
