//! Boundary-spanning deletions, complex delins, and MNV annotation.

use super::helpers::{
    add_splice_spill, build_base_result, complete_codons_end, compute_cdna_position_for_cds,
    ends_in_stop_codon, fetch_cds_sequence, fetch_ref_codon, finalize_consequences,
    has_complete_start_codon, is_incomplete_terminal_codon,
};
use super::indel::{
    annotate_cds_frameshift, build_noncds_indel_result, build_splice_indel_result,
    indel_into_partial_codon_result,
};
use super::{Consequence, ConsequenceResult, Impact, PtcStatus};
use crate::codon::{
    format_amino_acids_indel, format_codons_indel, reverse_complement, translate_reference_codon,
    translate_reference_sequence, translate_sequence,
};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_c;
use crate::locate::helpers::{exonic_splice_spill, frameshift_intron_endpoint};
use crate::locate::{
    IndelLocation, IndelRegion, LocateIndex, format_exon_number, format_intron_number, locate_indel,
};
use crate::types::{Strand, TranscriptModel, TranscriptTier};

/// Compute the number of CDS bases deleted by a genomic range `[del_start, del_end)`.
///
/// Walks all CDS segments and sums the overlap with the deletion footprint.
pub(super) fn compute_cds_projected_length(
    del_start: u64,
    del_end: u64,
    transcript: &TranscriptModel,
) -> u32 {
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

/// Result for an indel or MNV with an end inside frameshift intron
/// `intron_index` (see [`frameshift_intron_endpoint`]).
///
/// VEP 116 gives such a variant only `coding_sequence_variant` (or
/// `non_coding_transcript_variant`), no splice term, the intron number, the
/// exon number when coding bases are removed, and the CDS positions of those
/// bases (ENST00000400999.7 `c.228_228+1del`: exon 2/6, intron 2/5, cDNA
/// 295, CDS 228, protein 76).
pub(super) fn frameshift_intron_indel_result(
    start: u64,
    end: u64,
    intron_index: u16,
    location: &IndelLocation,
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
) -> ConsequenceResult {
    let consequence = if transcript.cds_segments.is_empty() {
        Consequence::NonCodingTranscriptVariant
    } else {
        Consequence::CodingSequenceVariant
    };
    let mut consequences = vec![consequence];
    // Another intron's canonical site in the same footprint keeps its term:
    // VEP skips only the frameshift intron's windows, and the splice detail
    // already leaves them out for a deletion (an insertion's detail still
    // holds the gap's own windows, so it adds nothing).
    if start < end
        && let Some(detail) = location.splice_detail.as_ref()
    {
        if detail.overlaps_donor {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if detail.overlaps_acceptor {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }
    }
    let impact = finalize_consequences(&mut consequences);
    let mut result = build_base_result(transcript, consequences);
    result.impact = impact;
    if start < end {
        let (exon, intron) = footprint_exon_intron_numbers(start, end, transcript);
        result.exon = exon;
        result.intron = intron;
    } else {
        result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
    }
    // VEP 116 reports CDS positions only when the footprint's 5' end (in
    // transcript order) is a coding base: c.228_228+1del has CDS 228,
    // c.228+1_229del none.
    let five_prime_base = match transcript.strand {
        Strand::Plus => start,
        Strand::Minus => end.saturating_sub(1),
    };
    let five_prime_exonic = start < end
        && transcript
            .exons
            .iter()
            .any(|e| e.genomic_start <= five_prime_base && five_prime_base < e.genomic_end);
    let offsets = if five_prime_exonic {
        compute_cds_overlap_offsets(start, end, transcript, locate_index)
    } else {
        Vec::new()
    };
    if let (Some(first), Some(last)) = (offsets.first(), offsets.last()) {
        let (cds_start, cds_end) = (first.0, last.1 - 1);
        result.cds_position = Some(cds_start + 1);
        result.cds_position_end = Some(cds_end + 1);
        result.cdna_position = Some(compute_cdna_position_for_cds(cds_start, locate_index));
        result.cdna_position_end = Some(compute_cdna_position_for_cds(cds_end, locate_index));
        result.protein_start = Some(cds_start / 3 + 1);
        result.protein_end = Some(cds_end / 3 + 1);
    }
    result
}

/// `true` when `[start, end)` crosses at least one intron and every intron it
/// overlaps is a frameshift intron.
fn spans_only_frameshift_introns(start: u64, end: u64, transcript: &TranscriptModel) -> bool {
    let mut crossed = false;
    for pair in transcript.exons.windows(2) {
        let (lo, hi) = if pair[0].genomic_end <= pair[1].genomic_start {
            (pair[0].genomic_end, pair[1].genomic_start)
        } else {
            (pair[1].genomic_end, pair[0].genomic_start)
        };
        if start < hi && lo < end {
            if crate::locate::helpers::frameshift_gap(&pair[0], &pair[1]).is_none() {
                return false;
            }
            crossed = true;
        }
    }
    crossed
}

/// VEP's EXON and INTRON fields for a deletion footprint `[start, end)`:
/// the exons and introns it overlaps, as `N/total` or `first-last/total`
/// (VEP 116, ENST00000400999.7 `c.226_229del` across a one-base gap: exon
/// `2-3/6`, intron `2/5`).
fn footprint_exon_intron_numbers(
    start: u64,
    end: u64,
    transcript: &TranscriptModel,
) -> (Option<String>, Option<String>) {
    let span = |indices: &mut dyn Iterator<Item = usize>, total: usize| {
        let indices: Vec<usize> = indices.collect();
        match (indices.first(), indices.last()) {
            (Some(&first), Some(&last)) if first == last => Some(format!("{}/{total}", first + 1)),
            (Some(&first), Some(&last)) => Some(format!("{}-{}/{total}", first + 1, last + 1)),
            _ => None,
        }
    };
    let exons = &transcript.exons;
    let exon = span(
        &mut exons
            .iter()
            .enumerate()
            .filter(|(_, e)| start < e.genomic_end && e.genomic_start < end)
            .map(|(i, _)| i),
        exons.len(),
    );
    let intron = span(
        &mut exons
            .windows(2)
            .enumerate()
            .filter(|(_, pair)| {
                let (lo, hi) = if pair[0].genomic_end <= pair[1].genomic_start {
                    (pair[0].genomic_end, pair[1].genomic_start)
                } else {
                    (pair[1].genomic_end, pair[0].genomic_start)
                };
                start < hi && lo < end
            })
            .map(|(i, _)| i),
        exons.len().saturating_sub(1),
    );
    (exon, intron)
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

    // A non-coding transcript's deletion that crosses only frameshift introns
    // removes exon bases and no real intron: VEP 116 ENST00000609389.5
    // n.281_282del across a 12-base gap is
    // non_coding_transcript_exon_variant, exon 3-4/6, intron 3/5.
    if transcript.cds_segments.is_empty() && spans_only_frameshift_introns(start, end, transcript) {
        let mut consequences = vec![Consequence::NonCodingTranscriptExonVariant];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        let first_exonic = (start..end).find_map(|pos| {
            crate::consequence::helpers::compute_cdna_position_exonic(pos, transcript)
        });
        let last_exonic = (start..end).rev().find_map(|pos| {
            crate::consequence::helpers::compute_cdna_position_exonic(pos, transcript)
        });
        let (cdna_start, cdna_end) = match transcript.strand {
            Strand::Plus => (first_exonic, last_exonic),
            Strand::Minus => (last_exonic, first_exonic),
        };
        result.cdna_position = cdna_start;
        result.cdna_position_end = cdna_end;
        let (exon, intron) = footprint_exon_intron_numbers(start, end, transcript);
        result.exon = exon;
        result.intron = intron;
        return Ok(result);
    }

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

    // Entire CDS deleted -> transcript ablation. `cds_del_len` counts genomic
    // bases, so it is compared with the unpadded coding length.
    if cds_del_len >= locate_index.coding_length() {
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

    // Both UTR-spanning branches below need the CDS genomic bounds.
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

    // Start-codon removal, the 5' mirror of the stop-codon branch below. The
    // deletion reaches the first CDS base, so the initiator methionine cannot
    // survive. When the footprint also extends into the 5'UTR, VEP emits
    // `{5_prime_UTR_variant, start_lost}` and **no** `frameshift_variant`,
    // however many CDS bases are removed: its `frameshift` predicate
    // short-circuits on the undefined `cds_start` exactly as it does on
    // `cds_end` at the 3' end (`ensembl-variation`,
    // `VariationEffect.pm::frameshift`). Ground truth for a frameshifting
    // case: `ALPL` NM_000478.6 chr1:21554077 `GCACCATGATTT>G` deletes 7 CDS
    // bases and VEP still reports only that pair.
    //
    // A deletion starting at CDS offset 1 or 2 instead leaves the first base
    // intact, so VEP's `cds_start` *is* defined and both `start_lost` and
    // `frameshift_variant` fire. That case belongs to the frameshift branch
    // below, which handles it.
    let extends_into_5utr = match transcript.strand {
        Strand::Plus => start < cds_genomic_start,
        Strand::Minus => end > cds_genomic_end,
    };
    let start_codon_complete = has_complete_start_codon(transcript, locate_index);
    // A 5'-incomplete CDS has no start codon to lose, and VEP's frame
    // predicates short-circuit on the undefined `cds_start` the same way, so
    // a deletion running in from the 5'UTR is only a coding change.
    if extends_into_5utr && !start_codon_complete {
        let mut consequences = vec![
            Consequence::CodingSequenceVariant,
            Consequence::FivePrimeUtrVariant,
        ];
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }
    if start_codon_complete && global_cds_start == 0 {
        let mut consequences = vec![Consequence::StartLost];
        if extends_into_5utr {
            consequences.push(Consequence::FivePrimeUtrVariant);
        }
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
    let touches_stop = total_cds >= 3 && global_cds_end + 3 > total_cds;
    let extends_into_3utr = match transcript.strand {
        Strand::Plus => end > cds_genomic_end,
        Strand::Minus => start < cds_genomic_start,
    };
    if touches_stop
        && extends_into_3utr
        && ends_in_stop_codon(chrom, transcript, locate_index, fasta)?
    {
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

    // A deletion from the CDS into the 3'UTR of a CDS that does not end in a
    // stop codon: there is no stop to lose and no reading frame left to
    // shift, so it is only a coding and UTR change (VEP 116, NM_001424184.1
    // c.658_*2del: incomplete_terminal_codon_variant, coding_sequence_variant,
    // 3_prime_UTR_variant). VEP's `stop_retained_variant` for the same
    // deletion starting one codon earlier (c.656_*2del) names a stop the
    // transcript does not have (VEP_DIVERGENCES.md).
    if extends_into_3utr && !ends_in_stop_codon(chrom, transcript, locate_index, fasta)? {
        let mut consequences = vec![
            Consequence::CodingSequenceVariant,
            Consequence::ThreePrimeUtrVariant,
        ];
        if !total_cds.is_multiple_of(3) {
            consequences.push(Consequence::IncompleteTerminalCodonVariant);
        }
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);
        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        let protein = global_cds_start / 3 + 1;
        result.protein_start = Some(protein);
        result.protein_end = Some(protein);
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        result.cdna_position = Some(compute_cdna_position_for_cds(
            global_cds_start,
            locate_index,
        ));
        // VEP 116 gives this deletion no HGVS at all (NM_001424184.1).
        result.hgvs_p = None;
        if let Some(exon_idx) = location.exon_index {
            result.exon = Some(format_exon_number(exon_idx, transcript.exon_count));
        }
        return Ok(result);
    }

    if !cds_del_len.is_multiple_of(3) {
        // Frameshift
        let first_codon = global_cds_start / 3;
        // Offset 0 was taken by the start-codon branch above, so this would be
        // a deletion beginning at CDS offset 1 or 2: the first CDS base
        // survives, VEP's `cds_start` is therefore defined, and `start_lost`
        // and `frameshift` both fire (`ensembl-variation`,
        // `VariationEffect.pm::start_lost` via `_ins_del_start_altered`) -- the
        // pairing VEP reports for `chr1:152114079 AT>A`. Currently unreachable:
        // arriving here means spanning a CDS/UTR boundary, which starts at
        // offset 0, or an exon boundary, which covers the canonical splice
        // dinucleotides and is diverted above. Handle it rather than assume a
        // surviving reading frame.
        let start_lost = start_codon_complete && global_cds_start < 3;
        let mut consequences = vec![Consequence::FrameshiftVariant];
        if start_lost {
            consequences.push(Consequence::StartLost);
        }
        if location.overlaps_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let impact = finalize_consequences(&mut consequences);

        let codon_start_offset = first_codon * 3;
        let ref_aa = if codon_start_offset + 3 <= total_cds {
            let ref_codon =
                fetch_ref_codon(codon_start_offset, chrom, transcript, locate_index, fasta)?;
            translate_reference_codon(
                &ref_codon,
                first_codon + 1,
                locate_index.selenocysteine_codons(),
                is_mito,
            )
        } else {
            b'X'
        };

        // With codon 1 destroyed there is no defined reading frame: HGVS names
        // that `p.Met1?`, and no termination codon may be claimed against a
        // frame that does not exist.
        let (hgvs_p, ptc) = if start_lost {
            (Some("p.Met1?".to_string()), PtcStatus::Indeterminate)
        } else {
            let frameshift = crate::hgvs_p::format_hgvs_p_frameshift(
                global_cds_start,
                global_cds_end,
                &[], // pure deletion, no inserted bases
                chrom,
                transcript,
                locate_index,
                fasta,
                0, // no 3' shift for boundary-spanning deletions
            )?;
            let ptc = super::nmd::resolve_ptc(frameshift.ptc, &frameshift.edit, locate_index);
            (frameshift.notation, ptc)
        };

        // Codons and amino acids over the whole codon window, as for an
        // exonic frameshift (VEP 116, ENST00000400999.7 c.226_229del across a
        // one-base gap: `TCCGAg/ag`, `SE/X`).
        let window_end = (((global_cds_end - 1) / 3 + 1) * 3).min(total_cds);
        let ref_window = fetch_cds_sequence(
            codon_start_offset,
            window_end,
            chrom,
            transcript,
            locate_index,
            fasta,
        )?;
        let local_del_start = (global_cds_start - codon_start_offset) as usize;
        let local_del_end = (global_cds_end - codon_start_offset) as usize;
        let mut alt_window = ref_window[..local_del_start].to_vec();
        alt_window.extend_from_slice(&ref_window[local_del_end..]);
        let ref_window_aas = if ref_window.len().is_multiple_of(3) {
            translate_reference_sequence(
                &ref_window,
                first_codon + 1,
                locate_index.selenocysteine_codons(),
                is_mito,
            )?
        } else {
            vec![ref_aa]
        };
        let ref_aa_str: String = ref_window_aas.iter().map(|&aa| char::from(aa)).collect();

        let mut result = build_base_result(transcript, consequences);
        result.impact = impact;
        result.protein_start = Some(first_codon + 1);
        result.protein_end = Some(first_codon + 1);
        result.codons = Some(format_codons_indel(
            &ref_window,
            &alt_window,
            local_del_start,
            local_del_end,
        ));
        result.amino_acids = Some(format!("{ref_aa_str}/X"));
        result.cds_position = Some(global_cds_start + 1);
        result.cds_position_end = Some(global_cds_end);
        result.cdna_position = Some(compute_cdna_position_for_cds(
            global_cds_start,
            locate_index,
        ));
        result.cdna_position_end = Some(compute_cdna_position_for_cds(
            global_cds_end.saturating_sub(1),
            locate_index,
        ));
        result.hgvs_p = hgvs_p;
        result.predicts_nmd =
            !start_lost && super::nmd::predicts_nmd(global_cds_start + 1, locate_index);
        result.ptc = ptc;
        let (exon, intron) = footprint_exon_intron_numbers(start, end, transcript);
        result.exon = exon;
        result.intron = intron;
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

    let ref_aas = translate_reference_sequence(
        &ref_seq,
        codon_start / 3 + 1,
        locate_index.selenocysteine_codons(),
        is_mito,
    )?;
    let alt_aas = if alt_seq.is_empty() || !alt_seq.len().is_multiple_of(3) {
        Vec::new()
    } else {
        translate_sequence(&alt_seq, is_mito)?
    };

    let first_codon = codon_start / 3;
    let mut consequences = Vec::new();

    if first_codon == 0 && start_codon_complete && alt_aas.first() != Some(&b'M') {
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
    let ext_codon_end = (codon_end + 90).min(complete_codons_end(locate_index));
    let (ref_aas_ext, alt_aas_ext) = if ext_codon_end > codon_end {
        let ext_seq = fetch_cds_sequence(
            codon_end,
            ext_codon_end,
            chrom,
            transcript,
            locate_index,
            fasta,
        )?;
        let ext_aas = translate_reference_sequence(
            &ext_seq,
            codon_end / 3 + 1,
            locate_index.selenocysteine_codons(),
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
            locate_index.has_selenocysteine_before(protein_start),
            fasta,
        )?);
    }

    let predicts_nmd = consequences.contains(&Consequence::StopGained)
        && super::nmd::predicts_nmd(global_cds_start + 1, locate_index);

    let ptc = super::nmd::ptc_for_window(
        &consequences,
        &alt_aas,
        &ref_aas,
        first_codon,
        &crate::hgvs_p::CdsEdit {
            start: global_cds_start,
            del_len: global_cds_end - global_cds_start,
            ins_len: 0,
        },
        locate_index,
    );
    let footprint_numbers = footprint_exon_intron_numbers(start, end, transcript);

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
        exon: footprint_numbers.0,
        intron: footprint_numbers.1,
        cds_position: Some(global_cds_start + 1),
        cds_position_end: Some(global_cds_end),
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
    let mut result = annotate_complex_delins_core(
        chrom,
        start,
        trimmed_ref,
        trimmed_alt,
        transcript,
        locate_index,
        fasta,
    )?;
    add_splice_spill(
        &mut result,
        exonic_splice_spill(
            start,
            start + trimmed_ref.len() as u64,
            transcript,
            locate_index,
        ),
    );
    Ok(result)
}

/// [`annotate_complex_delins`] without the splice terms a one-base frameshift intron spills
/// onto exon bases.
#[allow(clippy::too_many_arguments)]
fn annotate_complex_delins_core(
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

    if let Some(intron) = frameshift_intron_endpoint(start, del_end, transcript, locate_index) {
        let mut result = frameshift_intron_indel_result(
            start,
            del_end,
            intron,
            &location,
            transcript,
            locate_index,
        );
        result.hgvs_c = hgvs;
        return Ok(result);
    }

    if location.overlaps_splice_canonical
        && exonic_splice_spill(start, del_end, transcript, locate_index).is_none()
    {
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
                // Same terms VEP 116 gives an SNV in a 3' partial codon; the
                // residue is unknown, so the protein effect is `p.?`.
                let consequences = vec![
                    Consequence::IncompleteTerminalCodonVariant,
                    Consequence::CodingSequenceVariant,
                ];
                let mut result = build_base_result(transcript, consequences);
                result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
                result.cds_position = Some(cds_offset_start + 1);
                result.hgvs_c = hgvs;
                result.hgvs_p = Some("p.?".to_string());
                return Ok(result);
            }

            if is_incomplete_terminal_codon(cds_offset_end - 1, locate_index) {
                let mut result = indel_into_partial_codon_result(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_end,
                    trimmed_alt,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                )?;
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

                let ref_aas = translate_reference_sequence(
                    &ref_seq,
                    codon_start / 3 + 1,
                    locate_index.selenocysteine_codons(),
                    is_mito,
                )?;
                let alt_aas = if alt_seq.len().is_multiple_of(3) {
                    translate_sequence(&alt_seq, is_mito)?
                } else {
                    Vec::new()
                };

                let first_codon = codon_start / 3;
                let mut consequences = Vec::new();

                if first_codon == 0
                    && has_complete_start_codon(transcript, locate_index)
                    && alt_aas.first() != Some(&b'M')
                {
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
                // A delins can carry its own stop codon, so the termination
                // codon may lie inside the inserted sequence.
                let ptc = super::nmd::ptc_for_window(
                    &consequences,
                    &alt_aas,
                    &ref_aas,
                    first_codon,
                    &crate::hgvs_p::CdsEdit {
                        start: *cds_offset_start,
                        del_len: cds_offset_end - *cds_offset_start,
                        ins_len: coding_alt.len() as u32,
                    },
                    locate_index,
                );
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
                    is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
                    is_mane_plus_clinical: transcript.tier
                        == Some(TranscriptTier::ManePlusClinical),
                    is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
                    hgvs_c: None,
                    hgvs_p,
                    predicts_nmd: nmd,
                    ptc,
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
    let mut result = annotate_mnv_core(
        chrom,
        start,
        ref_bases,
        alt_bases,
        transcript,
        locate_index,
        fasta,
    )?;
    add_splice_spill(
        &mut result,
        exonic_splice_spill(
            start,
            start + ref_bases.len() as u64,
            transcript,
            locate_index,
        ),
    );
    Ok(result)
}

/// [`annotate_mnv`] without the splice terms a one-base frameshift intron spills
/// onto exon bases.
#[allow(clippy::too_many_arguments)]
fn annotate_mnv_core(
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

    if let Some(intron) = frameshift_intron_endpoint(start, end, transcript, locate_index) {
        let mut result =
            frameshift_intron_indel_result(start, end, intron, &location, transcript, locate_index);
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
                // Same terms VEP 116 gives an SNV in a 3' partial codon; the
                // residue is unknown, so the protein effect is `p.?`.
                let consequences = vec![
                    Consequence::IncompleteTerminalCodonVariant,
                    Consequence::CodingSequenceVariant,
                ];
                let mut result = build_base_result(transcript, consequences);
                result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
                result.cds_position = Some(cds_offset_start + 1);
                result.hgvs_c = hgvs;
                result.hgvs_p = Some("p.?".to_string());
                return Ok(result);
            }

            if is_incomplete_terminal_codon(cds_offset_end - 1, locate_index) {
                let mut result = indel_into_partial_codon_result(
                    chrom,
                    *cds_offset_start,
                    *cds_offset_end,
                    alt_bases,
                    transcript,
                    locate_index,
                    fasta,
                    exon_index,
                    location.overlaps_splice_region,
                )?;
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

            let ref_aas = translate_reference_sequence(
                &ref_seq,
                codon_start / 3 + 1,
                locate_index.selenocysteine_codons(),
                is_mito,
            )?;
            let alt_aas = translate_sequence(&alt_seq, is_mito)?;

            let first_codon = codon_start / 3;
            let mut consequences = Vec::new();

            let all_same = ref_aas.iter().zip(alt_aas.iter()).all(|(r, a)| r == a);

            // The padded first codon of a 5'-incomplete CDS translates to X.
            // VEP calls a change there `coding_sequence_variant`
            // (`coding_unknown`) and never synonymous; a change to another,
            // known residue keeps its term and adds `coding_sequence_variant`
            // (VEP 116, ENST00000455979.1 c.1_2inv vs c.2_3delinsAA). An `N`
            // in the alleles elsewhere keeps today's behavior.
            let has_unknown = first_codon == 0 && locate_index.cds_start_pad() > 0;
            let known_change = ref_aas
                .iter()
                .zip(alt_aas.iter())
                .any(|(r, a)| r != a && *r != b'X' && *a != b'X');
            let is_start_codon =
                first_codon == 0 && has_complete_start_codon(transcript, locate_index);
            if has_unknown && !known_change {
                consequences.push(Consequence::CodingSequenceVariant);
            } else if all_same {
                if is_start_codon && ref_aas.first() == Some(&b'M') {
                    consequences.push(Consequence::StartRetainedVariant);
                } else if ref_aas.last() == Some(&b'*') {
                    consequences.push(Consequence::StopRetainedVariant);
                } else {
                    consequences.push(Consequence::SynonymousVariant);
                }
            } else if is_start_codon
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
            if has_unknown && known_change {
                consequences.push(Consequence::CodingSequenceVariant);
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

            let hgvs_p = if has_unknown && !known_change {
                Some("p.?".to_string())
            } else {
                crate::hgvs_p::format_hgvs_p_delins(
                    &ref_aas,
                    &alt_aas,
                    protein_start,
                    protein_end,
                    &consequences,
                )
            };

            let predicts_nmd = consequences.contains(&Consequence::StopGained)
                && super::nmd::predicts_nmd(cds_offset_start + 1, locate_index);
            // An MNV substitutes base-for-base, so the CDS length is unchanged.
            let span = cds_offset_end - *cds_offset_start;
            let ptc = super::nmd::ptc_for_window(
                &consequences,
                &alt_aas,
                &ref_aas,
                first_codon,
                &crate::hgvs_p::CdsEdit {
                    start: *cds_offset_start,
                    del_len: span,
                    ins_len: span,
                },
                locate_index,
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
                is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
                is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
                is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
                hgvs_c: hgvs,
                hgvs_p,
                predicts_nmd,
                ptc,
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
