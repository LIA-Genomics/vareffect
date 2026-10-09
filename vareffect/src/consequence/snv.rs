//! SNV (single-nucleotide variant) consequence annotation.

use super::helpers::{
    add_splice_spill, build_alt_codon, build_base_result, compute_cdna_position_exonic,
    compute_cdna_position_for_cds, fetch_cds_sequence, fetch_ref_codon, finalize_consequences,
    has_complete_start_codon, is_incomplete_terminal_codon,
};
use super::{Consequence, ConsequenceResult};
use crate::codon::{
    complement, format_amino_acids, format_codon_bases, format_codons,
    translate_codon_for_transcript, translate_reference_codon,
};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_c;
use crate::locate::helpers::exonic_splice_spill;
use crate::locate::{
    LocateIndex, VariantLocation, format_exon_number, format_intron_number, locate_variant,
};
use crate::types::{Strand, TranscriptModel, TranscriptTier};

/// Annotate a single-nucleotide variant against one transcript.
///
/// `ref_base` and `alt_base` are on the **plus strand** (VCF convention).
/// For minus-strand transcripts, internal complementing is performed to
/// derive the coding-strand alleles.
///
/// # Arguments
///
/// * `chrom` -- UCSC-style chromosome name (e.g., `"chr17"`)
/// * `pos` -- 0-based genomic position
/// * `ref_base` -- VCF REF base (uppercase ASCII)
/// * `alt_base` -- VCF ALT base (uppercase ASCII)
/// * `transcript` -- Transcript model to annotate against
/// * `locate_index` -- Precomputed locate index for the transcript
/// * `fasta` -- Reference FASTA reader for codon extraction and ref
///   verification
///
/// # Returns
///
/// A [`ConsequenceResult`] with consequence terms, amino acids, codons,
/// and IMPACT. HGVS fields (`hgvs_c`, `hgvs_p`) are left as `None`.
///
/// # Errors
///
/// - [`VarEffectError::RefMismatch`] if the VCF REF base does not match
///   the reference FASTA
/// - [`VarEffectError::ChromNotFound`] if the chromosome is not in the
///   FASTA index
/// - [`VarEffectError::CoordinateOutOfRange`] if the position exceeds
///   chromosome length
pub fn annotate_snv(
    chrom: &str,
    pos: u64,
    ref_base: u8,
    alt_base: u8,
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    if !fasta.verify_ref(chrom, pos, &[ref_base])? {
        let expected = fasta.fetch_base(chrom, pos)?;
        return Err(VarEffectError::RefMismatch {
            chrom: chrom.to_string(),
            pos,
            expected: String::from(expected as char),
            got: String::from(ref_base as char),
        });
    }
    let result = annotate_snv_verified(
        chrom,
        pos,
        ref_base,
        alt_base,
        transcript,
        locate_index,
        fasta,
    )?;
    super::helpers::without_unknown_alt_effects(
        &super::helpers::AltSite {
            alt: &[alt_base],
            is_indel: false,
            chrom,
            transcript,
            index: locate_index,
            fasta,
        },
        super::helpers::without_gapped_numbering(transcript, result),
    )
}

/// Internal SNV annotator that skips REF verification. Called by the
/// [`super::annotate`] dispatcher after it has already verified the full
/// VCF REF, avoiding a redundant FASTA seek per transcript.
pub(super) fn annotate_snv_verified(
    chrom: &str,
    pos: u64,
    ref_base: u8,
    alt_base: u8,
    transcript: &TranscriptModel,
    locate_index: &LocateIndex,
    fasta: &FastaReader,
) -> Result<ConsequenceResult, VarEffectError> {
    let location = locate_variant(chrom, pos, pos + 1, transcript, locate_index)?;
    let hgvs_c =
        hgvs_c::format_snv_hgvs(pos, ref_base, alt_base, &location, transcript, locate_index)?;

    let mut result = match location {
        VariantLocation::CdsExon {
            exon_index,
            cds_offset,
            codon_number,
            codon_position,
            is_splice_region,
            ..
        } => annotate_cds_snv(
            chrom,
            alt_base,
            transcript,
            locate_index,
            fasta,
            exon_index,
            cds_offset,
            codon_number,
            codon_position,
            is_splice_region,
        ),

        VariantLocation::SpliceDonor { intron_index, .. } => {
            let mut result = build_base_result(transcript, vec![Consequence::SpliceDonorVariant]);
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::SpliceAcceptor { intron_index, .. } => {
            let mut result =
                build_base_result(transcript, vec![Consequence::SpliceAcceptorVariant]);
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::SpliceRegion { intron_index, .. } => {
            let mut result = build_base_result(
                transcript,
                vec![Consequence::SpliceRegionVariant, Consequence::IntronVariant],
            );
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::Intron { intron_index, .. } => {
            let mut result = build_base_result(transcript, vec![Consequence::IntronVariant]);
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::FivePrimeUtr {
            exon_index,
            is_splice_region,
            ..
        } => {
            let mut consequences = vec![Consequence::FivePrimeUtrVariant];
            if is_splice_region {
                consequences.push(Consequence::SpliceRegionVariant);
            }
            let mut result = build_base_result(transcript, consequences);
            result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
            result.cdna_position = compute_cdna_position_exonic(pos, transcript);
            result.cdna_position_end = result.cdna_position;
            Ok(result)
        }

        VariantLocation::ThreePrimeUtr {
            exon_index,
            is_splice_region,
            ..
        } => {
            let mut consequences = vec![Consequence::ThreePrimeUtrVariant];
            if is_splice_region {
                consequences.push(Consequence::SpliceRegionVariant);
            }
            let mut result = build_base_result(transcript, consequences);
            result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
            result.cdna_position = compute_cdna_position_exonic(pos, transcript);
            result.cdna_position_end = result.cdna_position;
            Ok(result)
        }

        VariantLocation::Upstream { .. } => Ok(build_base_result(
            transcript,
            vec![Consequence::UpstreamGeneVariant],
        )),

        VariantLocation::Downstream { .. } => Ok(build_base_result(
            transcript,
            vec![Consequence::DownstreamGeneVariant],
        )),

        VariantLocation::NonCodingExon {
            exon_index,
            is_splice_region,
        } => {
            let mut consequences = vec![Consequence::NonCodingTranscriptExonVariant];
            if is_splice_region {
                consequences.push(Consequence::SpliceRegionVariant);
            }
            let mut result = build_base_result(transcript, consequences);
            result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
            result.cdna_position = compute_cdna_position_exonic(pos, transcript);
            result.cdna_position_end = result.cdna_position;
            Ok(result)
        }

        // VEP 116 gives a variant inside a frameshift intron no splice or
        // intron term, and no exon number; the intron is still numbered
        // (ENST00000400999.7 c.228+1T>A: intron 2/5; ENST00000609389.5
        // n.281+6T>A: intron 3/5). Every such intron in a
        // coding transcript lies inside the CDS in GENCODE v50 and RefSeq
        // RS_2025_08.
        VariantLocation::FrameshiftIntron { intron_index, .. } => {
            let consequence = if transcript.cds_segments.is_empty() {
                Consequence::NonCodingTranscriptVariant
            } else {
                Consequence::CodingSequenceVariant
            };
            let mut result = build_base_result(transcript, vec![consequence]);
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::NonCodingIntron { intron_index, .. } => {
            let mut result = build_base_result(transcript, vec![Consequence::IntronVariant]);
            result.intron = Some(format_intron_number(intron_index, transcript.exon_count));
            Ok(result)
        }

        VariantLocation::Distal => Err(VarEffectError::Malformed(
            "variant is distal to transcript -- callers should filter by overlap first".to_string(),
        )),
    }?;

    result.hgvs_c = hgvs_c;
    add_splice_spill(
        &mut result,
        exonic_splice_spill(pos, pos + 1, transcript, locate_index),
    );

    // Stop_lost: compute the actual extension distance via a 3'UTR scan.
    if result.consequences.contains(&Consequence::StopLost) {
        if let Some(ref aa_str) = result.amino_acids
            && let Some((_, alt_str)) = aa_str.split_once('/')
            && let Some(&alt_byte) = alt_str.as_bytes().first()
        {
            result.hgvs_p = Some(crate::hgvs_p::format_hgvs_p_extension(
                alt_byte,
                result
                    .protein_start
                    .expect("protein_start always set for coding SNV"),
                chrom,
                transcript,
                locate_index,
                fasta,
            )?);
        }
    } else {
        result.hgvs_p = crate::hgvs_p::format_hgvs_p_snv(
            &result.consequences,
            result.amino_acids.as_deref(),
            result.protein_start,
        );
    }
    Ok(result)
}

/// Annotate a CDS exon SNV: codon extraction, translation, AA comparison.
#[allow(clippy::too_many_arguments)]
fn annotate_cds_snv(
    chrom: &str,
    alt_base: u8,
    transcript: &TranscriptModel,
    index: &LocateIndex,
    fasta: &FastaReader,
    exon_index: u16,
    cds_offset: u32,
    codon_number: u32,
    codon_position: u8,
    is_splice_region: bool,
) -> Result<ConsequenceResult, VarEffectError> {
    let codon_start_offset = cds_offset - codon_position as u32;

    // Complement the VCF ALT for minus-strand to get coding-strand base
    let coding_alt = match transcript.strand {
        Strand::Plus => alt_base,
        Strand::Minus => complement(alt_base),
    };

    if is_incomplete_terminal_codon(cds_offset, index) {
        // VEP 116 on a 3' partial codon (ENST00000379099.3 c.925C>A,
        // NM_001424184.1 c.658T>A): incomplete_terminal_codon_variant plus
        // coding_sequence_variant, the partial codon's own bases as codons,
        // amino acid X, protein position of the partial codon.
        let ref_bases = fetch_cds_sequence(
            codon_start_offset,
            index.total_cds_length(),
            chrom,
            transcript,
            index,
            fasta,
        )?;
        let mut alt_bases = ref_bases.clone();
        alt_bases[codon_position as usize] = coding_alt;
        let mut consequences = vec![
            Consequence::IncompleteTerminalCodonVariant,
            Consequence::CodingSequenceVariant,
        ];
        if is_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let mut result = build_base_result(transcript, consequences);
        result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
        result.protein_start = Some(codon_number);
        result.protein_end = Some(codon_number);
        result.codons = Some(format_codon_bases(&ref_bases, &alt_bases, codon_position));
        result.amino_acids = Some("X".to_string());
        result.cds_position = Some(cds_offset + 1);
        result.cds_position_end = Some(cds_offset + 1);
        result.cdna_position = Some(compute_cdna_position_for_cds(cds_offset, index));
        result.cdna_position_end = result.cdna_position;
        return Ok(result);
    }

    let ref_codon = fetch_ref_codon(codon_start_offset, chrom, transcript, index, fasta)?;
    let alt_codon = build_alt_codon(&ref_codon, codon_position, coding_alt);

    let is_mito = transcript.chrom == "chrM";
    let ref_aa = translate_reference_codon(
        &ref_codon,
        codon_number,
        index.selenocysteine_codons(),
        is_mito,
    );
    let alt_aa = translate_codon_for_transcript(&alt_codon, is_mito);

    // Ambiguous codons (ref contains N)
    if ref_aa == b'X' || alt_aa == b'X' {
        let mut consequences = vec![Consequence::CodingSequenceVariant];
        if is_splice_region {
            consequences.push(Consequence::SpliceRegionVariant);
        }
        let mut result = build_base_result(transcript, consequences);
        result.exon = Some(format_exon_number(exon_index, transcript.exon_count));
        result.cds_position = Some(cds_offset + 1);
        result.cds_position_end = Some(cds_offset + 1);
        result.cdna_position = Some(compute_cdna_position_for_cds(cds_offset, index));
        result.cdna_position_end = Some(compute_cdna_position_for_cds(cds_offset, index));
        // The partial first codon of a 5'-incomplete CDS: VEP 116 reports
        // its `N`-padded codons (`nnT/nnA`), amino acid X and protein
        // position 1 (ENST00000379099.3 c.3T>A). A reference `N` elsewhere
        // keeps the bare result. An unknown alternate base in a known codon
        // keeps the codons and protein position, without amino acids
        // (VEP 116 NM_003051.4 `Tct/Nct`, protein 33).
        if codon_number == 1 && index.cds_start_pad() > 0 {
            result.protein_start = Some(1);
            result.protein_end = Some(1);
            result.codons = Some(format_codons(&ref_codon, &alt_codon, codon_position));
            result.amino_acids = Some("X".to_string());
        } else if ref_aa != b'X' {
            result.protein_start = Some(codon_number);
            result.protein_end = Some(codon_number);
            result.codons = Some(format_codons(&ref_codon, &alt_codon, codon_position));
        }
        return Ok(result);
    }

    let is_start_codon = codon_number == 1 && has_complete_start_codon(transcript, index);
    let coding_consequence = if ref_aa == alt_aa {
        if is_start_codon && ref_aa == b'M' {
            Consequence::StartRetainedVariant
        } else if ref_aa == b'*' {
            Consequence::StopRetainedVariant
        } else {
            Consequence::SynonymousVariant
        }
    } else if alt_aa == b'*' {
        Consequence::StopGained
    } else if ref_aa == b'*' {
        Consequence::StopLost
    } else if is_start_codon && ref_aa == b'M' {
        Consequence::StartLost
    } else {
        Consequence::MissenseVariant
    };

    let mut consequences = vec![coding_consequence];
    if is_splice_region {
        consequences.push(Consequence::SpliceRegionVariant);
    }
    let impact = finalize_consequences(&mut consequences);
    let predicts_nmd = coding_consequence == Consequence::StopGained
        && super::nmd::predicts_nmd(cds_offset + 1, index);

    // An SNV replaces one base with one base, so the CDS length is unchanged
    // and the variant's own codon is the termination codon.
    let ptc = if coding_consequence == Consequence::StopGained {
        super::nmd::resolve_ptc(
            crate::hgvs_p::PtcSite::Residue(codon_number),
            &crate::hgvs_p::CdsEdit {
                start: cds_offset,
                del_len: 1,
                ins_len: 1,
            },
            index,
        )
    } else {
        super::PtcStatus::NotApplicable
    };

    Ok(ConsequenceResult {
        transcript: transcript.accession.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        protein_accession: transcript.protein_accession.clone(),
        consequences,
        impact,
        protein_start: Some(codon_number),
        protein_end: Some(codon_number),
        codons: Some(format_codons(&ref_codon, &alt_codon, codon_position)),
        amino_acids: Some(format_amino_acids(ref_aa, alt_aa)),
        exon: Some(format_exon_number(exon_index, transcript.exon_count)),
        intron: None,
        cds_position: Some(cds_offset + 1),
        cds_position_end: Some(cds_offset + 1),
        cdna_position: Some(compute_cdna_position_for_cds(cds_offset, index)),
        cdna_position_end: Some(compute_cdna_position_for_cds(cds_offset, index)),
        strand: transcript.strand,
        biotype: transcript.biotype.clone(),
        is_mane_select: transcript.tier == Some(TranscriptTier::ManeSelect),
        is_mane_plus_clinical: transcript.tier == Some(TranscriptTier::ManePlusClinical),
        is_refseq_select: transcript.tier == Some(TranscriptTier::RefSeqSelect),
        hgvs_c: None,
        hgvs_p: None,
        predicts_nmd,
        ptc,
    })
}
