//! Variant consequence assignment -- determines SO consequence terms for
//! SNVs, simple indels, boundary-spanning deletions, complex delins, and
//! MNVs against transcript models.
//!
//! This module takes a variant's genomic position and alleles, locates it
//! within one or more transcripts (via [`crate::locate`]), and for CDS
//! variants extracts the reference codons from the FASTA, translates both
//! ref and alt sequences, and assigns the appropriate SO consequence term(s).
//!
//! # Variant types handled
//!
//! - **SNVs:** `ref_allele.len() == 1 && alt_allele.len() == 1`
//! - **Simple indels:** pure insertions/deletions where the entire
//!   footprint falls within a single transcript region (CDS, UTR, intron,
//!   etc.). Includes frameshift, inframe insertion/deletion, and splice
//!   site overlap detection for indels.
//! - **Boundary-spanning deletions:** multi-exon deletions,
//!   exon-intron boundary spans.
//! - **Complex delins:** deletion + insertion with different lengths.
//! - **MNVs:** same-length multi-base substitutions.
//!
//! # VEP concordance
//!
//! This module replicates the logic in:
//! - `Bio::EnsEMBL::Variation::Utils::VariationEffect` -- consequence predicates
//! - `Bio::EnsEMBL::Variation::TranscriptVariationAllele` -- codon extraction,
//!   AA comparison
//!
//! # Thread safety
//!
//! [`annotate_snv`] calls [`FastaReader::fetch_base`]
//! up to 4 times per invocation (1 ref-verify + 3 codon bases). Indel
//! annotation may call it more (proportional to the codon-expanded region
//! size). For bulk annotation across multiple threads, create one
//! [`FastaReader`] per thread via
//! [`FastaReader::try_clone`] to avoid mutex
//! contention on the internal seek-based reader.

mod complex;
pub(crate) mod helpers;
mod indel;
mod nmd;
mod snv;
#[cfg(test)]
mod tests;

pub use indel::{annotate_deletion, annotate_insertion};
pub use snv::annotate_snv;

use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::transcript::TranscriptStore;
use crate::types::Biotype;

/// Top-level annotation output bundling per-transcript consequences with
/// any structured warnings the annotator surfaced.
///
/// Consequences and warnings are kept as separate vectors (rather than
/// pushing warnings onto each [`ConsequenceResult`]) so a divergent
/// transcript that overlaps multiple consequence rows produces exactly
/// one entry in [`Self::warnings`] per affected transcript. Downstream
/// consumers (CSQ formatter, VEP-JSON serializer) iterate the warnings
/// list once instead of de-duplicating per-row.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AnnotationResult {
    /// Per-transcript consequence rows. Empty only when no overlap was
    /// found and the intergenic-variant fallback was suppressed (which
    /// the public API does not do today, but the type permits).
    pub consequences: Vec<ConsequenceResult>,
    /// Structured warnings raised during annotation. Empty when nothing
    /// noteworthy happened — clinical-grade callers gate output on this
    /// being empty *or* on every entry being acceptable.
    pub warnings: Vec<Warning>,
}

impl AnnotationResult {
    /// Construct an [`AnnotationResult`] from a consequence vector with no
    /// warnings. Convenience for the common case of a non-divergent
    /// annotation.
    pub fn from_consequences(consequences: Vec<ConsequenceResult>) -> Self {
        Self {
            consequences,
            warnings: Vec::new(),
        }
    }
}

/// Structured warning surfaced by [`crate::VarEffect::annotate`] on top of
/// the consequence vector.
///
/// `#[non_exhaustive]` so future variants — `TranscriptVersionDrift`,
/// `AmbiguousReference`, `LowCoverageRegion`, etc. — can be added without
/// a SemVer break for downstream `match`-on-warning consumers.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Warning {
    /// The chosen transcript is flagged by NCBI as having sequence that
    /// differs from the reference assembly at one or more positions.
    /// HGVS positions emitted against this transcript may not map back to
    /// the same genomic position they would against the reference, so
    /// clinical callers should consider falling back to a non-divergent
    /// transcript for variant reporting.
    ///
    /// Sourced from the GFF3 `Note=` attribute on the mRNA / CDS rows;
    /// see [`crate::TranscriptModel::genome_transcript_divergent`]. About
    /// 5 % of NCBI RefSeq Select transcripts on GRCh37 carry this flag;
    /// MANE Select / Plus Clinical transcripts on GRCh38 are curated to
    /// exclude divergence by construction.
    DivergentTranscript {
        /// Affected transcript accession with version (e.g.
        /// `"NM_001134380.2"`).
        accession: String,
    },
}

/// VEP-compatible severity rating for consequence terms.
///
/// Variants are ordered by declaration so [`Ord`] gives
/// `Modifier < Low < Moderate < High`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Impact {
    /// Lowest severity -- variants in non-coding regions, upstream/downstream.
    Modifier,
    /// Low severity -- synonymous, splice region, start/stop retained.
    Low,
    /// Moderate severity -- missense.
    Moderate,
    /// Highest severity -- stop gained/lost, start lost, frameshift,
    /// splice donor/acceptor.
    High,
}

impl std::fmt::Display for Impact {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Modifier => "MODIFIER",
            Self::Low => "LOW",
            Self::Moderate => "MODERATE",
            Self::High => "HIGH",
        })
    }
}

/// Sequence Ontology consequence term with its VEP-assigned IMPACT.
///
/// The string representation (via [`Consequence::as_str`]) matches VEP's
/// output exactly (e.g., `"missense_variant"`). Covers SNV, indel,
/// boundary-spanning, complex delins, and MNV consequence terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Consequence {
    // --- HIGH impact ---
    /// Variant causes loss of the transcription unit (typically large SVs).
    TranscriptAblation,
    /// Variant in the canonical splice acceptor site (intronic -1 or -2).
    SpliceAcceptorVariant,
    /// Variant in the canonical splice donor site (intronic +1 or +2).
    SpliceDonorVariant,
    /// Premature stop codon introduced by the variant.
    StopGained,
    /// Insertion or deletion shifts the reading frame (length change not
    /// divisible by 3).
    FrameshiftVariant,
    /// Stop codon changed to a coding amino acid.
    StopLost,
    /// Start codon (ATG) disrupted -- translation initiation site lost.
    StartLost,

    // --- MODERATE impact ---
    /// In-frame insertion of one or more codons (length divisible by 3).
    InframeInsertion,
    /// In-frame deletion of one or more codons (length divisible by 3).
    InframeDeletion,
    /// Non-synonymous coding change -- different amino acid.
    MissenseVariant,
    /// Coding variant whose exact protein effect is ambiguous (e.g.
    /// complex delins that changes both sequence and length).
    ProteinAlteringVariant,

    // --- LOW impact ---
    /// Variant in the splice region (intronic +3..+8 / -3..-8, or exonic
    /// 1-3 bases from boundary). Always paired with the primary consequence.
    SpliceRegionVariant,
    /// Synonymous change at the start codon (ATG preserved). Only possible
    /// on chrM where ATA also codes for Met (NCBI table 2).
    StartRetainedVariant,
    /// Synonymous change within the stop codon.
    StopRetainedVariant,
    /// SO:0001626 -- variant in the final partial codon of an incompletely
    /// annotated CDS (total CDS length not divisible by 3).
    IncompleteTerminalCodonVariant,
    /// Synonymous coding change -- same amino acid.
    SynonymousVariant,

    // --- MODIFIER impact ---
    /// Variant in a coding exon but codon could not be determined (e.g.,
    /// reference contains N).
    CodingSequenceVariant,
    /// Variant in the 5' UTR.
    FivePrimeUtrVariant,
    /// Variant in the 3' UTR.
    ThreePrimeUtrVariant,
    /// Variant in an exon of a non-coding transcript.
    NonCodingTranscriptExonVariant,
    /// Variant in an intron.
    IntronVariant,
    /// Variant upstream of the transcript (within 5 kb, 5' direction).
    UpstreamGeneVariant,
    /// Variant downstream of the transcript (within 5 kb, 3' direction).
    DownstreamGeneVariant,
    /// Variant does not overlap any transcript's `[tx_start, tx_end)` region.
    IntergenicVariant,
}

impl Consequence {
    /// VEP's IMPACT rating for this consequence.
    pub fn impact(&self) -> Impact {
        match self {
            Self::TranscriptAblation
            | Self::SpliceAcceptorVariant
            | Self::SpliceDonorVariant
            | Self::StopGained
            | Self::FrameshiftVariant
            | Self::StopLost
            | Self::StartLost => Impact::High,

            Self::InframeInsertion
            | Self::InframeDeletion
            | Self::MissenseVariant
            | Self::ProteinAlteringVariant => Impact::Moderate,

            Self::SpliceRegionVariant
            | Self::StartRetainedVariant
            | Self::StopRetainedVariant
            | Self::IncompleteTerminalCodonVariant
            | Self::SynonymousVariant => Impact::Low,

            Self::CodingSequenceVariant
            | Self::FivePrimeUtrVariant
            | Self::ThreePrimeUtrVariant
            | Self::NonCodingTranscriptExonVariant
            | Self::IntronVariant
            | Self::UpstreamGeneVariant
            | Self::DownstreamGeneVariant
            | Self::IntergenicVariant => Impact::Modifier,
        }
    }

    /// Severity rank (lower = more severe). Matches VEP's ordering.
    pub fn severity_rank(&self) -> u8 {
        match self {
            Self::TranscriptAblation => 1,
            Self::SpliceAcceptorVariant => 2,
            Self::SpliceDonorVariant => 3,
            Self::StopGained => 4,
            Self::FrameshiftVariant => 5,
            Self::StopLost => 6,
            Self::StartLost => 7,
            Self::InframeInsertion => 8,
            Self::InframeDeletion => 9,
            Self::MissenseVariant => 10,
            Self::ProteinAlteringVariant => 11,
            Self::SpliceRegionVariant => 12,
            Self::StartRetainedVariant => 13,
            Self::StopRetainedVariant => 14,
            Self::IncompleteTerminalCodonVariant => 15,
            Self::SynonymousVariant => 16,
            Self::CodingSequenceVariant => 17,
            Self::FivePrimeUtrVariant => 18,
            Self::ThreePrimeUtrVariant => 19,
            Self::NonCodingTranscriptExonVariant => 20,
            Self::IntronVariant => 21,
            Self::UpstreamGeneVariant => 22,
            Self::DownstreamGeneVariant => 23,
            Self::IntergenicVariant => 24,
        }
    }

    /// SO term string as VEP would output it (e.g., `"missense_variant"`).
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::TranscriptAblation => "transcript_ablation",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::StopGained => "stop_gained",
            Self::FrameshiftVariant => "frameshift_variant",
            Self::StopLost => "stop_lost",
            Self::StartLost => "start_lost",
            Self::InframeInsertion => "inframe_insertion",
            Self::InframeDeletion => "inframe_deletion",
            Self::MissenseVariant => "missense_variant",
            Self::ProteinAlteringVariant => "protein_altering_variant",
            Self::SpliceRegionVariant => "splice_region_variant",
            Self::StartRetainedVariant => "start_retained_variant",
            Self::StopRetainedVariant => "stop_retained_variant",
            Self::IncompleteTerminalCodonVariant => "incomplete_terminal_codon_variant",
            Self::SynonymousVariant => "synonymous_variant",
            Self::CodingSequenceVariant => "coding_sequence_variant",
            Self::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Self::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Self::NonCodingTranscriptExonVariant => "non_coding_transcript_exon_variant",
            Self::IntronVariant => "intron_variant",
            Self::UpstreamGeneVariant => "upstream_gene_variant",
            Self::DownstreamGeneVariant => "downstream_gene_variant",
            Self::IntergenicVariant => "intergenic_variant",
        }
    }
}

/// Per-transcript consequence annotation for a variant.
///
/// Populated by [`annotate_snv`] for SNVs, [`annotate_deletion`] /
/// [`annotate_insertion`] for indels, and the internal `annotate`
/// dispatcher (called via [`VarEffect::annotate`](crate::VarEffect::annotate)).
/// `hgvs_c` and `hgvs_p` are `None` when the variant does not affect
/// the transcript or protein (UTR, intron, splice-region).
#[derive(Debug, Clone, PartialEq)]
pub struct ConsequenceResult {
    /// RefSeq transcript accession with version (e.g., `"NM_006772.2"`).
    pub transcript: String,
    /// HGNC gene symbol.
    pub gene_symbol: String,
    /// RefSeq protein accession (e.g., `"NP_006763.2"`), `None` for
    /// non-coding transcripts.
    pub protein_accession: Option<String>,
    /// SO consequence terms, ordered by severity (most severe first).
    pub consequences: Vec<Consequence>,
    /// Highest IMPACT among all consequences.
    pub impact: Impact,
    /// 1-based protein position of the affected residue.
    pub protein_start: Option<u32>,
    /// 1-based protein position end (equal to `protein_start` for SNVs).
    pub protein_end: Option<u32>,
    /// Ref/alt codons (e.g., `"cGt/cAt"` -- VEP capitalizes the changed base).
    pub codons: Option<String>,
    /// Ref/alt amino acids (e.g., `"R/W"`, `"R/*"` for stop, `"R"` for
    /// synonymous).
    pub amino_acids: Option<String>,
    /// Exon number `"N/total"` or `None` if intronic/intergenic.
    pub exon: Option<String>,
    /// Intron number `"N/total"` or `None` if exonic.
    pub intron: Option<String>,
    /// CDS position start (1-based, VEP convention). `None` for non-CDS
    /// locations.
    pub cds_position: Option<u32>,
    /// CDS position end (1-based, inclusive). For SNVs this equals
    /// `cds_position`. For indels this is the last affected CDS position.
    /// `None` for non-CDS locations.
    pub cds_position_end: Option<u32>,
    /// cDNA position start (1-based from transcript start, includes UTR).
    /// `None` for intronic, upstream, downstream, and splice variants.
    pub cdna_position: Option<u32>,
    /// cDNA position end (1-based, inclusive). For SNVs this equals
    /// `cdna_position`. For indels this is the last affected cDNA position.
    /// `None` for the same locations as `cdna_position`.
    pub cdna_position_end: Option<u32>,
    /// Transcript strand.
    pub strand: crate::types::Strand,
    /// Transcript biotype.
    pub biotype: Biotype,
    /// `true` if the transcript is MANE Select.
    pub is_mane_select: bool,
    /// `true` if the transcript is MANE Plus Clinical.
    pub is_mane_plus_clinical: bool,
    /// `true` if the transcript is RefSeq Select.
    pub is_refseq_select: bool,
    /// HGVS coding notation.
    pub hgvs_c: Option<String>,
    /// HGVS protein notation.
    pub hgvs_p: Option<String>,
    /// Whether the variant is predicted to trigger nonsense-mediated mRNA
    /// decay via the 50-nucleotide rule. `true` when a PTC (from
    /// `StopGained` or `FrameshiftVariant`) is >50 nt upstream of the last
    /// exon-exon junction in CDS coordinates. Used by downstream ACMG PVS1
    /// strength modulation -- NOT a consequence term (VEP's
    /// `NMD_transcript_variant` SO:0001621 is biotype-based, irrelevant for
    /// MANE/RefSeq Select transcripts).
    pub predicts_nmd: bool,
}

/// Annotate a variant against all overlapping transcripts.
///
/// Returns one [`ConsequenceResult`] per transcript. Handles SNVs,
/// pure insertions/deletions, boundary-spanning deletions, complex
/// delins, and MNVs.
///
/// The full VCF REF allele is verified against the FASTA before trimming,
/// catching anchor-base mismatches for insertions.
///
/// # Arguments
///
/// * `chrom` -- UCSC-style chromosome name
/// * `pos` -- 0-based genomic position
/// * `ref_allele` -- VCF REF allele bytes
/// * `alt_allele` -- VCF ALT allele bytes
/// * `store` -- Transcript store for overlap queries
/// * `fasta` -- Reference FASTA reader
///
/// # Errors
///
/// Returns [`VarEffectError::RefMismatch`] if the VCF REF does not match
/// the FASTA, or [`VarEffectError::Malformed`] for unsupported variant
/// types (REF == ALT).
pub(crate) fn annotate(
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    store: &TranscriptStore,
    fasta: &FastaReader,
) -> Result<Vec<ConsequenceResult>, VarEffectError> {
    // Verify the full VCF REF against FASTA before trimming. This catches
    // anchor-base mismatches for insertions (where trimmed_ref is empty).
    if !ref_allele.is_empty() && !fasta.verify_ref(chrom, pos, ref_allele)? {
        let expected = fasta.fetch_sequence(chrom, pos, pos + ref_allele.len() as u64)?;
        return Err(VarEffectError::RefMismatch {
            chrom: chrom.to_string(),
            pos,
            expected: String::from_utf8_lossy(&expected).into_owned(),
            got: String::from_utf8_lossy(ref_allele).into_owned(),
        });
    }

    let (trimmed_ref, trimmed_alt, pos_adj) = helpers::trim_alleles(ref_allele, alt_allele);
    let trimmed_pos = pos + pos_adj;

    let results = match (trimmed_ref.len(), trimmed_alt.len()) {
        // SNV
        (1, 1) => {
            let overlaps = store.query_overlap(chrom, trimmed_pos, trimmed_pos + 1);
            let mut results = Vec::with_capacity(overlaps.len());
            for (tx, idx) in overlaps {
                // REF already verified above; use the unchecked path to
                // avoid a redundant FASTA seek per transcript.
                results.push(snv::annotate_snv_verified(
                    chrom,
                    trimmed_pos,
                    trimmed_ref[0],
                    trimmed_alt[0],
                    tx,
                    idx,
                    fasta,
                )?);
            }
            results
        }

        // Pure insertion
        (0, _) => {
            let query_start = trimmed_pos.saturating_sub(1);
            let query_end = trimmed_pos + 1;
            let overlaps = store.query_overlap(chrom, query_start, query_end);
            let mut results = Vec::with_capacity(overlaps.len());
            for (tx, idx) in overlaps {
                results.push(annotate_insertion(
                    chrom,
                    trimmed_pos,
                    trimmed_alt,
                    tx,
                    idx,
                    fasta,
                )?);
            }
            results
        }

        // Pure deletion
        (n, 0) if n > 0 => {
            let del_end = trimmed_pos + n as u64;
            let overlaps = store.query_overlap(chrom, trimmed_pos, del_end);
            let mut results = Vec::with_capacity(overlaps.len());
            for (tx, idx) in overlaps {
                results.push(annotate_deletion(
                    chrom,
                    trimmed_pos,
                    del_end,
                    trimmed_ref,
                    tx,
                    idx,
                    fasta,
                )?);
            }
            results
        }

        // MNV: same-length multi-base substitution
        (r, a) if r == a && r > 1 => {
            let query_end = trimmed_pos + r as u64;
            let overlaps = store.query_overlap(chrom, trimmed_pos, query_end);
            let mut results = Vec::with_capacity(overlaps.len());
            for (tx, idx) in overlaps {
                results.push(complex::annotate_mnv(
                    chrom,
                    trimmed_pos,
                    trimmed_ref,
                    trimmed_alt,
                    tx,
                    idx,
                    fasta,
                )?);
            }
            results
        }

        // Complex delins: different-length substitution
        (r, a) if r > 0 && a > 0 => {
            let query_end = trimmed_pos + r as u64;
            let overlaps = store.query_overlap(chrom, trimmed_pos, query_end);
            let mut results = Vec::with_capacity(overlaps.len());
            for (tx, idx) in overlaps {
                results.push(complex::annotate_complex_delins(
                    chrom,
                    trimmed_pos,
                    trimmed_ref,
                    trimmed_alt,
                    tx,
                    idx,
                    fasta,
                )?);
            }
            results
        }

        // REF == ALT after trimming
        _ => {
            return Err(VarEffectError::Malformed(
                "REF and ALT alleles are identical after trimming".to_string(),
            ));
        }
    };

    // No overlapping transcripts -- intergenic variant.
    if results.is_empty() {
        return Ok(vec![ConsequenceResult {
            transcript: String::new(),
            gene_symbol: String::new(),
            protein_accession: None,
            consequences: vec![Consequence::IntergenicVariant],
            impact: Impact::Modifier,
            protein_start: None,
            protein_end: None,
            codons: None,
            amino_acids: None,
            exon: None,
            intron: None,
            cds_position: None,
            cds_position_end: None,
            cdna_position: None,
            cdna_position_end: None,
            strand: crate::types::Strand::Plus,
            biotype: Biotype::Unknown,
            is_mane_select: false,
            is_mane_plus_clinical: false,
            is_refseq_select: false,
            hgvs_c: None,
            hgvs_p: None,
            predicts_nmd: false,
        }]);
    }

    Ok(results)
}
