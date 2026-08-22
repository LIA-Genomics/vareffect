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

mod breakend;
mod complex;
pub(crate) mod helpers;
mod indel;
mod nmd;
mod snv;
mod sv;
#[cfg(test)]
mod tests;

pub use breakend::{Breakend, BreakendMate, BreakendOrientation};
pub use indel::{annotate_deletion, annotate_insertion};
pub use snv::annotate_snv;
pub use sv::{SvConsequenceResult, SvKind};
pub(crate) use sv::{annotate_breakend, annotate_interval, annotate_sv_insertion};

use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::transcript::TranscriptStore;
use crate::types::Biotype;

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
    /// A deletion that removes an entire transcription unit (SO:0001893).
    /// Emitted for SV-shaped deletions that span the whole transcript.
    TranscriptAblation,
    /// A duplication/amplification that spans an entire transcription unit,
    /// raising its copy number (SO:0001889). Emitted for SV-shaped
    /// duplications that cover the whole transcript.
    TranscriptAmplification,
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
    /// A deletion that removes part (but not all) of a transcription unit
    /// (SO:0001906). Emitted for SV-shaped deletions that partially overlap a
    /// transcript. Promoted MODIFIER -> HIGH by Ensembl (release 111+).
    FeatureTruncation,
    /// A duplication/insertion that lengthens part of a transcription unit
    /// (SO:0001907). Emitted for SV-shaped duplications that partially overlap
    /// a transcript. Promoted MODIFIER -> HIGH by Ensembl (release 111+).
    FeatureElongation,

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
            | Self::TranscriptAmplification
            | Self::SpliceAcceptorVariant
            | Self::SpliceDonorVariant
            | Self::StopGained
            | Self::FrameshiftVariant
            | Self::StopLost
            | Self::StartLost
            | Self::FeatureTruncation
            | Self::FeatureElongation => Impact::High,

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
    ///
    /// `transcript_amplification` sits in the top HIGH block (right after
    /// `transcript_ablation`); `feature_truncation` / `feature_elongation`
    /// rank just above `intergenic_variant`, matching VEP's canonical SO
    /// severity list even though both are now HIGH-impact.
    pub fn severity_rank(&self) -> u8 {
        match self {
            Self::TranscriptAblation => 1,
            Self::TranscriptAmplification => 2,
            Self::SpliceAcceptorVariant => 3,
            Self::SpliceDonorVariant => 4,
            Self::StopGained => 5,
            Self::FrameshiftVariant => 6,
            Self::StopLost => 7,
            Self::StartLost => 8,
            Self::InframeInsertion => 9,
            Self::InframeDeletion => 10,
            Self::MissenseVariant => 11,
            Self::ProteinAlteringVariant => 12,
            Self::SpliceRegionVariant => 13,
            Self::StartRetainedVariant => 14,
            Self::StopRetainedVariant => 15,
            Self::IncompleteTerminalCodonVariant => 16,
            Self::SynonymousVariant => 17,
            Self::CodingSequenceVariant => 18,
            Self::FivePrimeUtrVariant => 19,
            Self::ThreePrimeUtrVariant => 20,
            Self::NonCodingTranscriptExonVariant => 21,
            Self::IntronVariant => 22,
            Self::UpstreamGeneVariant => 23,
            Self::DownstreamGeneVariant => 24,
            Self::FeatureTruncation => 25,
            Self::FeatureElongation => 26,
            Self::IntergenicVariant => 27,
        }
    }

    /// SO term string as VEP would output it (e.g., `"missense_variant"`).
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::TranscriptAblation => "transcript_ablation",
            Self::TranscriptAmplification => "transcript_amplification",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::StopGained => "stop_gained",
            Self::FrameshiftVariant => "frameshift_variant",
            Self::StopLost => "stop_lost",
            Self::StartLost => "start_lost",
            Self::FeatureTruncation => "feature_truncation",
            Self::FeatureElongation => "feature_elongation",
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

/// Where a variant's premature termination codon (PTC) sits, and whether that
/// codon is predicted to trigger nonsense-mediated mRNA decay.
///
/// Every verdict here is measured **at the termination codon**, unlike
/// [`ConsequenceResult::predicts_nmd`], which measures at the variant site for
/// Ensembl VEP `NMD.pm` parity. The termination codon is what the ClinGen SVI
/// PVS1 decision tree asks for (Abou Tayoun et al. 2018, PMID 30192042,
/// doi:10.1002/humu.23626), so this is the field to use for ACMG PVS1.
///
/// `#[non_exhaustive]` leaves room for future states (a distinguished
/// read-through / stop-loss terminus, say) without a SemVer break.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PtcStatus {
    /// No termination codon was evaluated: the consequence is neither
    /// `stop_gained` nor `frameshift_variant`, or the variant falls outside
    /// the CDS.
    ///
    /// This states scope, not absence. Canonical splice variants land here,
    /// and aberrant splicing routinely creates a downstream termination codon
    /// that this crate does not model -- so `NotApplicable` must not be read
    /// as "this variant is not truncating".
    NotApplicable,
    /// The variant is truncating but the termination codon could not be
    /// located -- the alternate allele yields no complete codon, produces no
    /// observable change, the alternate protein is shorter than the first
    /// changed residue, or the start codon is destroyed so no reading frame is
    /// defined.
    ///
    /// This is missing data. It is **not** a statement that NMD is absent, and
    /// must never be collapsed into one.
    Indeterminate,
    /// The alternate reading frame runs to the end of the 3'UTR without
    /// reaching a stop codon (HGVS `fsTer?`). There is no termination codon,
    /// so the 50-nucleotide rule does not apply and NMD is not predicted.
    NoStopCodon,
    /// The termination codon was located.
    At {
        /// 1-based residue position of the termination codon in the
        /// **alternate** protein.
        ///
        /// For a frameshift this is the residue the HGVS `fsTer` count points
        /// at: `p.Gln1756ProfsTer74` gives `1756 + 74 - 1 = 1829`. For a
        /// `stop_gained` it is the residue of the codon the variant turned
        /// into a stop.
        protein_position: u32,
        /// Whether this termination codon lies more than 50 nucleotides
        /// upstream of the last exon-exon junction in the CDS.
        ///
        /// The distance is measured from the codon's **3'-most base**, so the
        /// codon escapes as soon as any of its three bases falls inside the
        /// 3'-most 50 nucleotides of the penultimate coding exon -- the
        /// reading of the ClinGen SVI rule that does not overcall NMD.
        ///
        /// The last-junction rule is the **only** escape rule applied, because
        /// it is the only one the ClinGen SVI PVS1 decision tree uses. Two
        /// further rules from the NMD literature -- start-proximal escape
        /// (within 150 nt of the coding start site) and long-exon escape
        /// (a termination codon inside an exon longer than 407 bp), both after
        /// Lindeboom et al., PMID 27618451 and PMID 31659324 -- are
        /// deliberately **not** applied. Either would withdraw PVS1's full
        /// strength from a large class of established loss-of-function
        /// alleles, which is a guideline decision rather than an annotation
        /// one; `protein_position` is everything a consumer needs to layer
        /// them itself.
        ///
        /// Deliberately **not** named `predicts_nmd`: it answers the same
        /// question at a different position, and for a frameshift the two can
        /// disagree. See [`ConsequenceResult::predicts_nmd`].
        nmd_at_ptc: bool,
    },
}

impl PtcStatus {
    /// 1-based alternate-protein residue of the termination codon.
    ///
    /// # Returns
    ///
    /// `None` when no termination codon was located or none exists.
    pub fn protein_position(&self) -> Option<u32> {
        match self {
            Self::At {
                protein_position, ..
            } => Some(*protein_position),
            _ => None,
        }
    }

    /// NMD verdict measured at the termination codon.
    ///
    /// # Returns
    ///
    /// `Some(false)` for [`PtcStatus::NoStopCodon`] -- the alternate frame was
    /// fully scanned and no stop exists, so NMD cannot be triggered.
    /// `Some(v)` for [`PtcStatus::At`]. `None` for both
    /// [`PtcStatus::NotApplicable`] and [`PtcStatus::Indeterminate`], which
    /// are absences of an answer rather than a negative answer -- treating
    /// either as `false` would let a missense and a genuine last-exon
    /// termination codon take the same branch.
    ///
    /// **This accessor is lossy, and PVS1 consumers should match on the enum
    /// instead.** A `Some(false)` from [`PtcStatus::NoStopCodon`] and one from
    /// a last-exon [`PtcStatus::At`] are not the same clinical object: the
    /// first has no termination codon to assess at all (the transcript is a
    /// non-stop-decay substrate) and yields no `protein_position`, while the
    /// second carries the position the decision tree's next step needs. How a
    /// stop-less transcript should be scored is a clinical question this crate
    /// does not answer.
    pub fn nmd_at_ptc(&self) -> Option<bool> {
        match self {
            Self::At { nmd_at_ptc, .. } => Some(*nmd_at_ptc),
            Self::NoStopCodon => Some(false),
            Self::NotApplicable | Self::Indeterminate => None,
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
///
/// `#[non_exhaustive]` so future annotation fields can be added without a
/// SemVer break for downstream construction sites.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
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
    /// Whether the **variant site** is predicted to trigger nonsense-mediated
    /// mRNA decay under the 50-nucleotide rule.
    ///
    /// Measures the distance from the variant's own CDS position to the last
    /// exon-exon junction in the CDS, for `StopGained` and `FrameshiftVariant`
    /// consequences.
    ///
    /// Measuring at the variant site follows Ensembl VEP's `NMD.pm`, but the
    /// agreement is partial and worth stating exactly. Matching: the variant
    /// site as the anchor, the last-junction distance rule, and the intronless
    /// case. Diverging: `NMD.pm` also escapes a variant in the first ~100
    /// coding bases (its rule 3, which this crate does not implement), it
    /// walks all transcript exons where this crate uses CDS segments only, and
    /// its penultimate-exon window is 51 bases inclusive against the 50 used
    /// here. Both of the first two make this field predict NMD where `NMD.pm`
    /// escapes.
    ///
    /// For a `StopGained` variant the variant site is the termination codon,
    /// so this field is correct at codon granularity and is a sound PVS1
    /// input. **For a `FrameshiftVariant` they differ substantively**: the
    /// termination codon lies downstream of the variant, often in a later
    /// exon, and this flag ignores that displacement -- so it can report NMD
    /// for a frameshift whose actual termination codon escapes it. Use
    /// [`ConsequenceResult::ptc`] instead.
    ///
    /// The two anchors also differ by up to two bases even for `StopGained`:
    /// this field measures from the variant's own base, while
    /// [`PtcStatus::At::nmd_at_ptc`] measures from the 3'-most base of the
    /// termination codon. Where the 50-nucleotide boundary falls inside a
    /// codon, the two can disagree.
    ///
    /// Not a consequence term: VEP's `NMD_transcript_variant` (SO:0001621) is
    /// biotype-based and irrelevant for MANE / RefSeq Select transcripts.
    pub predicts_nmd: bool,
    /// The premature termination codon's position and its own NMD verdict.
    ///
    /// Populated for `stop_gained` and `frameshift_variant`;
    /// [`PtcStatus::NotApplicable`] otherwise. Prefer this over
    /// [`ConsequenceResult::predicts_nmd`] for ACMG PVS1 -- it measures at the
    /// termination codon, which is what the decision tree asks for.
    pub ptc: PtcStatus,
}

impl ConsequenceResult {
    /// A minimal result for one transcript, with no positional annotation.
    ///
    /// Every optional field is `None`, `consequences` is empty, `impact` is
    /// [`Impact::Modifier`], and both NMD fields carry no information
    /// (`predicts_nmd` is `false`, `ptc` is [`PtcStatus::NotApplicable`]).
    /// Callers fill in what applies.
    ///
    /// `ConsequenceResult` is `#[non_exhaustive]`, so this is the only way to
    /// build one outside this crate.
    ///
    /// # Arguments
    ///
    /// * `transcript` -- RefSeq accession with version.
    /// * `gene_symbol` -- HGNC gene symbol.
    /// * `strand` -- transcript strand.
    /// * `biotype` -- transcript biotype.
    pub fn new(
        transcript: impl Into<String>,
        gene_symbol: impl Into<String>,
        strand: crate::types::Strand,
        biotype: Biotype,
    ) -> Self {
        Self {
            transcript: transcript.into(),
            gene_symbol: gene_symbol.into(),
            protein_accession: None,
            consequences: Vec::new(),
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
            strand,
            biotype,
            is_mane_select: false,
            is_mane_plus_clinical: false,
            is_refseq_select: false,
            hgvs_c: None,
            hgvs_p: None,
            predicts_nmd: false,
            ptc: PtcStatus::NotApplicable,
        }
    }
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
            ptc: PtcStatus::NotApplicable,
        }]);
    }

    Ok(results)
}
