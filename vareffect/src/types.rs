//! Runtime types for the transcript model store.
//!
//! All coordinates use the **0-based, half-open** convention (BED/UCSC style).
//! The GFF3 parser in `vareffect-cli` converts NCBI's 1-based fully-closed
//! coordinates on ingest, so consumers of this crate never see 1-based indices.
//!
//! Most types in this module round-trip through MessagePack via
//! [`serde::Serialize`]/[`serde::Deserialize`]. [`Biotype`] has a hand-written
//! `Serialize`/`Deserialize` so the on-disk format stays a single flat string
//! and unknown upstream labels (`vault_RNA`, future biotypes) survive as
//! [`Biotype::Other`] without schema changes. The exception is
//! [`GenomicVariant`], a non-serialized in-memory value type (the resolver
//! output shared by the HGVS `c.`/`g.` reverse mappers).

use serde::{Deserialize, Serialize};

/// A genomic variant in VCF-style coordinates.
///
/// The shared output type of the HGVS reverse resolvers
/// ([`crate::VarEffect::resolve_hgvs_c`] and [`crate::VarEffect::resolve_hgvs_g`]).
/// All coordinates are 0-based (matching the internal convention used
/// throughout `vareffect`). Alleles are plus-strand, uppercase ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicVariant {
    /// UCSC-style chromosome name (e.g., `"chr17"`).
    pub chrom: String,
    /// 0-based genomic position of the first REF base.
    pub pos: u64,
    /// Reference allele (plus-strand, uppercase ASCII).
    pub ref_allele: Vec<u8>,
    /// Alternate allele (plus-strand, uppercase ASCII).
    pub alt_allele: Vec<u8>,
}

/// Transcript strand orientation relative to the reference genome.
///
/// Marked `#[non_exhaustive]` so future assemblies with unknown / ambiguous
/// strand annotations can extend this enum without a SemVer break.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Strand {
    /// Plus strand (5'→3' runs in the direction of increasing genomic coordinate).
    Plus,
    /// Minus strand (5'→3' runs in the direction of decreasing genomic coordinate).
    Minus,
}

/// Curation tier that produced a transcript model.
///
/// `ManeSelect` and `ManePlusClinical` are the primary clinical tiers;
/// `RefSeqSelect` is provided as a fallback for genes without MANE coverage.
/// `#[non_exhaustive]` leaves room for future tiers (Ensembl Canonical,
/// CCDS-only, …) without breaking downstream matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TranscriptTier {
    /// NCBI/Ensembl jointly-curated MANE Select transcript (the default
    /// clinical reference isoform for each gene).
    ManeSelect,
    /// MANE Plus Clinical — a second isoform curated for clinically
    /// actionable variants not captured on the MANE Select isoform.
    ManePlusClinical,
    /// RefSeq Select — NCBI's fallback canonical transcript for genes
    /// without MANE coverage.
    RefSeqSelect,
}

/// A single exon within a [`TranscriptModel`].
///
/// Exons inside a `TranscriptModel` are ordered 5'→3' on the *transcript*,
/// which for minus-strand genes is the reverse of genomic order.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Exon {
    /// 1-based exon number, counting from the 5' end of the transcript.
    /// `u16` is sufficient — the largest known human transcript (TTN) has
    /// ~363 exons, far below the 65,535 ceiling.
    pub exon_number: u16,
    /// Genomic start coordinate, 0-based inclusive.
    pub genomic_start: u64,
    /// Genomic end coordinate, 0-based exclusive (half-open).
    pub genomic_end: u64,
}

/// Where a transcript's genome-built numbering, outside its CDS, stops
/// matching the transcript's own (see
/// [`TranscriptModel::utr_numbering_breaks`]). Each bound is the first
/// wrong number; `None` means every number of that kind is right.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct UtrNumberingBreaks {
    /// cDNA positions (HGVS n. on a non-coding transcript) from this one on.
    pub cdna_from: Option<u32>,
    /// HGVS c. 5'UTR positions `c.-N` (and their intronic offsets) with `N`
    /// at least this.
    pub five_prime_from: Option<u32>,
    /// HGVS c. 3'UTR positions `c.*N` with `N` at least this.
    pub three_prime_from: Option<u32>,
}

impl UtrNumberingBreaks {
    /// Breaks with the given bounds.
    ///
    /// # Arguments
    ///
    /// * `cdna_from`, `five_prime_from`, `three_prime_from` - See the fields
    pub fn new(
        cdna_from: Option<u32>,
        five_prime_from: Option<u32>,
        three_prime_from: Option<u32>,
    ) -> Self {
        Self {
            cdna_from,
            five_prime_from,
            three_prime_from,
        }
    }

    /// `true` when some numbering is affected.
    pub fn any(&self) -> bool {
        self.cdna_from.is_some()
            || self.five_prime_from.is_some()
            || self.three_prime_from.is_some()
    }
}

/// A single CDS segment within a [`TranscriptModel`].
///
/// One `CdsSegment` corresponds to one GFF3 `CDS` row. Segments are ordered
/// 5'→3' on the *transcript* (reversed for minus-strand genes), matching the
/// `exons` vector on [`TranscriptModel`]. The per-segment [`phase`] captures
/// the reading-frame offset that VEP needs for frameshift detection and p. HGVS
/// notation across exon boundaries.
///
/// `exon_index` is the 0-based index into `TranscriptModel::exons` of the exon
/// that contains this CDS segment — this lets downstream code walk codons
/// without re-scanning the exon vector for every CDS row.
///
/// [`phase`]: Self::phase
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CdsSegment {
    /// 0-based index into `TranscriptModel::exons` of the containing exon.
    pub exon_index: u16,
    /// Genomic start coordinate, 0-based inclusive.
    pub genomic_start: u64,
    /// Genomic end coordinate, 0-based exclusive (half-open).
    pub genomic_end: u64,
    /// GFF3 column-8 phase: `0`, `1`, or `2` — the number of bases at the
    /// *transcript-5'* end of this CDS segment that belong to the final codon
    /// of the previous segment. A missing GFF3 phase (`.`) is normalized to
    /// `0` at build time. Any value > 2 is rejected as malformed input.
    pub phase: u8,
}

/// Transcript biotype.
///
/// Known biotypes are enum variants for type safety; unrecognized upstream
/// labels (e.g. a future `Y_RNA` added to MANE) survive verbatim in
/// [`Biotype::Other`] and round-trip through the MessagePack store without
/// code changes.
///
/// The on-disk representation is a single flat string — known variants use
/// their canonical NCBI/Ensembl label, and `Other` stores the raw upstream
/// label. This is implemented via a hand-written `Serialize`/`Deserialize`
/// below rather than `#[serde(untagged)]` to guarantee the flat wire format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Biotype {
    /// `NM_*` — protein-coding mRNA.
    ProteinCoding,
    /// Generic non-coding RNA — the default for `NR_*` accessions without a
    /// more specific gene-level biotype.
    NonCodingRna,
    /// Long non-coding RNA.
    LncRna,
    /// Antisense RNA.
    AntisenseRna,
    /// Small nucleolar RNA.
    SnoRna,
    /// Small nuclear RNA.
    SnRna,
    /// RNase MRP RNA.
    RnaseMrpRna,
    /// Telomerase RNA component.
    TelomeraseRna,
    /// Vault RNA.
    VaultRna,
    /// Upstream biotype label not recognized by this crate. The raw string is
    /// preserved verbatim so future gene types survive without a schema change.
    Other(String),
    /// No biotype signal was available at build time (neither the accession
    /// prefix nor a gene-level `gene_biotype` attribute).
    Unknown,
}

impl Biotype {
    /// Return the canonical on-disk label for this biotype.
    pub fn as_str(&self) -> &str {
        match self {
            Self::ProteinCoding => "protein_coding",
            Self::NonCodingRna => "non_coding_rna",
            Self::LncRna => "lncRNA",
            Self::AntisenseRna => "antisense_RNA",
            Self::SnoRna => "snoRNA",
            Self::SnRna => "snRNA",
            Self::RnaseMrpRna => "RNase_MRP_RNA",
            Self::TelomeraseRna => "telomerase_RNA",
            Self::VaultRna => "vault_RNA",
            Self::Other(s) => s.as_str(),
            Self::Unknown => "unknown",
        }
    }

    /// Parse a biotype label from its canonical string form.
    ///
    /// Unknown labels are returned as [`Biotype::Other`] (preserving the raw
    /// string) rather than an error — the transcript model store is meant to
    /// round-trip whatever MANE/RefSeq ships today or tomorrow.
    pub fn from_label(label: &str) -> Self {
        match label {
            "protein_coding" => Self::ProteinCoding,
            "non_coding_rna" => Self::NonCodingRna,
            "lncRNA" => Self::LncRna,
            "antisense_RNA" => Self::AntisenseRna,
            "snoRNA" => Self::SnoRna,
            "snRNA" => Self::SnRna,
            "RNase_MRP_RNA" => Self::RnaseMrpRna,
            "telomerase_RNA" => Self::TelomeraseRna,
            "vault_RNA" => Self::VaultRna,
            "unknown" | "" => Self::Unknown,
            other => Self::Other(other.to_string()),
        }
    }

    /// `true` if the biotype is `ProteinCoding`. Convenience helper for
    /// downstream filters.
    pub fn is_protein_coding(&self) -> bool {
        matches!(self, Self::ProteinCoding)
    }
}

impl Serialize for Biotype {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for Biotype {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        // Deserialize as an owned `String` then map into the enum. Avoids the
        // lifetime dance of `&str` for formats that might not hand out borrows
        // (rmp-serde does, but JSON with escapes does not).
        let raw = String::deserialize(deserializer)?;
        Ok(Self::from_label(&raw))
    }
}

/// Annotation source a transcript store was built from.
///
/// Recorded once per store (see [`crate::TranscriptStore::source`]); the
/// per-source meaning of the [`TranscriptModel`] identifier fields is
/// documented there. Serialized as the lowercase label returned by
/// [`TranscriptSource::as_str`]. `#[non_exhaustive]` so further sources can
/// be added without a SemVer break.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[non_exhaustive]
pub enum TranscriptSource {
    /// MANE Select and MANE Plus Clinical (the NCBI RefSeq MANE GFF3).
    Mane,
    /// NCBI RefSeq, curated `NM_` / `NR_` transcripts on the primary assembly.
    RefSeq,
    /// GENCODE comprehensive gene annotation.
    Gencode,
}

impl TranscriptSource {
    /// Return the lowercase label used on disk and in output headers.
    ///
    /// # Returns
    ///
    /// `"mane"`, `"refseq"` or `"gencode"`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Mane => "mane",
            Self::RefSeq => "refseq",
            Self::Gencode => "gencode",
        }
    }
}

impl std::fmt::Display for TranscriptSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One transcript model from MANE, RefSeq or GENCODE.
///
/// Every coordinate is 0-based half-open. Accessions include the version
/// suffix (e.g., `"NM_006772.2"`) — version-less lookup is intentionally out
/// of scope; resolving an unversioned HGVS string to "the latest" belongs
/// with the consumer, not the raw store.
///
/// The identifier fields mean the same thing for every
/// [`TranscriptSource`]:
///
/// | Field | MANE / RefSeq | GENCODE |
/// |---|---|---|
/// | `accession` | `NM_` / `NR_` | `ENST` |
/// | `protein_accession` | `NP_` | `ENSP` |
/// | `ensembl_accession` | the `ENST` cross-reference | `None` |
/// | `hgnc_id` | `HGNC:n` | `HGNC:n` |
///
/// `#[non_exhaustive]`: construct with [`TranscriptModel::new`] and the
/// `with_*` setters. The fields stay public for reading.
///
/// # Examples
///
/// ```
/// use vareffect::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};
///
/// let exons = vec![Exon { exon_number: 1, genomic_start: 100, genomic_end: 400 }];
/// let cds = vec![CdsSegment { exon_index: 0, genomic_start: 150, genomic_end: 300, phase: 0 }];
/// let model = TranscriptModel::new("NM_000001.1", "GENE1", "chr1", Strand::Plus, exons)
///     .with_cds_segments(cds)
///     .with_biotype(Biotype::ProteinCoding)
///     .with_tier(TranscriptTier::ManeSelect);
/// assert_eq!((model.tx_start, model.tx_end), (100, 400));
/// assert_eq!(model.cds_genomic_start, Some(150));
/// assert!(model.cds_start_complete && model.cds_end_complete);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct TranscriptModel {
    /// Transcript accession including version: `NM_` / `NR_` for MANE and
    /// RefSeq (e.g. `"NM_006772.2"`), `ENST` for GENCODE.
    pub accession: String,
    /// Protein accession including version: `NP_` for MANE and RefSeq,
    /// `ENSP` for GENCODE. `None` for non-coding transcripts.
    pub protein_accession: Option<String>,
    /// HGNC-approved gene symbol, e.g. `"SYNGAP1"`.
    pub gene_symbol: String,
    /// HGNC identifier including the `HGNC:` prefix, e.g. `"HGNC:11497"`.
    /// `None` if the source did not carry an HGNC cross-reference.
    pub hgnc_id: Option<String>,
    /// Ensembl transcript cross-reference with version, e.g.
    /// `"ENST00000418600.6"`, for RefSeq-based models. `None` if the Dbxref
    /// did not list one, and always `None` for GENCODE (whose `accession` is
    /// the Ensembl ID).
    pub ensembl_accession: Option<String>,
    /// Chromosome name in UCSC style (`"chr1"`, …, `"chrX"`, `"chrY"`, `"chrM"`).
    /// For transcripts on GRCh38 patch sequences, this field holds the UCSC
    /// contig name exactly as published by MANE GFF3 column 1 (e.g.
    /// `"chr9_KN196479v1_fix"`, `"chr22_KI270879v1_alt"`). Against an NCBI
    /// RefSeq FASTA, patch lookups work only when
    /// [`crate::FastaReader::open_with_patch_aliases`] is supplied a
    /// `patch_chrom_aliases.csv` that maps the UCSC form back to the
    /// matching `NW_*`/`NT_*` RefSeq accession.
    pub chrom: String,
    /// Transcript strand orientation.
    pub strand: Strand,
    /// Genomic start of the entire transcript (0-based, inclusive).
    /// Includes 5' and 3' UTRs.
    pub tx_start: u64,
    /// Genomic end of the entire transcript (0-based, exclusive).
    /// Includes 5' and 3' UTRs.
    pub tx_end: u64,
    /// Lowest genomic coordinate of any CDS segment (0-based, inclusive).
    ///
    /// **This is a genomic interval span, not a transcript-relative start.**
    /// For a minus-strand transcript, this is the 3' end of the protein in
    /// transcript order. Use [`cds_segments`](Self::cds_segments) when you
    /// need the true 5' coding start or per-exon CDS bounds. `None` for
    /// non-coding transcripts.
    pub cds_genomic_start: Option<u64>,
    /// Highest genomic coordinate of any CDS segment (0-based, exclusive).
    /// See [`cds_genomic_start`](Self::cds_genomic_start) for the interval
    /// semantics caveat. `None` for non-coding transcripts.
    pub cds_genomic_end: Option<u64>,
    /// Exons ordered 5'→3' on the *transcript*.
    ///
    /// Invariant: for plus-strand transcripts, `exons[i].genomic_start <
    /// exons[i+1].genomic_start`; for minus-strand transcripts, the inverse
    /// holds (exon 1 has the highest genomic coordinates).
    pub exons: Vec<Exon>,
    /// CDS segments ordered 5'→3' on the *transcript*, matching
    /// [`exons`](Self::exons). Empty for non-coding transcripts.
    ///
    /// Each segment carries the containing exon index (for O(1) exon lookup
    /// without re-scanning) and the GFF3 column-8 phase (needed for codon
    /// walking across exon boundaries and frameshift detection).
    pub cds_segments: Vec<CdsSegment>,
    /// The highest curation tier the transcript carries (MANE Select >
    /// MANE Plus Clinical > RefSeq Select), `None` when it has none. Every
    /// MANE model has a tier. The store requires the key (nil allowed).
    #[serde(deserialize_with = "Option::deserialize")]
    pub tier: Option<TranscriptTier>,
    /// Biotype, either a known variant or a raw upstream label preserved
    /// verbatim via [`Biotype::Other`].
    pub biotype: Biotype,
    /// Total number of exons — always equal to `exons.len()`.
    ///
    /// Stored for O(1) access without traversing the vec. `u16` suffices
    /// (TTN, the longest known, is ~363).
    pub exon_count: u16,
    /// `false` when the source marks the CDS 5' end as incomplete (GENCODE
    /// `cds_start_NF`, a RefSeq `start_range` at the 5' end). Without a
    /// complete start, VEP calls no `start_lost` or `start_retained_variant`.
    /// A CDS whose first phase is non-zero is always incomplete (the store
    /// refuses the contradiction). `true` for non-coding transcripts.
    pub cds_start_complete: bool,
    /// `false` when the source marks the CDS 3' end as incomplete (GENCODE
    /// `cds_end_NF`, a RefSeq `end_range` at the 3' end). Without a complete
    /// end, the CDS has no stop codon to lose or retain. `true` for
    /// non-coding transcripts.
    pub cds_end_complete: bool,
    /// `true` when the transcript sequence differs from the reference genome
    /// (a RefSeq `exception` or a gapped genome alignment). The model is
    /// built from the genome.
    pub genome_mismatch: bool,
    /// `true` when, inside the CDS, the genome-built numbering is not the
    /// transcript's own: a RefSeq `cDNA_match` `I` gap, transcript bases
    /// aligned nowhere (between blocks, at an interior CDS range marker, or
    /// past a CDS end left open on the genome), or a `D` gap inside an exon.
    /// A `D` gap where RefSeq splits the exon (a frameshift intron) does not
    /// set it. Annotation then withholds everything numbered along the
    /// transcript and HGVS c. is not resolved on it (`VEP_DIVERGENCES.md` 3.3,
    /// [`crate::ConsequenceResult`]).
    pub cds_alignment_gap: bool,
    /// Where, outside the CDS, the genome-built numbering stops being the
    /// transcript's own (RefSeq alignment gaps in a UTR, or transcript bases
    /// before the alignment's first block). Numbers past a break are not
    /// reported nor resolved; CDS-relative numbering is unaffected. Defaults to
    /// empty when absent (the JSON test fixtures).
    #[serde(default)]
    pub utr_numbering_breaks: UtrNumberingBreaks,
    /// 0-based genomic start of each selenocysteine codon (the lowest base of
    /// the codon on either strand), in ascending order. Empty when the
    /// source names none. `TGA` there reads as `U` in the reference peptide;
    /// a base outside the CDS fails [`crate::LocateIndex::build`].
    pub selenocysteine_codons: Vec<u64>,
}

impl TranscriptModel {
    /// Create a model from its required parts.
    ///
    /// `tx_start` / `tx_end` are the exon span and `exon_count` is
    /// `exons.len()`. Everything else starts empty: no CDS, protein or
    /// cross-references, no tier, [`Biotype::Unknown`], a complete CDS at
    /// both ends and no genome mismatch. Set the rest with the `with_*`
    /// setters.
    ///
    /// # Arguments
    ///
    /// * `accession` - Versioned transcript accession
    /// * `gene_symbol` - Gene symbol
    /// * `chrom` - UCSC-style chromosome name
    /// * `strand` - Transcript strand
    /// * `exons` - Exons in transcript order (0-based half-open)
    ///
    /// # Returns
    ///
    /// The model. With no exons the span is `0..0`, which
    /// [`crate::TranscriptStore`] refuses; more than `u16::MAX` exons
    /// saturate `exon_count`, which the store also refuses.
    pub fn new(
        accession: impl Into<String>,
        gene_symbol: impl Into<String>,
        chrom: impl Into<String>,
        strand: Strand,
        exons: Vec<Exon>,
    ) -> Self {
        let tx_start = exons.iter().map(|e| e.genomic_start).min().unwrap_or(0);
        let tx_end = exons.iter().map(|e| e.genomic_end).max().unwrap_or(0);
        let exon_count = u16::try_from(exons.len()).unwrap_or(u16::MAX);
        Self {
            accession: accession.into(),
            protein_accession: None,
            gene_symbol: gene_symbol.into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: chrom.into(),
            strand,
            tx_start,
            tx_end,
            cds_genomic_start: None,
            cds_genomic_end: None,
            exons,
            cds_segments: Vec::new(),
            tier: None,
            biotype: Biotype::Unknown,
            exon_count,
            cds_start_complete: true,
            cds_end_complete: true,
            genome_mismatch: false,
            cds_alignment_gap: false,
            utr_numbering_breaks: UtrNumberingBreaks::default(),
            selenocysteine_codons: Vec::new(),
        }
    }

    /// Set the CDS segments and derive `cds_genomic_start` / `cds_genomic_end`
    /// from them (`None` for an empty list).
    ///
    /// # Arguments
    ///
    /// * `segments` - CDS segments in transcript order
    pub fn with_cds_segments(mut self, segments: Vec<CdsSegment>) -> Self {
        self.cds_genomic_start = segments.iter().map(|s| s.genomic_start).min();
        self.cds_genomic_end = segments.iter().map(|s| s.genomic_end).max();
        self.cds_segments = segments;
        self
    }

    /// Set the protein accession.
    ///
    /// # Arguments
    ///
    /// * `protein_accession` - Versioned `NP_` / `ENSP` accession, or `None`
    pub fn with_protein_accession(mut self, protein_accession: Option<String>) -> Self {
        self.protein_accession = protein_accession;
        self
    }

    /// Set the HGNC identifier.
    ///
    /// # Arguments
    ///
    /// * `hgnc_id` - `HGNC:n`, or `None`
    pub fn with_hgnc_id(mut self, hgnc_id: Option<String>) -> Self {
        self.hgnc_id = hgnc_id;
        self
    }

    /// Set the Ensembl transcript cross-reference.
    ///
    /// # Arguments
    ///
    /// * `ensembl_accession` - Versioned `ENST` accession, or `None`
    pub fn with_ensembl_accession(mut self, ensembl_accession: Option<String>) -> Self {
        self.ensembl_accession = ensembl_accession;
        self
    }

    /// Set the curation tier.
    ///
    /// # Arguments
    ///
    /// * `tier` - The highest tier the transcript carries
    pub fn with_tier(mut self, tier: TranscriptTier) -> Self {
        self.tier = Some(tier);
        self
    }

    /// Set the biotype.
    ///
    /// # Arguments
    ///
    /// * `biotype` - Transcript biotype
    pub fn with_biotype(mut self, biotype: Biotype) -> Self {
        self.biotype = biotype;
        self
    }

    /// Set whether the CDS is complete at its 5' and 3' ends.
    ///
    /// # Arguments
    ///
    /// * `start_complete` - `false` for a 5'-incomplete CDS (`cds_start_NF`)
    /// * `end_complete` - `false` for a 3'-incomplete CDS (`cds_end_NF`)
    pub fn with_cds_completeness(mut self, start_complete: bool, end_complete: bool) -> Self {
        self.cds_start_complete = start_complete;
        self.cds_end_complete = end_complete;
        self
    }

    /// Set the genome-mismatch flags.
    ///
    /// # Arguments
    ///
    /// * `genome_mismatch` - The transcript sequence differs from the genome
    /// * `cds_alignment_gap` - The genome alignment has a gap inside the CDS
    pub fn with_genome_mismatch(mut self, genome_mismatch: bool, cds_alignment_gap: bool) -> Self {
        self.genome_mismatch = genome_mismatch;
        self.cds_alignment_gap = cds_alignment_gap;
        self
    }

    /// Set [`Self::utr_numbering_breaks`].
    ///
    /// # Arguments
    ///
    /// * `breaks` - Where the cDNA and UTR numbering built from the genome
    ///   stops being the transcript's own
    pub fn with_utr_numbering_breaks(mut self, breaks: UtrNumberingBreaks) -> Self {
        self.utr_numbering_breaks = breaks;
        self
    }

    /// Set the selenocysteine codon positions.
    ///
    /// # Arguments
    ///
    /// * `codons` - 0-based genomic start of each Sec codon; sorted here
    pub fn with_selenocysteine_codons(mut self, mut codons: Vec<u64>) -> Self {
        codons.sort_unstable();
        self.selenocysteine_codons = codons;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn biotype_roundtrips_known_variants() {
        for variant in [
            Biotype::ProteinCoding,
            Biotype::NonCodingRna,
            Biotype::LncRna,
            Biotype::AntisenseRna,
            Biotype::SnoRna,
            Biotype::SnRna,
            Biotype::RnaseMrpRna,
            Biotype::TelomeraseRna,
            Biotype::VaultRna,
            Biotype::Unknown,
        ] {
            let encoded = rmp_serde::to_vec_named(&variant).unwrap();
            let decoded: Biotype = rmp_serde::from_slice(&encoded).unwrap();
            assert_eq!(decoded, variant);
        }
    }

    #[test]
    fn biotype_preserves_unknown_upstream_label() {
        let custom = Biotype::from_label("misc_RNA");
        assert!(matches!(&custom, Biotype::Other(s) if s == "misc_RNA"));
        let encoded = rmp_serde::to_vec_named(&custom).unwrap();
        let decoded: Biotype = rmp_serde::from_slice(&encoded).unwrap();
        assert_eq!(decoded, custom);
        assert_eq!(decoded.as_str(), "misc_RNA");
    }

    #[test]
    fn transcript_source_labels_match_serde() {
        for (source, label) in [
            (TranscriptSource::Mane, "mane"),
            (TranscriptSource::RefSeq, "refseq"),
            (TranscriptSource::Gencode, "gencode"),
        ] {
            assert_eq!(source.as_str(), label);
            assert_eq!(source.to_string(), label);
            assert_eq!(
                serde_json::to_value(source).expect("json"),
                serde_json::json!(label)
            );
        }
    }

    #[test]
    fn new_derives_span_and_defaults() {
        let exons = vec![
            Exon {
                exon_number: 1,
                genomic_start: 900,
                genomic_end: 1_000,
            },
            Exon {
                exon_number: 2,
                genomic_start: 100,
                genomic_end: 200,
            },
        ];
        let model = TranscriptModel::new("NM_1.1", "G", "chr1", Strand::Minus, exons);
        assert_eq!(
            (model.tx_start, model.tx_end, model.exon_count),
            (100, 1_000, 2)
        );
        assert_eq!(model.tier, None);
        assert_eq!(model.biotype, Biotype::Unknown);
        assert!(model.cds_start_complete && model.cds_end_complete);
        assert!(!model.genome_mismatch && !model.cds_alignment_gap);
        assert!(model.cds_segments.is_empty() && model.cds_genomic_start.is_none());
    }

    #[test]
    fn new_without_exons_has_an_empty_span() {
        let model = TranscriptModel::new("NM_1.1", "G", "chr1", Strand::Plus, Vec::new());
        assert_eq!((model.tx_start, model.tx_end, model.exon_count), (0, 0, 0));
    }

    #[test]
    fn setters_fill_every_optional_field() {
        let exon = Exon {
            exon_number: 1,
            genomic_start: 0,
            genomic_end: 500,
        };
        let segments = vec![CdsSegment {
            exon_index: 0,
            genomic_start: 100,
            genomic_end: 400,
            phase: 0,
        }];
        let model = TranscriptModel::new("ENST1.1", "G", "chr1", Strand::Plus, vec![exon])
            .with_cds_segments(segments)
            .with_protein_accession(Some("ENSP1.1".into()))
            .with_hgnc_id(Some("HGNC:1".into()))
            .with_ensembl_accession(None)
            .with_tier(TranscriptTier::ManePlusClinical)
            .with_biotype(Biotype::ProteinCoding)
            .with_cds_completeness(false, true)
            .with_genome_mismatch(true, true)
            .with_selenocysteine_codons(vec![300, 120]);
        assert_eq!(
            (model.cds_genomic_start, model.cds_genomic_end),
            (Some(100), Some(400))
        );
        assert_eq!(model.protein_accession.as_deref(), Some("ENSP1.1"));
        assert_eq!(model.hgnc_id.as_deref(), Some("HGNC:1"));
        assert_eq!(model.tier, Some(TranscriptTier::ManePlusClinical));
        assert!(!model.cds_start_complete && model.cds_end_complete);
        assert!(model.genome_mismatch && model.cds_alignment_gap);
        assert_eq!(model.selenocysteine_codons, vec![120, 300]);
        let cleared = model.with_cds_segments(Vec::new());
        assert_eq!(
            (cleared.cds_genomic_start, cleared.cds_genomic_end),
            (None, None)
        );
    }

    #[test]
    fn biotype_from_label_normalizes_empty_to_unknown() {
        assert_eq!(Biotype::from_label(""), Biotype::Unknown);
        assert_eq!(Biotype::from_label("unknown"), Biotype::Unknown);
    }
}
