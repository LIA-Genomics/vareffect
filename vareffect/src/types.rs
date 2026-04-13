//! Runtime types for the transcript model store.
//!
//! All coordinates use the **0-based, half-open** convention (BED/UCSC style).
//! The GFF3 parser in `vareffect-cli` converts NCBI's 1-based fully-closed
//! coordinates on ingest, so consumers of this crate never see 1-based indices.
//!
//! Every type in this module round-trips through MessagePack via
//! [`serde::Serialize`]/[`serde::Deserialize`]. [`Biotype`] has a hand-written
//! `Serialize`/`Deserialize` so the on-disk format stays a single flat string
//! and unknown upstream labels (`vault_RNA`, future biotypes) survive as
//! [`Biotype::Other`] without schema changes.

use serde::{Deserialize, Serialize};

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

/// A canonical transcript model sourced from MANE or RefSeq Select.
///
/// Every coordinate is 0-based half-open. Accessions include the version
/// suffix (e.g., `"NM_006772.2"`) — version-less lookup is intentionally out
/// of scope; resolving an unversioned HGVS string to "the latest" belongs
/// with the consumer, not the raw store.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptModel {
    /// RefSeq transcript accession including version, e.g. `"NM_006772.2"`.
    pub accession: String,
    /// RefSeq protein accession including version, e.g. `"NP_006763.2"`.
    /// `None` for non-coding transcripts (`NR_*` accessions).
    pub protein_accession: Option<String>,
    /// HGNC-approved gene symbol, e.g. `"SYNGAP1"`.
    pub gene_symbol: String,
    /// HGNC identifier including the `HGNC:` prefix, e.g. `"HGNC:11497"`.
    /// `None` if the parent gene row did not carry an HGNC cross-reference.
    pub hgnc_id: Option<String>,
    /// Ensembl transcript accession with version, e.g. `"ENST00000418600.6"`.
    /// `None` if the Dbxref did not list an Ensembl cross-reference.
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
    /// MANE / RefSeq Select curation tier.
    pub tier: TranscriptTier,
    /// Biotype, either a known variant or a raw upstream label preserved
    /// verbatim via [`Biotype::Other`].
    pub biotype: Biotype,
    /// Total number of exons — always equal to `exons.len()`.
    ///
    /// Stored for O(1) access without traversing the vec. `u16` suffices
    /// (TTN, the longest known, is ~363).
    pub exon_count: u16,
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
    fn biotype_from_label_normalizes_empty_to_unknown() {
        assert_eq!(Biotype::from_label(""), Biotype::Unknown);
        assert_eq!(Biotype::from_label("unknown"), Biotype::Unknown);
    }
}
