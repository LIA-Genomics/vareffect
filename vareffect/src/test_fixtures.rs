//! Shared test transcript model builders.
//!
//! These fixtures are used by both `locate` and `consequence` test suites.

use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};

/// Plus-strand coding transcript (3 exons, chr1).
///
/// ```text
/// Exon 0: [1000, 2000)  5'UTR [1000,1500) + CDS [1500,2000)
/// Intron 0: [2000, 3000)
/// Exon 1: [3000, 3500)  CDS only
/// Intron 1: [3500, 4000)
/// Exon 2: [4000, 5000)  CDS [4000,4500) + 3'UTR [4500,5000)
/// ```
///
/// CDS total: 500 + 500 + 500 = 1500 bases.
pub fn plus_strand_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_PLUS.1".into(),
        protein_accession: Some("NP_TEST_PLUS.1".into()),
        gene_symbol: "TESTPLUS".into(),
        hgnc_id: Some("HGNC:99901".into()),
        ensembl_accession: None,
        chrom: "chr1".into(),
        strand: Strand::Plus,
        tx_start: 1_000,
        tx_end: 5_000,
        cds_genomic_start: Some(1_500),
        cds_genomic_end: Some(4_500),
        exons: vec![
            Exon {
                exon_number: 1,
                genomic_start: 1_000,
                genomic_end: 2_000,
            },
            Exon {
                exon_number: 2,
                genomic_start: 3_000,
                genomic_end: 3_500,
            },
            Exon {
                exon_number: 3,
                genomic_start: 4_000,
                genomic_end: 5_000,
            },
        ],
        cds_segments: vec![
            CdsSegment {
                exon_index: 0,
                genomic_start: 1_500,
                genomic_end: 2_000,
                phase: 0,
            },
            CdsSegment {
                exon_index: 1,
                genomic_start: 3_000,
                genomic_end: 3_500,
                phase: 2,
            },
            CdsSegment {
                exon_index: 2,
                genomic_start: 4_000,
                genomic_end: 4_500,
                phase: 0,
            },
        ],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 3,
    }
}

/// Minus-strand coding transcript (3 exons, chr17).
///
/// ```text
/// Exon 0 (5'-most): [18000, 20000)  5'UTR [19500,20000) + CDS [18000,19500)
/// Intron 0: genomic [16000, 18000)
/// Exon 1:            [14000, 16000)  CDS only
/// Intron 1: genomic [12000, 14000)
/// Exon 2 (3'-most): [10000, 12000)  CDS [11000,12000) + 3'UTR [10000,11000)
/// ```
///
/// CDS total: 1500 + 2000 + 1000 = 4500 bases.
pub fn minus_strand_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_MINUS.1".into(),
        protein_accession: Some("NP_TEST_MINUS.1".into()),
        gene_symbol: "TESTMINUS".into(),
        hgnc_id: Some("HGNC:99902".into()),
        ensembl_accession: None,
        chrom: "chr17".into(),
        strand: Strand::Minus,
        tx_start: 10_000,
        tx_end: 20_000,
        cds_genomic_start: Some(11_000),
        cds_genomic_end: Some(19_500),
        exons: vec![
            Exon {
                exon_number: 1,
                genomic_start: 18_000,
                genomic_end: 20_000,
            },
            Exon {
                exon_number: 2,
                genomic_start: 14_000,
                genomic_end: 16_000,
            },
            Exon {
                exon_number: 3,
                genomic_start: 10_000,
                genomic_end: 12_000,
            },
        ],
        cds_segments: vec![
            CdsSegment {
                exon_index: 0,
                genomic_start: 18_000,
                genomic_end: 19_500,
                phase: 0,
            },
            CdsSegment {
                exon_index: 1,
                genomic_start: 14_000,
                genomic_end: 16_000,
                phase: 0,
            },
            CdsSegment {
                exon_index: 2,
                genomic_start: 11_000,
                genomic_end: 12_000,
                phase: 0,
            },
        ],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 3,
    }
}

/// Non-coding transcript (2 exons, chr2, no CDS).
pub fn noncoding_2_exon() -> TranscriptModel {
    TranscriptModel {
        accession: "NR_TEST_NC.1".into(),
        protein_accession: None,
        gene_symbol: "TESTNC".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr2".into(),
        strand: Strand::Plus,
        tx_start: 100,
        tx_end: 600,
        cds_genomic_start: None,
        cds_genomic_end: None,
        exons: vec![
            Exon {
                exon_number: 1,
                genomic_start: 100,
                genomic_end: 300,
            },
            Exon {
                exon_number: 2,
                genomic_start: 400,
                genomic_end: 600,
            },
        ],
        cds_segments: vec![],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::NonCodingRna,
        exon_count: 2,
    }
}
