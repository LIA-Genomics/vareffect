//! Shared test transcript model builders.
//!
//! These fixtures are used by both `locate` and `consequence` test suites.

use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};

/// A model built with [`TranscriptModel::new`] from exon spans `(start,
/// end)` and CDS segments `(start, end, phase)`, both 0-based half-open and
/// in transcript order; each segment lies in the exon containing it.
pub fn model(
    accession: &str,
    gene_symbol: &str,
    chrom: &str,
    strand: Strand,
    exons: &[(u64, u64)],
    cds: &[(u64, u64, u8)],
) -> TranscriptModel {
    let exons: Vec<Exon> = exons
        .iter()
        .zip(1..)
        .map(|(&(genomic_start, genomic_end), exon_number)| Exon {
            exon_number,
            genomic_start,
            genomic_end,
        })
        .collect();
    let cds_segments = cds
        .iter()
        .map(|&(genomic_start, genomic_end, phase)| CdsSegment {
            exon_index: exons
                .iter()
                .position(|e| e.genomic_start <= genomic_start && genomic_end <= e.genomic_end)
                .and_then(|i| u16::try_from(i).ok())
                .expect("CDS segment inside an exon"),
            genomic_start,
            genomic_end,
            phase,
        })
        .collect();
    TranscriptModel::new(accession, gene_symbol, chrom, strand, exons)
        .with_cds_segments(cds_segments)
}

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
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 3,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        utr_numbering_breaks: Default::default(),
        selenocysteine_codons: Vec::new(),
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
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 3,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        utr_numbering_breaks: Default::default(),
        selenocysteine_codons: Vec::new(),
    }
}

/// Single-exon plus-strand coding transcript (1 exon, chr3).
///
/// The lone exon is simultaneously the most-5' (index 0) and most-3' exon, so
/// any overlap flags both `overlaps_first_exon` and `overlaps_last_exon`. Used
/// to pin the intragenic predicate, which must rely on the exact bound test
/// rather than the exon flags.
///
/// ```text
/// Exon 0: [2000, 5000)  5'UTR [2000,2200) + CDS [2200,4800) + 3'UTR [4800,5000)
/// ```
pub fn single_exon_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_SINGLE.1".into(),
        protein_accession: Some("NP_TEST_SINGLE.1".into()),
        gene_symbol: "TESTSINGLE".into(),
        hgnc_id: Some("HGNC:99903".into()),
        ensembl_accession: None,
        chrom: "chr3".into(),
        strand: Strand::Plus,
        tx_start: 2_000,
        tx_end: 5_000,
        cds_genomic_start: Some(2_200),
        cds_genomic_end: Some(4_800),
        exons: vec![Exon {
            exon_number: 1,
            genomic_start: 2_000,
            genomic_end: 5_000,
        }],
        cds_segments: vec![CdsSegment {
            exon_index: 0,
            genomic_start: 2_200,
            genomic_end: 4_800,
            phase: 0,
        }],
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        utr_numbering_breaks: Default::default(),
        selenocysteine_codons: Vec::new(),
    }
}

/// Non-coding transcript (2 exons, chr2, no CDS).
pub fn noncoding_2_exon() -> TranscriptModel {
    model(
        "NR_TEST_NC.1",
        "TESTNC",
        "chr2",
        Strand::Plus,
        &[(100, 300), (400, 600)],
        &[],
    )
    .with_biotype(Biotype::NonCodingRna)
    .with_tier(TranscriptTier::ManeSelect)
}
