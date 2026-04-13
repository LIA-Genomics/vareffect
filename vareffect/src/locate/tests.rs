use super::*;
use crate::test_fixtures::{minus_strand_coding, noncoding_2_exon, plus_strand_coding};
use crate::types::{Biotype, CdsSegment, Exon, TranscriptTier};

// -- Test fixtures (module-specific) -----------------------------------------

/// Single-exon coding transcript (no introns, no splice sites).
fn single_exon_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_SINGLE.1".into(),
        protein_accession: Some("NP_TEST_SINGLE.1".into()),
        gene_symbol: "TESTSINGLE".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr3".into(),
        strand: Strand::Plus,
        tx_start: 500,
        tx_end: 1_500,
        cds_genomic_start: Some(600),
        cds_genomic_end: Some(1_400),
        exons: vec![Exon {
            exon_number: 1,
            genomic_start: 500,
            genomic_end: 1_500,
        }],
        cds_segments: vec![CdsSegment {
            exon_index: 0,
            genomic_start: 600,
            genomic_end: 1_400,
            phase: 0,
        }],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
    }
}

fn locate_plus(pos: u64) -> VariantLocation {
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    locate_variant("chr1", pos, pos + 1, &tx, &idx).unwrap()
}

fn locate_minus(pos: u64) -> VariantLocation {
    let tx = minus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    locate_variant("chr17", pos, pos + 1, &tx, &idx).unwrap()
}

// -- Upstream / Downstream / Distal -----------------------------------------

#[test]
fn plus_upstream() {
    assert_eq!(
        locate_plus(500),
        VariantLocation::Upstream { distance: 500 },
    );
}

#[test]
fn plus_downstream() {
    assert_eq!(
        locate_plus(5_500),
        VariantLocation::Downstream { distance: 501 },
    );
}

#[test]
fn plus_distal() {
    assert_eq!(locate_plus(15_000), VariantLocation::Distal);
}

// -- UTR (plus strand) ------------------------------------------------------

#[test]
fn plus_5prime_utr() {
    assert_eq!(
        locate_plus(1_200),
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -300,
            is_splice_region: false,
        },
    );
}

#[test]
fn plus_3prime_utr() {
    assert_eq!(
        locate_plus(4_600),
        VariantLocation::ThreePrimeUtr {
            exon_index: 2,
            offset_from_cds_end: 101,
            is_splice_region: false,
        },
    );
}

// -- CDS exon positions (plus strand) ---------------------------------------

#[test]
fn plus_cds_first_base() {
    assert_eq!(
        locate_plus(1_500),
        VariantLocation::CdsExon {
            exon_index: 0,
            cds_segment_index: 0,
            cds_offset: 0,
            codon_number: 1,
            codon_position: 0,
            is_splice_region: false,
        },
    );
}

#[test]
fn plus_cds_codon_position() {
    assert_eq!(
        locate_plus(1_502),
        VariantLocation::CdsExon {
            exon_index: 0,
            cds_segment_index: 0,
            cds_offset: 2,
            codon_number: 1,
            codon_position: 2,
            is_splice_region: false,
        },
    );
}

#[test]
fn plus_cds_second_exon() {
    // Segment 0 = 500 bases, intra = 50. cds_offset = 550, codon 184, pos 1.
    assert_eq!(
        locate_plus(3_050),
        VariantLocation::CdsExon {
            exon_index: 1,
            cds_segment_index: 1,
            cds_offset: 550,
            codon_number: 184,
            codon_position: 1,
            is_splice_region: false,
        },
    );
}

#[test]
fn plus_cds_last_base() {
    assert_eq!(
        locate_plus(4_499),
        VariantLocation::CdsExon {
            exon_index: 2,
            cds_segment_index: 2,
            cds_offset: 1_499,
            codon_number: 500,
            codon_position: 2,
            is_splice_region: false,
        },
    );
}

// -- Splice donor/acceptor (plus strand) ------------------------------------

#[test]
fn plus_splice_donor_1() {
    assert_eq!(
        locate_plus(2_000),
        VariantLocation::SpliceDonor {
            intron_index: 0,
            offset: 1,
        },
    );
}

#[test]
fn plus_splice_donor_2() {
    assert_eq!(
        locate_plus(2_001),
        VariantLocation::SpliceDonor {
            intron_index: 0,
            offset: 2,
        },
    );
}

#[test]
fn plus_splice_acceptor_1() {
    assert_eq!(
        locate_plus(2_999),
        VariantLocation::SpliceAcceptor {
            intron_index: 0,
            offset: 1,
        },
    );
}

#[test]
fn plus_splice_acceptor_2() {
    assert_eq!(
        locate_plus(2_998),
        VariantLocation::SpliceAcceptor {
            intron_index: 0,
            offset: 2,
        },
    );
}

// -- Splice region (plus strand) --------------------------------------------

#[test]
fn plus_splice_region_donor() {
    assert_eq!(
        locate_plus(2_005),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Donor,
            distance: 6,
        },
    );
}

#[test]
fn plus_splice_region_acceptor() {
    assert_eq!(
        locate_plus(2_995),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Acceptor,
            distance: -5,
        },
    );
}

// -- Deep intron (plus strand) ----------------------------------------------

#[test]
fn plus_deep_intron() {
    // Acceptor closer (500) than donor (501) -> negative
    assert_eq!(
        locate_plus(2_500),
        VariantLocation::Intron {
            intron_index: 0,
            distance_to_nearest_exon: -500,
        },
    );
}

// -- Exonic splice region (plus strand) -------------------------------------

#[test]
fn plus_exonic_splice_region_donor() {
    let loc = locate_plus(1_998);
    match loc {
        VariantLocation::CdsExon {
            is_splice_region, ..
        } => assert!(is_splice_region),
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

#[test]
fn plus_exonic_splice_region_acceptor() {
    let loc = locate_plus(3_001);
    match loc {
        VariantLocation::CdsExon {
            is_splice_region, ..
        } => assert!(is_splice_region),
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

#[test]
fn plus_no_splice_region_first_exon_start() {
    // No intron before exon 0 -> no splice region at transcript start
    let loc = locate_plus(1_000);
    match loc {
        VariantLocation::FivePrimeUtr {
            is_splice_region, ..
        } => assert!(!is_splice_region),
        other => panic!("Expected FivePrimeUtr, got {:?}", other),
    }
}

// -- Minus-strand UTR and CDS -----------------------------------------------

#[test]
fn minus_5prime_utr() {
    assert_eq!(
        locate_minus(19_800),
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -301,
            is_splice_region: false,
        },
    );
}

#[test]
fn minus_3prime_utr() {
    assert_eq!(
        locate_minus(10_500),
        VariantLocation::ThreePrimeUtr {
            exon_index: 2,
            offset_from_cds_end: 500,
            is_splice_region: false,
        },
    );
}

#[test]
fn minus_cds_offset() {
    assert_eq!(
        locate_minus(19_000),
        VariantLocation::CdsExon {
            exon_index: 0,
            cds_segment_index: 0,
            cds_offset: 499,
            codon_number: 167,
            codon_position: 1,
            is_splice_region: false,
        },
    );
}

// -- Minus-strand splice sites ----------------------------------------------

#[test]
fn minus_splice_donor() {
    assert_eq!(
        locate_minus(17_999),
        VariantLocation::SpliceDonor {
            intron_index: 0,
            offset: 1,
        },
    );
}

#[test]
fn minus_splice_acceptor() {
    assert_eq!(
        locate_minus(16_000),
        VariantLocation::SpliceAcceptor {
            intron_index: 0,
            offset: 1,
        },
    );
}

// -- Non-coding transcript --------------------------------------------------

#[test]
fn non_coding_exon() {
    let tx = noncoding_2_exon();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr2", 150, 151, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::NonCodingExon {
            exon_index: 0,
            is_splice_region: false,
        },
    );
}

#[test]
fn non_coding_intron() {
    let tx = noncoding_2_exon();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr2", 350, 351, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::NonCodingIntron {
            intron_index: 0,
            distance_to_nearest_exon: -50,
        },
    );
}

// -- Single-exon gene -------------------------------------------------------

#[test]
fn single_exon_no_splice() {
    let tx = single_exon_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr3", 700, 701, &tx, &idx).unwrap();
    match loc {
        VariantLocation::CdsExon {
            is_splice_region, ..
        } => assert!(!is_splice_region),
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

// -- Exon boundary semantics (half-open) ------------------------------------

#[test]
fn at_exon_boundary_inclusive() {
    let loc = locate_plus(3_000);
    match loc {
        VariantLocation::CdsExon { exon_index, .. } => {
            assert_eq!(exon_index, 1);
        }
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

#[test]
fn at_exon_boundary_exclusive() {
    // pos == exon.genomic_end -> NOT in exon (half-open) -> intronic
    let loc = locate_plus(2_000);
    match loc {
        VariantLocation::SpliceDonor { .. } => {}
        other => panic!("Expected SpliceDonor (intronic), got {:?}", other),
    }
}

// -- Splice region boundaries -----------------------------------------------

#[test]
fn splice_region_boundary_at_9() {
    // donor_dist = 9 > 8 -> NOT splice region -> deep intron
    let loc = locate_plus(2_008);
    match loc {
        VariantLocation::Intron { .. } => {}
        other => panic!("Expected Intron, got {:?}", other),
    }
}

#[test]
fn splice_region_boundary_at_8() {
    // donor_dist = 8 <= 8 -> SpliceRegion
    assert_eq!(
        locate_plus(2_007),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Donor,
            distance: 8,
        },
    );
}

#[test]
fn splice_region_at_3() {
    // donor_dist = 3 -> SpliceRegion
    assert_eq!(
        locate_plus(2_002),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Donor,
            distance: 3,
        },
    );
}

// -- Acceptor splice region (polypyrimidine tract, -3..-17) -----------------
//
// These tests cover VEP's `splice_polypyrimidine_tract_variant` (SO:0002169),
// which extends the acceptor-side splice region from -8 down to -17. Positions
// -9..-17 were previously deep intron in vareffect; they are now tagged as
// `SpliceRegion` (side = Acceptor) to match VEP.
//
// Fixture: `plus_strand_coding` has intron 0 at genomic [2000, 3000). The
// acceptor boundary is at pos 3000 (first base of exon 1), so pos 2983 is
// at acceptor distance 17 and pos 2982 is at acceptor distance 18.

#[test]
fn acceptor_splice_region_polypyrimidine_boundary_at_17() {
    // acceptor_dist = 3000 - 2983 = 17 <= 17 -> SpliceRegion Acceptor
    assert_eq!(
        locate_plus(2_983),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Acceptor,
            distance: -17,
        },
    );
}

#[test]
fn acceptor_splice_region_polypyrimidine_boundary_at_18() {
    // acceptor_dist = 18 > 17 -> deep intron
    let loc = locate_plus(2_982);
    match loc {
        VariantLocation::Intron { .. } => {}
        other => panic!("Expected Intron (deep), got {other:?}"),
    }
}

#[test]
fn acceptor_splice_region_mid_polypyrimidine_at_10() {
    // acceptor_dist = 10 (in the -9..-17 extended polypyrimidine range)
    assert_eq!(
        locate_plus(2_990),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Acceptor,
            distance: -10,
        },
    );
}

#[test]
fn minus_acceptor_splice_region_polypyrimidine_at_17() {
    // Fixture: `minus_strand_coding` intron 0 at genomic [16000, 18000),
    // minus-strand downstream exon is exons[1] = [14000, 16000). The
    // acceptor boundary (transcript-order) is at down_exon.genomic_end =
    // 16000, so pos 16016 -> acceptor_dist = 16016 - 16000 + 1 = 17.
    assert_eq!(
        locate_minus(16_016),
        VariantLocation::SpliceRegion {
            intron_index: 0,
            side: SpliceSide::Acceptor,
            distance: -17,
        },
    );
}

// -- CDS boundary exact -----------------------------------------------------

#[test]
fn cds_boundary_exact() {
    let loc = locate_plus(1_500);
    match loc {
        VariantLocation::CdsExon { cds_offset, .. } => {
            assert_eq!(cds_offset, 0);
        }
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

// -- Second intron ----------------------------------------------------------

#[test]
fn plus_intron_1_splice_donor() {
    assert_eq!(
        locate_plus(3_500),
        VariantLocation::SpliceDonor {
            intron_index: 1,
            offset: 1,
        },
    );
}

#[test]
fn plus_intron_1_splice_acceptor() {
    assert_eq!(
        locate_plus(3_999),
        VariantLocation::SpliceAcceptor {
            intron_index: 1,
            offset: 1,
        },
    );
}

// -- Minus-strand upstream/downstream ---------------------------------------

#[test]
fn minus_upstream() {
    assert_eq!(
        locate_minus(20_500),
        VariantLocation::Upstream { distance: 501 },
    );
}

#[test]
fn minus_downstream() {
    assert_eq!(
        locate_minus(9_500),
        VariantLocation::Downstream { distance: 500 },
    );
}

// -- Minus-strand CDS in second exon ----------------------------------------

#[test]
fn minus_cds_second_exon() {
    // Accumulated = 1500, intra = 999, cds_offset = 2499, codon 834, pos 0.
    assert_eq!(
        locate_minus(15_000),
        VariantLocation::CdsExon {
            exon_index: 1,
            cds_segment_index: 1,
            cds_offset: 2_499,
            codon_number: 834,
            codon_position: 0,
            is_splice_region: false,
        },
    );
}

// -- Minus-strand deep intron -----------------------------------------------

#[test]
fn minus_deep_intron() {
    // donor closer (1000) than acceptor (1001) -> positive
    assert_eq!(
        locate_minus(17_000),
        VariantLocation::Intron {
            intron_index: 0,
            distance_to_nearest_exon: 1_000,
        },
    );
}

// -- Non-coding splice sites ------------------------------------------------

#[test]
fn non_coding_splice_donor() {
    let tx = noncoding_2_exon();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr2", 300, 301, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::SpliceDonor {
            intron_index: 0,
            offset: 1,
        },
    );
}

#[test]
fn non_coding_splice_acceptor() {
    let tx = noncoding_2_exon();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr2", 399, 400, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::SpliceAcceptor {
            intron_index: 0,
            offset: 1,
        },
    );
}

#[test]
fn non_coding_exonic_splice_region() {
    let tx = noncoding_2_exon();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr2", 298, 299, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::NonCodingExon {
            exon_index: 0,
            is_splice_region: true,
        },
    );
}

// -- Display helpers --------------------------------------------------------

#[test]
fn format_exon_number_works() {
    assert_eq!(format_exon_number(0, 3), "1/3");
    assert_eq!(format_exon_number(4, 11), "5/11");
}

#[test]
fn format_intron_number_works() {
    assert_eq!(format_intron_number(0, 3), "1/2");
    assert_eq!(format_intron_number(3, 11), "4/10");
}

// -- Single-exon UTR checks -------------------------------------------------

#[test]
fn single_exon_5prime_utr() {
    let tx = single_exon_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr3", 550, 551, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -50,
            is_splice_region: false,
        },
    );
}

#[test]
fn single_exon_3prime_utr() {
    let tx = single_exon_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr3", 1_450, 1_451, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::ThreePrimeUtr {
            exon_index: 0,
            offset_from_cds_end: 51,
            is_splice_region: false,
        },
    );
}

// -- Minus-strand exonic splice region --------------------------------------

#[test]
fn minus_exonic_splice_region_donor() {
    // Donor side = low genomic (18000). pos 18001: dist = 2 <= 3 -> true.
    let loc = locate_minus(18_001);
    match loc {
        VariantLocation::CdsExon {
            is_splice_region, ..
        } => assert!(is_splice_region),
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

#[test]
fn minus_exonic_splice_region_acceptor() {
    // Acceptor side = high genomic (16000). pos 15999: dist = 1 <= 3 -> true.
    let loc = locate_minus(15_999);
    match loc {
        VariantLocation::CdsExon {
            is_splice_region, ..
        } => assert!(is_splice_region),
        other => panic!("Expected CdsExon, got {:?}", other),
    }
}

// -- Last exon 3' boundary: no splice region --------------------------------

#[test]
fn plus_last_exon_3prime_no_splice() {
    // Last exon -> no intron after -> no splice region
    let loc = locate_plus(4_999);
    match loc {
        VariantLocation::ThreePrimeUtr {
            is_splice_region, ..
        } => assert!(!is_splice_region),
        other => panic!("Expected ThreePrimeUtr, got {:?}", other),
    }
}

// -- Minus-strand distal ----------------------------------------------------

#[test]
fn minus_distal_upstream() {
    assert_eq!(locate_minus(30_000), VariantLocation::Distal);
}

#[test]
fn minus_distal_downstream() {
    assert_eq!(locate_minus(1_000), VariantLocation::Distal);
}

// -- Minus-strand: 5'UTR at exact CDS boundary ------------------------------

#[test]
fn minus_5prime_utr_at_cds_boundary() {
    // pos 19500 = cds_genomic_end -> first UTR base, offset = -1
    assert_eq!(
        locate_minus(19_500),
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -1,
            is_splice_region: false,
        },
    );
}

// -- Minus-strand: first CDS base -------------------------------------------

#[test]
fn minus_cds_first_base() {
    // First CDS base = cds_genomic_end - 1 = 19499. cds_offset = 0.
    assert_eq!(
        locate_minus(19_499),
        VariantLocation::CdsExon {
            exon_index: 0,
            cds_segment_index: 0,
            cds_offset: 0,
            codon_number: 1,
            codon_position: 0,
            is_splice_region: false,
        },
    );
}

// -- Multi-exon UTR tests ---------------------------------------------------

/// Plus-strand transcript with 5'UTR and 3'UTR each spanning 2 exons.
///
/// ```text
/// Exon 0: [1000, 1500)  entirely 5'UTR
/// Intron 0: [1500, 2000)
/// Exon 1: [2000, 3000)  5'UTR [2000, 2200) + CDS [2200, 3000)
/// Intron 1: [3000, 4000)
/// Exon 2: [4000, 5000)  CDS [4000, 4800) + 3'UTR [4800, 5000)
/// Intron 2: [5000, 5500)
/// Exon 3: [5500, 6000)  entirely 3'UTR
/// ```
fn plus_strand_multi_exon_utr() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_PLUS_MUTR.1".into(),
        protein_accession: Some("NP_TEST_PLUS_MUTR.1".into()),
        gene_symbol: "TESTPMUTR".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr1".into(),
        strand: Strand::Plus,
        tx_start: 1_000,
        tx_end: 6_000,
        cds_genomic_start: Some(2_200),
        cds_genomic_end: Some(4_800),
        exons: vec![
            Exon {
                exon_number: 1,
                genomic_start: 1_000,
                genomic_end: 1_500,
            },
            Exon {
                exon_number: 2,
                genomic_start: 2_000,
                genomic_end: 3_000,
            },
            Exon {
                exon_number: 3,
                genomic_start: 4_000,
                genomic_end: 5_000,
            },
            Exon {
                exon_number: 4,
                genomic_start: 5_500,
                genomic_end: 6_000,
            },
        ],
        cds_segments: vec![
            CdsSegment {
                exon_index: 1,
                genomic_start: 2_200,
                genomic_end: 3_000,
                phase: 0,
            },
            CdsSegment {
                exon_index: 2,
                genomic_start: 4_000,
                genomic_end: 4_800,
                phase: 0,
            },
        ],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 4,
    }
}

/// Minus-strand transcript with 5'UTR and 3'UTR each spanning 2 exons.
fn minus_strand_multi_exon_utr() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_MINUS_MUTR.1".into(),
        protein_accession: Some("NP_TEST_MINUS_MUTR.1".into()),
        gene_symbol: "TESTMMUTR".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr17".into(),
        strand: Strand::Minus,
        tx_start: 13_000,
        tx_end: 20_000,
        cds_genomic_start: Some(15_200),
        cds_genomic_end: Some(17_500),
        exons: vec![
            Exon {
                exon_number: 1,
                genomic_start: 19_000,
                genomic_end: 20_000,
            },
            Exon {
                exon_number: 2,
                genomic_start: 17_000,
                genomic_end: 18_000,
            },
            Exon {
                exon_number: 3,
                genomic_start: 15_000,
                genomic_end: 16_000,
            },
            Exon {
                exon_number: 4,
                genomic_start: 13_000,
                genomic_end: 14_000,
            },
        ],
        cds_segments: vec![
            CdsSegment {
                exon_index: 1,
                genomic_start: 17_000,
                genomic_end: 17_500,
                phase: 0,
            },
            CdsSegment {
                exon_index: 2,
                genomic_start: 15_200,
                genomic_end: 16_000,
                phase: 0,
            },
        ],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 4,
    }
}

#[test]
fn plus_multi_exon_5prime_utr() {
    let tx = plus_strand_multi_exon_utr();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr1", 1_200, 1_201, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -500,
            is_splice_region: false,
        },
    );
}

#[test]
fn plus_multi_exon_3prime_utr() {
    let tx = plus_strand_multi_exon_utr();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr1", 5_700, 5_701, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::ThreePrimeUtr {
            exon_index: 3,
            offset_from_cds_end: 401,
            is_splice_region: false,
        },
    );
}

#[test]
fn minus_multi_exon_5prime_utr() {
    let tx = minus_strand_multi_exon_utr();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr17", 19_500, 19_501, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::FivePrimeUtr {
            exon_index: 0,
            offset_from_cds_start: -1001,
            is_splice_region: false,
        },
    );
}

#[test]
fn minus_multi_exon_3prime_utr() {
    let tx = minus_strand_multi_exon_utr();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_variant("chr17", 13_500, 13_501, &tx, &idx).unwrap();
    assert_eq!(
        loc,
        VariantLocation::ThreePrimeUtr {
            exon_index: 3,
            offset_from_cds_end: 700,
            is_splice_region: false,
        },
    );
}

// -- locate_indel tests -----------------------------------------------------

#[test]
fn locate_indel_cds_deletion() {
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_503, 1_506, &tx, &idx).unwrap();
    assert_eq!(
        loc.region,
        IndelRegion::Cds {
            cds_offset_start: 3,
            cds_offset_end: 6,
        },
    );
    assert!(!loc.overlaps_splice_canonical);
    assert!(!loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(0));
}

#[test]
fn locate_indel_splice_donor_overlap() {
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_998, 2_002, &tx, &idx).unwrap();
    assert!(loc.overlaps_splice_canonical);
    assert!(loc.crosses_exon_boundary);
}

#[test]
fn locate_indel_exon_boundary_cross() {
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_999, 2_001, &tx, &idx).unwrap();
    assert!(loc.crosses_exon_boundary);
    assert_eq!(loc.region, IndelRegion::BoundarySpanning);
}

#[test]
fn locate_indel_insertion_in_cds() {
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_510, 1_510, &tx, &idx).unwrap();
    assert_eq!(
        loc.region,
        IndelRegion::Cds {
            cds_offset_start: 10,
            cds_offset_end: 10,
        },
    );
    assert!(!loc.overlaps_splice_canonical);
    assert!(!loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(0));
}

// -- Partial-overlap (edge-straddling) regression tests --------------------
//
// Regressions: deletions/delins whose genomic footprint
// extends past one transcript edge but retains at least one base inside
// the transcript used to surface as `Err(Malformed)` from
// `find_exon_or_intron`. They must now classify as `BoundarySpanning`.

#[test]
fn locate_indel_plus_straddles_5p_edge() {
    // `[995, 1005)` starts 5 bp upstream of `tx_start=1000` and ends
    // 5 bp into exon 0 `[1000, 2000)`.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 995, 1_005, &tx, &idx).unwrap();
    assert_eq!(loc.region, IndelRegion::BoundarySpanning);
    assert!(loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(0));
    assert!(loc.intron_index.is_none());
}

#[test]
fn locate_indel_plus_straddles_3p_edge() {
    // `[4995, 5005)` ends 5 bp past `tx_end=5000`; the in-range portion
    // lies in the last exon `[4000, 5000)`.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 4_995, 5_005, &tx, &idx).unwrap();
    assert_eq!(loc.region, IndelRegion::BoundarySpanning);
    assert!(loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(2));
    assert!(loc.intron_index.is_none());
}

#[test]
fn locate_indel_minus_straddles_low_genomic_edge() {
    // Minus-strand fixture: `tx_start=10000` is the 3'-UTR end in
    // transcript order, exon `[10000, 12000)` is the last exon
    // (index 2). `[9995, 10005)` straddles `tx_start`, so the in-range
    // anchor falls inside exon 2. Mirrors the SDHB NM_003000.3 shape
    // from the large-scale concordance corpus.
    let tx = minus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr17", 9_995, 10_005, &tx, &idx).unwrap();
    assert_eq!(loc.region, IndelRegion::BoundarySpanning);
    assert!(loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(2));
    assert!(loc.intron_index.is_none());
}

#[test]
fn locate_indel_minus_straddles_high_genomic_edge() {
    // `[19995, 20005)` straddles `tx_end=20000`; on this minus-strand
    // fixture the in-range anchor falls inside the first exon in
    // transcript order `[18000, 20000)`.
    let tx = minus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr17", 19_995, 20_005, &tx, &idx).unwrap();
    assert_eq!(loc.region, IndelRegion::BoundarySpanning);
    assert!(loc.crosses_exon_boundary);
    assert_eq!(loc.exon_index, Some(0));
    assert!(loc.intron_index.is_none());
}

// -- Regression: canonical splice overlap must use half-open -----------------
//
// These tests guard the locate-layer fix for canonical splice overlap: a previous
// `check_splice_overlap_for_range` helper used a symmetric ±2 buffer around
// the indel footprint and tested the exon boundary point against that buffer,
// which over-claimed canonical overlap at intronic +3 (plus strand) and the
// last two exonic bases. The fix routes `locate_indel` through
// `check_splice_overlap_detailed`, which does proper half-open interval
// intersection against `[exon_end, exon_end + 2)`. These tests lock that
// behavior so the helper cannot be re-introduced.
//
// Fixture `plus_strand_coding`: exon 0 `[1000, 2000)`, intron 0
// `[2000, 3000)`. Donor canonical = intronic `[2000, 2002)`, donor splice
// region = `[2002, 2008)`.

#[test]
fn locate_indel_plus_intronic_plus3_not_canonical() {
    // 1-bp deletion at intronic +3 (genomic 2002). Must be SpliceRegion,
    // NOT SpliceDonor.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 2_002, 2_003, &tx, &idx).unwrap();
    assert!(
        !loc.overlaps_splice_canonical,
        "intronic +3 must NOT be canonical: {loc:?}",
    );
    assert!(
        loc.overlaps_splice_region,
        "intronic +3 must be splice region: {loc:?}",
    );
    assert_eq!(loc.region, IndelRegion::Intron);
    assert!(loc.splice_detail.is_none());
}

#[test]
fn locate_indel_plus_last_exonic_base_not_canonical() {
    // 1-bp deletion at last exonic base (genomic 1999). Must be splice
    // region (exonic, within 3 of boundary) but NOT canonical.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_999, 2_000, &tx, &idx).unwrap();
    assert!(
        !loc.overlaps_splice_canonical,
        "last exonic base must NOT be canonical: {loc:?}",
    );
    assert!(
        loc.overlaps_splice_region,
        "last exonic base must be splice region: {loc:?}",
    );
}

#[test]
fn locate_indel_plus_second_to_last_exonic_not_canonical() {
    // 1-bp deletion at exonic position exon_end - 2 (genomic 1998).
    // Same bug class as the last-base case. Must NOT be canonical.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 1_998, 1_999, &tx, &idx).unwrap();
    assert!(
        !loc.overlaps_splice_canonical,
        "exonic -2 must NOT be canonical: {loc:?}",
    );
    assert!(
        loc.overlaps_splice_region,
        "exonic -2 must be splice region: {loc:?}",
    );
}

#[test]
fn locate_indel_plus_insertion_intronic_plus3_not_canonical() {
    // Insertion at intronic +3 (genomic 2002). Insertions with zero-width
    // footprints were also bitten by the symmetric-buffer bug.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 2_002, 2_002, &tx, &idx).unwrap();
    assert!(
        !loc.overlaps_splice_canonical,
        "insertion at intronic +3 must NOT be canonical: {loc:?}",
    );
    assert!(
        loc.overlaps_splice_region,
        "insertion at intronic +3 must be splice region: {loc:?}",
    );
}

#[test]
fn locate_indel_plus_canonical_still_fires() {
    // Positive control: 1-bp deletion at intronic +1 (genomic 2000). Must
    // stay canonical after the fix.
    let tx = plus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr1", 2_000, 2_001, &tx, &idx).unwrap();
    assert!(
        loc.overlaps_splice_canonical,
        "intronic +1 must be canonical: {loc:?}",
    );
    let detail = loc.splice_detail.as_ref().expect("detail populated");
    assert!(detail.overlaps_donor);
    assert!(!detail.overlaps_acceptor);
}

#[test]
fn locate_indel_minus_intronic_plus3_not_canonical() {
    // Minus-strand mirror: intron 0 is genomic `[16000, 18000)`; upstream
    // exon (transcript-order) is exon 0 at `[18000, 20000)`. Donor
    // canonical = `[17998, 18000)`. Intronic +3 in transcript terms =
    // genomic 17997.
    let tx = minus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr17", 17_997, 17_998, &tx, &idx).unwrap();
    assert!(
        !loc.overlaps_splice_canonical,
        "minus-strand intronic +3 must NOT be canonical: {loc:?}",
    );
    assert!(
        loc.overlaps_splice_region,
        "minus-strand intronic +3 must be splice region: {loc:?}",
    );
}

#[test]
fn locate_indel_minus_canonical_still_fires() {
    // Positive control minus strand: 1-bp deletion at intronic +1 (genomic
    // 17999). Donor canonical on minus strand = `[17998, 18000)`.
    let tx = minus_strand_coding();
    let idx = LocateIndex::build(&tx).unwrap();
    let loc = locate_indel("chr17", 17_999, 18_000, &tx, &idx).unwrap();
    assert!(
        loc.overlaps_splice_canonical,
        "minus-strand intronic +1 must be canonical: {loc:?}",
    );
    let detail = loc.splice_detail.as_ref().expect("detail populated");
    assert!(detail.overlaps_donor);
}
