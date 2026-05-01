use super::complex::{annotate_complex_delins, annotate_mnv};
use super::helpers::AnnotateCtx;
use super::helpers::trim_alleles;
use super::*;
use crate::fasta::{FastaReader, write_genome_binary};
use crate::locate::{LocateIndex, VariantLocation, locate_variant};
use crate::test_fixtures::{minus_strand_coding, noncoding_2_exon, plus_strand_coding};
use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};
use tempfile::TempDir;

// -- Synthetic FASTA fixture ------------------------------------------------

/// Build a synthetic FASTA covering the test transcript models.
///
/// **chr1** (plus_strand_coding): 6000 bp.
///   CDS segments [1500,2000) + [3000,3500) + [4000,4500) = 1500 bp.
///   Start codon ATG at 1500-1502, body CGT (Arg), stop TAA at 4497-4499.
///
/// **chr17** (minus_strand_coding): 21000 bp.
///   CDS segments [18000,19500) + [14000,16000) + [11000,12000) = 4500 bp.
///   Coding strand: ATG + CGT repeating + stop TAA.
///   FASTA stores plus-strand (reverse complement of coding strand).
///
/// **chr2** (noncoding): 1000 bp of A's.
///
/// **chr3** (single_exon_coding): 2000 bp.
///   CDS [600,1400) = 800 bp. ATG at 600-602, stop TAA at 797-799.
///
/// **chrM** (mitochondrial): 1000 bp.
///   CDS [100,400) = 300 bp. ATG at 100-102, body TGA (Trp on mito),
///   stop TAA at 297-299.
fn write_test_fasta() -> (TempDir, FastaReader) {
    let tmp = TempDir::new().unwrap();

    // chr1: 6000 bp
    let mut chr1 = vec![b'A'; 6000];
    chr1[1500] = b'A';
    chr1[1501] = b'T';
    chr1[1502] = b'G';
    let mut i = 1503;
    while i + 2 < 2000 {
        chr1[i] = b'C';
        chr1[i + 1] = b'G';
        chr1[i + 2] = b'T';
        i += 3;
    }
    // Split codon alignment: offset 498-500 spans segments 0 and 1
    chr1[1998] = b'C';
    chr1[1999] = b'G';
    chr1[3000] = b'T';
    i = 3001;
    while i + 2 < 3500 {
        chr1[i] = b'C';
        chr1[i + 1] = b'G';
        chr1[i + 2] = b'T';
        i += 3;
    }
    while i < 3500 {
        chr1[i] = b'C';
        i += 1;
    }
    i = 4000;
    while i + 2 < 4497 {
        chr1[i] = b'C';
        chr1[i + 1] = b'G';
        chr1[i + 2] = b'T';
        i += 3;
    }
    while i < 4497 {
        chr1[i] = b'C';
        i += 1;
    }
    chr1[4497] = b'T';
    chr1[4498] = b'A';
    chr1[4499] = b'A';

    // chr17: 21000 bp (minus-strand)
    let mut chr17 = vec![b'A'; 21000];
    let mut coding = Vec::with_capacity(4500);
    coding.extend_from_slice(b"ATG");
    while coding.len() + 3 <= 4497 {
        coding.extend_from_slice(b"CGT");
    }
    while coding.len() < 4497 {
        coding.push(b'C');
    }
    coding.extend_from_slice(b"TAA");
    assert_eq!(coding.len(), 4500);

    // Place coding strand into FASTA (plus-strand) per CDS segment.
    // Minus-strand: coding offset 0 -> genomic 19499, complemented.
    for offset in 0u32..1500 {
        let genomic = 19499 - offset as usize;
        chr17[genomic] = crate::codon::complement(coding[offset as usize]);
    }
    for offset in 0u32..2000 {
        let genomic = 15999 - offset as usize;
        chr17[genomic] = crate::codon::complement(coding[(1500 + offset) as usize]);
    }
    for offset in 0u32..1000 {
        let genomic = 11999 - offset as usize;
        chr17[genomic] = crate::codon::complement(coding[(3500 + offset) as usize]);
    }

    let chr2 = vec![b'A'; 1000];

    // chr3: 2000 bp (single-exon coding)
    let mut chr3 = vec![b'A'; 2000];
    chr3[600] = b'A';
    chr3[601] = b'T';
    chr3[602] = b'G';
    i = 603;
    while i + 2 < 797 {
        chr3[i] = b'C';
        chr3[i + 1] = b'G';
        chr3[i + 2] = b'T';
        i += 3;
    }
    while i < 797 {
        chr3[i] = b'C';
        i += 1;
    }
    chr3[797] = b'T';
    chr3[798] = b'A';
    chr3[799] = b'A';

    // chrM: 1000 bp (mitochondrial)
    let mut chrm = vec![b'A'; 1000];
    chrm[100] = b'A';
    chrm[101] = b'T';
    chrm[102] = b'G';
    i = 103;
    while i + 2 < 297 {
        chrm[i] = b'T';
        chrm[i + 1] = b'G';
        chrm[i + 2] = b'A';
        i += 3;
    }
    while i < 297 {
        chrm[i] = b'T';
        i += 1;
    }
    chrm[297] = b'T';
    chrm[298] = b'A';
    chrm[299] = b'A';

    let contigs: Vec<(&str, &[u8])> = vec![
        ("chr1", &chr1),
        ("chr17", &chr17),
        ("chr2", &chr2),
        ("chr3", &chr3),
        ("chrM", &chrm),
    ];
    let bin_path = tmp.path().join("test.bin");
    let idx_path = tmp.path().join("test.bin.idx");
    write_genome_binary(&contigs, "test", &bin_path, &idx_path).unwrap();
    let reader = FastaReader::open_with_assembly(&bin_path, crate::Assembly::GRCh38).unwrap();
    (tmp, reader)
}

// -- Test transcript builders (module-specific) -----------------------------

/// Mitochondrial plus-strand coding transcript.
/// CDS [100, 400) = 300 bp = 100 codons.
fn mito_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_MITO.1".into(),
        protein_accession: Some("NP_TEST_MITO.1".into()),
        gene_symbol: "TESTMITO".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chrM".into(),
        strand: Strand::Plus,
        tx_start: 0,
        tx_end: 1000,
        cds_genomic_start: Some(100),
        cds_genomic_end: Some(400),
        exons: vec![Exon {
            exon_number: 1,
            genomic_start: 0,
            genomic_end: 1000,
        }],
        cds_segments: vec![CdsSegment {
            exon_index: 0,
            genomic_start: 100,
            genomic_end: 400,
            phase: 0,
        }],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        genome_transcript_divergent: false,
        translational_exception: None,
    }
}

/// Single-exon plus-strand coding transcript on chr3.
/// CDS [600, 1400) = 800 bp.
fn single_exon_coding() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_TEST_SE.1".into(),
        protein_accession: Some("NP_TEST_SE.1".into()),
        gene_symbol: "TESTSE".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr3".into(),
        strand: Strand::Plus,
        tx_start: 0,
        tx_end: 2000,
        cds_genomic_start: Some(600),
        cds_genomic_end: Some(1400),
        exons: vec![Exon {
            exon_number: 1,
            genomic_start: 0,
            genomic_end: 2000,
        }],
        cds_segments: vec![CdsSegment {
            exon_index: 0,
            genomic_start: 600,
            genomic_end: 1400,
            phase: 0,
        }],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        genome_transcript_divergent: false,
        translational_exception: None,
    }
}

fn build_index(tx: &TranscriptModel) -> LocateIndex {
    LocateIndex::build(tx).unwrap()
}

// -- SNV consequence tests --------------------------------------------------

#[test]
fn missense_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Position 1503: CDS offset 3, codon 2 pos 0. CGT->TGT = Arg->Cys.
    let result = annotate_snv("chr1", 1503, b'C', b'T', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::MissenseVariant));
    assert_eq!(result.impact, Impact::Moderate);
    assert_eq!(result.protein_start, Some(2));
    assert_eq!(result.amino_acids.as_deref(), Some("R/C"));
    assert_eq!(result.codons.as_deref(), Some("Cgt/Tgt"));
    assert_eq!(result.cds_position, Some(4));
}

#[test]
fn synonymous_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Position 1505: CDS offset 5, codon 2 pos 2. CGT->CGA = Arg->Arg.
    let result = annotate_snv("chr1", 1505, b'T', b'A', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SynonymousVariant)
    );
    assert_eq!(result.impact, Impact::Low);
    assert_eq!(result.amino_acids.as_deref(), Some("R"));
    assert_eq!(result.codons.as_deref(), Some("cgT/cgA"));
}

#[test]
fn stop_gained() {
    let (_tmp, fasta) = write_stop_gained_fasta();
    let tx = stop_gained_transcript();
    let idx = build_index(&tx);
    // CGA at offset 3-5 (genomic 103-105). CGA->TGA: C->T at 103.
    let result = annotate_snv("chr1", 103, b'C', b'T', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::StopGained));
    assert_eq!(result.impact, Impact::High);
    assert_eq!(result.amino_acids.as_deref(), Some("R/*"));
}

/// Tiny transcript + FASTA for stop_gained testing. CDS has CGA codon.
fn stop_gained_transcript() -> TranscriptModel {
    TranscriptModel {
        accession: "NM_STOP_GAINED.1".into(),
        protein_accession: Some("NP_STOP_GAINED.1".into()),
        gene_symbol: "TESTSTOP".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: "chr1".into(),
        strand: Strand::Plus,
        tx_start: 50,
        tx_end: 200,
        cds_genomic_start: Some(100),
        cds_genomic_end: Some(112),
        exons: vec![Exon {
            exon_number: 1,
            genomic_start: 50,
            genomic_end: 200,
        }],
        cds_segments: vec![CdsSegment {
            exon_index: 0,
            genomic_start: 100,
            genomic_end: 112,
            phase: 0,
        }],
        tier: TranscriptTier::ManeSelect,
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        genome_transcript_divergent: false,
        translational_exception: None,
    }
}

fn write_stop_gained_fasta() -> (TempDir, FastaReader) {
    let tmp = TempDir::new().unwrap();
    // CDS at [100,112): ATG CGA CGT TAA
    let mut seq = vec![b'A'; 300];
    seq[100] = b'A';
    seq[101] = b'T';
    seq[102] = b'G';
    seq[103] = b'C';
    seq[104] = b'G';
    seq[105] = b'A';
    seq[106] = b'C';
    seq[107] = b'G';
    seq[108] = b'T';
    seq[109] = b'T';
    seq[110] = b'A';
    seq[111] = b'A';

    let bin_path = tmp.path().join("stop.bin");
    let idx_path = tmp.path().join("stop.bin.idx");
    write_genome_binary(&[("chr1", seq.as_slice())], "test", &bin_path, &idx_path).unwrap();
    let reader = FastaReader::open_with_assembly(&bin_path, crate::Assembly::GRCh38).unwrap();
    (tmp, reader)
}

#[test]
fn stop_lost() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Stop codon TAA at CDS offsets 1497-1499, genomic 4497-4499.
    // T->C at offset 1497: TAA->CAA = Gln. StopLost.
    let result = annotate_snv("chr1", 4497, b'T', b'C', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::StopLost));
    assert_eq!(result.impact, Impact::High);
    assert_eq!(result.amino_acids.as_deref(), Some("*/Q"));
    assert_eq!(result.protein_start, Some(500));
}

#[test]
fn start_lost() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Start codon ATG at CDS 0-2, genomic 1500-1502. A->C: ATG->CTG = Leu.
    let result = annotate_snv("chr1", 1500, b'A', b'C', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::StartLost));
    assert_eq!(result.impact, Impact::High);
    assert_eq!(result.protein_start, Some(1));
}

#[test]
fn stop_retained_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Stop TAA at offsets 1497-1499. A->G at offset 1499: TAA->TAG = stop.
    let result = annotate_snv("chr1", 4499, b'A', b'G', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::StopRetainedVariant)
    );
    assert_eq!(result.impact, Impact::Low);
}

#[test]
fn splice_donor_maps_correctly() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Intron 0 starts at genomic 2000 (after exon 0). Donor +1 = 2000.
    let result = annotate_snv("chr1", 2000, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
    );
    assert_eq!(result.impact, Impact::High);
    assert!(result.intron.is_some());
}

#[test]
fn splice_acceptor_maps_correctly() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Intron 0 ends before exon 1 at 3000. Acceptor -1 = 2999.
    let result = annotate_snv("chr1", 2999, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant)
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn intron_variant_maps_correctly() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_snv("chr1", 2500, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::IntronVariant));
    assert_eq!(result.impact, Impact::Modifier);
    assert!(result.intron.is_some());
    assert!(result.cds_position.is_none());
    assert!(result.cdna_position.is_none());
}

#[test]
fn five_prime_utr_maps_correctly() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_snv("chr1", 1200, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FivePrimeUtrVariant)
    );
    assert_eq!(result.impact, Impact::Modifier);
    assert!(result.exon.is_some());
}

#[test]
fn three_prime_utr_maps_correctly() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_snv("chr1", 4600, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant)
    );
    assert_eq!(result.impact, Impact::Modifier);
}

#[test]
fn compound_missense_splice_region() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Exon 1 starts at 3000. First 3 bases are exonic splice region.
    // CDS offset 500 at pos 3000: codon 167, position 2. CGT->CGA = Arg (synonymous + splice_region).
    let result = annotate_snv("chr1", 3000, b'T', b'A', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SynonymousVariant)
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceRegionVariant)
    );
    assert_eq!(result.impact, Impact::Low);
}

#[test]
fn upstream_gene_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_snv("chr1", 500, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::UpstreamGeneVariant)
    );
    assert_eq!(result.impact, Impact::Modifier);
}

#[test]
fn downstream_gene_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_snv("chr1", 5500, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::DownstreamGeneVariant)
    );
    assert_eq!(result.impact, Impact::Modifier);
}

#[test]
fn non_coding_exon_variant() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = noncoding_2_exon();
    let idx = build_index(&tx);
    let result = annotate_snv("chr2", 200, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant)
    );
    assert_eq!(result.impact, Impact::Modifier);
}

#[test]
fn minus_strand_complement() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    // CDS offset 3 = codon 2, pos 0, at genomic 19496.
    // VCF ref = G (plus-strand), alt = A. Coding-strand: C->T. CGT->TGT = R->C.
    let result = annotate_snv("chr17", 19496, b'G', b'A', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::MissenseVariant));
    assert_eq!(result.amino_acids.as_deref(), Some("R/C"));
    assert_eq!(result.protein_start, Some(2));
}

#[test]
fn codon_spanning_exon_boundary() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Codon at CDS offsets 498-500 spans segments 0 and 1.
    // Positions 1998=C, 1999=G, 3000=T. Change G->A at 1999: CGT->CAT = R->H.
    let result = annotate_snv("chr1", 1999, b'G', b'A', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::MissenseVariant));
    assert_eq!(result.amino_acids.as_deref(), Some("R/H"));
    assert_eq!(result.codons.as_deref(), Some("cGt/cAt"));
    assert_eq!(result.cds_position, Some(500));
}

#[test]
fn mitochondrial_tga_is_trp() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = mito_coding();
    let idx = build_index(&tx);
    // TGA on chrM = Trp (W), NOT stop.
    let loc = locate_variant("chrM", 103, 104, &tx, &idx).unwrap();
    match loc {
        VariantLocation::CdsExon { codon_number, .. } => {
            assert_eq!(codon_number, 2);
        }
        other => panic!("expected CdsExon, got {:?}", other),
    }
    // T->A at 103 (codon pos 0): TGA->AGA. On mito AGA = stop.
    let result = annotate_snv("chrM", 103, b'T', b'A', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::StopGained));

    // TGA->TCA: G->C at 104. TCA = Ser. TGA(Trp)->TCA(Ser) = missense.
    let result2 = annotate_snv("chrM", 104, b'G', b'C', &tx, &idx, &fasta).unwrap();
    assert!(result2.consequences.contains(&Consequence::MissenseVariant));
    assert_eq!(result2.amino_acids.as_deref(), Some("W/S"));
}

#[test]
fn ref_mismatch_returns_error() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let err = annotate_snv("chr1", 1500, b'G', b'T', &tx, &idx, &fasta).unwrap_err();
    assert!(
        matches!(err, crate::VarEffectError::RefMismatch { .. }),
        "expected RefMismatch, got: {err:?}",
    );
}

// -- Integration tests (real data, #[ignore]-gated) -------------------------

#[test]
#[ignore]
fn tp53_r248w() {
    let store_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("data/vareffect/transcript_models_grch38.bin");
    let store = crate::TranscriptStore::load_from_path(&store_path).unwrap_or_else(|e| {
        panic!(
            "failed to load GRCh38 transcript store from {}: {e}. \
             Run `vareffect setup --assembly grch38` first.",
            store_path.display(),
        )
    });

    let fasta_path = std::env::var("GRCH38_FASTA").expect("GRCH38_FASTA env var");
    let fasta =
        FastaReader::open_with_assembly(std::path::Path::new(&fasta_path), crate::Assembly::GRCh38)
            .unwrap();

    let (tx, idx) = store
        .get_by_accession("NM_000546.6")
        .expect("NM_000546.6 not found");
    let result = annotate_snv("chr17", 7_674_220, b'G', b'A', tx, idx, &fasta).unwrap();

    assert!(
        result.consequences.contains(&Consequence::MissenseVariant)
            || result.consequences.contains(&Consequence::StopGained),
        "unexpected consequences: {:?}",
        result.consequences,
    );
    assert_eq!(result.protein_start, Some(248));
}

// -- trim_alleles tests -----------------------------------------------------

#[test]
fn trim_deletion() {
    let (r, a, adj) = trim_alleles(b"CA", b"C");
    assert_eq!(r, b"A");
    assert!(a.is_empty());
    assert_eq!(adj, 1);
}

#[test]
fn trim_insertion() {
    let (r, a, adj) = trim_alleles(b"C", b"CA");
    assert!(r.is_empty());
    assert_eq!(a, b"A");
    assert_eq!(adj, 1);
}

#[test]
fn trim_complex() {
    let (r, a, adj) = trim_alleles(b"CATG", b"CG");
    assert_eq!(r, b"AT");
    assert!(a.is_empty());
    assert_eq!(adj, 1);
}

#[test]
fn trim_with_suffix() {
    let (r, a, adj) = trim_alleles(b"ACGT", b"AT");
    assert_eq!(r, b"CG");
    assert!(a.is_empty());
    assert_eq!(adj, 1);
}

#[test]
fn trim_snv() {
    let (r, a, adj) = trim_alleles(b"A", b"T");
    assert_eq!(r, b"A");
    assert_eq!(a, b"T");
    assert_eq!(adj, 0);
}

#[test]
fn trim_no_shared() {
    let (r, a, adj) = trim_alleles(b"AC", b"GT");
    assert_eq!(r, b"AC");
    assert_eq!(a, b"GT");
    assert_eq!(adj, 0);
}

// -- CDS deletion tests -----------------------------------------------------

#[test]
fn frameshift_1bp_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1503, 1504, b"C", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert_eq!(result.impact, Impact::High);
    assert_eq!(result.protein_start, Some(2));
}

#[test]
fn frameshift_2bp_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1503, 1505, b"CG", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn inframe_3bp_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1503, 1506, b"CGT", &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::InframeDeletion));
    assert_eq!(result.impact, Impact::Moderate);
    assert!(result.amino_acids.is_some());
}

#[test]
fn inframe_6bp_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1503, 1509, b"CGTCGT", &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::InframeDeletion));
    assert_eq!(result.impact, Impact::Moderate);
}

// -- CDS insertion tests ----------------------------------------------------

#[test]
fn frameshift_1bp_insertion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 1503, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn frameshift_2bp_insertion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 1503, b"AT", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn inframe_3bp_insertion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 1503, b"GGG", &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::InframeInsertion));
    assert_eq!(result.impact, Impact::Moderate);
    assert!(result.amino_acids.is_some());
}

// -- Non-CDS indel tests ----------------------------------------------------

#[test]
fn intronic_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 2500, 2503, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::IntronVariant));
    assert_eq!(result.impact, Impact::Modifier);
}

#[test]
fn utr5_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1200, 1203, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FivePrimeUtrVariant)
    );
}

#[test]
fn utr3_insertion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 4600, b"GGG", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant)
    );
}

// -- Splice overlap tests ---------------------------------------------------

#[test]
fn deletion_overlaps_splice_donor() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [1998, 2002) -- spans exon/intron boundary, overlaps donor
    let result = annotate_deletion("chr1", 1998, 2002, b"AAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
            || result
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant),
        "expected splice consequence, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn deletion_overlaps_splice_acceptor() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [2998, 3002) -- spans intron/exon boundary, overlaps acceptor
    let result = annotate_deletion("chr1", 2998, 3002, b"AAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
            || result
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant),
        "expected splice consequence, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn insertion_at_splice_donor() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 2001, b"GGG", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
            || result
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant)
            || result
                .consequences
                .contains(&Consequence::SpliceRegionVariant),
        "expected splice consequence, got {:?}",
        result.consequences,
    );
}

// -- Minus-strand indel tests -----------------------------------------------

#[test]
fn minus_strand_frameshift() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    // 1bp deletion at codon 1 destroys start -> StartLost + FrameshiftVariant.
    // Regression: VEP emits both terms, vareffect previously
    // suppressed FrameshiftVariant when StartLost fired.
    let result = annotate_deletion("chr17", 19498, 19499, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StartLost),
        "expected start_lost (1bp del at codon 1 destroys Met), got {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "expected frameshift_variant alongside start_lost (1bp del is frameshifting), got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);

    // Frameshift well inside CDS (codon 4)
    let result2 = annotate_deletion("chr17", 19490, 19491, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        result2
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "expected frameshift_variant for deletion at codon 4, got {:?}",
        result2.consequences,
    );
    assert_eq!(result2.impact, Impact::High);
}

#[test]
fn minus_strand_inframe_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr17", 19490, 19493, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::InframeDeletion),
        "expected inframe_deletion, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::Moderate);
}

// -- Edge case tests --------------------------------------------------------

#[test]
fn deletion_of_start_codon() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1500, 1503, b"ATG", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StartLost),
        "expected start_lost, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn plus_strand_frameshift_start_lost() {
    // Regression (plus-strand mirror of `minus_strand_frameshift`):
    // a 1-bp deletion at codon 1 destroys the ATG AND shifts the reading
    // frame. VEP emits both `start_lost` and `frameshift_variant`; vareffect
    // must match.
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // CDS starts at genomic 1500 (first base of ATG). Delete the A.
    let result = annotate_deletion("chr1", 1500, 1501, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StartLost),
        "expected start_lost, got {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "expected frameshift_variant alongside start_lost, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn deletion_of_stop_codon() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 4497, 4500, b"TAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StopLost),
        "expected stop_lost, got {:?}",
        result.consequences,
    );
}

#[test]
fn plus_strand_deletion_spans_stop_into_3utr_frameshift_length() {
    // Regression (plus strand). Corresponds to the ClinVar
    // examples `1:1336134 GCCACGAGTCACATGATGT>G` and `1:2306758
    // AGCCGTAGATTCCGTGCCT>A`: an 18-bp genomic deletion whose CDS-projected
    // length is NOT a multiple of 3 (because part of the deletion lies in
    // the 3'UTR), yet which removes the stop codon. VEP emits
    // `{stop_lost, 3_prime_UTR_variant}` without `frameshift_variant`.
    //
    // Fixture `plus_strand_coding`: stop codon at genomic [4497, 4500),
    // 3'UTR at [4500, 5000), all within exon 2 [4000, 5000). Delete
    // [4498, 4508) -- 10 bases: 2 CDS bases (part of the stop) + 8 UTR
    // bases. cds_del_len = 2, not multiple of 3 -> previously hit the
    // frameshift branch and emitted only `FrameshiftVariant`.
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 4498, 4508, b"AAAAAAAAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StopLost),
        "expected stop_lost, got {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "expected 3_prime_UTR_variant, got {:?}",
        result.consequences,
    );
    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "must NOT emit frameshift_variant (VEP suppresses it via undefined cds_end), got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn plus_strand_deletion_spans_stop_into_3utr_inframe_length() {
    // Sanity check: a deletion whose CDS-projected length IS a multiple
    // of 3 but still spans the stop into the 3'UTR should also hit the
    // new stop-lost early-return and emit
    // `{stop_lost, 3_prime_UTR_variant}`.
    //
    // Delete [4497, 4506) -- 9 bases: 3 CDS bases (the stop codon) + 6
    // UTR bases. cds_del_len = 3.
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 4497, 4506, b"AAAAAAAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StopLost),
        "expected stop_lost, got {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "expected 3_prime_UTR_variant, got {:?}",
        result.consequences,
    );
}

#[test]
fn minus_strand_deletion_spans_stop_into_3utr_frameshift_length() {
    // Regression (minus-strand mirror). Fixture
    // `minus_strand_coding`: stop codon at genomic [11000, 11003) (last
    // 3 CDS bases on minus strand = genomically lowest), 3'UTR at
    // [10000, 11000), all within exon 2 [10000, 12000). Delete
    // [10995, 11002) -- 7 bases: 2 CDS bases + 5 UTR bases. cds_del_len
    // = 2, not multiple of 3.
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr17", 10_995, 11_002, b"AAAAAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result.consequences.contains(&Consequence::StopLost),
        "expected stop_lost, got {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant),
        "expected 3_prime_UTR_variant, got {:?}",
        result.consequences,
    );
    assert!(
        !result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "must NOT emit frameshift_variant, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

#[test]
fn non_coding_transcript_deletion() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = noncoding_2_exon();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr2", 150, 153, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant)
    );
}

#[test]
fn exon_intron_boundary_deferred() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1999, 2002, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "expected splice_donor_variant, got {:?}",
        result.consequences,
    );
}

// -- Boundary-spanning, complex delins, MNV tests ---------------------------

#[test]
fn boundary_del_splice_donor_priority() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [1998, 2003): 2 CDS + 3 intronic. Donor +1/+2 at 2000-2001.
    let result = annotate_deletion("chr1", 1998, 2003, b"AAAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "splice donor should win over frame: {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::High);
}

/// CDS/UTR boundary deletion (within same exon, no exon-intron crossing).
#[test]
fn boundary_del_cds_utr() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [4498, 4502): 2 CDS bases (last 2 of stop TAA) + 2 UTR bases.
    let result = annotate_deletion("chr1", 4498, 4502, b"AAAA", &tx, &idx, &fasta).unwrap();
    let has_frame_or_stop = result
        .consequences
        .contains(&Consequence::FrameshiftVariant)
        || result.consequences.contains(&Consequence::StopLost);
    assert!(
        has_frame_or_stop,
        "expected frameshift or stop_lost for CDS/UTR boundary del: {:?}",
        result.consequences,
    );
}

#[test]
fn insertion_at_exon_intron_junction() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 2000, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "insertion at donor junction: {:?}",
        result.consequences,
    );
}

#[test]
fn multi_exon_del_full_exon() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [2500, 3600): spans intron + exon 1 + intron.
    let result = annotate_deletion("chr1", 2500, 3600, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
            || result
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant),
        "multi-exon del should have splice: {:?}",
        result.consequences,
    );
}

#[test]
fn delins_cds_frameshift() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete 4bp, insert 2bp at CDS offset 3. Net -2 -> frameshift.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_complex_delins(&ctx, 1503, b"CGTC", b"AT").unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "delins -2 should be frameshift: {:?}",
        result.consequences,
    );
}

#[test]
fn delins_cds_inframe() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete 3bp, insert 6bp. Net +3 -> inframe delins.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_complex_delins(&ctx, 1503, b"CGT", b"AAAAAA").unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "inframe delins should be protein_altering: {:?}",
        result.consequences,
    );
}

#[test]
fn delins_splice_overlap() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // Delete [1999, 2002), replace with 1bp. Overlaps donor +1/+2.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_complex_delins(&ctx, 1999, b"AAA", b"T").unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "delins at donor: {:?}",
        result.consequences,
    );
}

#[test]
fn mnv_single_codon_missense() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // 2-base MNV at CDS offset 3-4 (codon 2). REF=CG, ALT=TA.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_mnv(&ctx, 1503, b"CG", b"TA").unwrap();
    assert!(
        result.consequences.contains(&Consequence::MissenseVariant)
            || result.consequences.contains(&Consequence::StopGained),
        "MNV in single codon: {:?}",
        result.consequences,
    );
}

#[test]
fn mnv_two_codons() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // CDS offset 5-6 spans codon boundary between codon 1 (3-5) and codon 2 (6-8).
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_mnv(&ctx, 1505, b"TC", b"AA").unwrap();
    let has_coding = result.consequences.contains(&Consequence::MissenseVariant)
        || result.consequences.contains(&Consequence::StopGained)
        || result
            .consequences
            .contains(&Consequence::SynonymousVariant);
    assert!(
        has_coding,
        "MNV spanning codons should have coding consequence: {:?}",
        result.consequences,
    );
}

#[test]
fn mnv_creates_stop() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // CGT at 1503-1505 -> TAA (stop). 3-base MNV.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_mnv(&ctx, 1503, b"CGT", b"TAA").unwrap();
    assert!(
        result.consequences.contains(&Consequence::StopGained),
        "MNV creating stop: {:?}",
        result.consequences,
    );
}

#[test]
fn mnv_synonymous() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // CGT -> AGA: all 3 bases change, both Arg. Synonymous 3-base MNV.
    let ctx = AnnotateCtx {
        chrom: "chr1",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_mnv(&ctx, 1503, b"CGT", b"AGA").unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SynonymousVariant),
        "CGT->AGA (Arg->Arg) should be synonymous: {:?}",
        result.consequences,
    );
}

#[test]
fn minus_strand_boundary_del() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    // Delete [17998, 18002): 2 intronic donor bases + 2 exonic CDS bases.
    let result = annotate_deletion("chr17", 17998, 18002, b"AAAA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "minus strand boundary del at donor: {:?}",
        result.consequences,
    );
}

#[test]
fn minus_strand_mnv() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = minus_strand_coding();
    let idx = build_index(&tx);
    // 2-base MNV at 19498-19499 (CDS offset 0-1, codon 0)
    let ctx = AnnotateCtx {
        chrom: "chr17",
        transcript: &tx,
        index: &idx,
        fasta: &fasta,
    };
    let result = annotate_mnv(&ctx, 19498, b"AT", b"CC").unwrap();
    let has_coding = result.consequences.contains(&Consequence::MissenseVariant)
        || result.consequences.contains(&Consequence::StartLost)
        || result.consequences.contains(&Consequence::StopGained);
    assert!(
        has_coding,
        "minus strand MNV should have coding consequence: {:?}",
        result.consequences,
    );
}

#[test]
fn incomplete_terminal_codon() {
    let (_tmp, fasta) = write_test_fasta();
    let mut tx = single_exon_coding();
    // Make CDS length not divisible by 3 (802 bases, remainder 1)
    tx.cds_genomic_end = Some(1402);
    tx.cds_segments[0].genomic_end = 1402;
    let idx = build_index(&tx);
    assert_eq!(idx.total_cds_length() % 3, 1);
    // SNV at the last CDS base (in the incomplete terminal codon)
    let result = annotate_snv("chr3", 1401, b'A', b'T', &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::IncompleteTerminalCodonVariant),
        "variant in incomplete terminal codon: {:?}",
        result.consequences,
    );
}

#[test]
fn splice_donor_boundary_regression() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1999, 2002, b"AAA", &tx, &idx, &fasta).unwrap();
    assert!(
        !result
            .consequences
            .contains(&Consequence::ProteinAlteringVariant),
        "should not return old placeholder: {:?}",
        result.consequences,
    );
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant),
        "should be splice_donor_variant: {:?}",
        result.consequences,
    );
}

#[test]
fn trim_complex_dispatches_delins() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx]);
    let results = annotate("chr1", 1503, b"CGTC", b"TG", &store, &fasta).unwrap();
    assert!(
        !results.is_empty(),
        "should annotate against the transcript"
    );
}

#[test]
fn trim_mnv_dispatches() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx]);
    let results = annotate("chr1", 1503, b"CG", b"TA", &store, &fasta).unwrap();
    assert!(!results.is_empty(), "should annotate MNV");
}

// -- Intergenic variant tests -----------------------------------------------

#[test]
fn intergenic_no_overlap() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx]);
    // Position 50 is outside the transcript [1000, 5000).
    let results = annotate("chr1", 50, b"A", b"T", &store, &fasta).unwrap();
    assert_eq!(results.len(), 1);
    let r = &results[0];
    assert_eq!(r.consequences, vec![Consequence::IntergenicVariant]);
    assert_eq!(r.impact, Impact::Modifier);
    assert!(r.transcript.is_empty());
    assert!(r.gene_symbol.is_empty());
    assert!(r.protein_start.is_none());
    assert!(r.cds_position.is_none());
    assert!(!r.predicts_nmd);
}

#[test]
fn intergenic_different_chrom() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx]);
    // chr2 exists in the FASTA but has no transcripts in the store.
    let results = annotate("chr2", 150, b"A", b"T", &store, &fasta).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(
        results[0].consequences,
        vec![Consequence::IntergenicVariant]
    );
}

#[test]
fn non_intergenic_has_transcript() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(crate::Assembly::GRCh38, vec![tx]);
    // Position 1505 is inside the transcript CDS.
    let results = annotate("chr1", 1505, b"T", b"A", &store, &fasta).unwrap();
    assert!(!results.is_empty());
    assert!(
        !results[0]
            .consequences
            .contains(&Consequence::IntergenicVariant)
    );
}

// -- predicts_nmd integration tests -----------------------------------------

#[test]
fn stop_gained_single_exon_no_nmd() {
    let (_tmp, fasta) = write_stop_gained_fasta();
    let tx = stop_gained_transcript();
    let idx = build_index(&tx);
    // CGA->TGA: stop gained in a single-exon transcript.
    let result = annotate_snv("chr1", 103, b'C', b'T', &tx, &idx, &fasta).unwrap();
    assert!(result.consequences.contains(&Consequence::StopGained));
    assert!(!result.predicts_nmd);
}

#[test]
fn frameshift_multi_exon_nmd() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // 1bp deletion at early CDS position: CDS offset 3, 3 exons -> NMD.
    let result = annotate_deletion("chr1", 1503, 1504, b"C", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert!(result.predicts_nmd);
}

#[test]
fn frameshift_immediate_stop_still_gets_nmd() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    // 2bp deletion at early CDS: should be FrameshiftVariant (not StopGained)
    // even if the first changed AA is * (immediate nonsense). NMD should
    // still be predicted because the variant is far from the last junction.
    let result = annotate_deletion("chr1", 1503, 1505, b"CG", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant),
        "expected FrameshiftVariant, got: {:?}",
        result.consequences,
    );
    assert!(
        !result.consequences.contains(&Consequence::StopGained),
        "should not contain StopGained for a frameshift: {:?}",
        result.consequences,
    );
    assert!(result.predicts_nmd);
}
