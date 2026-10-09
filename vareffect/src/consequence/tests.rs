use super::complex::{annotate_complex_delins, annotate_mnv};
use super::helpers::trim_alleles;
use super::*;
use crate::fasta::{FastaReader, write_genome_binary};
use crate::locate::{LocateIndex, VariantLocation, locate_variant};
use crate::test_fixtures::{minus_strand_coding, noncoding_2_exon, plus_strand_coding};
use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};
use tempfile::TempDir;

// -- PTC fixture ------------------------------------------------------------

/// Transcript for the premature-termination-codon tests.
///
/// Same geometry as [`plus_strand_coding`] — 3 exons, CDS segments
/// `[1500,2000) + [3000,3500) + [4000,4500)`, last CDS junction at 1-based CDS
/// position 1001 — but on `chr4`, so its planted stop codons cannot perturb
/// any existing test.
fn ptc_transcript() -> TranscriptModel {
    let mut tx = plus_strand_coding();
    tx.accession = "NM_TEST_PTC.1".into();
    tx.protein_accession = Some("NP_TEST_PTC.1".into());
    tx.gene_symbol = "TESTPTC".into();
    tx.hgnc_id = Some("HGNC:99904".into());
    tx.chrom = "chr4".into();
    tx
}

/// Genome for [`ptc_transcript`], carrying stops we can actually reach.
///
/// The CDS body is `CGT` repeating, which reads `CGT` / `GTC` / `TCG` in the
/// three frames — no stop in any of them. Every stop below is therefore one
/// planted deliberately, at a CDS offset chosen so it is read in exactly one
/// shifted frame:
///
/// - offset % 3 == 1 → reached after a **1-base deletion**
/// - offset % 3 == 2 → reached after a **2-base deletion**
///
/// That separation is what lets one fixture serve both an NMD-predicted and an
/// NMD-escaping frameshift. `write_test_fasta`'s chr1 cannot do this at all:
/// its shifted frames read `GTC`/`TCG` forever and its 3'UTR is poly-A, so
/// every frameshift there ends in `fsTer?`.
fn write_ptc_fasta() -> (TempDir, FastaReader) {
    let tmp = TempDir::new().unwrap();
    let mut chr4 = vec![b'A'; 6000];

    let mut cds = Vec::with_capacity(1500);
    cds.extend_from_slice(b"ATG");
    while cds.len() < 1500 {
        cds.extend_from_slice(b"CGT");
    }
    cds.truncate(1500);

    // Reference terminator, preceded by a codon chosen so an in-frame
    // deletion can fuse across it into a *new* stop (`TAC|TAA`, delete
    // `AC`+`T` -> `TAA`). That is the one case where `alt_aas` holds a stop
    // while `ref_aas` still holds the reference one.
    cds[1494..1497].copy_from_slice(b"TAC");
    cds[1497..1500].copy_from_slice(b"TAA");
    // Premature stop reached only after a 2-base deletion (401 % 3 == 2).
    // Upstream of the junction, so the PTC is NMD-predicted.
    cds[401..404].copy_from_slice(b"TAA");
    // Premature stop reached only after a 1-base deletion (1150 % 3 == 1).
    // Inside the last coding exon, so the PTC escapes NMD even though the
    // variant itself sits far upstream — the divergence this field exists for.
    cds[1150..1153].copy_from_slice(b"TAA");

    for (offset, base) in cds.iter().enumerate() {
        let genomic = match offset {
            0..500 => 1500 + offset,
            500..1000 => 3000 + (offset - 500),
            _ => 4000 + (offset - 1000),
        };
        chr4[genomic] = *base;
    }

    // 3'UTR [4500,5000) is poly-A — no stop in any frame — except one planted
    // TAA that only the 2-base-deletion frame reads as a codon. The 1-base
    // frame reads it as `ATA`, so it still runs off the end as `fsTer?`.
    chr4[4502..4505].copy_from_slice(b"TAA");

    let contigs: Vec<(&str, &[u8])> = vec![("chr4", &chr4)];
    let bin_path = tmp.path().join("ptc.bin");
    let idx_path = tmp.path().join("ptc.bin.idx");
    write_genome_binary(&contigs, "ptc", &bin_path, &idx_path).unwrap();
    let reader = FastaReader::open(&bin_path).unwrap();
    (tmp, reader)
}

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
    let reader = FastaReader::open(&bin_path).unwrap();
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
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        selenocysteine_codons: Vec::new(),
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
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        selenocysteine_codons: Vec::new(),
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
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 1,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        selenocysteine_codons: Vec::new(),
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
    let reader = FastaReader::open(&bin_path).unwrap();
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
        .expect("workspace root")
        .join("data/vareffect/transcript_models.bin");
    let store = crate::TranscriptStore::load_from_path(&store_path).unwrap_or_else(|e| {
        panic!(
            "failed to load transcript store from {}: {e}",
            store_path.display()
        )
    });

    let fasta_path = std::env::var("FASTA_PATH").expect("FASTA_PATH env var");
    let fasta = FastaReader::open(std::path::Path::new(&fasta_path)).unwrap();

    let (tx, idx) = store
        .get_by_accession("NM_000546.6")
        .expect("NM_000546.6 not found");
    let result = annotate_snv("chr17", 7_674_220, b'G', b'A', tx, idx, &fasta).unwrap();

    // VEP: NM_000546.6:c.742C>T, Cgg/Tgg, missense only (vep_concordance_snv.rs #1).
    assert_eq!(result.consequences, vec![Consequence::MissenseVariant]);
    assert_eq!(result.hgvs_p.as_deref(), Some("p.Arg248Trp"));
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

#[test]
fn inframe_insertion_with_an_unknown_base_is_coding_sequence_variant() {
    // VEP 116 gives a non-ACGT allele no peptide: `coding_sequence_variant`
    // alone, no amino acids and no HGVS p. (ClinVar NM_000057.4 insertion
    // with `NNNNNNNNNN`, VEP REST).
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 1503, b"GNG", &tx, &idx, &fasta).unwrap();
    assert_eq!(
        result.consequences,
        vec![Consequence::CodingSequenceVariant]
    );
    assert_eq!(result.amino_acids, None);
    assert_eq!(result.hgvs_p, None);
    assert_eq!(result.codons.as_deref(), Some("-/GNG"));
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
    // Between +1 and +2: the canonical term alone, with the intron number
    // (VEP 116 BRCA2 NM_000059.4 c.67+1_67+2insA: splice_donor_variant).
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 2001, b"GGG", &tx, &idx, &fasta).unwrap();
    assert_eq!(result.consequences, vec![Consequence::SpliceDonorVariant]);
    assert_eq!(result.intron.as_deref(), Some("1/2"));
    assert_eq!(result.exon, None);
    let minus = minus_strand_coding();
    let minus_idx = build_index(&minus);
    let result = annotate_insertion("chr17", 17_999, b"C", &minus, &minus_idx, &fasta).unwrap();
    assert_eq!(result.consequences, vec![Consequence::SpliceDonorVariant]);
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
    // CDS starts at genomic 1500 (first base of ATG). Delete the T: `AGC`
    // is no start codon. (Deleting the A would only shift the start into
    // this fixture's all-`A` 5'UTR, which VEP's `_ins_del_start_altered`
    // calls retained; see the next assertion.)
    let retained = annotate_deletion("chr1", 1500, 1501, b"A", &tx, &idx, &fasta).unwrap();
    assert!(
        retained
            .consequences
            .contains(&Consequence::StartRetainedVariant)
            && retained.consequences.contains(&Consequence::StartLost)
    );
    let result = annotate_deletion("chr1", 1501, 1502, b"T", &tx, &idx, &fasta).unwrap();
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
    // With codon 1 destroyed there is no defined reading frame, so no
    // termination codon may be claimed — the formatter is skipped entirely.
    assert_eq!(result.ptc, PtcStatus::Indeterminate);
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
    // Between the last exon base and +1: an mRNA insertion before exon 2's
    // first coding base, plus splice_region_variant, with no exon or intron
    // number (VEP 116 BRCA2 NM_000059.4 c.67_67+1insA: frameshift_variant,
    // splice_region_variant, CDS 67-68).
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_insertion("chr1", 2000, b"A", &tx, &idx, &fasta).unwrap();
    assert_eq!(
        result.consequences,
        vec![
            Consequence::FrameshiftVariant,
            Consequence::SpliceRegionVariant
        ]
    );
    assert_eq!(
        (result.cds_position, result.cds_position_end),
        (Some(500), Some(501))
    );
    assert_eq!(
        (result.exon.as_deref(), result.intron.as_deref()),
        (None, None)
    );
    // The acceptor side of the same junction is the same mRNA insertion.
    let acceptor = annotate_insertion("chr1", 3000, b"A", &tx, &idx, &fasta).unwrap();
    assert_eq!(acceptor.cds_position, Some(500));
    assert!(
        acceptor
            .consequences
            .contains(&Consequence::SpliceRegionVariant)
    );
}

#[test]
fn insertion_at_a_non_coding_junction_is_a_transcript_variant() {
    // VEP 116 ENST00000424518.6 n.197_197+1insT: splice_region_variant,
    // non_coding_transcript_variant, cDNA at the mRNA point, no exon.
    let (_tmp, fasta) = write_test_fasta();
    let tx = noncoding_2_exon();
    let idx = build_index(&tx);
    let result = annotate_insertion(&tx.chrom.clone(), 300, b"A", &tx, &idx, &fasta).unwrap();
    assert_eq!(
        result.consequences,
        vec![
            Consequence::SpliceRegionVariant,
            Consequence::NonCodingTranscriptVariant
        ]
    );
    assert_eq!(
        (result.cdna_position, result.cdna_position_end),
        (Some(200), Some(201))
    );
    assert_eq!(result.exon, None);
}

#[test]
fn deletion_at_a_donor_site_keeps_its_terms() {
    // Deletions keep `build_splice_indel_result`: a deletion over the exon
    // end and +1 is splice_donor_variant with its coding and intron terms.
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);
    let result = annotate_deletion("chr1", 1_999, 2_001, b"AA", &tx, &idx, &fasta).unwrap();
    assert!(
        result
            .consequences
            .contains(&Consequence::SpliceDonorVariant)
    );
    assert!(result.consequences.contains(&Consequence::IntronVariant));
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
    let result = annotate_complex_delins("chr1", 1503, b"CGTC", b"AT", &tx, &idx, &fasta).unwrap();
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
    let result =
        annotate_complex_delins("chr1", 1503, b"CGT", b"AAAAAA", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_complex_delins("chr1", 1999, b"AAA", b"T", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_mnv("chr1", 1503, b"CG", b"TA", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_mnv("chr1", 1505, b"TC", b"AA", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_mnv("chr1", 1503, b"CGT", b"TAA", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_mnv("chr1", 1503, b"CGT", b"AGA", &tx, &idx, &fasta).unwrap();
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
    let result = annotate_mnv("chr17", 19498, b"AT", b"CC", &tx, &idx, &fasta).unwrap();
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
    let store = crate::TranscriptStore::from_transcripts(vec![tx]);
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
    let store = crate::TranscriptStore::from_transcripts(vec![tx]);
    let results = annotate("chr1", 1503, b"CG", b"TA", &store, &fasta).unwrap();
    assert!(!results.is_empty(), "should annotate MNV");
}

// -- Intergenic variant tests -----------------------------------------------

#[test]
fn intergenic_no_overlap() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let store = crate::TranscriptStore::from_transcripts(vec![tx]);
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
    let store = crate::TranscriptStore::from_transcripts(vec![tx]);
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
    let store = crate::TranscriptStore::from_transcripts(vec![tx]);
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

// -- PTC through the real annotate pipeline ---------------------------------

/// The divergence the `ptc` field exists for: the variant sits far upstream of
/// the last junction (so `predicts_nmd` reports NMD), while the termination
/// codon it actually creates lands inside the last coding exon and escapes.
#[test]
fn ptc_frameshift_escapes_nmd_while_variant_site_predicts_it() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    // 1-base deletion at CDS offset 100 (genomic 1600) -> reads the %3==1
    // frame -> first stop at CDS 1150, i.e. alt residue 384.
    let result = annotate_deletion("chr4", 1600, 1601, b"G", &tx, &idx, &fasta).unwrap();

    assert!(
        result
            .consequences
            .contains(&Consequence::FrameshiftVariant)
    );
    assert!(
        result.predicts_nmd,
        "the variant site is 900 nt upstream of the junction, so the \
         VEP-parity field must still report NMD"
    );
    assert_eq!(
        result.ptc,
        PtcStatus::At {
            protein_position: 384,
            nmd_at_ptc: false,
        },
        "the termination codon is in the last coding exon and escapes"
    );
}

/// The agreeing case: both fields report NMD when the termination codon is
/// also well upstream of the junction.
#[test]
fn ptc_frameshift_predicts_nmd_when_codon_is_upstream() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    // 2-base deletion at CDS offset 100 -> reads the %3==2 frame -> first stop
    // at CDS 401, i.e. alt residue 134.
    let result = annotate_deletion("chr4", 1600, 1602, b"GT", &tx, &idx, &fasta).unwrap();

    assert!(result.predicts_nmd);
    assert_eq!(
        result.ptc,
        PtcStatus::At {
            protein_position: 134,
            nmd_at_ptc: true,
        }
    );
}

/// A frameshift whose alternate frame reaches a stop only in the 3'UTR still
/// yields a located codon — the continuation scan is part of the PTC search.
#[test]
fn ptc_frameshift_finds_a_stop_in_the_3prime_utr() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    // 2-base deletion at CDS offset 1200 (genomic 4200): the %3==2 frame has
    // no remaining CDS stop, so the scan continues into the UTR.
    let result = annotate_deletion("chr4", 4200, 4202, b"GT", &tx, &idx, &fasta).unwrap();

    let position = result
        .ptc
        .protein_position()
        .expect("a UTR stop is still a located termination codon");
    assert!(
        position > 499,
        "a UTR-derived codon lies past the 499-residue reference protein, got {position}"
    );
    assert_eq!(result.ptc.nmd_at_ptc(), Some(false), "past the junction");
}

/// HGVS `fsTer?`: the alternate frame runs off the end of the 3'UTR with no
/// stop at all. Distinct from `Indeterminate`, and never a located position.
#[test]
fn ptc_frameshift_with_no_reachable_stop_is_no_stop_codon() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    // 1-base deletion at CDS offset 1200: the %3==1 frame reads the planted
    // UTR stop as `ATA`, so nothing terminates it.
    let result = annotate_deletion("chr4", 4200, 4201, b"G", &tx, &idx, &fasta).unwrap();

    assert!(
        result
            .hgvs_p
            .as_deref()
            .is_some_and(|p| p.ends_with("fsTer?")),
        "expected fsTer? notation, got {:?}",
        result.hgvs_p
    );
    assert_eq!(result.ptc, PtcStatus::NoStopCodon);
    assert_eq!(result.ptc.protein_position(), None);
}

/// Every located codon must agree with the `fsTer` count the same call emits.
#[test]
fn ptc_agrees_with_the_emitted_hgvs_notation() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    for (start, end, deleted) in [
        (1600u64, 1601u64, &b"G"[..]),
        (1600, 1602, &b"GT"[..]),
        (4200, 4202, &b"GT"[..]),
    ] {
        let result = annotate_deletion("chr4", start, end, deleted, &tx, &idx, &fasta).unwrap();
        let hgvs_p = result.hgvs_p.clone().expect("frameshift notation");
        let (first, count) = parse_fs_ter(&hgvs_p).expect("fsTer notation");
        assert_eq!(
            result.ptc.protein_position(),
            Some(first + count - 1),
            "ptc must equal N + M - 1 from {hgvs_p}"
        );
    }
}

/// Extract `(N, M)` from `p.{Aaa}{N}{Bbb}fsTer{M}`.
fn parse_fs_ter(hgvs_p: &str) -> Option<(u32, u32)> {
    let (head, tail) = hgvs_p.split_once("fsTer")?;
    let count: u32 = tail.parse().ok()?;
    let digits: String = head
        .chars()
        .skip_while(|c| !c.is_ascii_digit())
        .take_while(|c| c.is_ascii_digit())
        .collect();
    Some((digits.parse().ok()?, count))
}

/// An in-frame deletion fusing across the reference terminator leaves a stop
/// in the alternate window, but it is the normal terminator re-formed one
/// codon earlier, not a premature one. `annotate_cds_inframe_deletion` must
/// not call that `stop_gained`: VEP's own predicate is
/// `($alt_pep =~ /\*/) and ($ref_pep !~ /\*/)`, and a HIGH-impact nonsense
/// term here would be a truncation claim about a protein that terminates
/// where it always did.
#[test]
fn inframe_deletion_across_the_terminator_is_not_stop_gained() {
    let (_tmp, fasta) = write_ptc_fasta();
    let tx = ptc_transcript();
    let idx = build_index(&tx);

    // Delete CDS offsets 1495-1497 (genomic 4495-4497): `TAC|TAA` -> `TAA`.
    let result = annotate_deletion("chr4", 4495, 4498, b"ACT", &tx, &idx, &fasta).unwrap();

    assert!(
        !result.consequences.contains(&Consequence::StopGained),
        "the reference window still holds the normal terminator, so no \
         premature stop was gained; got {:?}",
        result.consequences,
    );
    assert!(
        result.consequences.contains(&Consequence::InframeDeletion),
        "expected the inframe-deletion term to stand alone, got {:?}",
        result.consequences,
    );
    assert_eq!(result.impact, Impact::Moderate);
    assert!(!result.predicts_nmd);
    assert_eq!(
        result.ptc,
        PtcStatus::NotApplicable,
        "nothing truncating happened, so there is no termination codon to \
         report",
    );
}

/// A deletion spanning the 5'UTR/CDS boundary reaches the start-codon branch
/// whatever its length, and must carry the `5_prime_UTR_variant` term without
/// `frameshift_variant`.
///
/// VEP reports exactly `{5_prime_UTR_variant, start_lost}` for this geometry
/// even when the CDS bases removed are not a multiple of three, because its
/// `frameshift` predicate short-circuits on the undefined `cds_start` — the
/// same short-circuit this file already replicates at the 3' end. Recorded
/// ground truth: `ALPL` NM_000478.6 chr1:21554077 `GCACCATGATTT>G`, 7 CDS
/// bases deleted.
///
/// Both lengths are asserted together because the branch boundary used to sit
/// between them, sending a 1-2 base deletion down the frameshift path.
#[test]
fn utr_spanning_deletion_into_the_start_codon_never_frameshifts() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let idx = build_index(&tx);

    // CDS begins at genomic 1500, `ATG` at 1500-1502, body `CGT` from 1503.
    let one_cds_base = annotate_deletion("chr1", 1498, 1501, b"AAA", &tx, &idx, &fasta).unwrap();
    let four_cds_bases =
        annotate_deletion("chr1", 1497, 1504, b"AAAATGC", &tx, &idx, &fasta).unwrap();

    for (label, result) in [
        ("1 CDS base", one_cds_base),
        ("4 CDS bases", four_cds_bases),
    ] {
        assert!(
            result.consequences.contains(&Consequence::StartLost),
            "{label}: the deletion reaches the initiator codon, got {:?}",
            result.consequences,
        );
        assert!(
            result
                .consequences
                .contains(&Consequence::FivePrimeUtrVariant),
            "{label}: the footprint extends into the 5'UTR, got {:?}",
            result.consequences,
        );
        assert!(
            !result
                .consequences
                .contains(&Consequence::FrameshiftVariant),
            "{label}: VEP withholds the frameshift term on this geometry, got {:?}",
            result.consequences,
        );
        assert_eq!(result.impact, Impact::High, "{label}");
        assert!(
            !result.predicts_nmd,
            "{label}: with no reading frame there is no termination codon to place",
        );
        assert_eq!(result.ptc, PtcStatus::NotApplicable, "{label}");
    }
}

// -- Incomplete CDS framing and frameshift introns ---------------------------
//
// Expected values follow VEP: the padded numbering and partial-codon output
// come from the Ensembl VEP REST release 116 rows pinned in
// `tests/incomplete_cds_concordance.rs`; the start-codon rule from
// `VariationEffect.pm::_overlaps_start_codon` (ensembl-variation release/115),
// which returns 0 for a `cds_start_NF` transcript.
mod incomplete_cds {
    use super::*;
    use crate::locate::helpers::{exonic_splice_spill, frameshift_intron_endpoint, intron_gap};

    /// [`plus_strand_coding`] with GFF3 phase `phase` on the first segment,
    /// marked 5'-incomplete when the phase pads codon 1. On
    /// `write_test_fasta`'s chr1 the first genomic coding bases are `ATG`.
    fn padded_plus(phase: u8) -> (TranscriptModel, LocateIndex) {
        let mut tx = plus_strand_coding();
        tx.cds_segments[0].phase = phase;
        tx.cds_start_complete = phase == 0;
        let index = LocateIndex::build(&tx).unwrap();
        (tx, index)
    }

    #[test]
    fn first_segment_phase_sets_the_pad() {
        for (phase, pad) in [(0u8, 0u32), (1, 2), (2, 1)] {
            let (_, index) = padded_plus(phase);
            assert_eq!(index.cds_start_pad(), pad, "phase {phase}");
            assert_eq!(index.total_cds_length(), 1500 + pad, "phase {phase}");
            assert_eq!(index.coding_length(), 1500, "phase {phase}");
        }
    }

    #[test]
    fn padded_cds_offsets_follow_vep_numbering() {
        // Phase 1 -> pad 2: the first genomic coding base is c.3 in codon 1,
        // the next one c.4 opens codon 2 (VEP: ENST00000379099.3 c.3 / c.4).
        let (tx, index) = padded_plus(1);
        let first = locate_variant("chr1", 1500, 1501, &tx, &index).unwrap();
        let second = locate_variant("chr1", 1501, 1502, &tx, &index).unwrap();
        assert!(matches!(
            first,
            VariantLocation::CdsExon {
                cds_offset: 2,
                codon_number: 1,
                codon_position: 2,
                ..
            }
        ));
        assert!(matches!(
            second,
            VariantLocation::CdsExon {
                cds_offset: 3,
                codon_number: 2,
                codon_position: 0,
                ..
            }
        ));
    }

    #[test]
    fn fetch_cds_sequence_writes_the_pad_as_n() {
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(1);
        let seq =
            super::super::helpers::fetch_cds_sequence(0, 5, "chr1", &tx, &index, &fasta).unwrap();
        assert_eq!(seq, b"NNATG");
        let pad_only =
            super::super::helpers::fetch_cds_sequence(0, 2, "chr1", &tx, &index, &fasta).unwrap();
        assert_eq!(pad_only, b"NN");
    }

    #[test]
    fn snv_in_the_padded_codon_is_coding_sequence_variant_with_x() {
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(1);
        // chr1:1500 is `A`, the third base of the padded codon `NNA`.
        let result = annotate_snv("chr1", 1500, b'A', b'G', &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(result.codons.as_deref(), Some("nnA/nnG"));
        assert_eq!(result.amino_acids.as_deref(), Some("X"));
        assert_eq!(result.protein_start, Some(1));
        assert_eq!(result.cds_position, Some(3));
        // cDNA is unpadded: 500 bases of 5'UTR, then the first coding base.
        assert_eq!(result.cdna_position, Some(501));
        assert_eq!(result.hgvs_p.as_deref(), Some("p.?"));
        assert_eq!(result.hgvs_c.as_deref(), Some("NM_TEST_PLUS.1:c.3A>G"));
    }

    #[test]
    fn indels_in_a_padded_codon_never_call_start_lost() {
        let (_tmp, fasta) = write_test_fasta();
        // Phase 2 -> pad 1: codon 1 is `N` + chr1:1500-1501 (`AT`).
        let (tx, index) = padded_plus(2);

        let deletion = annotate_deletion("chr1", 1500, 1501, b"A", &tx, &index, &fasta).unwrap();
        assert!(
            deletion
                .consequences
                .contains(&Consequence::FrameshiftVariant)
        );
        assert!(!deletion.consequences.contains(&Consequence::StartLost));
        assert_ne!(deletion.hgvs_p.as_deref(), Some("p.Met1?"));

        let insertion = annotate_insertion("chr1", 1501, b"C", &tx, &index, &fasta).unwrap();
        assert!(
            insertion
                .consequences
                .contains(&Consequence::FrameshiftVariant)
        );
        assert!(!insertion.consequences.contains(&Consequence::StartLost));

        let inframe = annotate_deletion("chr1", 1500, 1503, b"ATG", &tx, &index, &fasta).unwrap();
        assert!(inframe.consequences.contains(&Consequence::InframeDeletion));
        assert!(!inframe.consequences.contains(&Consequence::StartLost));
    }

    #[test]
    fn complete_start_codon_keeps_start_lost() {
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(0);
        let result = annotate_snv("chr1", 1501, b'T', b'C', &tx, &index, &fasta).unwrap();
        assert_eq!(result.consequences, vec![Consequence::StartLost]);
    }

    #[test]
    fn ends_in_stop_codon_reads_the_last_codon() {
        let (_tmp, fasta) = write_test_fasta();
        let tx = plus_strand_coding();
        let index = LocateIndex::build(&tx).unwrap();
        assert!(super::super::helpers::ends_in_stop_codon("chr1", &tx, &index, &fasta).unwrap());

        // One base short: the CDS no longer ends on a codon boundary, as for a
        // stop completed only by polyadenylation.
        let mut short = plus_strand_coding();
        short.cds_segments[2].genomic_end = 4_499;
        short.cds_genomic_end = Some(4_499);
        let short_index = LocateIndex::build(&short).unwrap();
        assert!(
            !super::super::helpers::ends_in_stop_codon("chr1", &short, &short_index, &fasta)
                .unwrap()
        );

        // `cds_end_NF`: no stop codon, whatever the last codon reads.
        let mut end_nf = plus_strand_coding();
        end_nf.cds_end_complete = false;
        let end_nf_index = LocateIndex::build(&end_nf).unwrap();
        assert!(
            !super::super::helpers::ends_in_stop_codon("chr1", &end_nf, &end_nf_index, &fasta)
                .unwrap()
        );
    }

    /// VEP 116 REST ENST00000616125.5 (`cds_start_NF`, phase 0, codon 1
    /// `ATG`): `c.1A>G` is `missense_variant` `p.Met1Val` and `c.1_3del`
    /// `inframe_deletion` `p.Met1del`, never `start_lost`
    /// (`tests/data/tools/vep_rest/vep_rest_wp2.tsv`).
    #[test]
    fn whole_atg_on_a_start_incomplete_cds_is_not_a_start_codon() {
        let (_tmp, fasta) = write_test_fasta();
        let mut tx = plus_strand_coding();
        tx.cds_start_complete = false;
        let index = LocateIndex::build(&tx).unwrap();
        assert!(!super::super::helpers::has_complete_start_codon(
            &tx, &index
        ));

        let snv = annotate_snv("chr1", 1500, b'A', b'G', &tx, &index, &fasta).unwrap();
        assert_eq!(snv.consequences, vec![Consequence::MissenseVariant]);
        assert_eq!(snv.hgvs_p.as_deref(), Some("p.Met1Val"));

        let deletion = annotate_deletion("chr1", 1500, 1503, b"ATG", &tx, &index, &fasta).unwrap();
        assert!(
            deletion
                .consequences
                .contains(&Consequence::InframeDeletion)
        );
        assert!(!deletion.consequences.contains(&Consequence::StartLost));
        assert_ne!(deletion.hgvs_p.as_deref(), Some("p.Met1?"));
    }

    #[test]
    fn gapped_cds_alignment_drops_transcript_numbering() {
        let (_tmp, fasta) = write_test_fasta();
        let mut tx = plus_strand_coding();
        tx.cds_alignment_gap = true;
        let index = LocateIndex::build(&tx).unwrap();
        let snv = annotate_snv("chr1", 1503, b'C', b'A', &tx, &index, &fasta).unwrap();
        assert!(!snv.consequences.is_empty());
        let unnumbered = |r: &ConsequenceResult| {
            r.hgvs_c.is_none()
                && r.hgvs_p.is_none()
                && r.cdna_position.is_none()
                && r.cdna_position_end.is_none()
                && r.cds_position.is_none()
                && r.cds_position_end.is_none()
                && r.protein_start.is_none()
                && r.protein_end.is_none()
                && r.codons.is_none()
                && r.amino_acids.is_none()
                && !r.predicts_nmd
                && matches!(r.ptc, PtcStatus::NotApplicable | PtcStatus::Indeterminate)
                && r.exon.is_some()
        };
        assert!(unnumbered(&snv), "{snv:?}");
        // A one-base deletion is a frameshift: its termination codon is
        // undetermined, not located.
        let deletion = annotate_deletion("chr1", 1503, 1504, b"C", &tx, &index, &fasta).unwrap();
        assert!(unnumbered(&deletion), "{deletion:?}");
        assert_eq!(deletion.ptc, PtcStatus::Indeterminate);
        let insertion = annotate_insertion("chr1", 1504, b"T", &tx, &index, &fasta).unwrap();
        assert!(unnumbered(&insertion), "{insertion:?}");
        // The dispatcher's MNV and delins paths.
        let store = crate::TranscriptStore::from_transcripts(vec![tx.clone()]);
        for (reference, alternate) in [(&b"CG"[..], &b"TA"[..]), (&b"CGTC"[..], &b"TG"[..])] {
            let results = annotate("chr1", 1503, reference, alternate, &store, &fasta).unwrap();
            assert!(!results.is_empty());
            assert!(results.iter().all(unnumbered), "{results:?}");
        }

        tx.cds_alignment_gap = false;
        let index = LocateIndex::build(&tx).unwrap();
        let snv = annotate_snv("chr1", 1503, b'C', b'A', &tx, &index, &fasta).unwrap();
        assert!(snv.hgvs_c.is_some() && snv.cds_position.is_some() && snv.codons.is_some());
    }

    #[test]
    fn slippage_cds_is_refused_at_index_build() {
        let mut tx = plus_strand_coding();
        tx.cds_segments.insert(
            1,
            CdsSegment {
                exon_index: 0,
                genomic_start: 1_900,
                genomic_end: 2_000,
                phase: 0,
            },
        );
        tx.cds_segments[0].genomic_end = 1_901;
        let err = LocateIndex::build(&tx).unwrap_err();
        assert!(
            err.to_string()
                .contains("more than one CDS segment in one exon"),
            "{err}"
        );
    }

    /// Two-exon coding transcript on chr1 whose single intron is `gap` bases.
    fn transcript_with_gap(gap: u64) -> TranscriptModel {
        let mut tx = plus_strand_coding();
        tx.exons.truncate(2);
        tx.exons[1].genomic_start = 2_000 + gap;
        tx.exons[1].genomic_end = 3_500;
        tx.cds_segments.truncate(2);
        tx.cds_segments[1].genomic_start = 2_000 + gap;
        tx.cds_genomic_end = Some(3_500);
        tx.tx_end = 3_500;
        tx.exon_count = 2;
        tx
    }

    #[test]
    fn frameshift_intron_bound_is_a_thirteen_base_gap() {
        for (gap, frameshift) in [(1u64, true), (12, true), (13, true), (14, false)] {
            let tx = transcript_with_gap(gap);
            let index = LocateIndex::build(&tx).unwrap();
            let location = locate_variant("chr1", 2_000, 2_001, &tx, &index).unwrap();
            assert_eq!(
                matches!(location, VariantLocation::FrameshiftIntron { .. }),
                frameshift,
                "gap {gap}: {location:?}"
            );
        }
    }

    #[test]
    fn frameshift_intron_snv_is_coding_sequence_variant_with_intronic_hgvs() {
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_snv("chr1", 2_000, b'A', b'G', &tx, &index, &fasta).unwrap();
        // VEP 116 (numbers=1): ENST00000400999.7 c.228+1T>A has intron 2/5
        // and no exon number.
        assert_eq!(
            result.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(result.intron.as_deref(), Some("1/1"));
        assert_eq!(result.exon, None);
        assert_eq!(result.hgvs_c.as_deref(), Some("NM_TEST_PLUS.1:c.500+1A>G"));
    }

    #[test]
    fn frameshift_intron_endpoint_follows_vep_overlap_rules() {
        // Gap [2000, 2002). VEP: an end inside the gap makes the variant
        // `within_frameshift_intron`; an insertion needs both flanking bases
        // inside it; a deletion spanning the gap from exon to exon does not
        // count (ENST00000400999.7 c.228_228+1del vs c.226_229del,
        // c.228_228+1insA).
        let tx = transcript_with_gap(2);
        let index = LocateIndex::build(&tx).unwrap();
        assert_eq!(intron_gap(&tx.exons[0], &tx.exons[1]), 2);
        assert_eq!(
            frameshift_intron_endpoint(1_999, 2_001, &tx, &index),
            Some(0)
        );
        assert_eq!(frameshift_intron_endpoint(1_998, 2_004, &tx, &index), None);
        assert_eq!(
            frameshift_intron_endpoint(2_001, 2_001, &tx, &index),
            Some(0)
        );
        assert_eq!(frameshift_intron_endpoint(2_000, 2_000, &tx, &index), None);
        assert_eq!(frameshift_intron_endpoint(2_002, 2_002, &tx, &index), None);
        let normal = transcript_with_gap(14);
        let normal_index = LocateIndex::build(&normal).unwrap();
        assert_eq!(
            frameshift_intron_endpoint(1_999, 2_001, &normal, &normal_index),
            None
        );
    }

    #[test]
    fn one_base_gap_windows_spill_onto_flanking_exon_bases() {
        // Plus strand, gap [2000, 2001): donor window [2000, 2002) holds the
        // first exon base after the gap, acceptor window [1999, 2001) the
        // last base before it (VEP 116 ENST00000400999.7 c.229G>A
        // splice_donor_variant, c.228C>A splice_acceptor_variant).
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let after = exonic_splice_spill(2_001, 2_002, &tx, &index).unwrap();
        assert!(after.donor && !after.acceptor);
        let before = exonic_splice_spill(1_999, 2_000, &tx, &index).unwrap();
        assert!(before.acceptor && !before.donor);
        assert_eq!(exonic_splice_spill(1_997, 1_998, &tx, &index), None);
        // A footprint that touches the gap is not a spill.
        assert_eq!(exonic_splice_spill(1_999, 2_001, &tx, &index), None);
        // Two-base gaps keep their windows inside the gap.
        let two = transcript_with_gap(2);
        let two_index = LocateIndex::build(&two).unwrap();
        assert_eq!(exonic_splice_spill(2_002, 2_003, &two, &two_index), None);
    }

    #[test]
    fn cds_phase_above_two_is_refused() {
        let mut tx = plus_strand_coding();
        tx.cds_segments[1].phase = 3;
        let err = LocateIndex::build(&tx).unwrap_err();
        assert!(err.to_string().contains("CDS phase outside 0..=2"), "{err}");
    }

    #[test]
    fn deletion_spanning_a_frameshift_intron_is_a_frameshift_without_splice_terms() {
        // VEP 116 ENST00000400999.7 c.226_229del across a one-base gap:
        // frameshift_variant, no splice term, exon 2-3/6, intron 2/5.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_deletion("chr1", 1_999, 2_002, b"GAA", &tx, &index, &fasta).unwrap();
        assert!(
            result
                .consequences
                .contains(&Consequence::FrameshiftVariant)
        );
        assert!(!result.consequences.iter().any(|c| matches!(
            c,
            Consequence::SpliceDonorVariant | Consequence::SpliceAcceptorVariant
        )));
        assert_eq!(result.exon.as_deref(), Some("1-2/2"));
        assert_eq!(result.intron.as_deref(), Some("1/1"));
    }

    #[test]
    fn deletion_ending_inside_a_frameshift_intron_is_coding_sequence_variant() {
        // VEP 116 ENST00000400999.7 c.228_228+1del: coding_sequence_variant,
        // exon 2/6, intron 2/5, CDS 228.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_deletion("chr1", 1_999, 2_001, b"GA", &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(result.exon.as_deref(), Some("1/2"));
        assert_eq!(result.intron.as_deref(), Some("1/1"));
        assert_eq!(result.cds_position, Some(500));
    }

    #[test]
    fn snvs_beside_a_one_base_gap_get_the_spilled_splice_term() {
        // VEP 116 ENST00000400999.7: c.228C>A splice_acceptor_variant +
        // synonymous_variant, c.229G>A splice_donor_variant + missense, no
        // splice_region_variant on either.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let before = annotate_snv("chr1", 1_999, b'G', b'A', &tx, &index, &fasta).unwrap();
        assert!(
            before
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant)
        );
        assert!(
            !before
                .consequences
                .contains(&Consequence::SpliceRegionVariant)
        );
        assert_eq!(before.impact, Impact::High);
        let after = annotate_snv("chr1", 2_001, b'A', b'G', &tx, &index, &fasta).unwrap();
        assert!(
            after
                .consequences
                .contains(&Consequence::SpliceDonorVariant)
        );
        assert!(
            !after
                .consequences
                .contains(&Consequence::SpliceRegionVariant)
        );
    }

    #[test]
    fn insertion_at_a_frameshift_intron_junction_is_a_frameshift_with_splice_term() {
        // VEP 116 ENST00000400999.7 c.228_228+1insA: splice_acceptor_variant
        // + frameshift_variant, no exon number.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_insertion("chr1", 2_000, b"C", &tx, &index, &fasta).unwrap();
        assert!(
            result
                .consequences
                .contains(&Consequence::FrameshiftVariant)
        );
        assert!(
            result
                .consequences
                .contains(&Consequence::SpliceAcceptorVariant)
        );
        assert_eq!(result.exon, None);
    }

    /// [`plus_strand_coding`] with the last CDS base dropped: 1499 coding
    /// bases, so the last codon is the partial `TA` at chr1:4497-4498.
    fn partial_last_codon() -> (TranscriptModel, LocateIndex) {
        let mut tx = plus_strand_coding();
        tx.cds_segments[2].genomic_end = 4_499;
        tx.cds_genomic_end = Some(4_499);
        let index = LocateIndex::build(&tx).unwrap();
        (tx, index)
    }

    #[test]
    fn deletions_in_a_partial_last_codon_follow_vep() {
        // VEP 116 NM_001424184.1: c.658del gives incomplete_terminal_codon +
        // coding_sequence_variant (`Ta/a`, X); c.658_659del removes the whole
        // partial codon: protein_altering + incomplete_terminal (`TA/-`, X/-).
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = partial_last_codon();
        let part = annotate_deletion("chr1", 4_497, 4_498, b"T", &tx, &index, &fasta).unwrap();
        assert_eq!(
            part.consequences,
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant
            ]
        );
        assert_eq!(part.codons.as_deref(), Some("Ta/a"));
        assert_eq!(part.amino_acids.as_deref(), Some("X"));
        assert_eq!(part.hgvs_p.as_deref(), Some("p.?"));
        let whole = annotate_deletion("chr1", 4_497, 4_499, b"TA", &tx, &index, &fasta).unwrap();
        assert_eq!(
            whole.consequences,
            vec![
                Consequence::ProteinAlteringVariant,
                Consequence::IncompleteTerminalCodonVariant
            ]
        );
        assert_eq!(whole.codons.as_deref(), Some("TA/-"));
        assert_eq!(whole.amino_acids.as_deref(), Some("X/-"));
    }

    #[test]
    fn insertion_completing_a_stop_in_a_partial_codon_is_stop_gained() {
        // VEP 116 NM_001424184.1 c.658_659insG: `ta/tGa`, X/*, stop_gained.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = partial_last_codon();
        let result = annotate_insertion("chr1", 4_498, b"G", &tx, &index, &fasta).unwrap();
        assert!(result.consequences.contains(&Consequence::StopGained));
        assert_eq!(result.codons.as_deref(), Some("ta/tGa"));
        assert_eq!(result.amino_acids.as_deref(), Some("X/*"));
        assert_eq!(result.cds_position, Some(1_498));
    }

    #[test]
    fn mnv_inside_the_padded_codon_is_coding_sequence_variant() {
        // VEP 116 ENST00000455979.1 (pad 1) c.1_2inv: coding_sequence_variant
        // only, never synonymous; X-only change gets `p.?`.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(2);
        let result = annotate_mnv("chr1", 1_500, b"AT", b"GC", &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(result.hgvs_p.as_deref(), Some("p.?"));
    }

    #[test]
    fn padded_transcript_numbers_indels_on_genomic_bases() {
        // VEP 116 ENST00000455979.1 (pad 1): the SNV at CDS position 5 is
        // c.5, the deletion of the same base is c.4del.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(2);
        let snv = annotate_snv("chr1", 1_503, b'C', b'A', &tx, &index, &fasta).unwrap();
        assert_eq!(snv.hgvs_c.as_deref(), Some("NM_TEST_PLUS.1:c.5C>A"));
        let deletion = annotate_deletion("chr1", 1_503, 1_504, b"C", &tx, &index, &fasta).unwrap();
        assert_eq!(deletion.cds_position, Some(5));
        assert_eq!(deletion.hgvs_c.as_deref(), Some("NM_TEST_PLUS.1:c.4del"));
    }

    #[test]
    fn fetch_cds_sequence_refuses_offsets_past_the_cds() {
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = partial_last_codon();
        let err =
            super::super::helpers::fetch_cds_sequence(1_497, 1_500, "chr1", &tx, &index, &fasta)
                .unwrap_err();
        assert!(matches!(err, VarEffectError::Malformed(_)), "{err}");
    }

    #[test]
    fn reverse_c_of_an_indel_on_a_padded_transcript_counts_genomic_bases() {
        // VEP 116 numbers indels without the pad (ENST00000455979.1 c.4del
        // is the base at CDS position 5), so `c.4del` with a one-base pad is
        // the fourth genomic coding base, chr1:1503.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, _) = padded_plus(2);
        let engine =
            crate::VarEffect::new(crate::TranscriptStore::from_transcripts(vec![tx]), fasta);
        let variant = engine.resolve_hgvs_c("NM_TEST_PLUS.1:c.4del").unwrap();
        assert_eq!(variant.pos, 1_502);
        assert_eq!(variant.ref_allele.len(), 2);
    }

    #[test]
    fn deletion_from_the_5utr_into_a_padded_cds_is_only_a_coding_change() {
        // No start codon to lose, and no defined cds_start for VEP's frame
        // predicates.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = padded_plus(1);
        let result = annotate_deletion("chr1", 1_498, 1_502, b"AAAT", &tx, &index, &fasta).unwrap();
        let mut terms = result.consequences.clone();
        terms.sort();
        assert_eq!(
            terms,
            vec![
                Consequence::CodingSequenceVariant,
                Consequence::FivePrimeUtrVariant
            ]
        );
    }

    #[test]
    fn insertion_at_a_two_base_gap_edge_adds_splice_region_and_intron_terms() {
        // VEP 116 ENST00000436261.7 c.1475_1475+1insT (two-base gap):
        // frameshift_variant + splice_region_variant + intron_variant.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(2);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_insertion("chr1", 2_000, b"C", &tx, &index, &fasta).unwrap();
        for term in [
            Consequence::FrameshiftVariant,
            Consequence::SpliceRegionVariant,
            Consequence::IntronVariant,
        ] {
            assert!(
                result.consequences.contains(&term),
                "{:?}",
                result.consequences
            );
        }
        assert_eq!(result.exon, None);
    }

    #[test]
    fn deletion_starting_in_a_gap_has_no_cds_position() {
        // VEP 116 ENST00000400999.7 c.228+1_229del: coding_sequence_variant,
        // exon 3/6, intron 2/5, no CDS positions.
        let (_tmp, fasta) = write_test_fasta();
        let tx = transcript_with_gap(1);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_deletion("chr1", 2_000, 2_002, b"AA", &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(result.exon.as_deref(), Some("2/2"));
        assert_eq!(result.intron.as_deref(), Some("1/1"));
        assert_eq!(result.cds_position, None);
    }

    #[test]
    fn non_coding_deletion_across_a_gap_is_an_exon_variant() {
        // VEP 116 ENST00000609389.5 n.281_282del across a 12-base gap:
        // non_coding_transcript_exon_variant, exon 3-4/6, intron 3/5.
        let (_tmp, fasta) = write_test_fasta();
        let mut tx = noncoding_2_exon();
        tx.exons[1].genomic_start = 305;
        let index = LocateIndex::build(&tx).unwrap();
        let result =
            annotate_deletion("chr2", 298, 308, b"AAAAAAAAAA", &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::NonCodingTranscriptExonVariant]
        );
        assert_eq!(result.exon.as_deref(), Some("1-2/2"));
        assert_eq!(result.intron.as_deref(), Some("1/1"));
    }

    #[test]
    fn deletions_reaching_the_partial_codon_from_an_earlier_codon() {
        // `partial_last_codon`: complete codons end at CDS offset 1497, the
        // partial codon is chr1:4497-4498.
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = partial_last_codon();
        // Stays in the CDS: a frameshift of the complete codon it starts in
        // (VEP 116 NM_001424184.1 c.657_658del), protein effect unknown.
        let shift = annotate_deletion("chr1", 4_496, 4_498, b"GT", &tx, &index, &fasta).unwrap();
        assert!(shift.consequences.contains(&Consequence::FrameshiftVariant));
        assert_eq!(shift.hgvs_p.as_deref(), Some("p.?"));
        assert_eq!(shift.ptc, PtcStatus::Indeterminate);
        // Removes the CDS end: no stop term (maintainer decision; VEP gives
        // frameshift_variant + stop_lost for c.656_659del).
        let to_end = annotate_deletion("chr1", 4_496, 4_499, b"GTA", &tx, &index, &fasta).unwrap();
        assert_eq!(
            to_end.consequences,
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant
            ]
        );
        assert_eq!(to_end.hgvs_p, None);
        // Runs on into the 3'UTR: coding and UTR change, still no stop term
        // (VEP gives stop_retained_variant for c.656_*2del).
        let into_utr =
            annotate_deletion("chr1", 4_495, 4_502, b"CGTAAAA", &tx, &index, &fasta).unwrap();
        assert!(
            into_utr
                .consequences
                .contains(&Consequence::ThreePrimeUtrVariant)
        );
        assert!(
            into_utr
                .consequences
                .contains(&Consequence::IncompleteTerminalCodonVariant)
        );
        assert!(!into_utr.consequences.iter().any(|c| matches!(
            c,
            Consequence::StopLost
                | Consequence::StopRetainedVariant
                | Consequence::FrameshiftVariant
        )));
    }

    #[test]
    fn minus_strand_insertion_in_a_partial_codon_goes_after_the_classified_base() {
        // VEP 116 ENST00000518351.2 c.544_545insG: inserted between the two
        // partial-codon bases. Here the partial codon is coding `TA`
        // (chr17:11002, 11001); a coding `G` between them completes `TGA`.
        let (_tmp, fasta) = write_test_fasta();
        let mut tx = minus_strand_coding();
        tx.cds_segments[2].genomic_start = 11_001;
        tx.cds_genomic_start = Some(11_001);
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_insertion("chr17", 11_002, b"C", &tx, &index, &fasta).unwrap();
        assert!(
            result.consequences.contains(&Consequence::StopGained),
            "{result:?}"
        );
        assert_eq!(result.codons.as_deref(), Some("ta/tGa"));
        assert_eq!(result.cds_position, Some(4_498));
        assert!(
            !matches!(result.ptc, PtcStatus::NotApplicable),
            "{:?}",
            result.ptc
        );
    }

    /// [`minus_strand_coding`] with the last CDS base dropped: 4499 coding
    /// bases, partial codon `TA` at chr17:11002 (coding offset 4497) and
    /// chr17:11001 (4498), 3'UTR below.
    fn minus_partial_last_codon() -> (TranscriptModel, LocateIndex) {
        let mut tx = minus_strand_coding();
        tx.cds_segments[2].genomic_start = 11_001;
        tx.cds_genomic_start = Some(11_001);
        let index = LocateIndex::build(&tx).unwrap();
        (tx, index)
    }

    #[test]
    fn insertion_at_the_partial_codon_edge_follows_vep_on_both_strands() {
        // VEP 116 NM_001424184.1 c.657_658insG (plus) and ENST00000518351.2
        // c.543dup (minus): inframe_insertion + incomplete_terminal_codon +
        // coding_sequence_variant, MODERATE, `-/X`, no HGVS p.
        let (_tmp, fasta) = write_test_fasta();
        let expected = vec![
            Consequence::InframeInsertion,
            Consequence::IncompleteTerminalCodonVariant,
            Consequence::CodingSequenceVariant,
        ];
        let (plus, plus_index) = partial_last_codon();
        let plus_edge =
            annotate_insertion("chr1", 4_497, b"G", &plus, &plus_index, &fasta).unwrap();
        let mut plus_terms = plus_edge.consequences.clone();
        plus_terms.sort();
        let mut want = expected.clone();
        want.sort();
        assert_eq!(plus_terms, want);
        assert_eq!(plus_edge.impact, Impact::Moderate);
        assert_eq!(plus_edge.amino_acids.as_deref(), Some("-/X"));
        assert_eq!(plus_edge.hgvs_p, None);
        assert_eq!(plus_edge.cds_position, Some(1_497));

        let (minus, minus_index) = minus_partial_last_codon();
        let minus_edge =
            annotate_insertion("chr17", 11_003, b"C", &minus, &minus_index, &fasta).unwrap();
        let mut minus_terms = minus_edge.consequences.clone();
        minus_terms.sort();
        assert_eq!(minus_terms, want);
        assert_eq!(minus_edge.cds_position, Some(4_497));
    }

    #[test]
    fn edge_insertions_translate_the_inserted_bases() {
        // VEP 116 NM_001424184.1 after c.657: insGG `-/X` with
        // coding_sequence_variant, insGGG `-/G` without it, insTGA and
        // insTGAG `-/*` with stop_gained (HIGH); ENST00000518351.2 insTAG
        // (VCF `CTA`) the same on the minus strand.
        let (_tmp, fasta) = write_test_fasta();
        let (plus, plus_index) = partial_last_codon();
        let edge = |bases: &[u8]| {
            let mut result =
                annotate_insertion("chr1", 4_497, bases, &plus, &plus_index, &fasta).unwrap();
            result.consequences.sort();
            result
        };
        let mut base_terms = vec![
            Consequence::InframeInsertion,
            Consequence::IncompleteTerminalCodonVariant,
        ];
        base_terms.sort();
        let with = |extra: Consequence| {
            let mut terms = base_terms.clone();
            terms.push(extra);
            terms.sort();
            terms
        };

        let two = edge(b"GG");
        assert_eq!(two.consequences, with(Consequence::CodingSequenceVariant));
        assert_eq!(two.amino_acids.as_deref(), Some("-/X"));

        let three = edge(b"GGG");
        assert_eq!(three.consequences, base_terms);
        assert_eq!(three.amino_acids.as_deref(), Some("-/G"));
        assert_eq!(three.hgvs_p, None);
        assert_eq!(three.ptc, PtcStatus::NotApplicable);

        // Codons after a stop are kept; a left-over partial codon is `X`
        // unless the complete codons are a lone stop (VEP 116 insGGGTGAG).
        assert_eq!(edge(b"TGAGGG").amino_acids.as_deref(), Some("-/*G"));
        let stop_then_partial = edge(b"GGGTGAG");
        assert_eq!(stop_then_partial.amino_acids.as_deref(), Some("-/G*X"));
        assert!(
            stop_then_partial
                .consequences
                .contains(&Consequence::CodingSequenceVariant)
        );

        for bases in [&b"TGA"[..], b"TGAG"] {
            let stop = edge(bases);
            assert_eq!(stop.consequences, with(Consequence::StopGained));
            assert_eq!(stop.impact, Impact::High);
            assert_eq!(stop.amino_acids.as_deref(), Some("-/*"));
            assert_eq!(stop.hgvs_p.as_deref(), Some("p.?"));
            assert!(
                matches!(
                    stop.ptc,
                    PtcStatus::At {
                        protein_position: 500,
                        ..
                    }
                ),
                "{:?}",
                stop.ptc
            );
        }

        let (minus, minus_index) = minus_partial_last_codon();
        let mut minus_stop =
            annotate_insertion("chr17", 11_003, b"CTA", &minus, &minus_index, &fasta).unwrap();
        minus_stop.consequences.sort();
        assert_eq!(minus_stop.consequences, with(Consequence::StopGained));
        assert_eq!(minus_stop.codons.as_deref(), Some("-/TAG"));
    }

    #[test]
    fn insertion_inside_the_partial_codon_translates_every_completed_codon() {
        // VEP 116 NM_001424184.1 c.658_659insGGTGA: `ta/tGGTGAa`, `X/W*X`,
        // stop_gained (HIGH). Here the partial codon is `TA` (chr1:4497-4498).
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = partial_last_codon();
        let result = annotate_insertion("chr1", 4_498, b"GGTGA", &tx, &index, &fasta).unwrap();
        assert!(result.consequences.contains(&Consequence::StopGained));
        assert_eq!(result.impact, Impact::High);
        assert_eq!(result.amino_acids.as_deref(), Some("X/W*X"));
        assert!(
            matches!(
                result.ptc,
                PtcStatus::At {
                    protein_position: 501,
                    ..
                }
            ),
            "{:?}",
            result.ptc
        );
    }

    #[test]
    fn a_cds_shorter_than_a_codon_has_no_partial_codon_edge() {
        // With no complete codon there is no 5' edge: an insertion at offset
        // 0 is inside the only (partial) codon.
        let mut tx = plus_strand_coding();
        tx.cds_segments.truncate(1);
        tx.cds_segments[0].genomic_end = 1_502;
        tx.cds_genomic_end = Some(1_502);
        let index = LocateIndex::build(&tx).unwrap();
        assert_eq!(
            super::indel::partial_codon_insertion_site(0, &index),
            Some(super::indel::PartialCodonInsertion::Inside)
        );
    }

    #[test]
    fn minus_strand_insertion_after_the_cds_end_without_a_utr_is_downstream() {
        // VEP 116 ENST00000518351.2, whose CDS ends at the transcript end:
        // downstream_gene_variant, no exon, cDNA or HGVS.
        let (_tmp, fasta) = write_test_fasta();
        let (mut tx, _) = minus_partial_last_codon();
        let last = tx.exons.len() - 1;
        tx.exons[last].genomic_start = 11_001;
        tx.tx_start = 11_001;
        let index = LocateIndex::build(&tx).unwrap();
        let result = annotate_insertion("chr17", 11_001, b"C", &tx, &index, &fasta).unwrap();
        assert_eq!(
            result.consequences,
            vec![Consequence::DownstreamGeneVariant]
        );
        assert_eq!(result.exon, None);
        assert_eq!(result.cdna_position, None);
        assert_eq!(result.hgvs_c, None);
    }

    #[test]
    fn minus_strand_insertion_after_the_last_cds_base_is_utr() {
        // Between the last CDS base (chr17:11001) and the 3'UTR: outside the
        // CDS (VEP 116 ENST00000654733.1 c.*1dup: 3_prime_UTR_variant).
        let (_tmp, fasta) = write_test_fasta();
        let (tx, index) = minus_partial_last_codon();
        let result = annotate_insertion("chr17", 11_001, b"C", &tx, &index, &fasta).unwrap();
        assert!(
            result
                .consequences
                .contains(&Consequence::ThreePrimeUtrVariant),
            "{:?}",
            result.consequences
        );
        assert!(!result.consequences.contains(&Consequence::StopGained));
    }

    #[test]
    fn frameshift_gap_is_not_the_last_exon_junction_for_ptc() {
        // Exons 2 and 3 of `plus_strand_coding` separated by a one-base gap:
        // the last real exon-exon junction is the one before exon 2 (CDS
        // position 501), so a termination codon at residue 200 (bases
        // 598-600) lies in the last real exon and escapes. `predicts_nmd`
        // keeps VEP's NMD plugin behavior and still counts the gap.
        let mut tx = plus_strand_coding();
        tx.exons[2].genomic_start = 3_501;
        tx.cds_segments[2].genomic_start = 3_501;
        let index = LocateIndex::build(&tx).unwrap();
        let edit = crate::hgvs_p::CdsEdit {
            start: 599,
            del_len: 1,
            ins_len: 1,
        };
        let ptc =
            super::super::nmd::resolve_ptc(crate::hgvs_p::PtcSite::Residue(200), &edit, &index);
        assert!(
            matches!(
                ptc,
                PtcStatus::At {
                    nmd_at_ptc: false,
                    ..
                }
            ),
            "{ptc:?}"
        );
        assert!(super::super::nmd::predicts_nmd(600, &index));
        let normal = LocateIndex::build(&plus_strand_coding()).unwrap();
        let normal_ptc =
            super::super::nmd::resolve_ptc(crate::hgvs_p::PtcSite::Residue(200), &edit, &normal);
        assert!(
            matches!(
                normal_ptc,
                PtcStatus::At {
                    nmd_at_ptc: true,
                    ..
                }
            ),
            "{normal_ptc:?}"
        );
    }

    #[test]
    fn non_coding_transcript_variant_is_a_modifier() {
        let term = Consequence::NonCodingTranscriptVariant;
        assert_eq!(term.as_str(), "non_coding_transcript_variant");
        assert_eq!(term.impact(), Impact::Modifier);
        assert!(term.severity_rank() > Consequence::IntronVariant.severity_rank());
        assert!(term.severity_rank() < Consequence::UpstreamGeneVariant.severity_rank());
    }

    #[test]
    fn reverse_c_inside_the_pad_has_no_genomic_position() {
        let (_tmp, fasta) = write_test_fasta();
        let (tx, _) = padded_plus(1);
        let store = crate::TranscriptStore::from_transcripts(vec![tx]);
        let engine = crate::VarEffect::new(store, fasta);
        let err = engine.resolve_hgvs_c("NM_TEST_PLUS.1:c.2A>G").unwrap_err();
        assert!(
            matches!(err, VarEffectError::PositionOutOfRange { .. }),
            "{err}"
        );
        let variant = engine.resolve_hgvs_c("NM_TEST_PLUS.1:c.3A>G").unwrap();
        assert_eq!(variant.pos, 1_500);
    }
}

// -- Insertion HGVS 3' shift ------------------------------------------------

/// Plus-strand genome for the shift tests: exon 1 `[10, 40)` (CDS offsets
/// 0..30, `ATG GAT CAG CCC.. CAG`), intron `[40, 60)` starting `GTAAG`, exon 2
/// `[60, 90)` (CDS offsets 30..57 ending `TAA`, then 3 UTR bases).
fn shift_plus_genome() -> Vec<u8> {
    let mut genome = b"CCCCCCCCCC".to_vec();
    genome.extend_from_slice(b"ATGGATCAG");
    genome.extend_from_slice(&b"CCC".repeat(6));
    genome.extend_from_slice(b"CAG");
    genome.extend_from_slice(b"GTAAG");
    genome.extend_from_slice(&b"C".repeat(13));
    genome.extend_from_slice(b"AG");
    genome.extend_from_slice(&b"CCC".repeat(8));
    genome.extend_from_slice(b"TAA");
    genome.extend_from_slice(&b"C".repeat(13));
    assert_eq!(genome.len(), 100);
    genome
}

/// The shift-test genome on chr5 and its reverse complement on chr6.
fn write_shift_fasta() -> (TempDir, FastaReader) {
    let tmp = TempDir::new().unwrap();
    let plus = shift_plus_genome();
    let minus = crate::codon::reverse_complement(&plus);
    let bin_path = tmp.path().join("shift.bin");
    let idx_path = tmp.path().join("shift.bin.idx");
    let contigs: &[(&str, &[u8])] = &[("chr5", plus.as_slice()), ("chr6", minus.as_slice())];
    write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();
    (tmp, FastaReader::open(&bin_path).unwrap())
}

/// Two-exon transcript over [`shift_plus_genome`]: on chr5 for the plus
/// strand, mirrored onto chr6 (`x -> 99 - x`) for the minus strand.
fn shift_transcript(strand: Strand) -> TranscriptModel {
    let mirror = |start: u64, end: u64| match strand {
        Strand::Plus => (start, end),
        Strand::Minus => (100 - end, 100 - start),
    };
    let exon = |number: u16, start, end| {
        let (genomic_start, genomic_end) = mirror(start, end);
        Exon {
            exon_number: number,
            genomic_start,
            genomic_end,
        }
    };
    let segment = |exon_index: u16, start, end| {
        let (genomic_start, genomic_end) = mirror(start, end);
        CdsSegment {
            exon_index,
            genomic_start,
            genomic_end,
            phase: 0,
        }
    };
    let (tx_start, tx_end) = mirror(10, 90);
    let (cds_start, cds_end) = mirror(10, 87);
    TranscriptModel {
        accession: "NM_SHIFT.1".into(),
        protein_accession: Some("NP_SHIFT.1".into()),
        gene_symbol: "SHIFT".into(),
        hgnc_id: None,
        ensembl_accession: None,
        chrom: match strand {
            Strand::Plus => "chr5".into(),
            Strand::Minus => "chr6".into(),
        },
        strand,
        tx_start,
        tx_end,
        cds_genomic_start: Some(cds_start),
        cds_genomic_end: Some(cds_end),
        exons: vec![exon(1, 10, 40), exon(2, 60, 90)],
        cds_segments: vec![segment(0, 10, 40), segment(1, 60, 87)],
        tier: Some(TranscriptTier::ManeSelect),
        biotype: Biotype::ProteinCoding,
        exon_count: 2,
        cds_start_complete: true,
        cds_end_complete: true,
        genome_mismatch: false,
        cds_alignment_gap: false,
        selenocysteine_codons: Vec::new(),
    }
}

/// Annotate an insertion given as a plus-strand point and bases on the
/// plus-strand genome, mirrored for the minus strand.
fn annotate_shift_insertion(strand: Strand, pos: u64, bases: &[u8]) -> ConsequenceResult {
    let (_tmp, fasta) = write_shift_fasta();
    let tx = shift_transcript(strand);
    let idx = build_index(&tx);
    let (pos, bases) = match strand {
        Strand::Plus => (pos, bases.to_vec()),
        Strand::Minus => (100 - pos, crate::codon::reverse_complement(bases)),
    };
    annotate_insertion(&tx.chrom, pos, &bases, &tx, &idx, &fasta).unwrap()
}

#[test]
fn shifted_inframe_insertion_rotates_bases() {
    // ATA inserted after c.4 (GA|TCAG) shifts two bases to c.6_7 as AAT
    // (VEP 116 NM_001378454.1 c.7539_7540insAAT p.His2513_Ala2514insAsn,
    // where the unrotated bases gave insIle).
    for strand in [Strand::Plus, Strand::Minus] {
        let result = annotate_shift_insertion(strand, 14, b"ATA");
        assert_eq!(
            result.hgvs_c.as_deref(),
            Some("NM_SHIFT.1:c.6_7insAAT"),
            "{strand:?}"
        );
        assert_eq!(
            result.hgvs_p.as_deref(),
            Some("p.Asp2_Gln3insAsn"),
            "{strand:?}"
        );
        assert_eq!(result.consequences, vec![Consequence::InframeInsertion]);
    }
}

#[test]
fn junction_insertion_shifts_into_intron() {
    // GT inserted after the last exon base c.30 shifts into the intron, a
    // duplication of +1/+2 with no HGVS p.; the consequences stay on the
    // junction insertion (VEP 116 NM_000548.5 c.336+1_336+2dup:
    // frameshift_variant, splice_region_variant, no HGVS p.).
    for strand in [Strand::Plus, Strand::Minus] {
        let result = annotate_shift_insertion(strand, 40, b"GT");
        assert_eq!(
            result.hgvs_c.as_deref(),
            Some("NM_SHIFT.1:c.30+1_30+2dup"),
            "{strand:?}"
        );
        assert_eq!(result.hgvs_p, None, "{strand:?}");
        assert_eq!(
            result.consequences,
            vec![
                Consequence::FrameshiftVariant,
                Consequence::SpliceRegionVariant
            ],
            "{strand:?}"
        );
    }
}

#[test]
fn intronic_insertion_shifts_into_exon() {
    // GCC inserted between -2 and -1 (A|G CCC) shifts three bases into
    // exon 2, a duplication of -1 and c.31_32 (VEP 116 NM_001005361.3
    // c.691_699dup where the edge-bounded shift gave
    // c.689-1_689insGGCTACATT).
    for strand in [Strand::Plus, Strand::Minus] {
        let result = annotate_shift_insertion(strand, 59, b"GCC");
        assert_eq!(
            result.hgvs_c.as_deref(),
            Some("NM_SHIFT.1:c.31-1_32dup"),
            "{strand:?}"
        );
        assert_eq!(result.hgvs_p, None, "{strand:?}");
    }
}

#[test]
fn insertion_duplicating_following_bases_is_a_dup_of_them() {
    // CCC before c.13 inside a CCC run (a shift capped short of the run's
    // end) duplicates the 3' bases c.13_15, which VEP checks before the 5'
    // ones (`hgvs_variant_notation`).
    let (_tmp, fasta) = write_shift_fasta();
    for strand in [Strand::Plus, Strand::Minus] {
        let tx = shift_transcript(strand);
        let idx = build_index(&tx);
        let (pos, bases) = match strand {
            Strand::Plus => (22, b"CCC".to_vec()),
            Strand::Minus => (78, b"GGG".to_vec()),
        };
        let hgvs = crate::hgvs_c::format_insertion_hgvs(pos, &bases, &tx.chrom, &tx, &idx, &fasta)
            .unwrap();
        assert_eq!(hgvs.as_deref(), Some("NM_SHIFT.1:c.13_15dup"), "{strand:?}");
    }
}

mod inframe_insertion_terms {
    //! VEP 116 terms of in-frame insertions (`VariationEffect.pm`), on codon
    //! windows taken from VEP REST rows of NM_000546.6, NM_000374.5,
    //! NM_001625.4 and NM_003036.4.

    use super::super::indel::{CdsProtein, InsertionCodons, inframe_insertion_terms};
    use super::*;
    use crate::hgvs_p::ReferenceProtein;

    /// An in-memory reference protein without its terminal stop.
    struct Protein(&'static [u8]);

    impl ReferenceProtein for Protein {
        fn residue_count(&self) -> u32 {
            self.0.len() as u32
        }

        fn residues(&self, start: u32, end: u32) -> Result<Vec<u8>, VarEffectError> {
            Ok(self.0[start as usize..end as usize].to_vec())
        }
    }

    fn terms(
        ref_codon: &[u8],
        inserted: &[u8],
        local: usize,
        codon_index: u32,
        has_utr5: bool,
        protein: &'static [u8],
    ) -> Vec<Consequence> {
        let codons =
            InsertionCodons::from_codon(ref_codon.to_vec(), inserted, local, codon_index, false)
                .expect("whole codons");
        inframe_insertion_terms(&codons, inserted.len(), true, has_utr5, &Protein(protein))
            .expect("in-memory protein")
    }

    const PROTEIN: &[u8] = b"MEEPQSDPSVEPPLSQETFSDLWKLLPEN";

    #[test]
    fn between_codons_is_inframe_insertion() {
        assert_eq!(
            terms(b"", b"TTA", 0, 14, true, PROTEIN),
            vec![Consequence::InframeInsertion]
        );
    }

    #[test]
    fn inside_a_codon_without_the_reference_residue_is_protein_altering() {
        // NM_000546.6 c.523_524insTTT: `cgc/cTTTgc`, `R/LC`.
        assert_eq!(
            terms(b"CGC", b"TTT", 1, 10, true, PROTEIN),
            vec![Consequence::ProteinAlteringVariant]
        );
    }

    #[test]
    fn stop_insertions_inside_a_codon() {
        // NM_000374.5 c.385_386insAGT: `tct/tAGTct`, `S/*S`.
        assert_eq!(
            terms(b"TCT", b"AGT", 1, 10, true, PROTEIN),
            vec![Consequence::StopGained]
        );
        // NM_001625.4 c.597_599dup: `gag/gaGTAg`, `E/E*`.
        let mut got = terms(b"GAG", b"GTA", 2, 10, true, PROTEIN);
        got.sort_by_key(|c| c.severity_rank());
        assert_eq!(
            got,
            vec![Consequence::StopGained, Consequence::InframeInsertion]
        );
    }

    #[test]
    fn insertions_at_the_stop_codon() {
        let stop = PROTEIN.len() as u32;
        // NM_000546.6 c.1180_1181insAAG (`*/*G`) and insGGG (`*/WG`).
        assert_eq!(
            terms(b"TGA", b"AAG", 1, stop, true, PROTEIN),
            vec![
                Consequence::StopRetainedVariant,
                Consequence::InframeInsertion
            ]
        );
        assert_eq!(
            terms(b"TGA", b"GGG", 1, stop, true, PROTEIN),
            vec![Consequence::StopLost]
        );
        // NM_000546.6 c.1179_1180insCTGTAA (`-/L*`): the stop is not next to
        // the protein's last residue, so it is gained.
        assert_eq!(
            terms(b"", b"CTGTAA", 0, stop, true, PROTEIN),
            vec![Consequence::StopGained, Consequence::InframeInsertion]
        );
    }

    #[test]
    fn stop_inserted_after_a_copy_of_the_protein_end_is_retained() {
        // NM_003036.4 c.2152_*3dup shape: the spliced protein is the
        // reference followed by a stop.
        let protein: &'static [u8] = b"MAG";
        assert_eq!(
            terms(b"", b"GCAGGGTAGATA", 0, 1, true, protein),
            vec![
                Consequence::StopRetainedVariant,
                Consequence::InframeInsertion
            ]
        );
    }

    #[test]
    fn start_codon_insertions() {
        // NM_000546.6 c.1_2insCGA (`atg/aCGAtg`) and c.2_3insGAT
        // (`atg/atGATg`).
        assert_eq!(
            terms(b"ATG", b"CGA", 1, 0, true, PROTEIN),
            vec![Consequence::StartLost, Consequence::StartRetainedVariant]
        );
        assert_eq!(
            terms(b"ATG", b"GAT", 2, 0, true, PROTEIN),
            vec![Consequence::StartRetainedVariant]
        );
        // Without a 5'UTR the start-position check does not run
        // (`_inv_start_altered` needs one): derived from the source.
        assert_eq!(
            terms(b"ATG", b"CGA", 1, 0, false, PROTEIN),
            vec![Consequence::StartRetainedVariant]
        );
    }

    #[test]
    fn insertion_codons_between_and_inside_codons() {
        let between = InsertionCodons::from_codon(Vec::new(), b"TTA", 0, 14, false).unwrap();
        assert_eq!(between.peptides.ref_aas, b"");
        assert_eq!(between.peptides.alt_aas, b"L");
        assert_eq!(
            (between.peptides.tl_start, between.peptides.tl_end),
            (15, 14)
        );
        let inside = InsertionCodons::from_codon(b"CGC".to_vec(), b"TTT", 1, 174, false).unwrap();
        assert_eq!(inside.peptides.alt_aas, b"LC");
        assert_eq!(
            (inside.peptides.tl_start, inside.peptides.tl_end),
            (175, 175)
        );
    }

    #[test]
    fn cds_protein_reads_residues_from_the_genome() {
        let (_tmp, fasta) = write_test_fasta();
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        let protein = CdsProtein::new("chr1", &tx, &idx, &fasta).unwrap();
        assert_eq!(protein.residues(0, 3).unwrap(), b"MRR");
        assert!(protein.residues(2, 2).unwrap().is_empty());
    }
}

// -- Unknown alternate bases (VEP 116 `coding_unknown`) ----------------------

#[test]
fn unknown_alt_bases_are_detected_case_insensitively() {
    assert!(super::helpers::has_unknown_bases(b"N"));
    assert!(super::helpers::has_unknown_bases(b"ACNT"));
    assert!(super::helpers::has_unknown_bases(b"R"));
    assert!(!super::helpers::has_unknown_bases(b"acgtACGT"));
    assert!(!super::helpers::has_unknown_bases(b""));
}

#[test]
fn unknown_alt_snv_is_coding_unknown_with_codons_and_no_hgvs() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    let known = annotate_snv("chr1", 1503, b'C', b'A', &tx, &index, &fasta).unwrap();
    let unknown = annotate_snv("chr1", 1503, b'C', b'N', &tx, &index, &fasta).unwrap();
    assert_eq!(
        unknown.consequences,
        vec![Consequence::CodingSequenceVariant]
    );
    assert_eq!(unknown.impact, Impact::Modifier);
    assert_eq!((&unknown.hgvs_c, &unknown.hgvs_p), (&None, &None));
    assert!(unknown.amino_acids.is_none());
    // Codons and the protein position stay, as VEP's `Tct/Nct`.
    assert_eq!(unknown.protein_start, known.protein_start);
    assert!(unknown.codons.as_deref().is_some_and(|c| c.contains('N')));
}

#[test]
fn unknown_alt_mnv_and_delins_lose_peptide_terms() {
    let (_tmp, fasta) = write_test_fasta();
    let store = crate::TranscriptStore::from_transcripts(vec![plus_strand_coding()]);
    for (reference, alternate) in [(&b"CG"[..], &b"NA"[..]), (&b"CGTC"[..], &b"NG"[..])] {
        let results = annotate("chr1", 1503, reference, alternate, &store, &fasta).unwrap();
        for result in &results {
            assert!(
                !result.consequences.iter().any(|c| matches!(
                    c,
                    Consequence::MissenseVariant
                        | Consequence::SynonymousVariant
                        | Consequence::StopGained
                        | Consequence::InframeDeletion
                        | Consequence::ProteinAlteringVariant
                )),
                "{result:?}"
            );
            assert!(
                result.hgvs_c.is_none() && result.hgvs_p.is_none() && result.amino_acids.is_none()
            );
            // CGTC>NG shortens the CDS by two bases: a frameshift whose stop
            // is undetermined.
            let expected = if result
                .consequences
                .contains(&Consequence::FrameshiftVariant)
            {
                PtcStatus::Indeterminate
            } else {
                PtcStatus::NotApplicable
            };
            assert_eq!(result.ptc, expected);
        }
    }
}

#[test]
fn unknown_alt_frameshift_insertion_keeps_frameshift_without_a_located_stop() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    let result = annotate_insertion("chr1", 1504, b"N", &tx, &index, &fasta).unwrap();
    assert_eq!(result.consequences, vec![Consequence::FrameshiftVariant]);
    assert_eq!(result.impact, Impact::High);
    assert!(result.hgvs_c.is_none() && result.hgvs_p.is_none() && result.amino_acids.is_none());
    assert_eq!(result.ptc, PtcStatus::Indeterminate);
}

#[test]
fn unknown_alt_indels_at_the_start_and_stop_codons_follow_vep() {
    // VEP 116 NM_000546.6: c.1_2insNNN is `start_retained_variant`, c.1_2insN
    // `frameshift_variant,start_retained_variant`, and an N insertion or
    // delins in the stop codon `stop_retained_variant` alone.
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    let inframe = annotate_insertion("chr1", 1501, b"NNN", &tx, &index, &fasta).unwrap();
    assert_eq!(
        inframe.consequences,
        vec![Consequence::StartRetainedVariant]
    );
    let frameshift = annotate_insertion("chr1", 1501, b"N", &tx, &index, &fasta).unwrap();
    assert_eq!(
        frameshift.consequences,
        vec![
            Consequence::FrameshiftVariant,
            Consequence::StartRetainedVariant
        ]
    );
    for bases in [&b"N"[..], &b"NNN"[..]] {
        let at_stop = annotate_insertion("chr1", 4498, bases, &tx, &index, &fasta).unwrap();
        assert_eq!(
            at_stop.consequences,
            vec![Consequence::StopRetainedVariant],
            "{bases:?}"
        );
        assert_eq!(at_stop.ptc, PtcStatus::NotApplicable);
    }
    let store = crate::TranscriptStore::from_transcripts(vec![tx.clone()]);
    let delins = annotate("chr1", 4497, b"TA", b"N", &store, &fasta).unwrap();
    assert_eq!(
        delins[0].consequences,
        vec![Consequence::StopRetainedVariant]
    );
}

// -- Frameshift codon window, peptide and terms (VEP 116) -------------------

#[test]
fn frameshift_insertion_shows_the_codon_window() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    // Inside codon 2 (`CGT`), after its first base: `cgt/cAgt`, `R/QX`.
    let inside = annotate_insertion("chr1", 1504, b"A", &tx, &index, &fasta).unwrap();
    assert_eq!(inside.codons.as_deref(), Some("cgt/cAgt"));
    assert_eq!(inside.amino_acids.as_deref(), Some("R/QX"));
    // Between codons 1 and 2: no reference codon (`-/T`, `-/X`).
    let between = annotate_insertion("chr1", 1503, b"T", &tx, &index, &fasta).unwrap();
    assert_eq!(between.codons.as_deref(), Some("-/T"));
    assert_eq!(between.amino_acids.as_deref(), Some("-/X"));
    assert_eq!(between.consequences, vec![Consequence::FrameshiftVariant]);
}

#[test]
fn frameshift_reading_a_stop_in_its_window_is_also_stop_gained() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    // `cgt/cAATAgt`: CAA TAG T reads `Q*X`.
    let result = annotate_insertion("chr1", 1504, b"AATA", &tx, &index, &fasta).unwrap();
    assert_eq!(
        result.consequences,
        vec![Consequence::StopGained, Consequence::FrameshiftVariant]
    );
    assert_eq!(result.amino_acids.as_deref(), Some("R/Q*X"));
}

#[test]
fn frameshift_peptide_drops_x_after_a_lone_stop() {
    // VEP 116 `TranscriptVariationAllele::peptide`: `CE/*` (NM_017635.5).
    assert_eq!(frameshift_peptide_for_test(b"TGAC"), b"*");
    assert_eq!(frameshift_peptide_for_test(b"TGGC"), b"WX");
    assert_eq!(frameshift_peptide_for_test(b"TGATGGC"), b"*WX");
    assert_eq!(frameshift_peptide_for_test(b"T"), b"X");
}

fn frameshift_peptide_for_test(seq: &[u8]) -> Vec<u8> {
    super::indel::translate_with_partial_codon(seq, false).unwrap()
}

/// A minus-strand transcript on `chrT` whose 5'UTR spans two exons:
/// exon 1 `[150, 160)` (UTR only) and exon 2 `[80, 130)` with the CDS on
/// `[91, 127)`. Coding strand: UTR `...{last}` + `TGC`, CDS `ATGC...`.
fn two_exon_utr_minus(last_utr_base: u8) -> (TempDir, FastaReader, TranscriptModel) {
    let tmp = TempDir::new().unwrap();
    let mut chrt = vec![b'C'; 200];
    // Plus-strand bases are the complements of the coding bases, read 3'->5'.
    chrt[150] = crate::codon::complement(last_utr_base);
    chrt[127..130].copy_from_slice(b"GCA"); // coding `TGC`
    chrt[123..127].copy_from_slice(b"GCAT"); // coding `ATGC`
    let contigs: Vec<(&str, &[u8])> = vec![("chrT", &chrt)];
    let bin_path = tmp.path().join("t.bin");
    let idx_path = tmp.path().join("t.bin.idx");
    write_genome_binary(&contigs, "t", &bin_path, &idx_path).unwrap();
    let fasta = FastaReader::open(&bin_path).unwrap();
    let exons = vec![
        Exon {
            exon_number: 1,
            genomic_start: 150,
            genomic_end: 160,
        },
        Exon {
            exon_number: 2,
            genomic_start: 80,
            genomic_end: 130,
        },
    ];
    let tx = TranscriptModel::new("NM_UTR2.1", "UTR2", "chrT", Strand::Minus, exons)
        .with_cds_segments(vec![CdsSegment {
            exon_index: 1,
            genomic_start: 91,
            genomic_end: 127,
            phase: 0,
        }])
        .with_biotype(Biotype::ProteinCoding);
    (tmp, fasta, tx)
}

#[test]
fn utr5_tail_reads_across_the_utr_intron_on_the_minus_strand() {
    let (_tmp, fasta, tx) = two_exon_utr_minus(b'A');
    let index = LocateIndex::build(&tx).unwrap();
    let tail = |n| super::indel::utr5_tail(n, &tx, &index, &fasta).unwrap();
    assert_eq!(tail(2), b"GC");
    assert_eq!(tail(4), b"ATGC");
    // Exon 1 is all `G` on the coding strand except its last base.
    assert_eq!(tail(13), b"GGGGGGGGGATGC");
    assert_eq!(tail(50).len(), 13);
}

#[test]
fn deletion_shifting_the_start_into_a_two_exon_utr_is_start_retained() {
    // Deleting c.1_4 (`ATGC`) leaves `...A` + `TGC` + `...`: the same
    // sequence as deleting the UTR copy, so the start is kept (and moved).
    let (_tmp, fasta, tx) = two_exon_utr_minus(b'A');
    let index = LocateIndex::build(&tx).unwrap();
    let kept = annotate_deletion("chrT", 123, 127, b"GCAT", &tx, &index, &fasta).unwrap();
    assert!(
        kept.consequences
            .contains(&Consequence::StartRetainedVariant),
        "{kept:?}"
    );
    assert!(kept.consequences.contains(&Consequence::StartLost));
    // With a `G` before `TGC` the UTR no longer copies the start codon.
    let (_tmp, fasta, tx) = two_exon_utr_minus(b'G');
    let index = LocateIndex::build(&tx).unwrap();
    let lost = annotate_deletion("chrT", 123, 127, b"GCAT", &tx, &index, &fasta).unwrap();
    assert!(
        lost.consequences.contains(&Consequence::StartLost),
        "{lost:?}"
    );
    assert!(
        !lost
            .consequences
            .contains(&Consequence::StartRetainedVariant)
    );
}

#[test]
fn frameshift_at_the_start_codon_follows_vep_start_rules() {
    let (_tmp, fasta) = write_test_fasta();
    let tx = plus_strand_coding();
    let index = LocateIndex::build(&tx).unwrap();
    // A duplicated first base only shifts the start into the 5'UTR (`A`s):
    // VEP 116 NM_000546.6 c.1dup is `start_lost,start_retained_variant`.
    let dup = annotate_insertion("chr1", 1501, b"A", &tx, &index, &fasta).unwrap();
    assert_eq!(
        dup.consequences,
        vec![
            Consequence::FrameshiftVariant,
            Consequence::StartLost,
            Consequence::StartRetainedVariant
        ]
    );
    // Deleting the `G` leaves `ATC`: start lost, not retained.
    let deletion = annotate_deletion("chr1", 1502, 1503, b"G", &tx, &index, &fasta).unwrap();
    assert!(deletion.consequences.contains(&Consequence::StartLost));
    assert!(
        !deletion
            .consequences
            .contains(&Consequence::StartRetainedVariant)
    );
    // An insertion right after the start codon does not touch it.
    let after = annotate_insertion("chr1", 1503, b"T", &tx, &index, &fasta).unwrap();
    assert!(!after.consequences.iter().any(|c| matches!(
        c,
        Consequence::StartLost | Consequence::StartRetainedVariant
    )));
}
