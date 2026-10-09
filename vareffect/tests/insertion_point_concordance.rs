//! VEP concordance for where an insertion is classified: CDS on both
//! strands, exon/intron and UTR/CDS junctions, the canonical and
//! splice-region windows, frameshift introns, non-coding transcripts.
//!
//! The transcripts in `data/insertion_point_transcripts.jsonl` are built from
//! their MANE 1.5 and GENCODE v50 GFF3 rows (`data/tools/extract_fixture.py`).
//! The expected values in `data/insertion_point_vep116.tsv` are Ensembl VEP
//! REST release 116 output (GRCh38, fetched 2026-10-09 with
//! `data/tools/vep_insertion_rows.py`). Each inserted base matches neither
//! flanking base, so no row is a duplication or a shiftable repeat, except
//! two ClinVar duplications that pin HGVS p. after a 3' shift: a minus-strand
//! in-frame duplication splitting a codon (NM_006231.4 `c.2465_2467dup`) and
//! a junction insertion in a repeat (NM_000548.5 `c.336+1_336+2dup`).
//!
//! `#[ignore]`-gated: needs the reference genome. Run with:
//! ```bash
//! FASTA_PATH="$PWD/data/vareffect/GRCh38.bin" \
//!   cargo test -p vareffect --release --test insertion_point_concordance -- --ignored
//! ```
//!
//! # Normalization of VEP's output
//!
//! 1. `splice_donor_region_variant`, `splice_donor_5th_base_variant` and
//!    `splice_polypyrimidine_tract_variant` fold into `splice_region_variant`,
//!    as everywhere in vareffect.
//! 2. Codons and amino acids of a frameshift insertion are not compared
//!    (a pre-existing display difference on every transcript).
//! 3. Known pre-existing difference, not part of this suite's subject: an
//!    in-frame insertion inside a codon is VEP's `protein_altering_variant`
//!    and vareffect's `inframe_insertion`. Terms and impact are not compared
//!    on those rows (`INFRAME_NAMING_ROWS`).
//! 4. Known pre-existing difference, not part of this suite's subject: where
//!    a CDS ends at an exon end, VEP writes the HGVS c. of an insertion at
//!    the following junction at its mRNA point (ENST00000372299.4
//!    `c.1170_*1insA`) and vareffect at the genomic one (`c.1170_1170+1insA`).
//!    HGVS c. is not compared on those rows (`HGVS_C_JUNCTION_ROWS`).
//! 5. Known pre-existing difference, not part of this suite's subject: VEP's
//!    HGVS 3' shift carries an insertion in a repeat across an exon/intron
//!    junction and builds HGVS p. there (NM_000548.5 `c.336+1_336+2dup`, no
//!    HGVS p.); vareffect's stops at the exon edge (`c.336_336+1insGT`,
//!    `p.Gly113ValfsTer70`). HGVS c. and p. are not compared on that row
//!    (`HGVS_SHIFT_ROWS`).

use std::path::{Path, PathBuf};

use vareffect::{ConsequenceResult, FastaReader, TranscriptModel, TranscriptStore, VarEffect};

/// Rows hit by the pre-existing in-frame naming difference (normalization 3).
const INFRAME_NAMING_ROWS: &[&str] = &["tp53_mid_if", "brca2_mid_if"];

/// Rows hit by the pre-existing HGVS c. junction difference (normalization 4).
const HGVS_C_JUNCTION_ROWS: &[&str] = &["cds_exon_end_plus", "cds_exon_end_minus"];

/// Rows hit by the pre-existing HGVS shift limit at a junction
/// (normalization 5).
const HGVS_SHIFT_ROWS: &[&str] = &["tsc2_junction_dup_shifted"];

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn open_fasta() -> FastaReader {
    let path = std::env::var("FASTA_PATH")
        .expect("FASTA_PATH must point to the GRCh38 genome binary (data/vareffect/GRCh38.bin)");
    FastaReader::open(Path::new(&path)).unwrap_or_else(|e| panic!("open {path}: {e}"))
}

fn load_engine() -> VarEffect {
    let text = std::fs::read_to_string(data_dir().join("insertion_point_transcripts.jsonl"))
        .expect("read transcript fixture");
    let models: Vec<TranscriptModel> = text
        .lines()
        .map(|line| serde_json::from_str(line).expect("transcript fixture line"))
        .collect();
    let store = TranscriptStore::try_from_transcripts(models).expect("fixture store");
    VarEffect::new(store, open_fasta())
}

/// One recorded VEP row.
struct VepRow {
    label: String,
    chrom: String,
    /// 1-based position of the base before the insertion (the VCF anchor).
    left_flank: u64,
    inserted: String,
    transcript: String,
    consequences: Vec<String>,
    impact: String,
    hgvsc: Option<String>,
    hgvsp: Option<String>,
    cdna_start: Option<u32>,
    cdna_end: Option<u32>,
    cds_start: Option<u32>,
    cds_end: Option<u32>,
    protein_start: Option<u32>,
    codons: Option<String>,
    amino_acids: Option<String>,
    exon: Option<String>,
    intron: Option<String>,
}

fn field(value: &str) -> Option<String> {
    (value != "None").then(|| value.to_string())
}

fn number(value: &str) -> Option<u32> {
    field(value).map(|v| v.parse().expect("numeric VEP field"))
}

fn load_rows() -> Vec<VepRow> {
    let text = std::fs::read_to_string(data_dir().join("insertion_point_vep116.tsv"))
        .expect("read VEP fixture");
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("label\t"))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 18, "fixture line: {line}");
            VepRow {
                label: f[0].to_string(),
                chrom: f[1].to_string(),
                left_flank: f[2].parse().expect("left flank"),
                inserted: f[3].to_string(),
                transcript: f[4].to_string(),
                consequences: f[5].split(',').map(str::to_string).collect(),
                impact: f[6].to_string(),
                hgvsc: field(f[7]),
                hgvsp: field(f[8]),
                cdna_start: number(f[9]),
                cdna_end: number(f[10]),
                cds_start: number(f[11]),
                cds_end: number(f[12]),
                protein_start: number(f[13]),
                codons: field(f[14]),
                amino_acids: field(f[15]),
                exon: field(f[16]),
                intron: field(f[17]),
            }
        })
        .collect()
}

/// VEP's terms with normalization 1 applied.
fn expected_terms(row: &VepRow) -> Vec<String> {
    let mut terms: Vec<String> = row
        .consequences
        .iter()
        .map(|term| match term.as_str() {
            "splice_donor_region_variant"
            | "splice_donor_5th_base_variant"
            | "splice_polypyrimidine_tract_variant" => "splice_region_variant".to_string(),
            other => other.to_string(),
        })
        .collect();
    terms.sort();
    terms.dedup();
    terms
}

fn strip_accession(value: Option<&str>) -> Option<String> {
    value.map(|v| v.split_once(':').map_or(v, |(_, rest)| rest).to_string())
}

/// The accession without its version (VEP reports Ensembl IDs unversioned).
fn unversioned(accession: &str) -> &str {
    accession.split_once('.').map_or(accession, |(id, _)| id)
}

fn annotate_row(engine: &VarEffect, fasta: &FastaReader, row: &VepRow) -> ConsequenceResult {
    let anchor_pos = row.left_flank - 1;
    let anchor = fasta
        .fetch_base(&row.chrom, anchor_pos)
        .expect("anchor base");
    let mut alt = vec![anchor];
    alt.extend_from_slice(row.inserted.as_bytes());
    engine
        .annotate(&row.chrom, anchor_pos, &[anchor], &alt)
        .unwrap_or_else(|e| panic!("{} {}:{}: {e}", row.label, row.chrom, row.left_flank))
        .into_iter()
        .find(|r| unversioned(&r.transcript) == unversioned(&row.transcript))
        .unwrap_or_else(|| panic!("{}: no result for {}", row.label, row.transcript))
}

#[test]
#[ignore = "requires FASTA_PATH (GRCh38 genome binary)"]
fn insertions_are_classified_at_their_insertion_point_as_vep_116() {
    let engine = load_engine();
    let fasta = open_fasta();
    let rows = load_rows();
    assert_eq!(rows.len(), 46, "fixture rows");

    let mut failures = Vec::new();
    for row in &rows {
        let got = annotate_row(&engine, &fasta, row);
        let mut check = |what: &str, got: String, want: String| {
            if got != want {
                failures.push(format!("{} {what}: {got} != VEP {want}", row.label));
            }
        };
        if !INFRAME_NAMING_ROWS.contains(&row.label.as_str()) {
            let mut got_terms: Vec<String> = got
                .consequences
                .iter()
                .map(|c| c.as_str().to_string())
                .collect();
            got_terms.sort();
            check(
                "consequences",
                format!("{got_terms:?}"),
                format!("{:?}", expected_terms(row)),
            );
            check("impact", got.impact.to_string(), row.impact.clone());
        }
        // With its accession, so the fixture model's version must be VEP's.
        let hgvs_shift = HGVS_SHIFT_ROWS.contains(&row.label.as_str());
        if !hgvs_shift && !HGVS_C_JUNCTION_ROWS.contains(&row.label.as_str()) {
            check(
                "hgvs_c",
                format!("{:?}", got.hgvs_c),
                format!("{:?}", row.hgvsc),
            );
        }
        if !hgvs_shift {
            check(
                "hgvs_p",
                format!("{:?}", strip_accession(got.hgvs_p.as_deref())),
                format!("{:?}", strip_accession(row.hgvsp.as_deref())),
            );
        }
        check(
            "cdna",
            format!("{:?}", (got.cdna_position, got.cdna_position_end)),
            format!("{:?}", (row.cdna_start, row.cdna_end)),
        );
        check(
            "cds",
            format!("{:?}", (got.cds_position, got.cds_position_end)),
            format!("{:?}", (row.cds_start, row.cds_end)),
        );
        check(
            "protein_start",
            format!("{:?}", got.protein_start),
            format!("{:?}", row.protein_start),
        );
        let frameshift = row.consequences.iter().any(|t| t == "frameshift_variant");
        if !frameshift {
            check(
                "codons",
                format!("{:?}", got.codons),
                format!("{:?}", row.codons),
            );
            check(
                "amino_acids",
                format!("{:?}", got.amino_acids),
                format!("{:?}", row.amino_acids),
            );
        }
        check("exon", format!("{:?}", got.exon), format!("{:?}", row.exon));
        check(
            "intron",
            format!("{:?}", got.intron),
            format!("{:?}", row.intron),
        );
    }
    assert!(
        failures.is_empty(),
        "{} differences from VEP 116 over {} rows:\n{}",
        failures.len(),
        rows.len(),
        failures.join("\n")
    );
}
