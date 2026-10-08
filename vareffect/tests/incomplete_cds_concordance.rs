//! VEP concordance for incomplete CDS framing and frameshift introns.
//!
//! The transcripts in `data/incomplete_cds_transcripts.jsonl` are built from
//! their GENCODE v50 / RefSeq RS_2025_08 / MANE 1.5 GFF3 rows. The expected
//! values in `data/incomplete_cds_vep116.tsv` are Ensembl VEP REST release
//! 116 output (GRCh38, fetched 2026-10-07; how both files are made is in
//! `data/tools/README.md`). They cover:
//!
//! - a 5'-incomplete CDS (GFF3 phase 1 and 2, both strands): VEP's padded
//!   `c.` numbering for SNVs, unpadded numbering for indels, MNVs and
//!   intronic anchors, unpadded cDNA, the `N`-padded first codon;
//! - a 3' partial codon: a cds_end_NF transcript, and stops completed only by
//!   polyadenylation (MT-ND4 ENST00000361381.2, TMEM247 NM_001424184.1), for
//!   SNVs, deletions and insertions (of one to four bases at the partial
//!   codon's 5' edge, and after the last CDS base with and without a 3'UTR);
//! - frameshift introns (exon gaps of 1, 2 and 12 bases) next to a normal
//!   14-base intron: variants inside the gap, deletions spanning it, indels
//!   at its junctions, and the splice windows of a one-base gap spilling
//!   onto the flanking exon bases.
//!
//! `#[ignore]`-gated: needs the reference genome. Run with:
//! ```bash
//! FASTA_PATH="$PWD/data/vareffect/GRCh38.bin" \
//!   cargo test -p vareffect --release --test incomplete_cds_concordance -- --ignored
//! ```
//!
//! # Normalization of VEP's output (documented divergences)
//!
//! 1. `splice_donor_region_variant`, `splice_donor_5th_base_variant` and
//!    `splice_polypyrimidine_tract_variant` fold into `splice_region_variant`,
//!    as everywhere in vareffect.
//! 2. VEP's `p.Ter…` on a partial codon (amino acids with `X`: `p.Ter1=`,
//!    `p.Ter309del`, `p.TerLeu1=`) is `p.?` in vareffect: the residue is
//!    unknown, not a stop. So is a stop an insertion creates at the partial
//!    codon (`stop_gained` with `incomplete_terminal_codon_variant`), where
//!    VEP writes `p.Ter220=`, `p.Ter220delinsTer` or nothing.
//! 3. The 14-base intron rows check only that the intron is *not* treated as
//!    a frameshift intron (a splice or intron term is present, HGVS equal).
//!    VEP's overlapping donor-region, polypyrimidine and
//!    `non_coding_transcript_variant` terms on so short a normal intron are a
//!    separate, pre-existing difference.
//! 4. `VarEffect::annotate` reports only transcripts the variant overlaps
//!    (for an insertion, those its flanking bases touch), so a VEP
//!    `downstream_gene_variant` row has either no result or the same term.
//! 5. VEP describes a reverse-complement MNV as `inv` (`c.1_2inv`);
//!    vareffect writes `delins` (`c.1_2delinsCT`).
//! 6. VEP omits or clips HGVS for deletions that reach the end of a CDS with
//!    no genomic stop (NM_001424184.1 `c.658_*2del`, `c.656_659del`,
//!    `c.656_*2del` written as `c.656_659del`; ENST00000518351.2's whole
//!    partial codon). HGVS is not compared on those rows (`HGVS_NOT_COMPARED`).
//! 7. Codons and amino acids of a frameshift *insertion* are not compared:
//!    vareffect prints `-/A` and `A/X` where VEP prints the codon context
//!    (`gcg/gcAg`, `A/AX`) or `-/X`. This predates incomplete-CDS support and
//!    applies to every transcript.
//! 8. For a deletion spanning a frameshift intron VEP uppercases as many
//!    codon bases as the genomic footprint is long, intron base included
//!    (`tCCGAg/tag` for three coding bases); codons are compared ignoring
//!    case there.
//! 9. Maintainer decisions (`DECIDED`): a deletion that removes the end of a
//!    CDS with no genomic stop gets no stop term (VEP's
//!    `frameshift_variant,stop_lost` for `c.656_659del` and
//!    `stop_retained_variant,3_prime_UTR_variant` for `c.656_*2del`), and
//!    for an SNV, VEP's splice-region windows spilling further than three
//!    bases from a short gap are not reproduced (`c.232T>A`); for an indel,
//!    vareffect's acceptor window reaches up to 16 bases before a short gap
//!    and gives `splice_region_variant` where VEP gives none (deletions 9 and
//!    14 bases before it). The donor window spills as VEP's does
//!    (`c.233_235del`).
//! 10. Known pre-existing defect, not part of this suite's subject: on the
//!     minus strand a frameshift insertion is placed one coding base early
//!     (TP53 NM_000546.6 c.523_524insT: vareffect CDS 522 and
//!     p.Arg175SerfsTer6, VEP 116 CDS 523 and p.Arg175LeufsTer6). HGVS p. and
//!     CDS/cDNA positions are not compared for minus-strand frameshift
//!     insertions (`MINUS_STRAND_INSERTION_ROWS`).
//! 11. Known pre-existing defect, not part of this suite's subject: a
//!     non-duplication insertion that HGVS shifts 3' keeps its inserted
//!     bases unrotated (NM_001424184.1 insTGA after c.657: vareffect
//!     `c.658_659insTGA`, VEP 116 `c.658_659insGAT`). HGVS c. is not compared
//!     on those rows (`INSERTION_ROTATION_ROWS`).

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use vareffect::{ConsequenceResult, FastaReader, TranscriptModel, TranscriptStore, VarEffect};

/// Rows whose HGVS VEP omits or clips (normalization 6).
const HGVS_NOT_COMPARED: &[&str] = &[
    "tmem247_del_into_utr",
    "tmem247_del_c656_659",
    "tmem247_del_c656_utr",
    "minus_partial3_del2",
];

/// Rows hit by the pre-existing unrotated HGVS insertion (normalization 11).
const INSERTION_ROTATION_ROWS: &[&str] = &[
    "partial3_ins_5p_edge_tga",
    "partial3_ins_5p_edge_tgag",
    "partial3_ins_5p_edge_stop_then_codon",
];

/// Rows hit by the pre-existing minus-strand insertion offset (normalization
/// 10).
const MINUS_STRAND_INSERTION_ROWS: &[&str] = &["ins_2bp_gap_edge_5p", "ins_2bp_gap_edge_3p"];

/// Rows where vareffect deliberately differs from VEP (normalization 9):
/// label, expected terms, expected impact.
const DECIDED: &[(&str, &[&str], &str)] = &[
    (
        "tmem247_del_c656_659",
        &[
            "coding_sequence_variant",
            "incomplete_terminal_codon_variant",
        ],
        "LOW",
    ),
    ("region_spill_c232", &["missense_variant"], "MODERATE"),
    (
        "gap1_acceptor_window_del_9",
        &["inframe_deletion", "splice_region_variant"],
        "MODERATE",
    ),
    (
        "gap1_acceptor_window_del_14",
        &["inframe_deletion", "splice_region_variant"],
        "MODERATE",
    ),
    (
        "gap2_acceptor_window_del_9",
        &["inframe_deletion", "splice_region_variant"],
        "MODERATE",
    ),
    (
        "tmem247_del_c656_utr",
        &[
            "3_prime_UTR_variant",
            "coding_sequence_variant",
            "incomplete_terminal_codon_variant",
        ],
        "LOW",
    ),
];

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn open_fasta() -> FastaReader {
    let path = std::env::var("FASTA_PATH")
        .expect("FASTA_PATH must point to the GRCh38 genome binary (data/vareffect/GRCh38.bin)");
    FastaReader::open(Path::new(&path)).unwrap_or_else(|e| panic!("open {path}: {e}"))
}

fn load_engine() -> VarEffect {
    let text = std::fs::read_to_string(data_dir().join("incomplete_cds_transcripts.jsonl"))
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
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    transcript: String,
    consequences: Vec<String>,
    impact: Option<String>,
    hgvsc: Option<String>,
    hgvsp: Option<String>,
    cdna_start: Option<u32>,
    cds_start: Option<u32>,
    protein_start: Option<u32>,
    codons: Option<String>,
    amino_acids: Option<String>,
    /// `None` when the row was queried without `numbers=1`.
    exon: Option<Option<String>>,
    /// `None` when the row was queried without `numbers=1`.
    intron: Option<Option<String>>,
}

fn field(value: &str) -> Option<String> {
    (value != "-").then(|| value.to_string())
}

fn load_rows() -> Vec<VepRow> {
    let text = std::fs::read_to_string(data_dir().join("incomplete_cds_vep116.tsv"))
        .expect("read VEP fixture");
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("label\t"))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 17, "malformed fixture row: {line}");
            let number = |v: &str| field(v).map(|s| s.parse::<u32>().expect("number"));
            let numbered = |v: &str| (v != "?").then(|| field(v));
            VepRow {
                label: f[0].to_string(),
                chrom: f[1].to_string(),
                pos: f[2].parse().expect("pos"),
                ref_allele: f[3].to_string(),
                alt_allele: f[4].to_string(),
                transcript: f[5].to_string(),
                consequences: f[6].split(',').map(str::to_string).collect(),
                impact: field(f[7]),
                hgvsc: field(f[8]),
                hgvsp: field(f[9]),
                cdna_start: number(f[10]),
                cds_start: number(f[11]),
                protein_start: number(f[12]),
                codons: field(f[13]),
                amino_acids: field(f[14]),
                exon: numbered(f[15]),
                intron: numbered(f[16]),
            }
        })
        .collect()
}

/// VEP's consequence terms with normalizations 1 and 9 applied.
fn expected_terms(row: &VepRow) -> BTreeSet<String> {
    if let Some((_, terms, _)) = DECIDED.iter().find(|(label, _, _)| *label == row.label) {
        return terms.iter().map(|t| t.to_string()).collect();
    }
    row.consequences
        .iter()
        .map(|term| match term.as_str() {
            "splice_donor_region_variant"
            | "splice_donor_5th_base_variant"
            | "splice_polypyrimidine_tract_variant" => "splice_region_variant".to_string(),
            other => other.to_string(),
        })
        .collect()
}

/// VEP's HGVS p. with normalization 2 applied, without the accession.
fn expected_hgvsp(row: &VepRow) -> Option<String> {
    let has_term = |term: &str| row.consequences.iter().any(|t| t == term);
    if has_term("stop_gained") && has_term("incomplete_terminal_codon_variant") {
        return Some("p.?".to_string());
    }
    let p = row.hgvsp.as_deref()?.split_once(':').map(|(_, p)| p)?;
    let unknown_residue = row
        .amino_acids
        .as_deref()
        .is_some_and(|aa| aa.contains('X'));
    if unknown_residue && p.starts_with("p.Ter") {
        return Some("p.?".to_string());
    }
    Some(p.to_string())
}

/// VEP's HGVS c. with normalization 5 applied.
fn expected_hgvsc(row: &VepRow) -> Option<String> {
    let hgvsc = row.hgvsc.as_deref()?;
    match hgvsc.strip_suffix("inv") {
        Some(range) => Some(format!("{range}delins{}", row.alt_allele)),
        None => Some(hgvsc.to_string()),
    }
}

fn strip_accession(value: Option<&str>) -> Option<String> {
    value.map(|v| v.split_once(':').map_or(v, |(_, rest)| rest).to_string())
}

fn annotate_row(
    engine: &VarEffect,
    fasta: &FastaReader,
    row: &VepRow,
) -> Option<ConsequenceResult> {
    // VEP's allele strings drop the anchor base of an indel; VCF style needs
    // it. A deletion's region starts at its first deleted base and an
    // insertion's (`start = end + 1`) at the base after the insertion point,
    // so in both cases the anchor is the base before `pos` (1-based).
    let (pos0, ref_allele, alt_allele) = if row.ref_allele == "-" || row.alt_allele == "-" {
        let anchor_pos = row.pos - 2;
        let anchor = fasta
            .fetch_base(&row.chrom, anchor_pos)
            .expect("anchor base");
        let with_anchor = |a: &str| {
            let mut v = vec![anchor];
            if a != "-" {
                v.extend_from_slice(a.as_bytes());
            }
            v
        };
        (
            anchor_pos,
            with_anchor(&row.ref_allele),
            with_anchor(&row.alt_allele),
        )
    } else {
        (
            row.pos - 1,
            row.ref_allele.as_bytes().to_vec(),
            row.alt_allele.as_bytes().to_vec(),
        )
    };
    engine
        .annotate(&row.chrom, pos0, &ref_allele, &alt_allele)
        .unwrap_or_else(|e| panic!("{} {}:{}: {e}", row.label, row.chrom, row.pos))
        .into_iter()
        .find(|r| r.transcript == row.transcript)
}

#[test]
#[ignore = "requires FASTA_PATH (GRCh38 genome binary)"]
fn incomplete_cds_and_frameshift_introns_match_vep_116() {
    let engine = load_engine();
    let fasta = open_fasta();
    let rows = load_rows();
    assert!(rows.len() > 150, "fixture has {} rows", rows.len());

    let mut failures = Vec::new();
    for row in &rows {
        let id = format!(
            "{} {}:{} {}>{} {}",
            row.label, row.chrom, row.pos, row.ref_allele, row.alt_allele, row.transcript
        );
        let flanking_only = row
            .consequences
            .iter()
            .all(|t| t == "upstream_gene_variant" || t == "downstream_gene_variant");
        let Some(got) = annotate_row(&engine, &fasta, row) else {
            if !flanking_only {
                failures.push(format!("{id}: no result"));
            }
            continue;
        };
        if flanking_only {
            // An insertion next to a transcript end is queried with its
            // flanking bases, so the transcript is reported; the term must
            // then be VEP's.
            let got_terms: BTreeSet<String> = got
                .consequences
                .iter()
                .map(|c| c.as_str().to_string())
                .collect();
            if got_terms != expected_terms(row) {
                failures.push(format!("{id}: {got_terms:?} != VEP {:?}", row.consequences));
            }
            continue;
        }
        let got_terms: BTreeSet<String> = got
            .consequences
            .iter()
            .map(|c| c.as_str().to_string())
            .collect();

        if row.label.starts_with("intron_14bp") {
            let normal_intron = got_terms
                .iter()
                .any(|t| t.contains("splice") || t == "intron_variant")
                || got_terms.contains("non_coding_transcript_exon_variant");
            if !normal_intron {
                failures.push(format!(
                    "{id}: 14-base intron treated as frameshift: {got_terms:?}"
                ));
            }
            if got.hgvs_c != row.hgvsc {
                failures.push(format!(
                    "{id}: hgvs_c {:?} != VEP {:?}",
                    got.hgvs_c, row.hgvsc
                ));
            }
            continue;
        }

        let mut check = |what: &str, got: String, want: String| {
            if got != want {
                failures.push(format!("{id}: {what} {got} != VEP {want}"));
            }
        };
        check(
            "consequences",
            format!("{got_terms:?}"),
            format!("{:?}", expected_terms(row)),
        );
        let decided_impact = DECIDED
            .iter()
            .find(|(label, _, _)| *label == row.label)
            .map(|(_, _, impact)| impact.to_string());
        if let Some(impact) = decided_impact.or_else(|| row.impact.clone()) {
            check("impact", got.impact.to_string(), impact);
        }
        let hgvs_compared = !HGVS_NOT_COMPARED.contains(&row.label.as_str());
        let positions_compared = !MINUS_STRAND_INSERTION_ROWS.contains(&row.label.as_str());
        if hgvs_compared && !INSERTION_ROTATION_ROWS.contains(&row.label.as_str()) {
            check(
                "hgvs_c",
                format!("{:?}", got.hgvs_c),
                format!("{:?}", expected_hgvsc(row)),
            );
        }
        if hgvs_compared && positions_compared {
            check(
                "hgvs_p",
                format!("{:?}", strip_accession(got.hgvs_p.as_deref())),
                format!("{:?}", expected_hgvsp(row)),
            );
        }
        if positions_compared {
            check(
                "cdna_start",
                format!("{:?}", got.cdna_position),
                format!("{:?}", row.cdna_start),
            );
            check(
                "cds_start",
                format!("{:?}", got.cds_position),
                format!("{:?}", row.cds_start),
            );
        }
        check(
            "protein_start",
            format!("{:?}", got.protein_start),
            format!("{:?}", row.protein_start),
        );
        let frameshift_insertion =
            row.ref_allele == "-" && row.consequences.iter().any(|t| t == "frameshift_variant");
        if !frameshift_insertion {
            let spans_gap = row.label.starts_with("fs_del_spanning");
            let display = |codons: Option<&str>| {
                codons.map(|c| {
                    if spans_gap {
                        c.to_ascii_lowercase()
                    } else {
                        c.to_string()
                    }
                })
            };
            check(
                "codons",
                format!("{:?}", display(got.codons.as_deref())),
                format!("{:?}", display(row.codons.as_deref())),
            );
            check(
                "amino_acids",
                format!("{:?}", got.amino_acids),
                format!("{:?}", row.amino_acids),
            );
        }
        if let Some(exon) = &row.exon {
            check("exon", format!("{:?}", got.exon), format!("{exon:?}"));
        }
        if let Some(intron) = &row.intron {
            check("intron", format!("{:?}", got.intron), format!("{intron:?}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} differences from VEP 116 over {} rows:\n{}",
        failures.len(),
        rows.len(),
        failures.join("\n")
    );
}
