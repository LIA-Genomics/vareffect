//! CSQ INFO field formatting for VEP-compatible VCF output.
//!
//! Produces the pipe-delimited `CSQ` INFO value that VEP emits in `--vcf`
//! mode. Each [`vareffect::ConsequenceResult`] maps to one CSQ entry; entries
//! for multiple transcripts or ALT alleles are comma-separated within the
//! single `CSQ` INFO field.
//!
//! The 19 fields match VEP's standard `--vcf` output with `--fields`:
//!
//! ```text
//! Allele|Consequence|IMPACT|SYMBOL|Feature|Feature_type|BIOTYPE|EXON|INTRON|
//! HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|
//! Codons|STRAND|MANE_SELECT|MANE_PLUS_CLINICAL
//! ```

use std::fmt::Write;

use vareffect::consequence::ConsequenceResult;
use vareffect::types::Strand;

/// VCF `##INFO` header line declaring the CSQ field and its sub-fields.
///
/// Insert this line before the `#CHROM` line in the output VCF.
pub const CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,\
    Description=\"Consequence annotations from vareffect. \
    Format: Allele|Consequence|IMPACT|SYMBOL|Feature|Feature_type|BIOTYPE|\
    EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|\
    Amino_acids|Codons|STRAND|MANE_SELECT|MANE_PLUS_CLINICAL\">";

/// Trim shared prefix and suffix from VCF REF/ALT alleles (VEP `--minimal`
/// convention).
///
/// Returns the trimmed ALT as a `String`. Empty trimmed ALT (pure deletion)
/// is rendered as `"-"` per VEP CSQ convention.
///
/// Logic mirrors `consequence::helpers::trim_alleles` (which is
/// `pub(super)` and not accessible from outside the `consequence` module).
pub fn trimmed_csq_allele(ref_allele: &[u8], alt_allele: &[u8]) -> String {
    let mut r = ref_allele;
    let mut a = alt_allele;

    // Left-prefix strip (matches VEP's `--minimal` behavior).
    while !r.is_empty() && !a.is_empty() && r[0] == a[0] {
        r = &r[1..];
        a = &a[1..];
    }

    // Right-suffix strip.
    while !r.is_empty() && !a.is_empty() && r[r.len() - 1] == a[a.len() - 1] {
        r = &r[..r.len() - 1];
        a = &a[..a.len() - 1];
    }

    if a.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(a).into_owned()
    }
}

/// Format a single `ConsequenceResult` as a pipe-delimited CSQ entry.
///
/// # Arguments
///
/// * `trimmed_alt` — the minimized ALT allele for the CSQ `Allele` field
///   (produced by [`trimmed_csq_allele`]).
/// * `r` — the per-transcript annotation result.
///
/// # Returns
///
/// A `String` with 19 pipe-separated fields matching the VEP CSQ layout.
pub fn format_csq(trimmed_alt: &str, r: &ConsequenceResult) -> String {
    // Pre-allocate; typical CSQ entry is 100-180 bytes.
    let mut out = String::with_capacity(200);

    // 1. Allele
    out.push_str(trimmed_alt);
    out.push('|');

    // 2. Consequence — SO terms joined by "&"
    for (i, c) in r.consequences.iter().enumerate() {
        if i > 0 {
            out.push('&');
        }
        out.push_str(c.as_str());
    }
    out.push('|');

    // 3. IMPACT — Display impl already gives uppercase
    let _ = write!(out, "{}", r.impact);
    out.push('|');

    // 4. SYMBOL
    out.push_str(&r.gene_symbol);
    out.push('|');

    // 5. Feature (transcript accession)
    out.push_str(&r.transcript);
    out.push('|');

    // 6. Feature_type — always "Transcript" when transcript is populated
    if !r.transcript.is_empty() {
        out.push_str("Transcript");
    }
    out.push('|');

    // 7. BIOTYPE
    out.push_str(r.biotype.as_str());
    out.push('|');

    // 8. EXON
    push_opt_str(&mut out, r.exon.as_deref());
    out.push('|');

    // 9. INTRON
    push_opt_str(&mut out, r.intron.as_deref());
    out.push('|');

    // 10. HGVSc — already carries the full accession prefix
    //     (e.g. "NM_000546.6:c.742C>T"). Use directly.
    push_opt_str(&mut out, r.hgvs_c.as_deref());
    out.push('|');

    // 11. HGVSp — bare notation; prepend protein accession if both present.
    //     hgvs_p is stored as "p.Arg248Trp", protein_accession as "NP_000537.3".
    //     VEP emits "NP_000537.3:p.Arg248Trp".
    if let (Some(p), Some(acc)) = (r.hgvs_p.as_deref(), r.protein_accession.as_deref()) {
        out.push_str(acc);
        out.push(':');
        out.push_str(p);
    }
    out.push('|');

    // 12. cDNA_position
    push_position_range(&mut out, r.cdna_position, r.cdna_position_end);
    out.push('|');

    // 13. CDS_position
    push_position_range(&mut out, r.cds_position, r.cds_position_end);
    out.push('|');

    // 14. Protein_position
    push_position_range(&mut out, r.protein_start, r.protein_end);
    out.push('|');

    // 15. Amino_acids
    push_opt_str(&mut out, r.amino_acids.as_deref());
    out.push('|');

    // 16. Codons
    push_opt_str(&mut out, r.codons.as_deref());
    out.push('|');

    // 17. STRAND — "1" for plus, "-1" for minus (Strand is #[non_exhaustive])
    out.push_str(match r.strand {
        Strand::Plus => "1",
        Strand::Minus => "-1",
        _ => "",
    });
    out.push('|');

    // 18. MANE_SELECT — transcript ID if is_mane_select, else empty
    if r.is_mane_select {
        out.push_str(&r.transcript);
    }
    out.push('|');

    // 19. MANE_PLUS_CLINICAL — transcript ID if is_mane_plus_clinical, else empty
    if r.is_mane_plus_clinical {
        out.push_str(&r.transcript);
    }

    out
}

/// Format CSQ entries for all transcripts of a single ALT allele.
///
/// Returns the comma-separated CSQ value (one entry per transcript).
/// Skips intergenic results (empty transcript) since VEP does not emit
/// CSQ entries for intergenic variants.
///
/// # Arguments
///
/// * `trimmed_alt` — the minimized ALT allele (from [`trimmed_csq_allele`]).
/// * `results` — per-transcript annotations from [`vareffect::VarEffect::annotate`].
pub fn format_variant_csq(trimmed_alt: &str, results: &[ConsequenceResult]) -> String {
    let mut out = String::with_capacity(results.len() * 200);
    let mut first = true;

    for r in results {
        // Skip intergenic results — they have empty transcript/gene and
        // VEP's --vcf mode does not emit CSQ entries for them.
        if r.transcript.is_empty() {
            continue;
        }
        if !first {
            out.push(',');
        }
        out.push_str(&format_csq(trimmed_alt, r));
        first = false;
    }

    out
}

/// Push an optional string value, or nothing if `None`.
fn push_opt_str(out: &mut String, value: Option<&str>) {
    if let Some(v) = value {
        out.push_str(v);
    }
}

/// Push a position or position range (1-based, VEP convention).
///
/// - `Some(n), Some(n)` (start == end) → `"N"`
/// - `Some(n), Some(m)` (start != end) → `"N-M"`
/// - `Some(n), None` → `"N"`
/// - `None, _` → empty
fn push_position_range(out: &mut String, start: Option<u32>, end: Option<u32>) {
    if let Some(s) = start {
        let _ = write!(out, "{s}");
        if let Some(e) = end
            && e != s
        {
            let _ = write!(out, "-{e}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vareffect::consequence::{Consequence, Impact};
    use vareffect::types::Biotype;

    /// Minimal fixture with sensible defaults. Tests override relevant fields.
    fn make_result(transcript: &str, gene: &str) -> ConsequenceResult {
        ConsequenceResult {
            transcript: transcript.to_string(),
            gene_symbol: gene.to_string(),
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
            strand: Strand::Plus,
            biotype: Biotype::ProteinCoding,
            is_mane_select: false,
            is_mane_plus_clinical: false,
            is_refseq_select: false,
            hgvs_c: None,
            hgvs_p: None,
            predicts_nmd: false,
        }
    }

    #[test]
    fn trimmed_csq_allele_snv() {
        assert_eq!(trimmed_csq_allele(b"C", b"T"), "T");
    }

    #[test]
    fn trimmed_csq_allele_deletion() {
        assert_eq!(trimmed_csq_allele(b"CA", b"C"), "-");
    }

    #[test]
    fn trimmed_csq_allele_insertion() {
        assert_eq!(trimmed_csq_allele(b"C", b"CA"), "A");
    }

    #[test]
    fn trimmed_csq_allele_complex() {
        assert_eq!(trimmed_csq_allele(b"AT", b"GC"), "GC");
    }

    #[test]
    fn trimmed_csq_allele_multi_base_insertion() {
        assert_eq!(trimmed_csq_allele(b"T", b"TAA"), "AA");
    }

    #[test]
    fn trimmed_csq_allele_shared_suffix() {
        // REF=ACG, ALT=AG → trimmed ref=C, trimmed alt=empty → "-"
        assert_eq!(trimmed_csq_allele(b"ACG", b"AG"), "-");
    }

    #[test]
    fn format_csq_missense_full() {
        let r = ConsequenceResult {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            protein_start: Some(248),
            protein_end: Some(248),
            codons: Some("cGg/cTg".to_string()),
            amino_acids: Some("R/W".to_string()),
            exon: Some("7/11".to_string()),
            cds_position: Some(742),
            cds_position_end: Some(742),
            cdna_position: Some(884),
            cdna_position_end: Some(884),
            strand: Strand::Minus,
            is_mane_select: true,
            hgvs_c: Some("NM_000546.6:c.742C>T".to_string()),
            hgvs_p: Some("p.Arg248Trp".to_string()),
            protein_accession: Some("NP_000537.3".to_string()),
            ..make_result("NM_000546.6", "TP53")
        };

        let csq = format_csq("T", &r);
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields.len(), 19);
        assert_eq!(fields[0], "T"); // Allele
        assert_eq!(fields[1], "missense_variant"); // Consequence
        assert_eq!(fields[2], "MODERATE"); // IMPACT
        assert_eq!(fields[3], "TP53"); // SYMBOL
        assert_eq!(fields[4], "NM_000546.6"); // Feature
        assert_eq!(fields[5], "Transcript"); // Feature_type
        assert_eq!(fields[6], "protein_coding"); // BIOTYPE
        assert_eq!(fields[7], "7/11"); // EXON
        assert_eq!(fields[8], ""); // INTRON
        assert_eq!(fields[9], "NM_000546.6:c.742C>T"); // HGVSc
        assert_eq!(fields[10], "NP_000537.3:p.Arg248Trp"); // HGVSp
        assert_eq!(fields[11], "884"); // cDNA_position
        assert_eq!(fields[12], "742"); // CDS_position
        assert_eq!(fields[13], "248"); // Protein_position
        assert_eq!(fields[14], "R/W"); // Amino_acids
        assert_eq!(fields[15], "cGg/cTg"); // Codons
        assert_eq!(fields[16], "-1"); // STRAND
        assert_eq!(fields[17], "NM_000546.6"); // MANE_SELECT
        assert_eq!(fields[18], ""); // MANE_PLUS_CLINICAL
    }

    #[test]
    fn format_csq_intronic_many_empty() {
        let r = ConsequenceResult {
            consequences: vec![Consequence::IntronVariant],
            intron: Some("7/10".to_string()),
            hgvs_c: Some("NM_000546.6:c.782+131C>T".to_string()),
            ..make_result("NM_000546.6", "TP53")
        };

        let csq = format_csq("A", &r);
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[0], "A");
        assert_eq!(fields[1], "intron_variant");
        assert_eq!(fields[2], "MODIFIER");
        assert_eq!(fields[7], ""); // EXON empty
        assert_eq!(fields[8], "7/10"); // INTRON populated
        assert_eq!(fields[9], "NM_000546.6:c.782+131C>T"); // HGVSc
        assert_eq!(fields[10], ""); // HGVSp empty
        assert_eq!(fields[11], ""); // cDNA empty
        assert_eq!(fields[12], ""); // CDS empty
        assert_eq!(fields[13], ""); // Protein empty
        assert_eq!(fields[14], ""); // Amino_acids empty
        assert_eq!(fields[15], ""); // Codons empty
    }

    #[test]
    fn format_csq_multiple_consequences() {
        let r = ConsequenceResult {
            consequences: vec![
                Consequence::SpliceRegionVariant,
                Consequence::MissenseVariant,
            ],
            impact: Impact::Moderate,
            ..make_result("NM_000546.6", "TP53")
        };

        let csq = format_csq("T", &r);
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[1], "splice_region_variant&missense_variant");
    }

    #[test]
    fn format_csq_hgvsp_omitted_without_protein_accession() {
        let r = ConsequenceResult {
            hgvs_p: Some("p.?".to_string()),
            protein_accession: None,
            ..make_result("NR_046018.2", "DDX11L1")
        };

        let csq = format_csq("A", &r);
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[10], ""); // HGVSp empty when protein_accession is None
    }

    #[test]
    fn format_csq_position_range_for_indel() {
        let r = ConsequenceResult {
            consequences: vec![Consequence::InframeDeletion],
            impact: Impact::Moderate,
            cds_position: Some(742),
            cds_position_end: Some(744),
            cdna_position: Some(884),
            cdna_position_end: Some(886),
            protein_start: Some(248),
            protein_end: Some(248),
            ..make_result("NM_000546.6", "TP53")
        };

        let csq = format_csq("-", &r);
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[11], "884-886"); // cDNA range
        assert_eq!(fields[12], "742-744"); // CDS range
        assert_eq!(fields[13], "248"); // Protein single (start == end)
    }

    #[test]
    fn format_variant_csq_multiple_transcripts() {
        let r1 = ConsequenceResult {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            ..make_result("NM_000546.6", "TP53")
        };
        let r2 = ConsequenceResult {
            consequences: vec![Consequence::SynonymousVariant],
            impact: Impact::Low,
            ..make_result("NM_001126112.3", "TP53")
        };

        let csq = format_variant_csq("T", &[r1, r2]);
        let entries: Vec<&str> = csq.split(',').collect();
        assert_eq!(entries.len(), 2);
        assert!(entries[0].starts_with("T|missense_variant|"));
        assert!(entries[1].starts_with("T|synonymous_variant|"));
    }

    #[test]
    fn format_variant_csq_skips_intergenic() {
        let intergenic = ConsequenceResult {
            consequences: vec![Consequence::IntergenicVariant],
            ..make_result("", "")
        };
        let coding = ConsequenceResult {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            ..make_result("NM_000546.6", "TP53")
        };

        let csq = format_variant_csq("T", &[intergenic, coding]);
        let entries: Vec<&str> = csq.split(',').collect();
        // Only the coding result should appear.
        assert_eq!(entries.len(), 1);
        assert!(entries[0].starts_with("T|missense_variant|"));
    }

    #[test]
    fn format_variant_csq_all_intergenic_returns_empty() {
        let intergenic = ConsequenceResult {
            consequences: vec![Consequence::IntergenicVariant],
            ..make_result("", "")
        };

        let csq = format_variant_csq("T", &[intergenic]);
        assert!(csq.is_empty());
    }
}
