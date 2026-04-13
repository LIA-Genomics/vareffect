//! VEP REST API-compatible JSON serialization for vareffect annotation results.
//!
//! This module lets callers use vareffect as a drop-in replacement for
//! Ensembl VEP REST for the annotation subset vareffect produces: consequence
//! terms, HGVS c./p. notation, exon/intron numbering, codon and amino acid
//! strings, CDS/cDNA positions, and MANE membership flags.
//!
//! The module is crate-private -- the only public surface is
//! [`crate::VarEffect::annotate_to_vep_json`], which wraps
//! [`crate::VarEffect::annotate`] and calls into the serializer here.
//!
//! # Output shape
//!
//! The returned `serde_json::Value` is a single-element JSON **array** whose
//! element matches the top-level object returned by
//!
//! ```text
//! GET https://rest.ensembl.org/vep/human/region/{chrom}:{pos}:{pos}/{alt}
//!     ?refseq=1&hgvs=1
//! ```
//!
//! An array rather than a bare object is emitted because VEP REST always
//! wraps even single-variant responses in an array, and downstream parsers
//! typically deserialize into `Vec<VepAnnotation>`.
//!
//! # Coordinate and allele conventions
//!
//! - vareffect input positions are 0-based (BED convention); the emitted
//!   JSON's `start` and `end` fields are 1-based (VEP convention).
//! - Alleles are plus-strand genomic bytes (matching VEP's `allele_string`).
//! - Empty alleles (pure insertions/deletions) are rendered as `"-"` in
//!   `allele_string`, matching VEP.
//!
//! # Field mapping (top-level)
//!
//! | JSON field                | Source                                                |
//! | ---                       | ---                                                   |
//! | `seq_region_name`         | `chrom` with any leading `"chr"` stripped             |
//! | `start`                   | `pos + 1`                                             |
//! | `end`                     | `pos + max(1, ref.len())` (1-based inclusive)         |
//! | `strand`                  | Always `1` -- top-level is plus-strand genomic        |
//! | `allele_string`           | `"{REF}/{ALT}"` with `"-"` for empty alleles          |
//! | `assembly_name`           | Caller-supplied `assembly` parameter                  |
//! | `most_severe_consequence` | Minimum [`Consequence::severity_rank`] across results |
//! | `transcript_consequences` | One element per [`ConsequenceResult`]                 |
//! | `colocated_variants`      | Always `[]` -- vareffect has no dbSNP knowledge       |
//!
//! # Field mapping (per-transcript)
//!
//! | JSON field           | Source                                                    |
//! | ---                  | ---                                                       |
//! | `transcript_id`      | `r.transcript`                                            |
//! | `gene_symbol`        | `r.gene_symbol`                                           |
//! | `biotype`            | `r.biotype.as_str()`                                      |
//! | `strand`             | `Plus -> 1`, `Minus -> -1`                                |
//! | `impact`             | `r.impact.to_string()` -- uppercase via `Display`         |
//! | `consequence_terms`  | `[c.as_str() for c in r.consequences]` in input order     |
//! | `hgvsc`              | `"{r.transcript}:{r.hgvs_c}"` (omitted if `None`)         |
//! | `hgvsp`              | `"{r.protein_accession}:{r.hgvs_p}"` (omitted unless both present) |
//! | `protein_start/end`  | `r.protein_start/protein_end` widened `u32 -> u64`        |
//! | `cds_start/end`      | `r.cds_position/cds_position_end`                         |
//! | `cdna_start/end`     | `r.cdna_position/cdna_position_end`                       |
//! | `amino_acids`        | `r.amino_acids` (e.g. `"R/W"`)                            |
//! | `codons`             | `r.codons` (e.g. `"Cgg/Tgg"`)                             |
//! | `exon`               | `r.exon` (e.g. `"7/11"`)                                  |
//! | `intron`             | `r.intron`                                                |
//! | `mane_select`        | `r.transcript` when `r.is_mane_select`                    |
//! | `mane_plus_clinical` | `r.transcript` when `r.is_mane_plus_clinical`             |
//!
//! Optional fields are **omitted** when their source is `None` rather than
//! emitted as `"field": null`, matching VEP's own behavior.
//!
//! # Fields NOT emitted
//!
//! vareffect has no knowledge of these, so they are deliberately absent and
//! downstream consumers must tolerate their absence:
//!
//! - `canonical` -- vareffect distinguishes MANE Select / MANE Plus Clinical /
//!   RefSeq Select instead of a single "canonical" boolean. Downstream
//!   transcript-selection cascades should prefer the MANE fields.
//! - `swissprot`, `uniparc`, `protein_id` -- no UniProt/Ensembl protein mapping.
//! - `sift_prediction`, `polyphen_prediction`, `revel_score`, `alphamissense`,
//!   `spliceai` -- these come from separate prediction stores, not from
//!   variant effect prediction itself.
//! - `colocated_variants[*].id` -- no dbSNP rsID lookup.
//! - `input`, `id` -- VEP echoes the query string; we don't carry one.

use serde_json::{Map, Value, json};

use crate::consequence::{Consequence, ConsequenceResult};
use crate::types::Strand;

/// Serialize a slice of [`ConsequenceResult`] as a VEP REST-compatible JSON
/// array.
///
/// This is the pure serialization helper -- no FASTA or transcript store is
/// required. Tests construct synthetic [`ConsequenceResult`] values and call
/// this directly to lock the JSON shape without going through the full
/// annotation pipeline.
///
/// The returned `Value` is always a single-element JSON array. If `results`
/// is empty (e.g. intergenic variant), the array still contains one element
/// whose `transcript_consequences` is `[]` and whose `most_severe_consequence`
/// is `"intergenic_variant"`.
///
/// # Arguments
///
/// * `chrom` -- Chromosome in vareffect format (`"chr17"`, `"chrX"`, or
///   plain `"17"`). Any leading `"chr"` is stripped for the emitted
///   `seq_region_name`.
/// * `pos` -- 0-based genomic start position (BED convention).
/// * `ref_allele` -- Plus-strand reference allele bytes. May be empty for
///   pure insertions.
/// * `alt_allele` -- Plus-strand alternate allele bytes. May be empty for
///   pure deletions.
/// * `assembly` -- Genome build label (e.g. `"GRCh38"`) written to the
///   top-level `assembly_name` field.
/// * `results` -- The per-transcript annotations returned by
///   [`crate::VarEffect::annotate`].
pub(crate) fn to_vep_json_array(
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    assembly: &str,
    results: &[ConsequenceResult],
) -> Value {
    Value::Array(vec![build_annotation(
        chrom, pos, ref_allele, alt_allele, assembly, results,
    )])
}

/// Build the single top-level annotation object.
///
/// Factored out of [`to_vep_json_array`] so the array-wrapping decision lives
/// in exactly one place; a future batch-endpoint caller that wants raw
/// objects can reuse this helper without re-wrapping.
fn build_annotation(
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    assembly: &str,
    results: &[ConsequenceResult],
) -> Value {
    let seq_region = chrom.strip_prefix("chr").unwrap_or(chrom);
    // 0-based -> 1-based.
    let start = pos + 1;
    // 1-based inclusive end. For pure insertions (ref empty) VEP uses
    // start > end, but vareffect callers typically pre-anchor indels via
    // `VarEffect::anchor_prepend_indel`, giving a non-empty ref, so this
    // branch rarely fires. The `max(1, ..)` floor keeps `end >= start`.
    let end = pos + (ref_allele.len().max(1) as u64);

    let allele_string = format!(
        "{}/{}",
        render_allele(ref_allele),
        render_allele(alt_allele),
    );

    // Smallest severity_rank is the most severe consequence (VEP's ordering:
    // rank 1 = transcript_ablation, rank 24 = intergenic_variant).
    let most_severe = results
        .iter()
        .flat_map(|r| r.consequences.iter())
        .min_by_key(|c| c.severity_rank())
        .map(Consequence::as_str)
        .unwrap_or("intergenic_variant");

    let transcript_consequences: Vec<Value> =
        results.iter().map(build_transcript_consequence).collect();

    json!({
        "seq_region_name": seq_region,
        "start": start,
        "end": end,
        "strand": 1,
        "allele_string": allele_string,
        "assembly_name": assembly,
        "most_severe_consequence": most_severe,
        "transcript_consequences": transcript_consequences,
        "colocated_variants": [],
    })
}

/// Build one VEP-compatible `transcript_consequences` entry.
///
/// Uses a `serde_json::Map` with conditional inserts rather than the `json!`
/// macro so optional fields can be omitted entirely instead of serialized as
/// `null` -- this matches VEP REST's own absence semantics.
fn build_transcript_consequence(r: &ConsequenceResult) -> Value {
    let mut map = Map::new();

    // --- Required / always-present fields ---
    map.insert("transcript_id".into(), Value::String(r.transcript.clone()));
    map.insert("gene_symbol".into(), Value::String(r.gene_symbol.clone()));
    map.insert(
        "biotype".into(),
        Value::String(r.biotype.as_str().to_string()),
    );
    map.insert(
        "strand".into(),
        // VEP uses `1` / `-1` for plus / minus. We emit them as signed
        // integers to match the schema the lia consumer deserializes into
        // (`strand: Option<i8>`).
        Value::from(match r.strand {
            Strand::Plus => 1_i8,
            Strand::Minus => -1_i8,
        }),
    );
    map.insert("impact".into(), Value::String(r.impact.to_string()));
    map.insert(
        "consequence_terms".into(),
        Value::Array(
            r.consequences
                .iter()
                .map(|c| Value::String(c.as_str().to_string()))
                .collect(),
        ),
    );

    // --- HGVS with accession prefix to match real VEP output format ---
    // `hgvs_c` carries the full accession prefix (e.g. "NM_000546.6:c.742C>T"),
    // matching VEP REST's `hgvsc` field format exactly. No wrapping needed.
    if let Some(c) = r.hgvs_c.as_ref() {
        map.insert("hgvsc".into(), Value::String(c.clone()));
    }
    // `hgvs_p` is stored bare (`"p.Arg248Trp"`). We require BOTH halves --
    // emitting `":p.Arg248Trp"` with an empty accession would be malformed
    // and indistinguishable from a broken parse upstream.
    if let (Some(p), Some(acc)) = (r.hgvs_p.as_ref(), r.protein_accession.as_ref()) {
        map.insert("hgvsp".into(), Value::String(format!("{acc}:{p}")));
    }

    // --- Optional positional fields (widened u32 -> u64 to match VEP) ---
    insert_opt_u64(&mut map, "protein_start", r.protein_start.map(u64::from));
    insert_opt_u64(&mut map, "protein_end", r.protein_end.map(u64::from));
    insert_opt_u64(&mut map, "cds_start", r.cds_position.map(u64::from));
    insert_opt_u64(&mut map, "cds_end", r.cds_position_end.map(u64::from));
    insert_opt_u64(&mut map, "cdna_start", r.cdna_position.map(u64::from));
    insert_opt_u64(&mut map, "cdna_end", r.cdna_position_end.map(u64::from));

    // --- Optional descriptive fields ---
    insert_opt_str(&mut map, "amino_acids", r.amino_acids.as_deref());
    insert_opt_str(&mut map, "codons", r.codons.as_deref());
    insert_opt_str(&mut map, "exon", r.exon.as_deref());
    insert_opt_str(&mut map, "intron", r.intron.as_deref());

    // --- MANE membership flags ---
    // VEP emits the transcript accession string (not a boolean) as the value
    // of these fields to indicate membership. Downstream consumers typically
    // read them as `Option<String>` and treat any non-None value as "this
    // transcript IS MANE Select / Plus Clinical".
    if r.is_mane_select {
        map.insert("mane_select".into(), Value::String(r.transcript.clone()));
    }
    if r.is_mane_plus_clinical {
        map.insert(
            "mane_plus_clinical".into(),
            Value::String(r.transcript.clone()),
        );
    }

    Value::Object(map)
}

/// Insert an optional `u64` only when present.
///
/// VEP's own JSON omits absent numeric fields rather than emitting `null`;
/// this helper mirrors that convention so vareffect output is byte-identical
/// to real VEP output for omitted fields.
fn insert_opt_u64(map: &mut Map<String, Value>, key: &str, value: Option<u64>) {
    if let Some(v) = value {
        map.insert(key.into(), Value::from(v));
    }
}

/// Insert an optional string slice only when present.
///
/// The value is cloned into an owned `String` because `serde_json::Value`
/// only holds owned strings. This is cheap compared to the surrounding JSON
/// allocation.
fn insert_opt_str(map: &mut Map<String, Value>, key: &str, value: Option<&str>) {
    if let Some(v) = value {
        map.insert(key.into(), Value::String(v.to_string()));
    }
}

/// Render a byte allele as its `allele_string` component.
///
/// Empty bytes become `"-"` (VEP's convention for the missing allele in an
/// indel). Non-empty bytes are rendered as lossy UTF-8; vareffect only
/// produces ACGTN bytes so the conversion is always lossless in practice,
/// but `from_utf8_lossy` is used defensively to avoid panics if a caller
/// supplies a non-ASCII byte.
fn render_allele(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(bytes).into_owned()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consequence::Impact;
    use crate::types::Biotype;

    /// Minimal [`ConsequenceResult`] fixture. Tests override the fields they
    /// care about, keeping each test scannable at a glance.
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

    /// TP53 R248W minus-strand SNV -- full VEP-compatible shape, every
    /// top-level field and every per-transcript field that should be present.
    #[test]
    fn vep_json_snv_missense_full_shape() {
        // Arrange
        let r = ConsequenceResult {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            protein_start: Some(248),
            protein_end: Some(248),
            codons: Some("Cgg/Tgg".to_string()),
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

        // Act
        let json = to_vep_json_array("chr17", 7_674_219, b"C", b"T", "GRCh38", &[r]);

        // Assert -- top-level
        let arr = json.as_array().expect("top-level should be an array");
        assert_eq!(arr.len(), 1);
        let top = &arr[0];
        assert_eq!(top["seq_region_name"], "17");
        assert_eq!(top["start"], 7_674_220);
        assert_eq!(top["end"], 7_674_220);
        assert_eq!(top["strand"], 1);
        assert_eq!(top["allele_string"], "C/T");
        assert_eq!(top["assembly_name"], "GRCh38");
        assert_eq!(top["most_severe_consequence"], "missense_variant");
        assert_eq!(top["colocated_variants"].as_array().unwrap().len(), 0);

        // Assert -- per-transcript
        let tc = &top["transcript_consequences"][0];
        assert_eq!(tc["transcript_id"], "NM_000546.6");
        assert_eq!(tc["gene_symbol"], "TP53");
        assert_eq!(tc["biotype"], "protein_coding");
        assert_eq!(tc["strand"], -1);
        assert_eq!(tc["impact"], "MODERATE");
        assert_eq!(
            tc["consequence_terms"],
            serde_json::json!(["missense_variant"])
        );
        assert_eq!(tc["hgvsc"], "NM_000546.6:c.742C>T");
        assert_eq!(tc["hgvsp"], "NP_000537.3:p.Arg248Trp");
        assert_eq!(tc["protein_start"], 248);
        assert_eq!(tc["protein_end"], 248);
        assert_eq!(tc["cds_start"], 742);
        assert_eq!(tc["cds_end"], 742);
        assert_eq!(tc["cdna_start"], 884);
        assert_eq!(tc["cdna_end"], 884);
        assert_eq!(tc["amino_acids"], "R/W");
        assert_eq!(tc["codons"], "Cgg/Tgg");
        assert_eq!(tc["exon"], "7/11");
        assert_eq!(tc["mane_select"], "NM_000546.6");
    }

    /// Non-coding transcript (no protein accession). `hgvsp` must be omitted
    /// entirely rather than emitted as `":p.?"` -- guards against the
    /// accession-less bug the original spec would have caused.
    #[test]
    fn vep_json_hgvsp_omitted_without_protein_accession() {
        // Arrange
        let r = ConsequenceResult {
            consequences: vec![Consequence::NonCodingTranscriptExonVariant],
            impact: Impact::Modifier,
            biotype: Biotype::LncRna,
            // Deliberately populate hgvs_p without a protein_accession to
            // exercise the "omit unless both halves present" branch.
            hgvs_p: Some("p.?".to_string()),
            protein_accession: None,
            ..make_result("NR_046018.2", "DDX11L1")
        };

        // Act
        let json = to_vep_json_array("chr1", 11_868, b"G", b"A", "GRCh38", &[r]);

        // Assert
        let tc = &json[0]["transcript_consequences"][0];
        assert!(
            tc.as_object().unwrap().get("hgvsp").is_none(),
            "hgvsp should be omitted when protein_accession is None; got: {tc:?}"
        );
    }

    /// Frameshift deletion `TG -> -`. Exercises the empty-alt allele_string
    /// branch and the HIGH impact category.
    #[test]
    fn vep_json_indel_frameshift_allele_string() {
        // Arrange
        let r = ConsequenceResult {
            consequences: vec![Consequence::FrameshiftVariant],
            impact: Impact::High,
            hgvs_c: Some("NM_006772.2:c.1861_1862del".to_string()),
            ..make_result("NM_006772.2", "SYNGAP1")
        };

        // Act
        let json = to_vep_json_array("chr6", 33_409_450, b"TG", b"", "GRCh38", &[r]);

        // Assert -- top-level
        let top = &json[0];
        assert_eq!(top["start"], 33_409_451);
        assert_eq!(top["end"], 33_409_452);
        assert_eq!(top["allele_string"], "TG/-");
        assert_eq!(top["most_severe_consequence"], "frameshift_variant");

        // Assert -- per-transcript
        let tc = &top["transcript_consequences"][0];
        assert_eq!(tc["impact"], "HIGH");
        assert_eq!(
            tc["consequence_terms"],
            serde_json::json!(["frameshift_variant"])
        );
        assert_eq!(tc["hgvsc"], "NM_006772.2:c.1861_1862del");
    }

    /// Intergenic variant: annotate returned an empty slice. We still emit
    /// a well-formed single-element array with an empty
    /// `transcript_consequences` and the `"intergenic_variant"` fallback.
    #[test]
    fn vep_json_intergenic_empty_results() {
        // Act
        let json = to_vep_json_array("chr1", 1_000_000, b"A", b"G", "GRCh38", &[]);

        // Assert
        let arr = json.as_array().expect("top-level should be an array");
        assert_eq!(arr.len(), 1);
        let top = &arr[0];
        assert_eq!(top["most_severe_consequence"], "intergenic_variant");
        assert_eq!(top["transcript_consequences"].as_array().unwrap().len(), 0);
    }

    /// Two consequences on one transcript, in the order
    /// [SpliceRegionVariant, MissenseVariant]. Rank 10 (missense) is more
    /// severe than rank 12 (splice_region), so the top-level
    /// `most_severe_consequence` must report missense, but the per-transcript
    /// `consequence_terms` must preserve the input order.
    #[test]
    fn vep_json_most_severe_from_min_severity_rank() {
        // Arrange
        let r = ConsequenceResult {
            consequences: vec![
                Consequence::SpliceRegionVariant,
                Consequence::MissenseVariant,
            ],
            impact: Impact::Moderate,
            ..make_result("NM_000546.6", "TP53")
        };

        // Act
        let json = to_vep_json_array("chr17", 7_674_219, b"C", b"T", "GRCh38", &[r]);

        // Assert
        let top = &json[0];
        assert_eq!(top["most_severe_consequence"], "missense_variant");
        let tc = &top["transcript_consequences"][0];
        assert_eq!(
            tc["consequence_terms"],
            serde_json::json!(["splice_region_variant", "missense_variant"])
        );
    }

    /// Intron variant. Every protein/CDS/exon field is `None`, so the
    /// serialized object must contain `intron` but have no keys for `exon`,
    /// `amino_acids`, `codons`, `protein_start`, `protein_end`, `hgvsp`,
    /// `cds_start`, `cds_end`.
    #[test]
    fn vep_json_optional_fields_omitted_for_intron_variant() {
        // Arrange
        let r = ConsequenceResult {
            consequences: vec![Consequence::IntronVariant],
            intron: Some("7/10".to_string()),
            ..make_result("NM_000546.6", "TP53")
        };

        // Act
        let json = to_vep_json_array("chr17", 7_675_000, b"A", b"G", "GRCh38", &[r]);

        // Assert
        let tc = &json[0]["transcript_consequences"][0];
        let obj = tc.as_object().unwrap();
        assert_eq!(tc["intron"], "7/10");
        for missing_key in [
            "exon",
            "amino_acids",
            "codons",
            "protein_start",
            "protein_end",
            "hgvsp",
            "hgvsc",
            "cds_start",
            "cds_end",
            "cdna_start",
            "cdna_end",
        ] {
            assert!(
                obj.get(missing_key).is_none(),
                "{missing_key} should be omitted for an intron variant; got: {tc:?}"
            );
        }
    }

    /// Two transcripts, one synonymous (rank 16) and one missense (rank 10).
    /// Top-level `most_severe_consequence` must report missense; the array
    /// must preserve input order.
    #[test]
    fn vep_json_multiple_transcripts_ordering() {
        // Arrange
        let synonymous = ConsequenceResult {
            consequences: vec![Consequence::SynonymousVariant],
            impact: Impact::Low,
            ..make_result("NM_000546.6", "TP53")
        };
        let missense = ConsequenceResult {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            ..make_result("NM_001126112.3", "TP53")
        };

        // Act
        let json = to_vep_json_array(
            "chr17",
            7_674_219,
            b"C",
            b"T",
            "GRCh38",
            &[synonymous, missense],
        );

        // Assert
        let top = &json[0];
        assert_eq!(top["most_severe_consequence"], "missense_variant");
        let tcs = top["transcript_consequences"].as_array().unwrap();
        assert_eq!(tcs.len(), 2);
        assert_eq!(tcs[0]["transcript_id"], "NM_000546.6");
        assert_eq!(tcs[0]["impact"], "LOW");
        assert_eq!(tcs[1]["transcript_id"], "NM_001126112.3");
        assert_eq!(tcs[1]["impact"], "MODERATE");
    }
}
