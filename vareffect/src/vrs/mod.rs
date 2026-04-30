//! GA4GH VRS computed-identifier emission.
//!
//! Produces `ga4gh:VA.<digest>` Allele identifiers under both the VRS
//! 1.3 and VRS 2.0 schemas, attached to [`crate::AnnotationResult`] for
//! cross-build evidence resolution against indexes that key on VRS ID.
//!
//! # Pipeline
//!
//! 1. **VOCA normalization** ([`normalize::voca_normalize`]) — fully-
//!    justified expand-and-trim using flanking reference reads from
//!    [`crate::FastaReader`]. Required for canonical IDs that match
//!    anyvar / ClinGen / ClinVar / MAVEDB on repeat-region indels.
//! 2. **SQ digest resolution** ([`sq_digest::sq_curie`]) — content-
//!    addressed digest of the chromosome sequence; computed lazily
//!    from FASTA bytes on first request and cached process-wide.
//! 3. **Recursive `ga4gh_serialize`** ([`serialize`]) — canonical JSON
//!    blob with sorted keys, no whitespace, and nested identifiable
//!    objects replaced by their bare 32-char digest strings.
//! 4. **Digest** ([`digest::ga4gh_identify`]) — SHA-512 → first 24
//!    bytes → URL-safe base64 (no padding) → `ga4gh:VA.` prefix.
//!
//! # Failure mode
//!
//! Per the implementation plan, VRS computation **never fails the
//! enclosing annotation**. Any internal error (chrom not in FASTA,
//! FASTA read past chromosome end during VOCA expansion, non-primary
//! contig with no canonical SQ digest) is swallowed and surfaces as
//! `(None, None)`. Callers see `vrs_id: None` / `vrs_id_v2: None` for
//! the affected variant and continue normally.

mod digest;
mod normalize;
mod serialize;
mod sq_digest;

use crate::chrom::Assembly;
use crate::fasta::FastaReader;

/// Compute the VRS 1.3 and VRS 2.0 Allele identifiers for a variant.
///
/// Returns `(vrs_1_3_id, vrs_2_0_id)` where each component is `Some`
/// only when the full pipeline succeeds for that schema. The two
/// schemas are computed independently; either can succeed while the
/// other fails (in practice they always co-succeed since they share
/// the VOCA + SQ steps and only differ in the final serialization).
///
/// # Arguments
///
/// * `assembly` — the genome build the variant is called against.
/// * `chrom` — UCSC-style chromosome name (`"chr17"`, `"chrM"`).
/// * `pos` — 0-based VCF position.
/// * `ref_allele` — VCF REF bytes (uppercase ASCII).
/// * `alt_allele` — VCF ALT bytes (uppercase ASCII).
/// * `fasta` — reference reader for VOCA flank reads and SQ digest.
pub(crate) fn compute_vrs_ids(
    assembly: Assembly,
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    fasta: &FastaReader,
) -> (Option<String>, Option<String>) {
    match compute_vrs_ids_inner(assembly, chrom, pos, ref_allele, alt_allele, fasta) {
        Some((v1, v2)) => (Some(v1), Some(v2)),
        None => (None, None),
    }
}

fn compute_vrs_ids_inner(
    assembly: Assembly,
    chrom: &str,
    pos: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    fasta: &FastaReader,
) -> Option<(String, String)> {
    // Short-circuit non-primary contigs (alt / random / unlocalized / patch)
    // before doing any FASTA reads in VOCA. These have no canonical
    // cross-pipeline SQ digest and would be discarded post-normalization.
    if !sq_digest::is_primary_contig(assembly, chrom) {
        return None;
    }

    let normalized =
        normalize::voca_normalize(chrom, pos, ref_allele, alt_allele, fasta).ok()??;
    let sq_curie = sq_digest::sq_curie(assembly, chrom, fasta)?;
    let alt_str = std::str::from_utf8(&normalized.alt).ok()?;

    let id_v1 = serialize::vrs1_allele_id(&sq_curie, normalized.start, normalized.end, alt_str);
    let id_v2 = serialize::vrs2_allele_id(&sq_curie, normalized.start, normalized.end, alt_str);
    Some((id_v1, id_v2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    fn synthetic_fasta() -> (TempDir, FastaReader) {
        let tmp = TempDir::new().unwrap();
        let bin = tmp.path().join("t.bin");
        let idx = tmp.path().join("t.bin.idx");
        let chr1: &[u8] = b"GGAAAAACCATATATGGATGATGTTTTTTT";
        write_genome_binary(&[("chr1", chr1)], "test", &bin, &idx).unwrap();
        let reader = FastaReader::open_with_assembly(&bin, Assembly::GRCh38).unwrap();
        (tmp, reader)
    }

    #[test]
    fn snv_produces_both_vrs_ids() {
        let (_tmp, fasta) = synthetic_fasta();
        let (v1, v2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 0, b"G", b"C", &fasta);
        let v1 = v1.expect("VRS 1.3 ID should be produced for SNV");
        let v2 = v2.expect("VRS 2.0 ID should be produced for SNV");
        assert!(v1.starts_with("ga4gh:VA."), "got {v1}");
        assert!(v2.starts_with("ga4gh:VA."), "got {v2}");
        assert_eq!(v1.len(), "ga4gh:VA.".len() + 32);
        assert_eq!(v2.len(), "ga4gh:VA.".len() + 32);
        assert_ne!(v1, v2, "1.3 and 2.0 schemas must produce different IDs");
    }

    #[test]
    fn snv_id_is_deterministic_across_calls() {
        // Same input must produce the same VRS ID — content addressing.
        let (_tmp, fasta) = synthetic_fasta();
        let (a1, a2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 0, b"G", b"C", &fasta);
        let (b1, b2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 0, b"G", b"C", &fasta);
        assert_eq!(a1, b1);
        assert_eq!(a2, b2);
    }

    #[test]
    fn ambiguous_indels_in_repeat_collapse_to_one_id() {
        // Three VCF representations of the same biological insertion
        // in the di-AT repeat should produce identical VRS IDs after
        // VOCA normalization. This is the entire point of VOCA.
        let (_tmp, fasta) = synthetic_fasta();

        // Three equivalent representations of "insert one AT unit into
        // the ATATAT repeat at [9,15)":
        // (a) anchored at the boundary G before the repeat (pos 8 'C').
        let (a1, a2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 8, b"C", b"CAT", &fasta);
        // (b) anchored mid-repeat at pos 10 'T' → GATATAT vs GTATATAT.
        //     Trim leaves AT inserted between repeat units.
        let (b1, b2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 10, b"T", b"TAT", &fasta);
        // (c) anchored at end of repeat at pos 14 'T' → adds an AT unit.
        let (c1, c2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 14, b"T", b"TAT", &fasta);

        let a1 = a1.unwrap();
        let b1 = b1.unwrap();
        let c1 = c1.unwrap();
        assert_eq!(a1, b1, "VRS 1.3 IDs must collide across equivalent calls");
        assert_eq!(a1, c1, "VRS 1.3 IDs must collide across equivalent calls");

        let a2 = a2.unwrap();
        let b2 = b2.unwrap();
        let c2 = c2.unwrap();
        assert_eq!(a2, b2, "VRS 2.0 IDs must collide across equivalent calls");
        assert_eq!(a2, c2, "VRS 2.0 IDs must collide across equivalent calls");
    }

    #[test]
    fn no_op_variant_yields_none() {
        let (_tmp, fasta) = synthetic_fasta();
        let (v1, v2) = compute_vrs_ids(Assembly::GRCh38, "chr1", 0, b"GG", b"GG", &fasta);
        assert!(v1.is_none());
        assert!(v2.is_none());
    }

    #[test]
    fn unknown_chrom_yields_none() {
        let (_tmp, fasta) = synthetic_fasta();
        let (v1, v2) = compute_vrs_ids(Assembly::GRCh38, "chrUNKNOWN", 0, b"A", b"G", &fasta);
        assert!(v1.is_none());
        assert!(v2.is_none());
    }
}
