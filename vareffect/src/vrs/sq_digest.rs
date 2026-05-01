//! Per-chromosome `ga4gh:SQ.<digest>` lookup for primary contigs.
//!
//! VRS Allele identifiers reference their underlying sequence by
//! content-addressed digest (`ga4gh:SQ.<sha512t24u(uppercase ASCII bases)>`).
//! For canonical references, these digests are stable across the GA4GH
//! ecosystem (anyvar, ClinGen Allele Registry, ClinVar, MAVEDB, seqrepo).
//!
//! # Sourcing strategy
//!
//! Each digest is computed once from the bytes the local [`FastaReader`]
//! actually serves and cached **on the reader itself**, so two readers
//! opened against different references in the same process produce
//! digests of their own bytes (not whichever populated a process-global
//! cache first). For a canonical GRCh38.p14 / GRCh37.p13 reference the
//! computed digests match the published seqrepo ones byte-for-byte; if
//! they ever diverge, that signals our shipped FASTA differs from the
//! canonical reference — the same divergence that would also break HGVS
//! notation, MANE matching, and every other downstream.
//!
//! # Cache
//!
//! Lives on [`FastaReader::vrs_sq_cached`]. Read-heavy, write-rare
//! `RwLock<HashMap<String, Arc<str>>>` — cache hits clone the `Arc` under
//! a shared read lock (no allocation, no contention with other readers).

use std::sync::Arc;

use rayon::prelude::*;

use crate::chrom::{Assembly, ucsc_to_refseq};
use crate::fasta::FastaReader;

/// Set of UCSC chromosome names for which a canonical SQ digest is
/// well-defined — the 25 primary contigs (1-22, X, Y, MT) for both
/// GRCh37 and GRCh38. Patch contigs, alt haplotypes, random and
/// unlocalized scaffolds are deliberately excluded — they do not have a
/// stable cross-pipeline SQ digest in the GA4GH ecosystem and are
/// therefore unindexable by VRS ID.
pub(super) const PRIMARY_UCSC_CHROMS: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM",
];

/// Return whether `chrom` has a canonical cross-pipeline SQ digest on
/// `assembly`. Primary contigs only — patch / alt / random / unlocalized
/// scaffolds return `false`.
pub(super) fn is_primary_contig(assembly: Assembly, chrom: &str) -> bool {
    if !PRIMARY_UCSC_CHROMS.contains(&chrom) {
        return false;
    }
    // `ucsc_to_refseq` returns the input unchanged for non-primary contigs;
    // the round-trip identity is the assembly-aware sanity check.
    ucsc_to_refseq(assembly, chrom) != chrom
}

/// Return the `ga4gh:SQ.<digest>` CURIE for the given chromosome on the
/// given assembly, computing it from the FASTA bytes on first request and
/// caching the result on the reader thereafter.
///
/// Returns `None` for non-primary contigs (patches, alt haplotypes,
/// random / unlocalized scaffolds) and for chromosomes not present in the
/// reader (e.g. trimmed test fixtures).
///
/// Special case: `chrM` resolves to the same digest on both GRCh37 and
/// GRCh38 because both NCBI assemblies use the rCRS mitochondrial
/// reference (`NC_012920.1`).
pub(super) fn sq_curie(assembly: Assembly, chrom: &str, fasta: &FastaReader) -> Option<Arc<str>> {
    if !is_primary_contig(assembly, chrom) {
        return None;
    }
    fasta.vrs_sq_cached(chrom, || {
        let length = fasta.chrom_length(chrom)?;
        let bytes = fasta.fetch_sequence(chrom, 0, length).ok()?;
        let digest = super::digest::sha512t24u(&bytes);
        Some(format!("ga4gh:SQ.{digest}"))
    })
}

/// Strip the `ga4gh:SQ.` prefix from a CURIE, returning just the bare
/// 32-char digest. Used by the recursive `ga4gh_serialize` step that
/// reduces existing CURIEs to their digest portion before hashing.
pub(super) fn bare_digest_from_curie(curie: &str) -> &str {
    curie.strip_prefix("ga4gh:SQ.").unwrap_or(curie)
}

/// Eagerly fill the per-`FastaReader` SQ-digest cache for every primary
/// contig in parallel, fanning the SHA-512 work across rayon's pool.
///
/// Concurrency is safe because:
/// - The compute closure inside [`FastaReader::vrs_sq_cached`] runs
///   outside the cache's write lock, so parallel SHA-512 work does not
///   serialize on lock acquisition.
/// - Concurrent inserts on the same key are deduplicated via
///   `entry().or_insert_with(...)`. Our fan-out gives each rayon task a
///   unique chrom, so this is belt-and-suspenders rather than load-
///   bearing.
///
/// Per-chrom failures (chrom absent from this reader's index, e.g. a
/// trimmed test fixture) are silently dropped — the same swallow-and-
/// continue semantics as the lazy fill path.
pub(super) fn warm_cache(fasta: &FastaReader) {
    let assembly = fasta.assembly();
    PRIMARY_UCSC_CHROMS.par_iter().for_each(|chrom| {
        let _ = sq_curie(assembly, chrom, fasta);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_primary_contigs_excluded() {
        assert!(!PRIMARY_UCSC_CHROMS.contains(&"chr9_KN196479v1_fix"));
        assert!(!PRIMARY_UCSC_CHROMS.contains(&"chr22_KI270879v1_alt"));
        assert!(!PRIMARY_UCSC_CHROMS.contains(&"chrUn_KI270302v1"));

        assert!(!is_primary_contig(Assembly::GRCh38, "chr9_KN196479v1_fix"));
        assert!(!is_primary_contig(Assembly::GRCh38, "totally_unknown"));
    }

    #[test]
    fn all_primary_chroms_have_refseq_accessions() {
        for &chrom in PRIMARY_UCSC_CHROMS {
            for assembly in [Assembly::GRCh38, Assembly::GRCh37] {
                let acc = ucsc_to_refseq(assembly, chrom);
                assert_ne!(
                    acc, chrom,
                    "primary chrom {chrom} on {assembly:?} has no RefSeq accession"
                );
                assert!(
                    acc.starts_with("NC_"),
                    "expected NC_ accession for {chrom} on {assembly:?}, got {acc}"
                );
                assert!(is_primary_contig(assembly, chrom));
            }
        }
    }

    #[test]
    fn bare_digest_strips_curie_prefix() {
        assert_eq!(
            bare_digest_from_curie("ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"),
            "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
        );
        assert_eq!(
            bare_digest_from_curie("Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"),
            "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
        );
    }
}
