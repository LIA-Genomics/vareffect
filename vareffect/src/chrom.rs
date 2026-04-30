//! Chromosome name conversion and genome-build identification.
//!
//! Maps NCBI RefSeq accessions (e.g. `"NC_000006.12"`) to UCSC-style
//! chromosome names (e.g. `"chr6"`). Patch-sequence accessions (`NW_*`,
//! `NT_*`) are returned unchanged so the transcripts on them remain
//! round-trippable, even though they can't be looked up by standard
//! chromosome name.
//!
//! The accession-version table is **assembly-specific** — GRCh37 and GRCh38
//! have different patch versions for every primary chromosome except chrM
//! (`NC_012920.1`, the rCRS, is identical between assemblies). All
//! chromosome-name conversion functions take an [`Assembly`] selector.
//!
//! Unknown inputs (anything not in the selected table and not recognized as
//! a patch accession) are returned as-is.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::error::VarEffectError;

/// Reference genome build that a [`crate::TranscriptStore`] /
/// [`crate::FastaReader`] / [`crate::VarEffect`] is bound to.
///
/// vareffect supports loading both assemblies into a single `VarEffect`
/// instance, so every public API that takes coordinates also takes an
/// `Assembly` to route the lookup to the correct chrom table and the
/// correct genome.
///
/// The enum has **no `Default` impl on purpose** — silent fallback to
/// GRCh38 in a multi-assembly process would corrupt every coordinate
/// lookup for GRCh37 callers. Construction must be explicit at every call
/// site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Assembly {
    /// NCBI GRCh38 / UCSC `hg38`. Patched up through GRCh38.p14.
    GRCh38,
    /// NCBI GRCh37 / UCSC `hg19`. Patched up through GRCh37.p13.
    GRCh37,
}

impl Assembly {
    /// Stable lowercase identifier for serialization, CLI flags, and
    /// error messages.
    pub fn as_str(&self) -> &'static str {
        match self {
            Assembly::GRCh38 => "GRCh38",
            Assembly::GRCh37 => "GRCh37",
        }
    }
}

impl fmt::Display for Assembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for Assembly {
    type Err = VarEffectError;

    /// Parse `"grch38"` / `"GRCh38"` (case-insensitive) into [`Assembly::GRCh38`]
    /// and `"grch37"` / `"GRCh37"` into [`Assembly::GRCh37`].
    ///
    /// `"hg19"` and `"hg38"` are **explicitly rejected**: UCSC `hg19` chrM
    /// (`NC_001807`) differs from GRCh37 chrMT (`NC_012920.1`, the rCRS) by
    /// ~10 bases and shifts coordinates throughout the mitochondrial genome.
    /// Silently accepting the alias would mis-annotate every chrM variant for
    /// users who genotyped against the actual UCSC distribution. The error
    /// message points the caller at the explicit `grch37` / `grch38` form.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "grch38" => Ok(Assembly::GRCh38),
            "grch37" => Ok(Assembly::GRCh37),
            "hg19" | "hg38" => Err(VarEffectError::UnsupportedAssemblyAlias {
                input: s.to_string(),
            }),
            _ => Err(VarEffectError::UnsupportedAssemblyAlias {
                input: s.to_string(),
            }),
        }
    }
}

/// GRCh38.p14 `NC_*` accession → UCSC chromosome name.
///
/// The table uses the versioned accessions from MANE v1.5 / GRCh38.p14.
/// If a future MANE release bumps the patch version (e.g. `NC_000001.12`),
/// add it here alongside the existing entry rather than replacing it —
/// older cached GFF3 files should continue to resolve correctly.
const NC_TO_UCSC_GRCH38: &[(&str, &str)] = &[
    ("NC_000001.11", "chr1"),
    ("NC_000002.12", "chr2"),
    ("NC_000003.12", "chr3"),
    ("NC_000004.12", "chr4"),
    ("NC_000005.10", "chr5"),
    ("NC_000006.12", "chr6"),
    ("NC_000007.14", "chr7"),
    ("NC_000008.11", "chr8"),
    ("NC_000009.12", "chr9"),
    ("NC_000010.11", "chr10"),
    ("NC_000011.10", "chr11"),
    ("NC_000012.12", "chr12"),
    ("NC_000013.11", "chr13"),
    ("NC_000014.9", "chr14"),
    ("NC_000015.10", "chr15"),
    ("NC_000016.10", "chr16"),
    ("NC_000017.11", "chr17"),
    ("NC_000018.10", "chr18"),
    ("NC_000019.10", "chr19"),
    ("NC_000020.11", "chr20"),
    ("NC_000021.9", "chr21"),
    ("NC_000022.11", "chr22"),
    ("NC_000023.11", "chrX"),
    ("NC_000024.10", "chrY"),
    ("NC_012920.1", "chrM"),
];

/// GRCh37.p13 `NC_*` accession → UCSC chromosome name.
///
/// Versioned to `GCF_000001405.25_GRCh37.p13`. Note that `NC_012920.1`
/// (chrM) is **identical to the GRCh38 entry**: both NCBI assemblies use
/// the rCRS mitochondrial reference. UCSC `hg19` chrM is a different
/// sequence (`NC_001807`) — vareffect does not support `hg19` directly;
/// see [`Assembly::from_str`] for the rejection rationale.
const NC_TO_UCSC_GRCH37: &[(&str, &str)] = &[
    ("NC_000001.10", "chr1"),
    ("NC_000002.11", "chr2"),
    ("NC_000003.11", "chr3"),
    ("NC_000004.11", "chr4"),
    ("NC_000005.9", "chr5"),
    ("NC_000006.11", "chr6"),
    ("NC_000007.13", "chr7"),
    ("NC_000008.10", "chr8"),
    ("NC_000009.11", "chr9"),
    ("NC_000010.10", "chr10"),
    ("NC_000011.9", "chr11"),
    ("NC_000012.11", "chr12"),
    ("NC_000013.10", "chr13"),
    ("NC_000014.8", "chr14"),
    ("NC_000015.9", "chr15"),
    ("NC_000016.9", "chr16"),
    ("NC_000017.10", "chr17"),
    ("NC_000018.9", "chr18"),
    ("NC_000019.9", "chr19"),
    ("NC_000020.10", "chr20"),
    ("NC_000021.8", "chr21"),
    ("NC_000022.10", "chr22"),
    ("NC_000023.10", "chrX"),
    ("NC_000024.9", "chrY"),
    ("NC_012920.1", "chrM"),
];

/// Pick the right `NC_*`-to-UCSC table for the requested assembly.
fn table_for(assembly: Assembly) -> &'static [(&'static str, &'static str)] {
    match assembly {
        Assembly::GRCh38 => NC_TO_UCSC_GRCH38,
        Assembly::GRCh37 => NC_TO_UCSC_GRCH37,
    }
}

/// Map a RefSeq chromosome accession to a UCSC-style name for the given
/// assembly.
///
/// Inputs not in the hardcoded standard-chromosome table (including patch
/// sequences with prefixes `NW_*` or `NT_*`) are returned unchanged.
///
/// # Examples
///
/// ```
/// use vareffect::chrom::{refseq_to_ucsc, Assembly};
/// assert_eq!(refseq_to_ucsc(Assembly::GRCh38, "NC_000006.12"), "chr6");
/// assert_eq!(refseq_to_ucsc(Assembly::GRCh37, "NC_000006.11"), "chr6");
/// assert_eq!(refseq_to_ucsc(Assembly::GRCh38, "NW_025791820.1"), "NW_025791820.1");
/// assert_eq!(refseq_to_ucsc(Assembly::GRCh38, "gibberish"), "gibberish");
/// ```
pub fn refseq_to_ucsc(assembly: Assembly, acc: &str) -> &str {
    // Linear scan over 25 entries — a HashMap would be overkill for such a
    // small static table and would require a one-time allocation per
    // process. The inner loop happens once per mRNA/gene row during GFF3
    // ingest; even for 19k transcripts the cumulative cost is < 1 ms.
    for (key, value) in table_for(assembly) {
        if *key == acc {
            return value;
        }
    }
    acc
}

/// Map a UCSC-style chromosome name to a RefSeq accession for the given
/// assembly.
///
/// Inverse of [`refseq_to_ucsc`] using the same 25-entry const table.
/// Inputs not in the table (patch sequences, UCSC-style patch names like
/// `chr9_KN196479v1_fix`, arbitrary strings) are returned unchanged — the
/// [`FastaReader`](crate::FastaReader) consults a runtime alias table for
/// patch-contig translation and relies on the pass-through for any
/// caller that already holds a RefSeq accession.
///
/// # Examples
///
/// ```
/// use vareffect::chrom::{ucsc_to_refseq, Assembly};
/// assert_eq!(ucsc_to_refseq(Assembly::GRCh38, "chr6"), "NC_000006.12");
/// assert_eq!(ucsc_to_refseq(Assembly::GRCh37, "chr6"), "NC_000006.11");
/// assert_eq!(ucsc_to_refseq(Assembly::GRCh38, "chrM"), "NC_012920.1");
/// assert_eq!(ucsc_to_refseq(Assembly::GRCh37, "chrM"), "NC_012920.1");
/// assert_eq!(ucsc_to_refseq(Assembly::GRCh38, "chr9_KN196479v1_fix"), "chr9_KN196479v1_fix");
/// ```
pub fn ucsc_to_refseq(assembly: Assembly, ucsc: &str) -> &str {
    // Symmetric linear scan over the same 25-entry table as `refseq_to_ucsc`.
    // Inverting via a HashMap would duplicate the table at startup for no
    // measurable benefit — the reverse lookup fires once per `FastaReader`
    // query (hot path is variant resolution, ~1k queries/s peak) and the
    // const-table scan is faster than a HashMap probe at that size.
    for (refseq, value) in table_for(assembly) {
        if *value == ucsc {
            return refseq;
        }
    }
    ucsc
}

/// Return `true` if `chrom` is a non-primary contig (GRCh37 or GRCh38) —
/// either in the RefSeq accession form (`NW_*`, `NT_*` — used by NCBI MANE
/// summary TSVs) or the UCSC alt/fix/random/unlocalized form (used by NCBI
/// MANE GFF3 column 1).
///
/// These are parsed and retained in the store, but transcripts on them
/// cannot be looked up by a standard chromosome name like `chr6`. The
/// prefix patterns are identical across GRCh37 and GRCh38, so this
/// function does not need an [`Assembly`] selector.
///
/// Recognized UCSC patterns:
/// - `chr*_*_fix` — GRCh38 fix patches (e.g. `chr9_KN196479v1_fix`)
/// - `chr*_*_alt` — GRCh38 alternate haplotypes (e.g. `chr22_KI270879v1_alt`)
/// - `chr*_*_random` — unlocalized contigs (e.g. `chr1_KI270706v1_random`,
///   hg19 `chr1_gl000191_random`)
/// - `chrUn_*` — unplaced scaffolds (e.g. `chrUn_KI270302v1`,
///   hg19 `chrUn_gl000211`)
/// - `chr*_*_hap[1-7]` — hg19 MHC alternate haplotypes (e.g.
///   `chr6_apd_hap1`, `chr6_cox_hap2`, `chr4_ctg9_hap1`,
///   `chr17_ctg5_hap1`). hg19's MHC region carries seven haplotypes
///   (`apd`/`cox`/`dbb`/`mann`/`mcf`/`qbl`/`ssto`); the suffix range
///   `_hap1`..`_hap7` covers the full set.
pub fn is_patch_sequence(chrom: &str) -> bool {
    // RefSeq forms (used by MANE summary TSV and other NCBI products).
    if chrom.starts_with("NW_") || chrom.starts_with("NT_") {
        return true;
    }
    // UCSC forms shared between hg19 and hg38.
    if chrom.ends_with("_alt")
        || chrom.ends_with("_fix")
        || chrom.ends_with("_random")
        || chrom.starts_with("chrUn_")
    {
        return true;
    }
    // hg19 MHC alt haplotypes: `_hap1` .. `_hap7`. Allocation-free
    // suffix check — `format!("_hap{h}")` would alloc seven Strings on
    // every call and this sits on the cross-validation hot path.
    if let Some((_, suffix)) = chrom.rsplit_once("_hap")
        && suffix.len() == 1
        && matches!(suffix.as_bytes()[0], b'1'..=b'7')
    {
        return true;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maps_all_standard_chromosomes_grch38() {
        // Exhaustive round-trip of the 25 standard GRCh38.p14 chromosomes.
        let expected = [
            ("NC_000001.11", "chr1"),
            ("NC_000002.12", "chr2"),
            ("NC_000003.12", "chr3"),
            ("NC_000004.12", "chr4"),
            ("NC_000005.10", "chr5"),
            ("NC_000006.12", "chr6"),
            ("NC_000007.14", "chr7"),
            ("NC_000008.11", "chr8"),
            ("NC_000009.12", "chr9"),
            ("NC_000010.11", "chr10"),
            ("NC_000011.10", "chr11"),
            ("NC_000012.12", "chr12"),
            ("NC_000013.11", "chr13"),
            ("NC_000014.9", "chr14"),
            ("NC_000015.10", "chr15"),
            ("NC_000016.10", "chr16"),
            ("NC_000017.11", "chr17"),
            ("NC_000018.10", "chr18"),
            ("NC_000019.10", "chr19"),
            ("NC_000020.11", "chr20"),
            ("NC_000021.9", "chr21"),
            ("NC_000022.11", "chr22"),
            ("NC_000023.11", "chrX"),
            ("NC_000024.10", "chrY"),
            ("NC_012920.1", "chrM"),
        ];
        for (acc, expected_ucsc) in expected {
            assert_eq!(
                refseq_to_ucsc(Assembly::GRCh38, acc),
                expected_ucsc,
                "GRCh38 mapping failed for {acc}"
            );
        }
    }

    #[test]
    fn maps_all_standard_chromosomes_grch37() {
        // Exhaustive round-trip of the 25 standard GRCh37.p13 chromosomes.
        let expected = [
            ("NC_000001.10", "chr1"),
            ("NC_000002.11", "chr2"),
            ("NC_000003.11", "chr3"),
            ("NC_000004.11", "chr4"),
            ("NC_000005.9", "chr5"),
            ("NC_000006.11", "chr6"),
            ("NC_000007.13", "chr7"),
            ("NC_000008.10", "chr8"),
            ("NC_000009.11", "chr9"),
            ("NC_000010.10", "chr10"),
            ("NC_000011.9", "chr11"),
            ("NC_000012.11", "chr12"),
            ("NC_000013.10", "chr13"),
            ("NC_000014.8", "chr14"),
            ("NC_000015.9", "chr15"),
            ("NC_000016.9", "chr16"),
            ("NC_000017.10", "chr17"),
            ("NC_000018.9", "chr18"),
            ("NC_000019.9", "chr19"),
            ("NC_000020.10", "chr20"),
            ("NC_000021.8", "chr21"),
            ("NC_000022.10", "chr22"),
            ("NC_000023.10", "chrX"),
            ("NC_000024.9", "chrY"),
            ("NC_012920.1", "chrM"),
        ];
        for (acc, expected_ucsc) in expected {
            assert_eq!(
                refseq_to_ucsc(Assembly::GRCh37, acc),
                expected_ucsc,
                "GRCh37 mapping failed for {acc}"
            );
        }
    }

    #[test]
    fn grch37_and_grch38_accessions_differ_except_chrm() {
        // Property: every primary chromosome must have a *different* RefSeq
        // accession version between assemblies (e.g. NC_000001.10 vs .11).
        // chrM is the documented exception — both NCBI assemblies use the
        // rCRS (NC_012920.1).
        assert_eq!(NC_TO_UCSC_GRCH37.len(), NC_TO_UCSC_GRCH38.len());
        for ((g37_acc, g37_ucsc), (g38_acc, g38_ucsc)) in
            NC_TO_UCSC_GRCH37.iter().zip(NC_TO_UCSC_GRCH38.iter())
        {
            assert_eq!(
                g37_ucsc, g38_ucsc,
                "table ordering must align between assemblies"
            );
            if *g37_ucsc == "chrM" {
                assert_eq!(
                    g37_acc, g38_acc,
                    "chrM must be NC_012920.1 in both assemblies (rCRS)"
                );
            } else {
                assert_ne!(
                    g37_acc, g38_acc,
                    "primary chromosome accessions must differ between GRCh37 and GRCh38: {g37_ucsc}"
                );
            }
        }
    }

    #[test]
    fn patch_sequences_round_trip_unchanged() {
        assert_eq!(
            refseq_to_ucsc(Assembly::GRCh38, "NW_025791820.1"),
            "NW_025791820.1"
        );
        assert_eq!(
            refseq_to_ucsc(Assembly::GRCh37, "NW_025791820.1"),
            "NW_025791820.1"
        );
        assert_eq!(
            refseq_to_ucsc(Assembly::GRCh38, "NT_187633.1"),
            "NT_187633.1"
        );
    }

    #[test]
    fn unknown_accession_returned_as_is() {
        assert_eq!(refseq_to_ucsc(Assembly::GRCh38, "gibberish"), "gibberish");
        assert_eq!(refseq_to_ucsc(Assembly::GRCh37, "gibberish"), "gibberish");
        assert_eq!(refseq_to_ucsc(Assembly::GRCh38, ""), "");
    }

    #[test]
    fn ucsc_to_refseq_maps_all_standard_chromosomes_grch38() {
        // Mirror of `maps_all_standard_chromosomes_grch38` for the inverse
        // direction. Exhaustive so a drift between the two mapping functions
        // is caught at CI time rather than surfacing as a silent FASTA lookup
        // miss.
        let expected = [
            ("chr1", "NC_000001.11"),
            ("chr2", "NC_000002.12"),
            ("chr3", "NC_000003.12"),
            ("chr4", "NC_000004.12"),
            ("chr5", "NC_000005.10"),
            ("chr6", "NC_000006.12"),
            ("chr7", "NC_000007.14"),
            ("chr8", "NC_000008.11"),
            ("chr9", "NC_000009.12"),
            ("chr10", "NC_000010.11"),
            ("chr11", "NC_000011.10"),
            ("chr12", "NC_000012.12"),
            ("chr13", "NC_000013.11"),
            ("chr14", "NC_000014.9"),
            ("chr15", "NC_000015.10"),
            ("chr16", "NC_000016.10"),
            ("chr17", "NC_000017.11"),
            ("chr18", "NC_000018.10"),
            ("chr19", "NC_000019.10"),
            ("chr20", "NC_000020.11"),
            ("chr21", "NC_000021.9"),
            ("chr22", "NC_000022.11"),
            ("chrX", "NC_000023.11"),
            ("chrY", "NC_000024.10"),
            ("chrM", "NC_012920.1"),
        ];
        for (ucsc, expected_refseq) in expected {
            assert_eq!(
                ucsc_to_refseq(Assembly::GRCh38, ucsc),
                expected_refseq,
                "GRCh38 inverse mapping failed for {ucsc}"
            );
        }
    }

    #[test]
    fn ucsc_to_refseq_maps_all_standard_chromosomes_grch37() {
        let expected = [
            ("chr1", "NC_000001.10"),
            ("chr17", "NC_000017.10"),
            ("chrX", "NC_000023.10"),
            ("chrY", "NC_000024.9"),
            ("chrM", "NC_012920.1"),
        ];
        for (ucsc, expected_refseq) in expected {
            assert_eq!(
                ucsc_to_refseq(Assembly::GRCh37, ucsc),
                expected_refseq,
                "GRCh37 inverse mapping failed for {ucsc}"
            );
        }
    }

    #[test]
    fn ucsc_to_refseq_passes_through_patches_and_unknowns() {
        // UCSC-style patch names must pass through unchanged — the reverse
        // FastaReader alias table handles them at runtime, not this const
        // table.
        assert_eq!(
            ucsc_to_refseq(Assembly::GRCh38, "chr9_KN196479v1_fix"),
            "chr9_KN196479v1_fix"
        );
        assert_eq!(
            ucsc_to_refseq(Assembly::GRCh38, "chr22_KI270879v1_alt"),
            "chr22_KI270879v1_alt"
        );
        assert_eq!(
            ucsc_to_refseq(Assembly::GRCh38, "chrUn_KI270302v1"),
            "chrUn_KI270302v1"
        );
        // A RefSeq accession should round-trip to itself (not re-translate).
        assert_eq!(
            ucsc_to_refseq(Assembly::GRCh38, "NC_000001.11"),
            "NC_000001.11"
        );
        assert_eq!(
            ucsc_to_refseq(Assembly::GRCh38, "NW_025791820.1"),
            "NW_025791820.1"
        );
        // Arbitrary strings — treated as pass-through like refseq_to_ucsc.
        assert_eq!(ucsc_to_refseq(Assembly::GRCh38, "gibberish"), "gibberish");
        assert_eq!(ucsc_to_refseq(Assembly::GRCh38, ""), "");
    }

    #[test]
    fn ucsc_refseq_round_trip_identity_on_primary_chromosomes() {
        // `refseq_to_ucsc(ucsc_to_refseq(x)) == x` for every UCSC primary
        // chrom on each assembly. Guard against either function silently
        // dropping an entry from its table.
        for assembly in [Assembly::GRCh38, Assembly::GRCh37] {
            for (_refseq, ucsc) in table_for(assembly) {
                assert_eq!(
                    refseq_to_ucsc(assembly, ucsc_to_refseq(assembly, ucsc)),
                    *ucsc,
                    "round-trip failed for {ucsc} on {assembly}"
                );
            }
        }
    }

    #[test]
    fn is_patch_sequence_detects_prefixes() {
        assert!(is_patch_sequence("NW_025791820.1"));
        assert!(is_patch_sequence("NT_187633.1"));
        assert!(!is_patch_sequence("NC_000001.11"));
        assert!(!is_patch_sequence("chr1"));
        assert!(!is_patch_sequence(""));
    }

    #[test]
    fn is_patch_sequence_detects_ucsc_alt_fix_contigs() {
        // UCSC-style contig names as they appear in column 1 of the real
        // MANE v1.5 GFF3 file.
        assert!(is_patch_sequence("chr9_KN196479v1_fix"));
        assert!(is_patch_sequence("chr22_KI270879v1_alt"));
        assert!(is_patch_sequence("chr1_KI270706v1_random"));
        assert!(is_patch_sequence("chrUn_KI270302v1"));
        // Primary chromosomes must NOT be flagged.
        assert!(!is_patch_sequence("chr1"));
        assert!(!is_patch_sequence("chrX"));
        assert!(!is_patch_sequence("chrY"));
        assert!(!is_patch_sequence("chrM"));
    }

    #[test]
    fn is_patch_sequence_detects_hg19_haplotype_contigs() {
        // hg19 MHC region — 7 alternate haplotypes named with the
        // `_hap1`..`_hap7` suffix.
        assert!(is_patch_sequence("chr6_apd_hap1"));
        assert!(is_patch_sequence("chr6_cox_hap2"));
        assert!(is_patch_sequence("chr6_dbb_hap3"));
        assert!(is_patch_sequence("chr6_mann_hap4"));
        assert!(is_patch_sequence("chr6_mcf_hap5"));
        assert!(is_patch_sequence("chr6_qbl_hap6"));
        assert!(is_patch_sequence("chr6_ssto_hap7"));
        // hg19 also carries `_hap1` for chr4 (ctg9) and chr17 (ctg5).
        assert!(is_patch_sequence("chr4_ctg9_hap1"));
        assert!(is_patch_sequence("chr17_ctg5_hap1"));
        // hg19 random contigs use the `gl\d+_random` form already covered
        // by the existing `_random` suffix; lock the behavior in.
        assert!(is_patch_sequence("chr1_gl000191_random"));
        assert!(is_patch_sequence("chr17_gl000204_random"));
        // hg19 unplaced scaffolds use `chrUn_gl\d+`, covered by the
        // `chrUn_` prefix.
        assert!(is_patch_sequence("chrUn_gl000211"));
        assert!(is_patch_sequence("chrUn_gl000247"));
        // `_hap8` and beyond don't exist on hg19; reject so future
        // accidental input doesn't slip through.
        assert!(!is_patch_sequence("chr6_hypothetical_hap8"));
        assert!(!is_patch_sequence("chr6_hypothetical_hap0"));
        // Primary chroms still negative.
        assert!(!is_patch_sequence("chr6"));
    }

    #[test]
    fn assembly_from_str_parses_canonical_and_lowercase() {
        assert_eq!("GRCh38".parse::<Assembly>().unwrap(), Assembly::GRCh38);
        assert_eq!("grch38".parse::<Assembly>().unwrap(), Assembly::GRCh38);
        assert_eq!("GRCh37".parse::<Assembly>().unwrap(), Assembly::GRCh37);
        assert_eq!("grch37".parse::<Assembly>().unwrap(), Assembly::GRCh37);
    }

    #[test]
    fn assembly_from_str_rejects_hg19_and_hg38_aliases() {
        // Rationale: UCSC `hg19` chrM (NC_001807) ≠ GRCh37 chrMT (NC_012920.1).
        // Silently accepting the alias would mis-annotate every chrM variant
        // for users who genotyped against the actual UCSC distribution.
        let err = "hg19".parse::<Assembly>().unwrap_err();
        match err {
            VarEffectError::UnsupportedAssemblyAlias { input } => assert_eq!(input, "hg19"),
            other => panic!("expected UnsupportedAssemblyAlias, got {other:?}"),
        }
        let err = "hg38".parse::<Assembly>().unwrap_err();
        match err {
            VarEffectError::UnsupportedAssemblyAlias { input } => assert_eq!(input, "hg38"),
            other => panic!("expected UnsupportedAssemblyAlias, got {other:?}"),
        }
    }

    #[test]
    fn assembly_from_str_rejects_garbage() {
        assert!("nonsense".parse::<Assembly>().is_err());
        assert!("".parse::<Assembly>().is_err());
    }

    #[test]
    fn assembly_display_round_trips_via_from_str() {
        for assembly in [Assembly::GRCh38, Assembly::GRCh37] {
            let s = assembly.to_string();
            assert_eq!(s.parse::<Assembly>().unwrap(), assembly);
        }
    }
}
