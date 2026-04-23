//! Error type for `vareffect`.
//!
//! Covers the runtime failure modes of every subsystem: transcript store
//! load, reference genome access, variant localization, HGVS parsing, and
//! consequence assignment.

use std::path::PathBuf;

/// Errors returned by `vareffect` APIs.
#[derive(Debug, thiserror::Error)]
pub enum VarEffectError {
    /// Failed to read a file from disk (transcript store, genome binary, or
    /// `.bin.idx` index). The path field identifies which file failed.
    #[error("failed to read file at {path}: {source}")]
    Io {
        /// Path the caller attempted to read.
        path: PathBuf,
        /// Underlying I/O error.
        #[source]
        source: std::io::Error,
    },

    /// Failed to deserialize the MessagePack payload into
    /// [`crate::types::TranscriptModel`] records.
    #[error("failed to deserialize transcript store: {0}")]
    Deserialize(#[from] rmp_serde::decode::Error),

    /// The deserialized data violated an invariant (e.g. `tx_end <= tx_start`
    /// or a coordinate that exceeds `i32::MAX`). The inner string describes
    /// the offending record.
    #[error("transcript store is malformed: {0}")]
    Malformed(String),

    /// A lookup referenced a chromosome not present in the genome index.
    ///
    /// The chromosome name is passed back verbatim (UCSC-style, e.g. `chrZZ`)
    /// so the caller can surface it in a diagnostic.
    #[error("chromosome '{chrom}' not found in genome index")]
    ChromNotFound {
        /// UCSC-style chromosome name that was requested.
        chrom: String,
    },

    /// A lookup used a coordinate outside the chromosome's valid range,
    /// either because `start >= end` or because `end` exceeded the chromosome
    /// length recorded in the genome index.
    #[error("coordinate out of range: {chrom}:{start}-{end} (chrom length: {chrom_len})")]
    CoordinateOutOfRange {
        /// UCSC-style chromosome name.
        chrom: String,
        /// Requested 0-based inclusive start.
        start: u64,
        /// Requested 0-based exclusive end.
        end: u64,
        /// Chromosome length as recorded in the `.fai` index.
        chrom_len: u64,
    },

    /// The `.bin.idx` index file was not found alongside the genome binary.
    /// Building the index is the responsibility of `vareffect-cli setup`,
    /// not the runtime reader — the reader fails fast rather than silently
    /// indexing.
    #[error("genome index not found at {path}")]
    IndexNotFound {
        /// Path where the `.bin.idx` was expected.
        path: String,
    },

    /// The VCF REF allele does not match the reference FASTA at the given
    /// position. This typically indicates a VCF built against a different
    /// assembly, a liftover error, or a corrupt input file.
    ///
    /// For SNVs the strings are single characters; for indels they may be
    /// multi-base (e.g. `"ACGT"` vs `"ACGG"`).
    ///
    /// For HGVS-input mismatches (where the REF is stated inside the HGVS
    /// notation itself) see [`VarEffectError::HgvsRefMismatch`], which
    /// includes the HGVS input string in its payload so callers can
    /// reproduce VEP-style diagnostics.
    #[error("REF allele mismatch at {chrom}:{pos}: genome has {expected}, VCF has {got}")]
    RefMismatch {
        /// UCSC-style chromosome name.
        chrom: String,
        /// 0-based genomic position.
        pos: u64,
        /// Sequence found in the reference genome (uppercase ASCII).
        expected: String,
        /// Sequence provided by the caller (from VCF REF column).
        got: String,
    },

    /// The HGVS-stated REF base disagrees with the reference FASTA at the
    /// resolved genomic position. Only raised from
    /// [`crate::VarEffect::resolve_hgvs_c`] /
    /// [`crate::VarEffect::resolve_hgvs_c_with_meta`]; for VCF-driven
    /// mismatches see [`VarEffectError::RefMismatch`].
    ///
    /// The `Display` format matches Ensembl VEP's wording, so
    /// `eprintln!("{err}")` produces a VEP-concordant diagnostic without
    /// further formatting.
    #[error(
        "ref allele mismatch at position {pos} for '{hgvs}': \
         genome has '{expected}', HGVS states '{got}'"
    )]
    HgvsRefMismatch {
        /// Full HGVS notation the caller passed in
        /// (e.g. `"NM_000546.6:c.738C>T"`).
        hgvs: String,
        /// UCSC-style chromosome name of the resolved position.
        chrom: String,
        /// 0-based genomic position where the disagreement was detected.
        pos: u64,
        /// Plus-strand reference base from the FASTA (uppercase ASCII).
        expected: String,
        /// HGVS-stated ref base projected to the plus strand (uppercase
        /// ASCII). For minus-strand transcripts this is the complement of
        /// the base as typed in the HGVS string.
        got: String,
    },

    /// The requested transcript accession was not found in the
    /// [`crate::TranscriptStore`].
    #[error("transcript not found: {accession}")]
    TranscriptNotFound {
        /// Full accession string as provided by the caller.
        accession: String,
    },

    /// An HGVS c. position exceeds the transcript boundaries (CDS length,
    /// UTR exonic extent, or intron bounds).
    #[error("HGVS position out of range for {accession}: {detail}")]
    PositionOutOfRange {
        /// Transcript accession for context.
        accession: String,
        /// Human-readable description of the violation.
        detail: String,
    },

    /// An HGVS c. notation string could not be parsed into structured
    /// components. The inner string describes what went wrong.
    #[error("failed to parse HGVS notation: {0}")]
    HgvsParseError(String),

    /// An allele string contained non-ASCII bytes after normalization.
    ///
    /// Should not occur with valid genomic input (ACGTN are ASCII).
    /// Surfaced instead of panicking per library error-handling policy.
    #[error("allele contains invalid (non-ASCII) bytes")]
    InvalidAllele,
}
