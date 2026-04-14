#![doc = include_str!("../README.md")]
//! `vareffect` — Variant consequence prediction and HGVS notation, targeting
//! near-100% concordance with Ensembl VEP (release 115/116).
//!
//! Consumers point the store loaders at whatever transcript and genome files
//! their build pipeline produces — `vareffect` ships no embedded reference
//! data and has no runtime dependency on an orchestrator CLI.
//!
//! # Transcript model store
//!
//! An in-memory store of MANE transcript models indexed by genomic interval
//! for O(log n + k) overlap queries. Each [`TranscriptModel`] carries
//! per-exon [`CdsSegment`]s with the GFF3 column-8 phase captured, so
//! downstream codon walks and frameshift detection don't have to re-derive
//! phase from scratch.
//!
//! ```no_run
//! use std::path::Path;
//! use vareffect::{Biotype, TranscriptStore};
//!
//! let store = TranscriptStore::load_from_path(
//!     Path::new("data/vareffect/transcript_models.bin"),
//! )?;
//!
//! // Overlap query: all transcripts whose tx_start..tx_end intersects the interval.
//! for (tx, _idx) in store.query_overlap("chr6", 33_409_450, 33_409_451) {
//!     println!(
//!         "{} ({}): cds [{:?}, {:?}), {} segments, biotype={:?}",
//!         tx.accession,
//!         tx.gene_symbol,
//!         tx.cds_genomic_start,
//!         tx.cds_genomic_end,
//!         tx.cds_segments.len(),
//!         tx.biotype,
//!     );
//!
//!     // Walk CDS segments in transcript 5'→3' order (reversed for minus strand):
//!     for seg in &tx.cds_segments {
//!         println!(
//!             "  segment in exon[{}], phase {}: [{}, {})",
//!             seg.exon_index, seg.phase, seg.genomic_start, seg.genomic_end
//!         );
//!     }
//! }
//!
//! // `biotype` is an enum with `Other(String)` for unknown upstream labels.
//! let total_protein_coding = store
//!     .transcripts()
//!     .iter()
//!     .filter(|t| matches!(t.biotype, Biotype::ProteinCoding))
//!     .count();
//! # let _ = total_protein_coding;
//! # Ok::<(), vareffect::VarEffectError>(())
//! ```
//!
//! # Reference genome reader
//!
//! Memory-mapped random access to the reference genome via [`FastaReader`].
//! Pair it with `TranscriptStore` to extract codons, verify REF alleles, and
//! walk downstream for frameshift termination. See the `fasta` module for
//! the on-disk format, coordinate conventions, and chromosome-name handling.
//!
//! The flat binary format stores uppercase IUPAC nucleotide codes, matching
//! GA4GH refget v2.0 conventions. Most bases are `A`/`C`/`G`/`T`/`N`; the
//! NCBI GRCh38.p14 assembly also uses ambiguity codes (`M`, `R`, `Y`, etc.)
//! in some patch-scaffold regions. Soft-mask information is not preserved.
//!
//! # Variant consequence assignment
//!
//! [`VarEffect::annotate`] takes a variant's position and alleles, locates it
//! within every overlapping transcript, extracts the reference codon(s) from
//! FASTA, translates ref and alt codons, and assigns SO consequence term(s)
//! with VEP-concordant IMPACT ratings. The [`codon`] module provides the
//! standard and mitochondrial genetic code translation tables.
//!
//! ```no_run
//! use std::path::Path;
//! use vareffect::VarEffect;
//!
//! let ve = VarEffect::open(
//!     Path::new("data/vareffect/transcript_models.bin"),
//!     Path::new("data/vareffect/GRCh38.bin"),
//! )?;
//!
//! // Annotate TP53 c.742C>T (p.R248W) — chr17, 0-based position 7,674,219.
//! let results = ve.annotate("chr17", 7_674_219, b"C", b"T")?;
//! for r in &results {
//!     for csq in &r.consequences {
//!         println!("{} ({})", csq.as_str(), r.impact);
//!     }
//! }
//! # Ok::<(), vareffect::VarEffectError>(())
//! ```
//!
//! For lower-level building blocks (per-transcript annotation when you
//! already hold a `&TranscriptModel`), see [`annotate_snv`],
//! [`annotate_deletion`], and [`annotate_insertion`].
//!
//! # Coordinate convention
//!
//! All coordinates in [`TranscriptModel`] are **0-based, half-open** (BED/UCSC
//! style). GFF3 input (1-based, fully-closed) is converted at build time by
//! `vareffect-cli`. See [`transcript`] for the interval-tree indexing details.
//!
//! `cds_genomic_start` / `cds_genomic_end` are the genomic `min` / `max`
//! coordinates across all CDS segments, **not** transcript-relative. For a
//! minus-strand gene, `cds_genomic_start` is biologically the 3' end of the
//! protein in transcript order. Walk `cds_segments` (ordered 5'→3' on the
//! transcript) when you need the true coding walk.
//!
//! # Thread safety
//!
//! Both [`TranscriptStore`] and [`FastaReader`] are `Send + Sync` (proven by
//! a compile-time assertion at the bottom of this file). `TranscriptStore`
//! is lock-free for reads. `FastaReader` is backed by a memory-mapped
//! `&[u8]` — inherently `Send + Sync` with zero contention. All threads
//! can read from the same `FastaReader` concurrently without cloning.

#![warn(missing_docs)]

pub mod chrom;
pub mod codon;
pub mod consequence;
pub mod error;
pub mod fasta;
pub mod hgvs_c;
pub mod hgvs_p;
pub mod hgvs_reverse;
pub mod locate;
mod normalize;
pub mod transcript;
pub mod types;
mod var_effect;
mod vep_json;

#[cfg(test)]
pub(crate) mod test_fixtures;

pub use consequence::{
    Consequence, ConsequenceResult, Impact, annotate_deletion, annotate_insertion, annotate_snv,
};
pub use error::VarEffectError;
pub use fasta::FastaReader;
pub use hgvs_reverse::GenomicVariant;
pub use locate::{
    IndelLocation, IndelRegion, LocateIndex, SpliceOverlapDetail, SpliceSide, VariantLocation,
    locate_indel, locate_variant,
};
pub use transcript::TranscriptStore;
pub use types::{Biotype, CdsSegment, Exon, Strand, TranscriptModel, TranscriptTier};
pub use var_effect::VarEffect;

// Compile-time proof that the runtime types are thread-safe to share across
// worker tasks. A field that silently breaks `Send + Sync` in a future
// refactor will fail this check at build time instead of at deploy time.
const _: fn() = || {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<TranscriptStore>();
    assert_send_sync::<FastaReader>();
    assert_send_sync::<VarEffect>();
};
