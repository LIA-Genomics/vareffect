//! Stateful entrypoint that bundles a [`TranscriptStore`] and [`FastaReader`].
//!
//! Most callers want to annotate variants without manually threading both
//! handles through every call. [`VarEffect`] holds an owned copy of each
//! store and exposes the high-level API as methods. Construct one at startup,
//! wrap it in `Arc`, and share it across worker tasks.
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
//! // Annotate TP53 c.742C>T (p.R248W) on chr17 (0-based position).
//! let results = ve.annotate("chr17", 7_674_219, b"C", b"T")?;
//! for csq in results.iter().flat_map(|r| r.consequences.iter()) {
//!     println!("{}", csq.as_str());
//! }
//! # Ok::<(), vareffect::VarEffectError>(())
//! ```

use std::path::Path;

use crate::consequence::{ConsequenceResult, SvConsequenceResult, SvKind};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_reverse::{GenomicVariant, ResolvedHgvsC};
use crate::locate::LocateIndex;
use crate::transcript::TranscriptStore;
use crate::types::TranscriptModel;

/// Stateful entrypoint to vareffect: bundles a [`TranscriptStore`] and a
/// [`FastaReader`] so callers don't have to thread both handles through
/// every annotation call.
///
/// `VarEffect` is `Send + Sync` (both inner stores are mmap- or `Arc`-backed
/// with no interior mutability), so a single `Arc<VarEffect>` can be shared
/// across all async tasks with zero contention. Cloning the underlying
/// stores is cheap (they're already `Arc`-backed internally), but the
/// recommended sharing pattern is `Arc<VarEffect>` for one ownership unit.
///
/// `VarEffect` does not derive `Clone` on purpose — sharing should go
/// through `Arc<VarEffect>`, not field-level cloning.
pub struct VarEffect {
    transcripts: TranscriptStore,
    fasta: FastaReader,
}

impl VarEffect {
    /// Construct a `VarEffect` from previously loaded stores.
    ///
    /// Use this when your build pipeline produces the stores in-memory or
    /// when you need to share existing handles. For the common path-based
    /// case, prefer [`VarEffect::open`].
    pub fn new(transcripts: TranscriptStore, fasta: FastaReader) -> Self {
        Self { transcripts, fasta }
    }

    /// Load both stores from disk and assemble a `VarEffect`.
    ///
    /// `transcript_models_path` should point at the MessagePack file produced
    /// by your transcript build pipeline (typically `transcript_models.bin`).
    /// `genome_path` should point at the flat-binary genome with its
    /// `.bin.idx` sidecar (typically `GRCh38.bin` + `GRCh38.bin.idx`).
    ///
    /// # Errors
    ///
    /// Propagates the union of errors from [`TranscriptStore::load_from_path`]
    /// and [`FastaReader::open`] — file I/O failures, malformed payloads, or
    /// version mismatches.
    pub fn open(transcript_models_path: &Path, genome_path: &Path) -> Result<Self, VarEffectError> {
        let transcripts = TranscriptStore::load_from_path(transcript_models_path)?;
        let fasta = FastaReader::open(genome_path)?;
        Ok(Self::new(transcripts, fasta))
    }

    /// Same as [`VarEffect::open`], plus an optional patch-contig alias CSV
    /// for NCBI RefSeq patch contig name resolution.
    ///
    /// `patch_aliases_csv` is a path to a `refseq,ucsc` CSV produced by
    /// `vareffect-cli setup`. When supplied *and* the genome binary uses
    /// NCBI naming, the reader loads it into a UCSC -> RefSeq map for patch
    /// contig lookups. Pass `None` if you only need primary chromosomes
    /// (`chr1`..`chrM`).
    ///
    /// # Errors
    ///
    /// Same as [`VarEffect::open`], plus errors from
    /// [`FastaReader::open_with_patch_aliases`] if the alias CSV is
    /// malformed.
    pub fn open_with_patch_aliases(
        transcript_models_path: &Path,
        genome_path: &Path,
        patch_aliases_csv: Option<&Path>,
    ) -> Result<Self, VarEffectError> {
        let transcripts = TranscriptStore::load_from_path(transcript_models_path)?;
        let fasta = FastaReader::open_with_patch_aliases(genome_path, patch_aliases_csv)?;
        Ok(Self::new(transcripts, fasta))
    }

    // -----------------------------------------------------------------
    // Variant annotation
    // -----------------------------------------------------------------

    /// Annotate a variant against every overlapping transcript.
    ///
    /// Returns one [`ConsequenceResult`] per overlapping transcript. An
    /// empty `Vec` means no transcript overlapped the variant locus
    /// (intergenic). The full VCF REF allele is verified against the FASTA
    /// before trimming, so a stated REF that doesn't match the genome
    /// produces [`VarEffectError::RefMismatch`].
    ///
    /// # Arguments
    ///
    /// * `chrom` — UCSC-style chromosome name (`"chr17"`, `"chrM"`).
    /// * `pos` — 0-based genomic position (BED convention).
    /// * `ref_allele` — VCF REF allele bytes (uppercase ASCII).
    /// * `alt_allele` — VCF ALT allele bytes (uppercase ASCII).
    ///
    /// # Errors
    ///
    /// See [`crate::consequence`] for the full error taxonomy:
    /// `RefMismatch`, `ChromNotFound`, `CoordinateOutOfRange`, `Malformed`.
    pub fn annotate(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> Result<Vec<ConsequenceResult>, VarEffectError> {
        crate::consequence::annotate(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            &self.transcripts,
            &self.fasta,
        )
    }

    /// Annotate a variant and return the result as VEP REST-compatible JSON.
    ///
    /// Convenience wrapper around [`Self::annotate`] that serializes the
    /// output as a single-element JSON array matching Ensembl VEP REST's
    /// response shape for
    /// `GET /vep/human/region/{chrom}:{pos}:{pos}/{alt}?refseq=1&hgvs=1`.
    ///
    /// Tools written against real VEP REST can consume the returned value
    /// directly -- the shape is a single-element `[{...}]` array whose
    /// element carries `seq_region_name`, `start`, `end`, `allele_string`,
    /// `assembly_name`, `most_severe_consequence`, and
    /// `transcript_consequences`.
    ///
    /// vareffect is a strict subset of VEP -- it does NOT populate
    /// `canonical`, `swissprot`, `protein_id`, `sift_prediction`,
    /// `polyphen_prediction`, `revel_score`, `alphamissense`, `spliceai`, or
    /// `colocated_variants[*].id`. Consumers must tolerate their absence.
    ///
    /// # Arguments
    ///
    /// * `chrom` -- Chromosome in vareffect format (`"chr17"`, `"chrX"`,
    ///   plain `"17"`, or any RefSeq/UCSC-style name the loaded FASTA
    ///   understands). Any leading `"chr"` is stripped from the emitted
    ///   `seq_region_name`.
    /// * `pos` -- 0-based genomic start position (BED convention). The
    ///   emitted `start` and `end` are converted to 1-based coordinates
    ///   (VEP convention).
    /// * `ref_allele` -- Plus-strand reference allele bytes. Verified
    ///   against the loaded FASTA before annotation.
    /// * `alt_allele` -- Plus-strand alternate allele bytes.
    /// * `assembly` -- Genome build label (e.g. `"GRCh38"`) written to the
    ///   top-level `assembly_name` field. vareffect has no build awareness;
    ///   the caller states which FASTA was loaded.
    ///
    /// # Returns
    ///
    /// `serde_json::Value` shaped as `[{ ... }]`. The array always has
    /// exactly one element. Intergenic variants (no overlapping transcripts)
    /// still produce a well-formed element whose `transcript_consequences`
    /// is `[]` and whose `most_severe_consequence` is `"intergenic_variant"`.
    ///
    /// # Errors
    ///
    /// Propagates any [`VarEffectError`] raised by [`Self::annotate`]:
    /// `RefMismatch`, `ChromNotFound`, `CoordinateOutOfRange`, `Malformed`.
    /// Serialization itself is infallible.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use std::path::Path;
    /// use vareffect::VarEffect;
    ///
    /// let ve = VarEffect::open(
    ///     Path::new("data/vareffect/transcript_models.bin"),
    ///     Path::new("data/vareffect/GRCh38.bin"),
    /// )?;
    ///
    /// // TP53 p.R248W -- 0-based chr17:7_674_219 C>T
    /// let json = ve.annotate_to_vep_json("chr17", 7_674_219, b"C", b"T", "GRCh38")?;
    ///
    /// // The returned value is a single-element array, same as VEP REST.
    /// let top = &json[0];
    /// assert_eq!(top["seq_region_name"], "17");
    /// assert_eq!(top["assembly_name"], "GRCh38");
    /// # Ok::<(), vareffect::VarEffectError>(())
    /// ```
    pub fn annotate_to_vep_json(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
        assembly: &str,
    ) -> Result<serde_json::Value, VarEffectError> {
        let results = self.annotate(chrom, pos, ref_allele, alt_allele)?;
        Ok(crate::vep_json::to_vep_json_array(
            chrom, pos, ref_allele, alt_allele, assembly, &results,
        ))
    }

    /// Annotate an interval structural variant (DEL / DUP / INV) against every
    /// overlapping transcript.
    ///
    /// Unlike [`Self::annotate`] (single ref/alt pair), this consumes a genomic
    /// *interval* plus a direction and emits the per-transcript SO term **set**
    /// VEP assigns to structural variants — the overlapped sub-regions plus, for
    /// unbalanced events, a copy-number headline term:
    ///
    /// | direction | spans whole transcript | partial overlap |
    /// |-----------|------------------------|-----------------|
    /// | deletion ([`SvKind::Deletion`]) | `transcript_ablation` (alone) | `feature_truncation` + sub-regions |
    /// | duplication ([`SvKind::Duplication`]) | `transcript_amplification` (alone) | `feature_elongation` + sub-regions |
    /// | inversion ([`SvKind::Inversion`]) | sub-regions only | sub-regions only |
    ///
    /// An inversion is copy-number-neutral, so VEP assigns **no** ablation /
    /// truncation headline — only the overlapped sub-region terms
    /// (`coding_sequence_variant`, `5_prime_UTR_variant`, `3_prime_UTR_variant`,
    /// `intron_variant`, `non_coding_transcript_exon_variant`). See
    /// Ensembl/ensembl-vep#79.
    ///
    /// Each [`SvConsequenceResult`] also carries per-transcript exon-overlap
    /// geometry (5'/3'/last-exon coverage, exon count, CDS bases affected,
    /// overlap fraction) so callers can drive region-centric ACMG/ClinGen CNV
    /// scoring without re-deriving exon math from the raw [`TranscriptModel`].
    ///
    /// # Arguments
    ///
    /// * `chrom` — UCSC-style chromosome name (`"chr17"`, `"chrX"`).
    /// * `start` — 0-based inclusive start of the SV footprint (BED convention).
    /// * `end` — 0-based exclusive end of the SV footprint.
    /// * `sv_kind` — interval SV shape ([`SvKind`]).
    ///
    /// # Returns
    ///
    /// One result per overlapping transcript. An **empty** vec is the valid
    /// "no transcript overlapped" answer — it is never conflated with an
    /// error, so a downstream CNV scorer can trust an empty result.
    ///
    /// # Errors
    ///
    /// Returns [`VarEffectError::Malformed`] only when `start >= end` (an
    /// inverted/zero-length interval is a caller bug, not bad reference data).
    pub fn annotate_interval(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        sv_kind: SvKind,
    ) -> Result<Vec<SvConsequenceResult>, VarEffectError> {
        crate::consequence::annotate_interval(chrom, start, end, sv_kind, &self.transcripts)
    }

    /// Annotate a single breakend (`<BND>`) breakpoint against every transcript
    /// it falls in.
    ///
    /// The breakpoint at 0-based genomic `pos` truncates any transcript it lands
    /// inside (`feature_truncation`), reported alongside the sub-region term at
    /// that base. Nearby transcripts the breakpoint does not enter are not
    /// returned (the store overlap query is exact, like [`Self::annotate`]).
    ///
    /// A translocation is two linked breakends; annotate each end with this
    /// method (or [`Self::annotate_breakend_pair`]). Mate pairing (MATEID /
    /// EVENT) is intentionally not performed — this matches VEP, which annotates
    /// each breakend independently. Use [`crate::Breakend::parse`] to extract a
    /// mate locus from a VCF breakend ALT string.
    ///
    /// Returns one [`SvConsequenceResult`] per overlapping transcript; an empty
    /// vec means the breakpoint fell in no transcript.
    pub fn annotate_breakend(&self, chrom: &str, pos: u64) -> Vec<SvConsequenceResult> {
        crate::consequence::annotate_breakend(chrom, pos, &self.transcripts)
    }

    /// Annotate both breakpoints of a breakend pair (a translocation or other
    /// two-locus rearrangement) and return the union of their transcript
    /// annotations.
    ///
    /// Each breakpoint is annotated independently via [`Self::annotate_breakend`];
    /// the two result lists are concatenated (not deduplicated across ends — the
    /// two loci are normally far apart or on different chromosomes). This mirrors
    /// VEP, which reports every transcript affected by either breakend.
    pub fn annotate_breakend_pair(
        &self,
        chrom_a: &str,
        pos_a: u64,
        chrom_b: &str,
        pos_b: u64,
    ) -> Vec<SvConsequenceResult> {
        let mut results = self.annotate_breakend(chrom_a, pos_a);
        results.extend(self.annotate_breakend(chrom_b, pos_b));
        results
    }

    /// Annotate a symbolic insertion (`<INS>`) at 0-based genomic `pos`.
    ///
    /// A symbolic insertion carries no copy-number headline — the inserted
    /// sequence is unknown to the annotator, so VEP reports only the sub-region
    /// the insertion point occupies (e.g. `intron_variant`,
    /// `coding_sequence_variant`). Purely positional.
    ///
    /// Returns one [`SvConsequenceResult`] per overlapping transcript.
    pub fn annotate_sv_insertion(&self, chrom: &str, pos: u64) -> Vec<SvConsequenceResult> {
        crate::consequence::annotate_sv_insertion(chrom, pos, &self.transcripts)
    }

    /// Resolve HGVS c. notation to plus-strand genomic coordinates.
    ///
    /// Returns a [`GenomicVariant`] with 0-based position, UCSC-style
    /// chromosome, and uppercase plus-strand alleles. Pair with
    /// [`VarEffect::annotate`] when you also want consequences:
    ///
    /// ```no_run
    /// # use vareffect::VarEffect;
    /// # let ve: VarEffect = unimplemented!();
    /// let gv = ve.resolve_hgvs_c("NM_000546.6:c.742C>T")?;
    /// let results = ve.annotate(&gv.chrom, gv.pos, &gv.ref_allele, &gv.alt_allele)?;
    /// # Ok::<(), vareffect::VarEffectError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// `HgvsParseError`, `TranscriptNotFound`, `PositionOutOfRange`,
    /// `HgvsRefMismatch`, `ChromNotFound`, `CoordinateOutOfRange`,
    /// `Malformed`.
    pub fn resolve_hgvs_c(&self, hgvs: &str) -> Result<GenomicVariant, VarEffectError> {
        crate::hgvs_reverse::resolve_hgvs_c(hgvs, &self.transcripts, &self.fasta)
    }

    /// Resolve HGVS c. notation and report which transcript version was
    /// actually used.
    ///
    /// Behaves like [`VarEffect::resolve_hgvs_c`] but returns a
    /// [`ResolvedHgvsC`] carrying both the genomic coordinates and the
    /// accession (with version) that was matched in the transcript store.
    /// When the store lacks the caller-specified version, the store's
    /// version-tolerant lookup falls through to the highest available
    /// version of the same base accession; comparing `resolved_accession`
    /// against the input lets callers surface a transcript-version-drift
    /// warning to end users.
    ///
    /// ```no_run
    /// # use vareffect::VarEffect;
    /// # let ve: VarEffect = unimplemented!();
    /// let r = ve.resolve_hgvs_c_with_meta("NM_000546.5:c.742C>T")?;
    /// // If the store only has `.6`, `r.resolved_accession == "NM_000546.6"`.
    /// # Ok::<(), vareffect::VarEffectError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Same as [`VarEffect::resolve_hgvs_c`].
    pub fn resolve_hgvs_c_with_meta(&self, hgvs: &str) -> Result<ResolvedHgvsC, VarEffectError> {
        crate::hgvs_reverse::resolve_hgvs_c_with_meta(hgvs, &self.transcripts, &self.fasta)
    }

    // -----------------------------------------------------------------
    // FastaReader forwarders
    // -----------------------------------------------------------------

    /// Fetch a single base at the given 0-based position. ~5 ns with the
    /// mmap backend. See [`FastaReader::fetch_base`].
    pub fn fetch_base(&self, chrom: &str, pos: u64) -> Result<u8, VarEffectError> {
        self.fasta.fetch_base(chrom, pos)
    }

    /// Fetch a genomic sequence as uppercase ASCII bytes. Coordinates are
    /// 0-based half-open `[start, end)`. See [`FastaReader::fetch_sequence`].
    pub fn fetch_sequence(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, VarEffectError> {
        self.fasta.fetch_sequence(chrom, start, end)
    }

    /// Verify that the reference allele at a position matches the genome.
    /// Zero-copy, case-insensitive. See [`FastaReader::verify_ref`].
    pub fn verify_ref(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
    ) -> Result<bool, VarEffectError> {
        self.fasta.verify_ref(chrom, pos, ref_allele)
    }

    /// Length of `chrom` in bases, or `None` if not present in the index.
    /// See [`FastaReader::chrom_length`].
    pub fn chrom_length(&self, chrom: &str) -> Option<u64> {
        self.fasta.chrom_length(chrom)
    }

    // -----------------------------------------------------------------
    // FASTA-driven helpers
    // -----------------------------------------------------------------

    /// Anchor-prepend HGVS-style indel alleles to VCF form using the loaded
    /// reference genome.
    ///
    /// VEP and other HGVS-first sources represent indels with a `"-"`
    /// placeholder. VCF convention requires the nucleotide immediately 5'
    /// of the event on both alleles, with the position shifted back by one:
    ///
    /// - deletion  (`ref="TG", alt="-"`): `pos -= 1`; `ref = anchor + "TG"`;
    ///   `alt = anchor`
    /// - insertion (`ref="-",  alt="C"`): `pos -= 1`; `ref = anchor`;
    ///   `alt = anchor + "C"`
    ///
    /// For VEP specifically: a pure insertion is reported as
    /// `start = N, end = N - 1`, with `pos` stored as `start`. The single
    /// formula `anchor_pos_0 = pos_1based - 2` applies uniformly to both
    /// deletion and insertion cases (verified against real fixtures for
    /// `NM_006772.2:c.1861_1862del` and `NM_007294.4:c.5266dupC`).
    ///
    /// SNVs, MNVs, and complex substitutions (neither allele `"-"`) pass
    /// through unchanged and return `Ok(None)`.
    ///
    /// # Arguments
    ///
    /// * `chrom` — Chromosome name accepted by the FASTA reader. The
    ///   reader's alias table handles UCSC (`"chr17"`), bare Ensembl
    ///   (`"17"`), and NCBI RefSeq (`"NC_000017.11"`) transparently.
    /// * `pos_1based` — 1-based position of the variant. Must be `>= 2`
    ///   for indels (there is no base 5' of position 1).
    /// * `ref_allele` / `alt_allele` — Allele strings; either may be `"-"`
    ///   for a pure deletion or insertion respectively.
    ///
    /// # Returns
    ///
    /// * `Ok(None)` — passthrough (neither allele is `"-"`).
    /// * `Ok(Some((new_pos, new_ref, new_alt)))` — anchor-prepended form,
    ///   `new_pos = pos_1based - 1`.
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::ChromNotFound`] — chromosome not in the index
    ///   (either from `fetch_base` or from the `pos_1based < 2` guard
    ///   fallback when `chrom_length` also returns `None`).
    /// * [`VarEffectError::CoordinateOutOfRange`] — valid chrom whose
    ///   length is exceeded; also emitted directly when `pos_1based < 2`.
    /// * [`VarEffectError::Io`] — propagated from the underlying mmap
    ///   reader.
    ///
    /// # Scope (known limitation)
    ///
    /// This helper performs *only* anchor prepend, not full VCF
    /// normalization. It does NOT left-align indels in repetitive regions
    /// (e.g. homopolymer runs). Downstream matchers (gnomAD, ClinVar) that
    /// have normalized to leftmost representation may still fail to match
    /// for repeat-context indels.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use vareffect::VarEffect;
    /// # let ve: VarEffect = unimplemented!();
    /// // Deletion: VEP reports pos=33409450, alleles "TG/-" on chr6.
    /// let (pos, r, a) = ve
    ///     .anchor_prepend_indel("chr6", 33_409_450, "TG", "-")?
    ///     .expect("deletion should normalize");
    /// assert_eq!((pos, r.as_str(), a.as_str()), (33_409_449, "ATG", "A"));
    /// # Ok::<(), vareffect::VarEffectError>(())
    /// ```
    pub fn anchor_prepend_indel(
        &self,
        chrom: &str,
        pos_1based: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<(u64, String, String)>, VarEffectError> {
        let is_deletion = alt_allele == "-";
        let is_insertion = ref_allele == "-";
        if !is_deletion && !is_insertion {
            return Ok(None);
        }

        // 1-based variant pos → 0-based anchor pos (base immediately 5').
        //   1-based anchor = pos_1based - 1
        //   0-based anchor = pos_1based - 2
        let anchor_pos_0 = match pos_1based.checked_sub(2) {
            Some(p) => p,
            None => {
                // pos_1based is 0 or 1: no base exists 5' of position 1.
                // Construct a CoordinateOutOfRange with the real chrom_len
                // so the error message is diagnostic. If the chrom is not
                // even in the index, surface the more specific ChromNotFound.
                let chrom_len =
                    self.chrom_length(chrom)
                        .ok_or_else(|| VarEffectError::ChromNotFound {
                            chrom: chrom.to_string(),
                        })?;
                return Err(VarEffectError::CoordinateOutOfRange {
                    chrom: chrom.to_string(),
                    start: 0,
                    end: 0,
                    chrom_len,
                });
            }
        };

        // `fetch_base` returns the byte from the mmap'd genome binary,
        // which is stored uppercase (soft-mask is discarded at build time
        // — see `fasta.rs`). No extra uppercase step is required.
        let anchor_byte = self.fetch_base(chrom, anchor_pos_0)?;
        let anchor = anchor_byte as char;
        let new_pos = pos_1based - 1;

        let (new_ref, new_alt) = if is_deletion {
            // ref="TG", alt="-" → ref="ATG", alt="A"
            (format!("{anchor}{ref_allele}"), anchor.to_string())
        } else {
            // ref="-", alt="TG" → ref="A", alt="ATG"
            (anchor.to_string(), format!("{anchor}{alt_allele}"))
        };
        Ok(Some((new_pos, new_ref, new_alt)))
    }

    /// Left-align a VCF-style variant to the leftmost equivalent position.
    ///
    /// Implements the Tan et al. 2015 normalization algorithm used by
    /// `vt normalize` and `bcftools norm`: shift-then-trim produces the
    /// unique left-aligned parsimonious representation.
    ///
    /// SNVs and MNVs pass through unchanged (the shift loop's
    /// rightmost-base comparison exits immediately when alleles differ).
    /// Complex inputs like `ref=ACT, alt=AT` are handled naturally: the
    /// matching rightmost `T` triggers the shift, exposing the underlying
    /// deletion.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name accepted by the FASTA reader (UCSC,
    ///   bare, or RefSeq naming -- the reader's alias table handles
    ///   translation).
    /// * `pos_1based` - 1-based VCF position of the variant.
    /// * `ref_allele` - Reference allele (uppercase ACGTN, no `"-"`
    ///   placeholders).
    /// * `alt_allele` - Alternate allele (uppercase ACGTN, no `"-"`
    ///   placeholders).
    ///
    /// # Returns
    ///
    /// * `Ok(None)` - No normalization needed (SNV, MNV, or already
    ///   leftmost).
    /// * `Ok(Some((pos, ref, alt)))` - Left-aligned representation with
    ///   new 1-based position and parsimonious alleles.
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::ChromNotFound`] - Chromosome not in the FASTA
    ///   index.
    /// * [`VarEffectError::CoordinateOutOfRange`] - Position exceeds
    ///   chromosome length during left-extension.
    /// * [`VarEffectError::InvalidAllele`] - Allele bytes are not valid
    ///   UTF-8 (should not occur with valid genomic input).
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use vareffect::VarEffect;
    /// # let ve: VarEffect = unimplemented!();
    /// // Right-shifted insertion in poly-A -> left-aligned
    /// let result = ve.left_align_indel("chr17", 43045705, "A", "AT")?;
    /// assert!(result.is_some()); // position shifted leftward
    /// # Ok::<(), vareffect::VarEffectError>(())
    /// ```
    pub fn left_align_indel(
        &self,
        chrom: &str,
        pos_1based: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<(u64, String, String)>, VarEffectError> {
        let orig_pos = pos_1based;
        let orig_ref = ref_allele;
        let orig_alt = alt_allele;

        let mut pos = pos_1based;
        let mut r: Vec<u8> = ref_allele.as_bytes().to_vec();
        let mut a: Vec<u8> = alt_allele.as_bytes().to_vec();

        // Step 1: Shift loop (Tan et al. 2015).
        //
        // Compare rightmost bases. While they match, drop them and shift left.
        // When trimming empties an allele, prepend the reference base immediately
        // 5' to BOTH alleles and step the position back by one. Extending both
        // alleles preserves the indel (the ref/alt length difference); extending
        // only the emptied allele would collapse an insertion/deletion into a
        // substitution.
        //
        // Coordinate convention: `fetch_base` takes a 0-based position. After
        // `pos -= 1`, `pos` is the new 1-based position, so the 0-based index for
        // the base AT `pos` is `pos - 1`.
        loop {
            if pos <= 1 || r.is_empty() || a.is_empty() {
                break;
            }
            let trimmed = if r.last() == a.last() {
                r.pop();
                a.pop();
                true
            } else {
                false
            };
            if r.is_empty() || a.is_empty() {
                pos -= 1;
                let prepend = self.fasta.fetch_base(chrom, pos - 1)?;
                r.insert(0, prepend);
                a.insert(0, prepend);
            } else if !trimmed {
                break;
            }
        }

        // Step 2: Left-trim for parsimony.
        //
        // Skips shared leading bases via an index (no O(n) `remove(0)`), keeping
        // at least one base per allele (the VCF anchor requirement).
        let mut prefix_skip = 0usize;
        while r.len() - prefix_skip > 1
            && a.len() - prefix_skip > 1
            && r[prefix_skip] == a[prefix_skip]
        {
            prefix_skip += 1;
            pos += 1;
        }

        // Step 3: Return None if nothing changed, Some if normalized.
        //
        // Alleles contain only ASCII bytes (ACGTN from the FASTA reader or the
        // original input validated upstream). `from_utf8` cannot fail on valid
        // ASCII but we propagate rather than panic.
        let new_ref = std::str::from_utf8(&r[prefix_skip..])
            .map_err(|_| VarEffectError::InvalidAllele)?
            .to_string();
        let new_alt = std::str::from_utf8(&a[prefix_skip..])
            .map_err(|_| VarEffectError::InvalidAllele)?
            .to_string();

        if pos == orig_pos && new_ref == orig_ref && new_alt == orig_alt {
            Ok(None)
        } else {
            Ok(Some((pos, new_ref, new_alt)))
        }
    }

    // -----------------------------------------------------------------
    // TranscriptStore forwarders
    // -----------------------------------------------------------------

    /// Direct lookup by full versioned accession (e.g. `"NM_000546.6"`).
    /// See [`TranscriptStore::get_by_accession`].
    pub fn get_by_accession(&self, accession: &str) -> Option<(&TranscriptModel, &LocateIndex)> {
        self.transcripts.get_by_accession(accession)
    }

    /// Return every transcript whose `tx_start..tx_end` overlaps the
    /// half-open interval `[start, end)` on `chrom`. See
    /// [`TranscriptStore::query_overlap`].
    pub fn query_overlap(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Vec<(&TranscriptModel, &LocateIndex)> {
        self.transcripts.query_overlap(chrom, start, end)
    }

    // -----------------------------------------------------------------
    // Field accessors (escape hatches)
    // -----------------------------------------------------------------

    /// Borrow the inner [`TranscriptStore`]. Useful when you need access
    /// to store-only methods like `len()`, `is_empty()`, or `transcripts()`.
    pub fn transcripts(&self) -> &TranscriptStore {
        &self.transcripts
    }

    /// Borrow the inner [`FastaReader`]. Useful when you need access to
    /// reader-only methods that aren't forwarded above.
    pub fn fasta(&self) -> &FastaReader {
        &self.fasta
    }
}

#[cfg(test)]
mod left_align_unit_tests {
    //! CI-runnable (no 3.1 GB FASTA) exact-output tests for `left_align_indel`,
    //! covering the shapes whose ref/alt share a rightmost base — duplications
    //! and homopolymer/tandem indels — which is where the per-allele-extend bug
    //! collapsed indels into substitutions. Run against a tiny synthetic genome.

    use super::*;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    /// Synthetic contig "1" (queried as "chr1"):
    ///   1-based pos: 1 2 3 4 5 6 7 8 9 10
    ///   base:        A A T G G G G T A A
    /// The `GGGG` run (pos 4-7) anchored by `T` at pos 3 exercises left-shifting.
    const CONTIG: &[u8] = b"AATGGGGTAA";

    /// Build a `VarEffect` over a synthetic in-memory genome (temp `.bin` +
    /// `.bin.idx`), with an empty transcript store — `left_align_indel` only
    /// touches the FASTA.
    fn ve_with_genome(contigs: &[(&str, &[u8])]) -> (TempDir, VarEffect) {
        let tmp = TempDir::new().expect("tempdir");
        let bin = tmp.path().join("g.bin");
        let idx = tmp.path().join("g.bin.idx");
        write_genome_binary(contigs, "test", &bin, &idx).expect("write synthetic genome");
        let fasta = FastaReader::open(&bin).expect("open synthetic genome");
        let ve = VarEffect::new(TranscriptStore::from_transcripts(Vec::new()), fasta);
        (tmp, ve)
    }

    /// A right-shifted 1 bp duplication inside the `GGGG` run must left-align to
    /// the `T` anchor at pos 3 and REMAIN an insertion (`T>TG`). Pre-fix the
    /// buggy per-allele extension collapsed this to a substitution.
    #[test]
    fn duplication_left_aligns_preserving_insertion() {
        let (_tmp, ve) = ve_with_genome(&[("1", CONTIG)]);
        let result = ve
            .left_align_indel("chr1", 6, "G", "GG")
            .expect("left_align_indel should not error");
        assert_eq!(result, Some((3, "T".to_string(), "TG".to_string())));
    }

    /// The symmetric 1 bp tandem deletion must left-align to the same anchor and
    /// REMAIN a deletion (`TG>T`).
    #[test]
    fn tandem_deletion_left_aligns_preserving_deletion() {
        let (_tmp, ve) = ve_with_genome(&[("1", CONTIG)]);
        let result = ve
            .left_align_indel("chr1", 6, "GG", "G")
            .expect("left_align_indel should not error");
        assert_eq!(result, Some((3, "TG".to_string(), "T".to_string())));
    }

    /// Feeding the already-left-aligned form back in is a no-op (idempotent).
    #[test]
    fn left_aligned_insertion_is_idempotent() {
        let (_tmp, ve) = ve_with_genome(&[("1", CONTIG)]);
        let result = ve
            .left_align_indel("chr1", 3, "T", "TG")
            .expect("left_align_indel should not error");
        assert_eq!(result, None);
    }

    /// An SNV (differing rightmost bases) never enters the shift loop → `None`.
    #[test]
    fn snv_passes_through() {
        let (_tmp, ve) = ve_with_genome(&[("1", CONTIG)]);
        let result = ve
            .left_align_indel("chr1", 8, "T", "A")
            .expect("left_align_indel should not error");
        assert_eq!(result, None);
    }
}
