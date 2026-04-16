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

use crate::consequence::ConsequenceResult;
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
    /// `RefMismatch`, `ChromNotFound`, `CoordinateOutOfRange`, `Malformed`.
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
        // Compare rightmost bases. If they match, drop them and shift
        // left. If an allele becomes empty after dropping, left-extend
        // from the reference genome.
        //
        // Coordinate convention: `fetch_base` takes 0-based position.
        // After `pos -= 1`, `pos` is the new 1-based position, so the
        // 0-based index for the base AT `pos` is `pos - 1`.
        loop {
            if pos <= 1 {
                break;
            }
            // Both alleles must be non-empty to compare rightmost
            // bases. On the first iteration both are always non-empty
            // (VCF input guarantees >= 1 base per allele).
            let (Some(&r_last), Some(&a_last)) = (r.last(), a.last()) else {
                break;
            };
            if r_last != a_last {
                break;
            }

            r.pop();
            a.pop();
            pos -= 1;

            // Left-extend empty alleles from the reference genome.
            // After decrement, `pos` is the new 1-based position;
            // 0-based index = `pos - 1`.
            if r.is_empty() {
                r.push(self.fasta.fetch_base(chrom, pos - 1)?);
            }
            if a.is_empty() {
                a.push(self.fasta.fetch_base(chrom, pos - 1)?);
            }
        }

        // Step 2: Left-trim for parsimony.
        //
        // Removes shared leading bases while preserving at least one
        // base per allele (VCF anchor requirement).
        while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
            r.remove(0);
            a.remove(0);
            pos += 1;
        }

        // Step 3: Return None if nothing changed, Some if normalized.
        //
        // Alleles contain only ASCII bytes (ACGTN from the FASTA reader
        // or the original input validated upstream). `from_utf8` cannot
        // fail on valid ASCII but we propagate rather than panic.
        let new_ref = String::from_utf8(r).map_err(|_| VarEffectError::InvalidAllele)?;
        let new_alt = String::from_utf8(a).map_err(|_| VarEffectError::InvalidAllele)?;

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
