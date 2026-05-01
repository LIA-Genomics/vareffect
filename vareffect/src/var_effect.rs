//! Stateful entrypoint that bundles up to two `(TranscriptStore, FastaReader)`
//! pairs — one per supported [`Assembly`] — and routes annotation calls to
//! the slot the caller selects.
//!
//! Most users want to annotate variants without manually threading FASTA and
//! transcript handles through every call. [`VarEffect`] owns the handles and
//! exposes the high-level API as methods. Construct one at startup with the
//! [`VarEffect::builder`], wrap it in `Arc`, and share across worker tasks.
//!
//! Loading both assemblies costs ~6 GB of mmap'd virtual address space (the
//! flat-binary genomes) plus ~50–150 MB resident per [`TranscriptStore`].
//! Pages of the FASTA only enter RSS when touched, so a typical sparse
//! annotation workload stays well under 1 GB of resident memory per assembly.
//!
//! ```no_run
//! use std::path::Path;
//! use vareffect::{Assembly, VarEffect};
//!
//! // Load only what you need — slots default to empty.
//! let ve = VarEffect::builder()
//!     .with_grch38(
//!         Path::new("data/vareffect/transcript_models_grch38.bin"),
//!         Path::new("data/vareffect/GRCh38.bin"),
//!     )?
//!     .build()?;
//!
//! // Annotate TP53 c.742C>T (p.R248W) on chr17 (0-based position).
//! let result = ve.annotate(Assembly::GRCh38, "chr17", 7_674_219, b"C", b"T")?;
//! for csq in result.consequences.iter().flat_map(|r| r.consequences.iter()) {
//!     println!("{}", csq.as_str());
//! }
//! # Ok::<(), vareffect::VarEffectError>(())
//! ```

use std::path::Path;

use crate::chrom::Assembly;
use crate::consequence::{AnnotationResult, Warning};
use crate::error::VarEffectError;
use crate::fasta::FastaReader;
use crate::hgvs_reverse::{GenomicVariant, ResolvedHgvsC};
use crate::locate::LocateIndex;
use crate::transcript::TranscriptStore;
use crate::types::TranscriptModel;

/// Bundled handles for a single assembly slot.
#[derive(Debug)]
struct AssemblyHandles {
    transcripts: TranscriptStore,
    fasta: FastaReader,
}

/// Stateful entrypoint to vareffect: holds up to one `(TranscriptStore,
/// FastaReader)` pair per [`Assembly`] and routes annotation calls to the
/// slot the caller selects.
///
/// `VarEffect` is `Send + Sync` (every inner field is mmap- or `Arc`-backed
/// with no interior mutability), so a single `Arc<VarEffect>` can be shared
/// across all async tasks with zero contention. Cloning the underlying
/// stores is cheap, but the recommended sharing pattern is
/// `Arc<VarEffect>` for one ownership unit.
///
/// `VarEffect` does not derive `Clone` on purpose — sharing should go
/// through `Arc<VarEffect>`, not field-level cloning.
pub struct VarEffect {
    grch38: Option<AssemblyHandles>,
    grch37: Option<AssemblyHandles>,
}

/// Caller-controlled toggles for [`VarEffect::annotate_with_options`].
///
/// Constructed via [`Default`] or struct-literal syntax with
/// `..Default::default()`. New fields may be added in minor releases —
/// the type is `#[non_exhaustive]` so call sites must use the default
/// path for forward compatibility.
///
/// # Examples
///
/// Because the type is `#[non_exhaustive]`, callers from outside this
/// crate cannot use struct-literal syntax. Default-construct, then set
/// the fields you want:
///
/// ```ignore
/// // Default: VRS off, behavior matches the simple `annotate(..)`.
/// let opts = AnnotateOptions::default();
///
/// // Opt in to the canonical VRS 2.0 Allele identifier.
/// let mut opts = AnnotateOptions::default();
/// opts.emit_vrs_v2 = true;
///
/// // Opt in to both schemas for callers indexing against legacy VRS 1.3
/// // consumers as well.
/// let mut opts = AnnotateOptions::default();
/// opts.emit_vrs_v1 = true;
/// opts.emit_vrs_v2 = true;
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct AnnotateOptions {
    /// If true, compute and attach the GA4GH VRS 1.3 Allele identifier
    /// (`vrs_id`) to the returned [`AnnotationResult`]. Off by default.
    ///
    /// VRS 1.3 is the legacy schema; most of the GA4GH ecosystem (anyvar,
    /// ClinGen Allele Registry, ClinVar) is now on 2.0. Enable this flag
    /// only when indexing against a consumer that still requires the 1.3
    /// form.
    ///
    /// Has no effect for variants on non-primary contigs (patches / alts
    /// / unlocalized scaffolds) — those have no canonical cross-pipeline
    /// SQ digest and always return `None` regardless of this flag.
    pub emit_vrs_v1: bool,

    /// If true, compute and attach the GA4GH VRS 2.0 Allele identifier
    /// (`vrs_id_v2`) to the returned [`AnnotationResult`]. Off by default.
    ///
    /// VRS 2.0 is the canonical schema used by anyvar, ClinGen Allele
    /// Registry, ClinVar, and MAVEDB. This is the flag to set when
    /// downstream consumers need content-addressed allele identifiers.
    ///
    /// Either flag (or both) triggers the shared VRS upstream — VOCA
    /// normalizer plus a one-time SHA-512 of the full sequence on each
    /// chromosome's first request — so per-variant cost is dominated by
    /// the schemas you actually request. Has no effect for variants on
    /// non-primary contigs.
    pub emit_vrs_v2: bool,
}

impl VarEffect {
    /// Begin building a multi-assembly `VarEffect`. Use the returned
    /// [`VarEffectBuilder`] to attach one or both assemblies, then call
    /// [`VarEffectBuilder::build`].
    pub fn builder() -> VarEffectBuilder {
        VarEffectBuilder::default()
    }

    fn handles(&self, assembly: Assembly) -> Result<&AssemblyHandles, VarEffectError> {
        let (slot, name) = match assembly {
            Assembly::GRCh38 => (self.grch38.as_ref(), "grch38"),
            Assembly::GRCh37 => (self.grch37.as_ref(), "grch37"),
        };
        slot.ok_or(VarEffectError::AssemblyNotLoaded {
            assembly,
            slot: name,
        })
    }

    /// Return `true` if the given assembly slot was populated by the
    /// builder. Useful for callers that want to gate behavior on which
    /// assemblies are actually available without catching an
    /// [`VarEffectError::AssemblyNotLoaded`] error.
    pub fn has_assembly(&self, assembly: Assembly) -> bool {
        match assembly {
            Assembly::GRCh38 => self.grch38.is_some(),
            Assembly::GRCh37 => self.grch37.is_some(),
        }
    }

    /// Borrow the [`TranscriptStore`] for the given assembly.
    pub fn transcripts(&self, assembly: Assembly) -> Result<&TranscriptStore, VarEffectError> {
        Ok(&self.handles(assembly)?.transcripts)
    }

    /// Borrow the [`FastaReader`] for the given assembly.
    pub fn fasta(&self, assembly: Assembly) -> Result<&FastaReader, VarEffectError> {
        Ok(&self.handles(assembly)?.fasta)
    }

    /// Eagerly fill the per-chromosome VRS SQ-digest cache for
    /// `assembly` in parallel.
    ///
    /// VRS ID emission needs the `ga4gh:SQ.<digest>` of every chromosome
    /// the variant lands on. Those digests are computed lazily on first
    /// use and cached on the matching [`FastaReader`]. The lazy fill
    /// pays a SHA-512 over the full chromosome (~150–180 ms each, ~5.7 s
    /// total for the 25 GRCh38 primary contigs) the first time a variant
    /// touches that chrom.
    ///
    /// For long-running jobs that lazy cost amortizes invisibly. For
    /// short jobs (CLI invocations on small VCFs, fresh process starts,
    /// latency-sensitive batches) it dominates. Calling
    /// `warm_vrs_cache` up front fans the SHA-512 work across rayon's
    /// pool and front-loads the cost out of the hot path.
    ///
    /// Has no effect on annotation correctness — the resulting digests
    /// are byte-identical to the lazy fill, since they are
    /// content-addressed over the same FASTA bytes. Skip this call if
    /// you do not enable [`AnnotateOptions::emit_vrs_v1`] or
    /// [`AnnotateOptions::emit_vrs_v2`]; you would pay the warm cost for
    /// a cache that is never read.
    ///
    /// # Errors
    ///
    /// [`VarEffectError::AssemblyNotLoaded`] if the matching slot was
    /// not populated by the builder.
    pub fn warm_vrs_cache(&self, assembly: Assembly) -> Result<(), VarEffectError> {
        let h = self.handles(assembly)?;
        crate::vrs::warm_cache(&h.fasta);
        Ok(())
    }

    /// Annotate a variant against every overlapping transcript in the
    /// requested assembly's transcript store.
    ///
    /// Returns an [`AnnotationResult`] bundling the per-transcript
    /// consequence rows with any structured warnings raised during
    /// annotation. Today the only warning variant is
    /// [`Warning::DivergentTranscript`], populated when a chosen
    /// transcript carries the
    /// [`crate::TranscriptModel::genome_transcript_divergent`] flag — the
    /// transcript's sequence disagrees with the reference at one or more
    /// positions, so emitted HGVS coordinates may not map back to the
    /// same genomic position they would against the reference.
    ///
    /// The full VCF REF allele is verified against the FASTA before
    /// trimming, so a stated REF that doesn't match the genome produces
    /// [`VarEffectError::RefMismatch`].
    ///
    /// # Arguments
    ///
    /// * `assembly` — Genome build to annotate against. Errors with
    ///   [`VarEffectError::AssemblyNotLoaded`] if the matching slot is
    ///   empty.
    /// * `chrom` — UCSC-style chromosome name (`"chr17"`, `"chrM"`).
    /// * `pos` — 0-based genomic position (BED convention).
    /// * `ref_allele` — VCF REF allele bytes (uppercase ASCII).
    /// * `alt_allele` — VCF ALT allele bytes (uppercase ASCII).
    pub fn annotate(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> Result<AnnotationResult, VarEffectError> {
        self.annotate_with_options(
            assembly,
            chrom,
            pos,
            ref_allele,
            alt_allele,
            &AnnotateOptions::default(),
        )
    }

    /// Annotate a variant with caller-controlled options.
    ///
    /// Same contract as [`Self::annotate`], but takes an
    /// [`AnnotateOptions`] reference for behavior toggles that the
    /// default-off path would otherwise pay for unconditionally. Use
    /// this when you need:
    ///
    /// * `emit_vrs_v1 = true` — compute and attach the VRS 1.3 Allele
    ///   identifier (`vrs_id`) to the returned [`AnnotationResult`].
    /// * `emit_vrs_v2 = true` — compute and attach the canonical VRS 2.0
    ///   Allele identifier (`vrs_id_v2`).
    ///
    /// Either flag triggers the shared VRS upstream: VOCA normalizer plus
    /// a one-time per-chromosome SHA-512 of the full sequence (~5.7 s for
    /// all 25 primary contigs of GRCh38) on first request per
    /// `FastaReader`. Per-variant cost on warm cache is roughly 6.5 µs
    /// per requested schema (~13 µs when both are on). For short jobs,
    /// call [`Self::warm_vrs_cache`] up front to amortize the
    /// per-chromosome SHA-512 cost out of the hot path.
    pub fn annotate_with_options(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
        options: &AnnotateOptions,
    ) -> Result<AnnotationResult, VarEffectError> {
        let h = self.handles(assembly)?;
        let consequences = crate::consequence::annotate(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            &h.transcripts,
            &h.fasta,
        )?;

        // Aggregate divergent-transcript warnings. A single divergent
        // transcript can show up in multiple consequence rows for the
        // same variant (e.g. one consequence per overlapping transcript
        // tier), so dedup by accession.
        let mut warnings: Vec<Warning> = Vec::new();
        let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for row in &consequences {
            if row.transcript.is_empty() {
                continue;
            }
            if !seen.insert(row.transcript.as_str()) {
                continue;
            }
            if let Some((tx, _)) = h.transcripts.get_by_accession(&row.transcript)
                && tx.genome_transcript_divergent
            {
                warnings.push(Warning::DivergentTranscript {
                    accession: row.transcript.clone(),
                });
            }
        }

        let (vrs_id, vrs_id_v2) = if options.emit_vrs_v1 || options.emit_vrs_v2 {
            crate::vrs::compute_vrs_ids(
                assembly, chrom, pos, ref_allele, alt_allele, &h.fasta, options,
            )
        } else {
            (None, None)
        };

        Ok(AnnotationResult {
            consequences,
            warnings,
            vrs_id,
            vrs_id_v2,
        })
    }

    /// Annotate a variant and return the result as VEP REST-compatible JSON.
    ///
    /// Convenience wrapper around [`Self::annotate`] that serializes the
    /// output as a single-element JSON array matching Ensembl VEP REST's
    /// response shape for
    /// `GET /vep/human/region/{chrom}:{pos}:{pos}/{alt}?refseq=1&hgvs=1`.
    ///
    /// `assembly` drives both the annotation routing and the JSON
    /// `assembly_name` field.
    ///
    /// vareffect is a strict subset of VEP -- it does NOT populate
    /// `canonical`, `swissprot`, `protein_id`, `sift_prediction`,
    /// `polyphen_prediction`, `revel_score`, `alphamissense`, `spliceai`, or
    /// `colocated_variants[*].id`. Consumers must tolerate their absence.
    ///
    /// # Errors
    ///
    /// Propagates any [`VarEffectError`] raised by [`Self::annotate`]:
    /// `AssemblyNotLoaded`, `RefMismatch`, `ChromNotFound`,
    /// `CoordinateOutOfRange`, `Malformed`. Serialization itself is
    /// infallible.
    pub fn annotate_to_vep_json(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> Result<serde_json::Value, VarEffectError> {
        let result = self.annotate(assembly, chrom, pos, ref_allele, alt_allele)?;
        Ok(crate::vep_json::to_vep_json_array(
            chrom,
            pos,
            ref_allele,
            alt_allele,
            assembly.as_str(),
            &result.consequences,
        ))
    }

    /// Resolve HGVS c. notation to plus-strand genomic coordinates against
    /// the given assembly's transcript store.
    ///
    /// Returns a [`GenomicVariant`] with 0-based position, UCSC-style
    /// chromosome, and uppercase plus-strand alleles. Pair with
    /// [`VarEffect::annotate`] when you also want consequences.
    ///
    /// # Errors
    ///
    /// `AssemblyNotLoaded`, `HgvsParseError`, `TranscriptNotFound`,
    /// `PositionOutOfRange`, `HgvsRefMismatch`, `ChromNotFound`,
    /// `CoordinateOutOfRange`, `Malformed`.
    pub fn resolve_hgvs_c(
        &self,
        assembly: Assembly,
        hgvs: &str,
    ) -> Result<GenomicVariant, VarEffectError> {
        let h = self.handles(assembly)?;
        crate::hgvs_reverse::resolve_hgvs_c(hgvs, &h.transcripts, &h.fasta)
    }

    /// Resolve HGVS c. notation and report which transcript version was
    /// actually used. See [`VarEffect::resolve_hgvs_c`] for the
    /// non-meta variant.
    pub fn resolve_hgvs_c_with_meta(
        &self,
        assembly: Assembly,
        hgvs: &str,
    ) -> Result<ResolvedHgvsC, VarEffectError> {
        let h = self.handles(assembly)?;
        crate::hgvs_reverse::resolve_hgvs_c_with_meta(hgvs, &h.transcripts, &h.fasta)
    }

    /// Fetch a single base at the given 0-based position from the chosen
    /// assembly's reference. ~5 ns with the mmap backend.
    pub fn fetch_base(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos: u64,
    ) -> Result<u8, VarEffectError> {
        self.handles(assembly)?.fasta.fetch_base(chrom, pos)
    }

    /// Fetch a genomic sequence as uppercase ASCII bytes from the chosen
    /// assembly. Coordinates are 0-based half-open `[start, end)`.
    pub fn fetch_sequence(
        &self,
        assembly: Assembly,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, VarEffectError> {
        self.handles(assembly)?
            .fasta
            .fetch_sequence(chrom, start, end)
    }

    /// Verify that the reference allele at a position matches the genome
    /// for the chosen assembly.
    pub fn verify_ref(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
    ) -> Result<bool, VarEffectError> {
        self.handles(assembly)?
            .fasta
            .verify_ref(chrom, pos, ref_allele)
    }

    /// Length of `chrom` in bases for the chosen assembly, or `None` if
    /// not present in that assembly's index.
    pub fn chrom_length(&self, assembly: Assembly, chrom: &str) -> Option<u64> {
        self.grch38
            .as_ref()
            .filter(|_| matches!(assembly, Assembly::GRCh38))
            .or_else(|| {
                self.grch37
                    .as_ref()
                    .filter(|_| matches!(assembly, Assembly::GRCh37))
            })
            .and_then(|h| h.fasta.chrom_length(chrom))
    }

    /// Anchor-prepend HGVS-style indel alleles to VCF form using the chosen
    /// assembly's reference genome.
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
    /// SNVs, MNVs, and complex substitutions (neither allele `"-"`) pass
    /// through unchanged and return `Ok(None)`.
    pub fn anchor_prepend_indel(
        &self,
        assembly: Assembly,
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

        let anchor_pos_0 = match pos_1based.checked_sub(2) {
            Some(p) => p,
            None => {
                let chrom_len = self.chrom_length(assembly, chrom).ok_or_else(|| {
                    VarEffectError::ChromNotFound {
                        chrom: chrom.to_string(),
                    }
                })?;
                return Err(VarEffectError::CoordinateOutOfRange {
                    chrom: chrom.to_string(),
                    start: 0,
                    end: 0,
                    chrom_len,
                });
            }
        };

        let anchor_byte = self.fetch_base(assembly, chrom, anchor_pos_0)?;
        let anchor = anchor_byte as char;
        let new_pos = pos_1based - 1;

        let (new_ref, new_alt) = if is_deletion {
            (format!("{anchor}{ref_allele}"), anchor.to_string())
        } else {
            (anchor.to_string(), format!("{anchor}{alt_allele}"))
        };
        Ok(Some((new_pos, new_ref, new_alt)))
    }

    /// Left-align a VCF-style variant to the leftmost equivalent position
    /// using the chosen assembly's reference.
    pub fn left_align_indel(
        &self,
        assembly: Assembly,
        chrom: &str,
        pos_1based: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<(u64, String, String)>, VarEffectError> {
        let h = self.handles(assembly)?;
        let orig_pos = pos_1based;
        let orig_ref = ref_allele;
        let orig_alt = alt_allele;

        let mut pos = pos_1based;
        let mut r: Vec<u8> = ref_allele.as_bytes().to_vec();
        let mut a: Vec<u8> = alt_allele.as_bytes().to_vec();

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
                let prepend = h.fasta.fetch_base(chrom, pos - 1)?;
                r.insert(0, prepend);
                a.insert(0, prepend);
            } else if !trimmed {
                break;
            }
        }

        let mut prefix_skip = 0usize;
        while r.len() - prefix_skip > 1
            && a.len() - prefix_skip > 1
            && r[prefix_skip] == a[prefix_skip]
        {
            prefix_skip += 1;
            pos += 1;
        }

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

    /// Direct lookup by full versioned accession in the chosen assembly's
    /// store (e.g. `"NM_000546.6"`).
    pub fn get_by_accession(
        &self,
        assembly: Assembly,
        accession: &str,
    ) -> Result<Option<(&TranscriptModel, &LocateIndex)>, VarEffectError> {
        Ok(self
            .handles(assembly)?
            .transcripts
            .get_by_accession(accession))
    }

    /// Return every transcript in the chosen assembly whose
    /// `tx_start..tx_end` overlaps the half-open interval `[start, end)`
    /// on `chrom`.
    pub fn query_overlap(
        &self,
        assembly: Assembly,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<(&TranscriptModel, &LocateIndex)>, VarEffectError> {
        Ok(self
            .handles(assembly)?
            .transcripts
            .query_overlap(chrom, start, end))
    }
}

/// Builder for [`VarEffect`]. Use [`VarEffect::builder`] to obtain one,
/// attach assemblies via [`VarEffectBuilder::with_grch38`] /
/// [`VarEffectBuilder::with_grch37`], then call
/// [`VarEffectBuilder::build`].
///
/// Builder semantics: **first-error short-circuit**. If `with_grch37(...)`
/// succeeds and a subsequent `with_grch38(...)` fails, the GRCh37 handle
/// is dropped (along with the partially-loaded builder) before the error
/// propagates. Constructor enforces
/// `transcript_store.assembly() == fasta_reader.assembly()` and errors
/// with [`VarEffectError::AssemblyMismatch`] if they disagree.
#[derive(Default, Debug)]
pub struct VarEffectBuilder {
    grch38: Option<AssemblyHandles>,
    grch37: Option<AssemblyHandles>,
}

impl VarEffectBuilder {
    /// Load the GRCh38 transcript store and reference genome.
    ///
    /// `transcripts_path` should point at the MessagePack file produced
    /// by `vareffect setup --assembly grch38` (typically
    /// `transcript_models_grch38.bin`); the sibling `.manifest.json`
    /// provides the assembly identifier. `genome_path` should point at
    /// the flat-binary genome with its `.bin.idx` sidecar.
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::Io`] / [`VarEffectError::Deserialize`] /
    ///   [`VarEffectError::AssemblyMissingFromManifest`] from
    ///   [`TranscriptStore::load_from_path`].
    /// * [`VarEffectError::IndexNotFound`] / [`VarEffectError::Io`] from
    ///   [`FastaReader::open_with_assembly`].
    /// * [`VarEffectError::AssemblyMismatch`] if the loaded transcript
    ///   store reports a different assembly than `Assembly::GRCh38`.
    pub fn with_grch38(
        mut self,
        transcripts_path: &Path,
        genome_path: &Path,
    ) -> Result<Self, VarEffectError> {
        self.grch38 = Some(load_handles(
            Assembly::GRCh38,
            transcripts_path,
            genome_path,
            None,
        )?);
        Ok(self)
    }

    /// Same as [`Self::with_grch38`] plus an optional patch-contig alias
    /// CSV. Use when annotating against patch contigs from the GRCh38.p14
    /// assembly report (`patch_chrom_aliases_grch38.csv`).
    pub fn with_grch38_and_patch_aliases(
        mut self,
        transcripts_path: &Path,
        genome_path: &Path,
        patch_aliases_csv: &Path,
    ) -> Result<Self, VarEffectError> {
        self.grch38 = Some(load_handles(
            Assembly::GRCh38,
            transcripts_path,
            genome_path,
            Some(patch_aliases_csv),
        )?);
        Ok(self)
    }

    /// Load the GRCh37 transcript store and reference genome.
    ///
    /// See [`Self::with_grch38`] for the on-disk layout requirements;
    /// the GRCh37 build produces `transcript_models_grch37.bin` /
    /// `GRCh37.bin` analogues via `vareffect setup --assembly grch37`.
    ///
    /// # Errors
    ///
    /// Same as [`Self::with_grch38`], but the manifest must report
    /// `Assembly::GRCh37`.
    pub fn with_grch37(
        mut self,
        transcripts_path: &Path,
        genome_path: &Path,
    ) -> Result<Self, VarEffectError> {
        self.grch37 = Some(load_handles(
            Assembly::GRCh37,
            transcripts_path,
            genome_path,
            None,
        )?);
        Ok(self)
    }

    /// Same as [`Self::with_grch37`] plus an optional patch-contig alias
    /// CSV (`patch_chrom_aliases_grch37.csv`).
    pub fn with_grch37_and_patch_aliases(
        mut self,
        transcripts_path: &Path,
        genome_path: &Path,
        patch_aliases_csv: &Path,
    ) -> Result<Self, VarEffectError> {
        self.grch37 = Some(load_handles(
            Assembly::GRCh37,
            transcripts_path,
            genome_path,
            Some(patch_aliases_csv),
        )?);
        Ok(self)
    }

    /// Attach a pair of pre-loaded handles for the given assembly.
    ///
    /// Preferred entry point for tests and build-time callers that
    /// produce stores in-memory rather than from disk. The constructor
    /// still enforces
    /// `transcripts.assembly() == fasta.assembly() == assembly` and
    /// errors with [`VarEffectError::AssemblyMismatch`] on disagreement.
    pub fn with_handles(
        mut self,
        assembly: Assembly,
        transcripts: TranscriptStore,
        fasta: FastaReader,
    ) -> Result<Self, VarEffectError> {
        if transcripts.assembly() != assembly || fasta.assembly() != assembly {
            return Err(VarEffectError::AssemblyMismatch {
                transcripts: transcripts.assembly(),
                fasta: fasta.assembly(),
            });
        }
        let handles = AssemblyHandles { transcripts, fasta };
        match assembly {
            Assembly::GRCh38 => self.grch38 = Some(handles),
            Assembly::GRCh37 => self.grch37 = Some(handles),
        }
        Ok(self)
    }

    /// Finalize the builder. Returns an error only if no assembly was
    /// attached — a `VarEffect` with both slots empty would error on
    /// every annotate call, so the build step fails fast instead.
    pub fn build(self) -> Result<VarEffect, VarEffectError> {
        if self.grch38.is_none() && self.grch37.is_none() {
            return Err(VarEffectError::Malformed(
                "VarEffectBuilder::build called without any assembly attached; \
                 use with_grch38 / with_grch37 / with_handles before build"
                    .to_string(),
            ));
        }
        Ok(VarEffect {
            grch38: self.grch38,
            grch37: self.grch37,
        })
    }
}

/// Load `(TranscriptStore, FastaReader)` for a single assembly slot and
/// validate that both stores agree on the assembly identifier.
fn load_handles(
    expected: Assembly,
    transcripts_path: &Path,
    genome_path: &Path,
    patch_aliases_csv: Option<&Path>,
) -> Result<AssemblyHandles, VarEffectError> {
    let transcripts = TranscriptStore::load_from_path(transcripts_path)?;
    let fasta = match patch_aliases_csv {
        Some(csv) => {
            FastaReader::open_with_patch_aliases_and_assembly(genome_path, Some(csv), expected)?
        }
        None => FastaReader::open_with_assembly(genome_path, expected)?,
    };
    if transcripts.assembly() != expected || fasta.assembly() != expected {
        return Err(VarEffectError::AssemblyMismatch {
            transcripts: transcripts.assembly(),
            fasta: fasta.assembly(),
        });
    }
    Ok(AssemblyHandles { transcripts, fasta })
}
