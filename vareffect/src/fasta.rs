//! Indexed reference genome reader backed by a flat memory-mapped binary.
//!
//! Thin wrapper around [`memmap2::Mmap`] that serves **0-based half-open**
//! sequences by **UCSC-style chromosome name** on the **plus strand**. All
//! three conventions match [`crate::types::TranscriptModel`] so downstream
//! consumers never need to translate between coordinate systems.
//!
//! # On-disk format
//!
//! The reader expects a pair of files produced by `vareffect-cli setup`:
//!
//! - **`GRCh38.bin`** — flat binary genome. Concatenated, newline-stripped,
//!   uppercased chromosome sequences. One byte per base, standard IUPAC
//!   nucleotide codes (`A`/`C`/`G`/`T`/`N` plus ambiguity codes `R`/`Y`/`S`/
//!   `W`/`K`/`M`/`B`/`D`/`H`/`V`). No headers, no line breaks, no padding
//!   between contigs. ~3.1 GB for GRCh38.p14 (primary + patches).
//! - **`GRCh38.bin.idx`** — MessagePack-serialized [`GenomeBinIndex`] mapping
//!   each contig name to its `(offset, length)` in the `.bin` file. ~10 KB.
//!
//! Use [`write_genome_binary`] to produce these files from raw contig data
//! (used by `vareffect-cli setup` and by unit tests).
//!
//! # Coordinate convention
//!
//! All coordinates exposed by [`FastaReader`] are 0-based half-open
//! `[start, end)` — identical to `TranscriptModel::tx_start` / `tx_end`.
//! Converting from a VCF 1-based position is `vcf_pos - 1`.
//!
//! # Chromosome name translation
//!
//! Callers always use UCSC-style names (`chr1`, `chr17`, `chrX`, `chrY`,
//! `chrM`, `chr9_KN196479v1_fix`, ...) to stay consistent with
//! `TranscriptModel::chrom`. The on-disk binary may use any of three naming
//! conventions (inherited from the source FASTA); the reader detects which
//! one at open time by scanning the index entries:
//!
//! - **NCBI RefSeq** (`NC_000001.11`, `NW_*`, ...) — produced by
//!   `vareffect-cli setup` from the NCBI GRCh38.p14 assembly. The reader
//!   translates primary chroms via [`crate::chrom::ucsc_to_refseq`] and
//!   patch contigs via an optional runtime alias CSV loaded through
//!   [`FastaReader::open_with_patch_aliases_and_assembly`].
//! - **UCSC-prefixed** (`chr1`, `chrM`, ...) — pass-through translation.
//! - **Ensembl bare** (`1`, `MT`, ...) — the reader strips the `chr` prefix
//!   and maps `chrM -> MT`.
//!
//! Patch contigs (`chr9_KN196479v1_fix`, `chr22_KI270879v1_alt`, ...) can
//! only be served against an NCBI-naming binary when a
//! `patch_chrom_aliases.csv` is supplied via
//! [`FastaReader::open_with_patch_aliases_and_assembly`].
//!
//! # Thread safety
//!
//! [`FastaReader`] is inherently `Send + Sync` — the underlying
//! [`memmap2::Mmap`] derefs to `&[u8]` with no Mutex required. All threads
//! can read from the same reader concurrently with zero contention.
//! [`FastaReader::try_clone`] is retained for API compatibility but simply
//! clones a handful of `Arc`s — it is no longer needed for parallel
//! workloads.
//!
//! # Soft-masking
//!
//! The flat binary stores uppercase IUPAC nucleotide codes. Soft-mask
//! information (Ensembl lowercase = repeat region) is destroyed at build
//! time. This matches VEP's internal behavior (which uppercases all fetched
//! bases) and GA4GH refget v2.0 (specifies uppercase IUPAC).
//! [`FastaReader::fetch_sequence_raw`] is retained for API compatibility but
//! returns the same uppercase bytes as [`FastaReader::fetch_sequence`].
//!
//! # Usage
//!
//! ```no_run
//! use std::path::Path;
//! use vareffect::{Assembly, FastaReader};
//!
//! // Open the flat binary genome produced by `vareffect-cli setup`.
//! // The assembly must be supplied explicitly — there is no default.
//! let reader = FastaReader::open_with_assembly(
//!     Path::new("data/vareffect/GRCh38.bin"),
//!     Assembly::GRCh38,
//! )?;
//!
//! // TP53 c.742C>T lives at chr17:7674221 (1-based VCF) = chr17:7674220 (0-based).
//! let base = reader.fetch_base("chr17", 7674220)?;
//! assert_eq!(base, b'C');
//!
//! // Patch-contig reads against an NCBI-source binary need the alias CSV.
//! let reader = FastaReader::open_with_patch_aliases_and_assembly(
//!     Path::new("data/vareffect/GRCh38.bin"),
//!     Some(Path::new("data/vareffect/patch_chrom_aliases_grch38.csv")),
//!     Assembly::GRCh38,
//! )?;
//! # Ok::<(), vareffect::VarEffectError>(())
//! ```

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use memmap2::Mmap;
use serde::{Deserialize, Serialize};

use crate::chrom::Assembly;
use crate::error::VarEffectError;

// ---------------------------------------------------------------------------
// Index types (serialized as MessagePack in .bin.idx)
// ---------------------------------------------------------------------------

/// Current format version for [`GenomeBinIndex`]. Increment on breaking
/// changes to the on-disk layout. Public so `vareffect-cli` uses the same
/// constant — there must be a single source of truth for reader/writer
/// version agreement.
pub const GENOME_BIN_INDEX_VERSION: u32 = 1;

#[cfg(not(target_pointer_width = "64"))]
compile_error!("vareffect requires a 64-bit target (genome files may exceed 4 GB)");

/// Flat binary genome index.
///
/// Serialized as MessagePack and stored alongside the `.bin` file with an
/// `.idx` suffix. Read once at [`FastaReader::open_with_assembly`] time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeBinIndex {
    /// Format version. Must equal [`GENOME_BIN_INDEX_VERSION`] at open time.
    pub version: u32,
    /// Assembly build label (e.g. `"GRCh38.p14"`). Informational — not
    /// currently validated at open time, but available for downstream
    /// logging and safety checks.
    pub build: String,
    /// Expected `.bin` file size in bytes. Verified against the mmap length
    /// at open time to catch truncation.
    pub expected_size: u64,
    /// Contigs in file order. Each entry maps a contig name to its byte
    /// offset and length in the `.bin` file.
    pub contigs: Vec<ContigEntry>,
}

/// One contig in the [`GenomeBinIndex`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigEntry {
    /// Raw contig name from the source FASTA (e.g. `"NC_000001.11"`).
    pub name: String,
    /// Byte offset in the `.bin` file where this contig's sequence starts.
    pub offset: u64,
    /// Sequence length in bytes (= number of bases).
    pub length: u64,
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

/// Build a flat binary genome from raw contig sequences.
///
/// Writes two files:
/// - `bin_path`: concatenated uppercase bases, one byte per base, no separators.
/// - `idx_path`: MessagePack-serialized [`GenomeBinIndex`].
///
/// Used by `vareffect-cli setup` for production builds and by unit tests
/// for synthetic genomes.
///
/// # Base validation
///
/// Every byte is uppercased and then checked against the standard IUPAC
/// nucleotide alphabet (`A`, `C`, `G`, `T`, `N`, `R`, `Y`, `S`, `W`, `K`,
/// `M`, `B`, `D`, `H`, `V`). Non-alphabetic or non-IUPAC bytes are rejected.
/// The NCBI GRCh38.p14 assembly uses ambiguity codes (e.g. `M`, `R`, `Y`)
/// in some patch-scaffold regions; these are preserved verbatim.
///
/// # Errors
///
/// * [`VarEffectError::Io`] on write failures.
/// * [`VarEffectError::Malformed`] if any byte is not a valid IUPAC
///   nucleotide code after uppercasing.
pub fn write_genome_binary(
    contigs: &[(&str, &[u8])],
    build: &str,
    bin_path: &Path,
    idx_path: &Path,
) -> Result<(), VarEffectError> {
    let bin_file = File::create(bin_path).map_err(|source| VarEffectError::Io {
        path: bin_path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(bin_file);

    let mut entries = Vec::with_capacity(contigs.len());
    let mut offset: u64 = 0;

    // 64 KB staging buffer for chunk-based uppercase + validation + write.
    // Avoids per-byte function-call overhead on multi-GB genomes.
    let mut buf = [0u8; 64 * 1024];

    for &(name, seq) in contigs {
        let entry_offset = offset;
        for chunk in seq.chunks(buf.len()) {
            let n = chunk.len();
            buf[..n].copy_from_slice(chunk);
            buf[..n].make_ascii_uppercase();
            for (i, &b) in buf[..n].iter().enumerate() {
                if !is_iupac_nucleotide(b) {
                    return Err(VarEffectError::Malformed(format!(
                        "non-IUPAC byte 0x{:02X} ('{}') in contig {name} at offset {}",
                        chunk[i],
                        chunk[i] as char,
                        offset + i as u64,
                    )));
                }
            }
            writer
                .write_all(&buf[..n])
                .map_err(|source| VarEffectError::Io {
                    path: bin_path.to_path_buf(),
                    source,
                })?;
            offset += n as u64;
        }
        entries.push(ContigEntry {
            name: name.to_string(),
            offset: entry_offset,
            length: seq.len() as u64,
        });
    }

    writer.flush().map_err(|source| VarEffectError::Io {
        path: bin_path.to_path_buf(),
        source,
    })?;
    // Ensure data hits disk before writing the index.
    writer
        .into_inner()
        .map_err(|e| VarEffectError::Io {
            path: bin_path.to_path_buf(),
            source: e.into_error(),
        })?
        .sync_all()
        .map_err(|source| VarEffectError::Io {
            path: bin_path.to_path_buf(),
            source,
        })?;

    let index = GenomeBinIndex {
        version: GENOME_BIN_INDEX_VERSION,
        build: build.to_string(),
        expected_size: offset,
        contigs: entries,
    };
    let idx_bytes = rmp_serde::to_vec(&index).map_err(|e| VarEffectError::Io {
        path: idx_path.to_path_buf(),
        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e),
    })?;
    // Atomic write: .tmp + rename so a crash never leaves a partial .bin.idx
    // that the AlreadyPresent fast-path would accept.
    let tmp_idx = idx_path.with_extension("tmp");
    std::fs::write(&tmp_idx, &idx_bytes).map_err(|source| VarEffectError::Io {
        path: idx_path.to_path_buf(),
        source,
    })?;
    std::fs::rename(&tmp_idx, idx_path).map_err(|source| VarEffectError::Io {
        path: idx_path.to_path_buf(),
        source,
    })?;

    Ok(())
}

/// Check whether a byte is a valid uppercase IUPAC nucleotide code.
///
/// Accepts the 15 standard codes: the four canonical bases (`A`, `C`, `G`,
/// `T`), the fully-ambiguous placeholder (`N`), and the ten two-/three-base
/// ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`).
/// The NCBI GRCh38.p14 assembly uses several of these in patch-scaffold
/// regions; rejecting them would prevent building a complete genome binary.
#[inline]
pub fn is_iupac_nucleotide(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'N'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
    )
}

// ---------------------------------------------------------------------------
// Contig naming detection
// ---------------------------------------------------------------------------

/// On-disk contig naming style. Detected by scanning the [`GenomeBinIndex`]
/// contig names at [`FastaReader::open_with_assembly`] time; drives the translation ladder
/// inside [`FastaReader::translate_chrom`].
///
/// Detection priority: `NC_*` -> `NcbiRefSeq` beats a later `chr*` ->
/// `UcscPrefixed` so a mixed-naming source (e.g. `NC_*` primaries + `NW_*`
/// patches) still resolves toward the translation-heavy branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ContigNaming {
    /// Contigs are NCBI RefSeq accessions (`NC_000001.11`, `NW_*`, `NT_*`).
    NcbiRefSeq,
    /// Contigs are UCSC-style (`chr1`, `chrM`, `chr9_KN196479v1_fix`).
    UcscPrefixed,
    /// Contigs are bare Ensembl names (`1`, `MT`, `X`, `Y`).
    EnsemblBare,
}

// ---------------------------------------------------------------------------
// FastaReader
// ---------------------------------------------------------------------------

/// Memory-mapped reference genome reader for random-access sequence retrieval.
///
/// Construct with [`FastaReader::open_with_assembly`] (or
/// [`FastaReader::open_with_patch_aliases_and_assembly`] when patch-contig lookups against
/// an NCBI-source binary are needed), then call
/// [`FastaReader::fetch_sequence`] to pull bases by `(chrom, start, end)`.
/// Coordinates are 0-based half-open; chromosome names are UCSC-style.
///
/// See the module docs for details on the on-disk format, coordinate
/// conventions, chromosome name translation, and thread safety.
pub struct FastaReader {
    /// Memory-mapped genome binary. `Mmap` derefs to `&[u8]` and is
    /// `Send + Sync`, so no Mutex is needed for concurrent reads.
    mmap: Arc<Mmap>,

    /// Per-contig `(byte_offset, length)` keyed by the **raw index name**
    /// (the exact name stored in the `.bin.idx`). The translation ladder
    /// in [`FastaReader::translate_chrom`] converts caller-side UCSC names
    /// *into* raw index keys before probing this map.
    contigs: Arc<HashMap<String, (u64, u64)>>,

    /// Detected on-disk naming convention. See [`ContigNaming`].
    naming: ContigNaming,

    /// Optional UCSC -> RefSeq patch alias table. Only populated when the
    /// caller supplies a `patch_chrom_aliases.csv` to
    /// [`FastaReader::open_with_patch_aliases_and_assembly`] *and* the
    /// binary uses [`ContigNaming::NcbiRefSeq`].
    patch_aliases: Option<Arc<HashMap<String, String>>>,

    /// Reference build the binary represents. Drives the chrom-table
    /// selection in [`crate::chrom::ucsc_to_refseq`].
    assembly: Assembly,

    /// Original `.bin` path, kept for diagnostic messages.
    path: PathBuf,
}

impl FastaReader {
    /// Open a flat binary genome file without a patch-contig alias table.
    ///
    /// This is the common-case constructor for callers that only need
    /// primary-chromosome reads (`chr1`..`chrM`). If the binary uses NCBI
    /// RefSeq naming and you need patch-contig reads as well, use
    /// [`FastaReader::open_with_patch_aliases_and_assembly`] instead.
    ///
    /// The flat-binary format has no embedded assembly tag, so the caller
    /// must declare the build (`Assembly::GRCh38` / `Assembly::GRCh37`).
    /// The reader uses it to select the right chrom-table for
    /// [`crate::chrom::ucsc_to_refseq`] when translating caller-side UCSC
    /// chromosome names into the raw index keys.
    ///
    /// The expected on-disk layout is:
    /// - `path` — the `.bin` file (flat binary genome)
    /// - `{path}.idx` — the MessagePack index (sibling file)
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::IndexNotFound`] if the `.idx` sidecar is missing.
    /// * [`VarEffectError::Io`] on I/O or deserialization errors.
    /// * [`VarEffectError::Malformed`] if the index version is unsupported
    ///   or the `.bin` file size doesn't match the expected size.
    pub fn open_with_assembly(path: &Path, assembly: Assembly) -> Result<Self, VarEffectError> {
        Self::open_with_patch_aliases_and_assembly(path, None, assembly)
    }

    /// Open a flat binary genome file, optionally with a patch-contig alias
    /// table for NCBI-source binaries.
    ///
    /// `patch_aliases_csv` is a path to a `refseq,ucsc` CSV produced by
    /// `vareffect-cli setup`. When supplied *and* the binary uses
    /// NCBI RefSeq contig naming, the reader loads the CSV into a
    /// UCSC -> RefSeq map for the second tier of the translation ladder.
    /// For any other combination, the argument is silently ignored.
    ///
    /// `assembly` declares which build the binary represents, driving
    /// the chrom-table selection used by [`crate::chrom::ucsc_to_refseq`].
    ///
    /// See [`FastaReader::open_with_assembly`] for the expected on-disk
    /// layout and error conditions.
    pub fn open_with_patch_aliases_and_assembly(
        path: &Path,
        patch_aliases_csv: Option<&Path>,
        assembly: Assembly,
    ) -> Result<Self, VarEffectError> {
        // 1. Locate and parse the .bin.idx sidecar.
        let idx_path = append_idx_extension(path);
        if !idx_path.exists() {
            return Err(VarEffectError::IndexNotFound {
                path: idx_path.display().to_string(),
            });
        }
        let idx_bytes = std::fs::read(&idx_path).map_err(|source| VarEffectError::Io {
            path: idx_path.clone(),
            source,
        })?;
        let index: GenomeBinIndex =
            rmp_serde::from_slice(&idx_bytes).map_err(|e| VarEffectError::Io {
                path: idx_path.clone(),
                source: std::io::Error::new(std::io::ErrorKind::InvalidData, e),
            })?;

        // 2. Validate index version.
        if index.version != GENOME_BIN_INDEX_VERSION {
            return Err(VarEffectError::Malformed(format!(
                "unsupported genome binary index version {} (expected {})",
                index.version, GENOME_BIN_INDEX_VERSION,
            )));
        }

        // 3. Memory-map the .bin file.
        let file = File::open(path).map_err(|source| VarEffectError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        // SAFETY: The .bin file is opened read-only and is write-once
        // (produced by `vareffect-cli setup`, never modified in place).
        // Hot-reload uses ArcSwap to atomically swap entire FastaReader
        // instances rather than modifying the underlying file.
        let mmap = unsafe {
            memmap2::MmapOptions::new()
                .map(&file)
                .map_err(|source| VarEffectError::Io {
                    path: path.to_path_buf(),
                    source,
                })?
        };

        // 4. Truncation guard: verify the mmap length matches the expected
        //    size stored in the index.
        if mmap.len() as u64 != index.expected_size {
            return Err(VarEffectError::Malformed(format!(
                "genome binary size mismatch: expected {} bytes, got {} — \
                 file may be truncated or corrupt",
                index.expected_size,
                mmap.len(),
            )));
        }

        // 5. Build the contigs lookup and detect naming convention.
        let mut contigs: HashMap<String, (u64, u64)> = HashMap::with_capacity(index.contigs.len());
        let mut naming = ContigNaming::EnsemblBare;

        for entry in &index.contigs {
            // Detection: NC_* wins outright; chr* flips to UcscPrefixed
            // only if we haven't already seen an NC_* record.
            if entry.name.starts_with("NC_") {
                naming = ContigNaming::NcbiRefSeq;
            } else if entry.name.starts_with("chr") && naming != ContigNaming::NcbiRefSeq {
                naming = ContigNaming::UcscPrefixed;
            }

            contigs.insert(entry.name.clone(), (entry.offset, entry.length));
        }

        // 6. Load the patch alias CSV only when it was provided *and* the
        //    binary actually uses NCBI RefSeq naming.
        let patch_aliases = match (patch_aliases_csv, naming) {
            (Some(csv_path), ContigNaming::NcbiRefSeq) => {
                Some(Arc::new(load_ucsc_to_refseq_aliases(csv_path)?))
            }
            _ => None,
        };

        Ok(Self {
            mmap: Arc::new(mmap),
            contigs: Arc::new(contigs),
            naming,
            patch_aliases,
            assembly,
            path: path.to_path_buf(),
        })
    }

    /// Reference build the binary was opened against.
    pub fn assembly(&self) -> Assembly {
        self.assembly
    }

    /// Mint an additional reader handle for the same genome binary.
    ///
    /// With the memory-mapped backend this is trivially cheap — it clones
    /// a handful of `Arc`s. The `Result` wrapper is retained for API
    /// compatibility with callers written for the previous seek-based
    /// reader; it always returns `Ok`.
    ///
    /// **Note:** With the mmap backend, `try_clone` is no longer needed for
    /// parallel workloads. All threads can read from the same `FastaReader`
    /// concurrently with zero contention.
    pub fn try_clone(&self) -> Result<Self, VarEffectError> {
        Ok(Self {
            mmap: Arc::clone(&self.mmap),
            contigs: Arc::clone(&self.contigs),
            naming: self.naming,
            patch_aliases: self.patch_aliases.as_ref().map(Arc::clone),
            assembly: self.assembly,
            path: self.path.clone(),
        })
    }

    /// Fetch a genomic sequence as uppercase ASCII bytes.
    ///
    /// Coordinates are 0-based half-open `[start, end)`. The returned
    /// sequence is the plus-strand reference — strand-aware complementing
    /// is the caller's responsibility (see `crate::types::Strand`).
    ///
    /// The flat binary stores uppercase bases only, so the returned bytes
    /// are always uppercase IUPAC nucleotide codes (typically `A`/`C`/`G`/
    /// `T`/`N`, occasionally ambiguity codes like `M`/`R`/`Y` in some
    /// GRCh38 patch regions).
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::ChromNotFound`] if `chrom` is not present in the
    ///   genome index.
    /// * [`VarEffectError::CoordinateOutOfRange`] if `start >= end` or `end`
    ///   exceeds the chromosome length.
    pub fn fetch_sequence(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, VarEffectError> {
        let translated = self.translate_chrom(chrom);

        let &(offset, length) =
            self.contigs
                .get(translated.as_ref())
                .ok_or_else(|| VarEffectError::ChromNotFound {
                    chrom: chrom.to_string(),
                })?;

        if start >= end || end > length {
            return Err(VarEffectError::CoordinateOutOfRange {
                chrom: chrom.to_string(),
                start,
                end,
                chrom_len: length,
            });
        }

        // Bounds validated — mmap indexing cannot panic.
        let slice_start = (offset + start) as usize;
        let slice_end = (offset + end) as usize;
        Ok(self.mmap[slice_start..slice_end].to_vec())
    }

    /// Fetch a genomic sequence as raw ASCII bytes.
    ///
    /// With the flat binary backend, this is functionally identical to
    /// [`FastaReader::fetch_sequence`] — the binary stores uppercase bases
    /// only, so soft-mask (lowercase) information is not preserved.
    ///
    /// Retained for API compatibility. Callers that need repeat-region
    /// annotations should use a dedicated RepeatMasker store rather than
    /// relying on FASTA soft-masking.
    pub fn fetch_sequence_raw(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, VarEffectError> {
        self.fetch_sequence(chrom, start, end)
    }

    /// Fetch a single base at the given 0-based position.
    ///
    /// Returns the uppercase ASCII byte. With the mmap backend, this is a
    /// single pointer dereference (~5 ns).
    ///
    /// # Errors
    ///
    /// Same as [`FastaReader::fetch_sequence`].
    pub fn fetch_base(&self, chrom: &str, pos: u64) -> Result<u8, VarEffectError> {
        let translated = self.translate_chrom(chrom);
        let &(offset, length) =
            self.contigs
                .get(translated.as_ref())
                .ok_or_else(|| VarEffectError::ChromNotFound {
                    chrom: chrom.to_string(),
                })?;

        if pos >= length {
            return Err(VarEffectError::CoordinateOutOfRange {
                chrom: chrom.to_string(),
                start: pos,
                end: pos + 1,
                chrom_len: length,
            });
        }

        Ok(self.mmap[(offset + pos) as usize])
    }

    /// Verify that the reference allele at a genomic position matches the
    /// genome binary. Returns `true` on match, `false` on mismatch.
    ///
    /// Both the stored bases and `ref_allele` are compared case-insensitively.
    /// I/O or coordinate errors propagate as `Err` — only a clean byte-level
    /// mismatch returns `Ok(false)`. This mirrors VCF `--check-ref` semantics.
    ///
    /// With the mmap backend, this is a zero-copy slice comparison — no
    /// `Vec<u8>` allocation.
    pub fn verify_ref(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &[u8],
    ) -> Result<bool, VarEffectError> {
        if ref_allele.is_empty() {
            return Ok(true);
        }

        let translated = self.translate_chrom(chrom);
        let &(offset, length) =
            self.contigs
                .get(translated.as_ref())
                .ok_or_else(|| VarEffectError::ChromNotFound {
                    chrom: chrom.to_string(),
                })?;

        let end_pos = pos + ref_allele.len() as u64;
        if end_pos > length {
            return Err(VarEffectError::CoordinateOutOfRange {
                chrom: chrom.to_string(),
                start: pos,
                end: end_pos,
                chrom_len: length,
            });
        }

        // Bounds validated — mmap indexing cannot panic.
        let start = (offset + pos) as usize;
        let end = start + ref_allele.len();
        let slice = &self.mmap[start..end];
        Ok(slice
            .iter()
            .zip(ref_allele.iter())
            .all(|(a, b)| a.eq_ignore_ascii_case(b)))
    }

    /// Return the length of `chrom` in bases, or `None` if the chromosome is
    /// not present in the genome index.
    ///
    /// O(1) — reads from the contig table built at [`FastaReader::open_with_assembly`] time.
    pub fn chrom_length(&self, chrom: &str) -> Option<u64> {
        let translated = self.translate_chrom(chrom);
        self.contigs
            .get(translated.as_ref())
            .map(|&(_, length)| length)
    }

    // ---------------------------------------------------------------
    // Internal helpers
    // ---------------------------------------------------------------

    /// Translate a caller-supplied UCSC chromosome name into the raw index
    /// name. The three-branch dispatcher covers every supported on-disk
    /// naming; see [`ContigNaming`] for the detection logic.
    ///
    /// Returns a `Cow` so the common-case branches (NCBI primary chrom,
    /// UCSC pass-through, Ensembl strip-`chr`) allocate nothing. The
    /// NcbiRefSeq patch-alias tier is the only path that allocates: the
    /// alias map stores owned `String` values, and returning a borrow of
    /// that string would require the returned `Cow`'s lifetime to extend
    /// beyond the caller-supplied `chrom`, which the function signature
    /// doesn't allow.
    fn translate_chrom<'a>(&'a self, chrom: &'a str) -> Cow<'a, str> {
        match self.naming {
            ContigNaming::NcbiRefSeq => {
                // Tier 1: primary chroms via the 25-entry const table.
                // `ucsc_to_refseq` returns `&'static str` on a hit, so the
                // `Cow::Borrowed` here has a 'static lifetime — always safe
                // to coerce into 'a.
                let primary = crate::chrom::ucsc_to_refseq(self.assembly, chrom);
                if !std::ptr::eq(primary.as_ptr(), chrom.as_ptr()) {
                    // `primary` points into `NC_TO_UCSC`, not into `chrom`
                    // — the pointer compare distinguishes "table hit" from
                    // "pass-through" without a second string compare.
                    return Cow::Borrowed(primary);
                }
                // Tier 2: patch contigs via the runtime alias map, if one
                // was supplied at open time.
                if let Some(aliases) = &self.patch_aliases
                    && let Some(refseq) = aliases.get(chrom)
                {
                    return Cow::Owned(refseq.clone());
                }
                // Tier 3: unknown name — pass through. The caller will get
                // a `ChromNotFound` from the contigs lookup, which is the
                // desired failure mode for typos or unknown patches.
                Cow::Borrowed(chrom)
            }
            ContigNaming::UcscPrefixed => {
                // Binary uses UCSC names — no translation needed.
                Cow::Borrowed(chrom)
            }
            ContigNaming::EnsemblBare => {
                if chrom == "chrM" {
                    // Mitochondrial contig is special: Ensembl calls it `MT`.
                    Cow::Owned("MT".to_string())
                } else if let Some(stripped) = chrom.strip_prefix("chr") {
                    Cow::Borrowed(stripped)
                } else {
                    // Caller passed a name that doesn't start with `chr`;
                    // hand it through unchanged so tests and ad-hoc
                    // callers can still query by the raw index name.
                    Cow::Borrowed(chrom)
                }
            }
        }
    }
}

/// Load `patch_chrom_aliases.csv` into a UCSC -> RefSeq lookup map.
///
/// The on-disk CSV format is `refseq,ucsc` (see
/// `vareffect-cli::builders::patch_chrom_aliases`), so this loader inverts
/// the direction at parse time. The same tolerances as the builder-side
/// loader apply: `#` comments, blank lines, the header row, and leading/
/// trailing whitespace are skipped.
///
/// We deliberately duplicate ~25 LOC of CSV parsing here rather than
/// depending on `vareffect-cli`: the `vareffect` crate is designed to be
/// publishable standalone on crates.io, and a dependency on the workspace
/// CLI would defeat that.
fn load_ucsc_to_refseq_aliases(path: &Path) -> Result<HashMap<String, String>, VarEffectError> {
    let file = File::open(path).map_err(|source| VarEffectError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::new(file);

    let mut map: HashMap<String, String> = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|source| VarEffectError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("refseq,") {
            continue;
        }
        let Some((refseq, ucsc)) = trimmed.split_once(',') else {
            continue;
        };
        let refseq = refseq.trim();
        let ucsc = ucsc.trim();
        if refseq.is_empty() || ucsc.is_empty() {
            continue;
        }
        // Inverted direction: key on UCSC, value is RefSeq.
        map.insert(ucsc.to_string(), refseq.to_string());
    }
    Ok(map)
}

impl std::fmt::Debug for FastaReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FastaReader")
            .field("path", &self.path)
            .field("naming", &self.naming)
            .field("patch_aliases", &self.patch_aliases.is_some())
            .field("contigs", &self.contigs.len())
            .finish()
    }
}

/// Derive the `.bin.idx` sidecar path from the `.bin` path by appending
/// `.idx` to the full filename.
fn append_idx_extension(path: &Path) -> PathBuf {
    let mut os = path.as_os_str().to_os_string();
    os.push(".idx");
    PathBuf::from(os)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// Build a flat binary genome from contigs and open it as a `FastaReader`.
    fn write_test_genome(contigs: &[(&str, &[u8])]) -> (TempDir, FastaReader) {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("test.bin");
        let idx_path = tmp.path().join("test.bin.idx");
        write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();
        let reader = FastaReader::open_with_assembly(&bin_path, Assembly::GRCh38).unwrap();
        (tmp, reader)
    }

    /// Write a `refseq,ucsc` alias CSV to a tempdir and return the path.
    fn write_patch_alias_csv(tmp: &TempDir) -> PathBuf {
        let csv_path = tmp.path().join("patch_chrom_aliases.csv");
        let contents = "\
# GRCh38 patch-contig aliases: RefSeq accession -> UCSC contig name.
refseq,ucsc
NW_009646194.1,chr1_KN196472v1_fix
NT_187633.1,chr22_KI270879v1_alt
";
        std::fs::write(&csv_path, contents).unwrap();
        csv_path
    }

    // -- Naming detection --------------------------------------------------

    #[test]
    fn chrom_name_translation_ensembl_style() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT"), ("MT", b"NNNNACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.naming, ContigNaming::EnsemblBare);

        // chr1 -> 1
        let seq = reader.fetch_sequence("chr1", 0, 4).unwrap();
        assert_eq!(seq, b"ACGT");

        // chrM -> MT, start at position 4 (the first `A`)
        let seq = reader.fetch_sequence("chrM", 4, 8).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn chrom_name_translation_ucsc_style() {
        let contigs: &[(&str, &[u8])] = &[("chr1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.naming, ContigNaming::UcscPrefixed);

        // UCSC name passes through unchanged.
        let seq = reader.fetch_sequence("chr1", 0, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn ucsc_naming_detected_even_without_chr1() {
        // A partial assembly without chr1 should still detect UCSC style.
        let contigs: &[(&str, &[u8])] = &[("chr22", b"ACGT"), ("chrY", b"TTTT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(
            reader.naming,
            ContigNaming::UcscPrefixed,
            "partial assembly should still be flagged UCSC-style",
        );
        assert_eq!(reader.fetch_base("chr22", 0).unwrap(), b'A');
        assert_eq!(reader.fetch_base("chrY", 0).unwrap(), b'T');
    }

    #[test]
    fn ncbi_naming_detected_with_nc_prefix() {
        let contigs: &[(&str, &[u8])] = &[
            ("NC_000001.11", b"ACGTACGTACGTACGT"),
            ("NW_009646194.1", b"GGGGCCCC"),
        ];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.naming, ContigNaming::NcbiRefSeq);
        // Caller supplies UCSC form — reader translates to NC_000001.11.
        assert_eq!(reader.fetch_base("chr1", 0).unwrap(), b'A');
        assert_eq!(reader.fetch_sequence("chr1", 0, 4).unwrap(), b"ACGT");
        assert_eq!(reader.chrom_length("chr1"), Some(16));
    }

    #[test]
    fn mixed_ncbi_naming_nc_plus_nw() {
        // Both NC_ and NW_ contigs — NW_ doesn't match NC_ or chr checks,
        // naming should stay NcbiRefSeq.
        let contigs: &[(&str, &[u8])] = &[
            ("NC_000001.11", b"AAAA"),
            ("NW_009646194.1", b"CCCC"),
            ("NT_187633.1", b"GGGG"),
        ];
        let (_tmp, reader) = write_test_genome(contigs);
        assert_eq!(reader.naming, ContigNaming::NcbiRefSeq);
    }

    // -- Coordinate validation ---------------------------------------------

    #[test]
    fn coordinate_validation_start_ge_end() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        let err = reader.fetch_sequence("chr1", 5, 5).unwrap_err();
        assert!(matches!(err, VarEffectError::CoordinateOutOfRange { .. }));

        let err = reader.fetch_sequence("chr1", 10, 5).unwrap_err();
        assert!(matches!(err, VarEffectError::CoordinateOutOfRange { .. }));
    }

    #[test]
    fn coordinate_validation_out_of_range() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        // Contig is 16 bases, so end > 16 is out of range.
        let err = reader.fetch_sequence("chr1", 0, 17).unwrap_err();
        assert!(matches!(err, VarEffectError::CoordinateOutOfRange { .. }));
    }

    #[test]
    fn fetch_base_rejects_position_past_end_of_chrom() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);
        // 16 bases — valid positions are 0..16.
        let err = reader.fetch_base("chr1", 16).unwrap_err();
        assert!(matches!(err, VarEffectError::CoordinateOutOfRange { .. }));
    }

    #[test]
    fn unknown_chrom_returns_err() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        let err = reader.fetch_sequence("chrZZ", 0, 4).unwrap_err();
        assert!(matches!(err, VarEffectError::ChromNotFound { .. }));
    }

    // -- Fetch operations --------------------------------------------------

    #[test]
    fn fetch_base_returns_single_byte() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.fetch_base("chr1", 0).unwrap(), b'A');
        assert_eq!(reader.fetch_base("chr1", 3).unwrap(), b'T');
    }

    #[test]
    fn fetch_sequence_returns_correct_range() {
        let contigs: &[(&str, &[u8])] = &[("chr1", b"ACGTACGTACGTACGT"), ("chrM", b"NNNNACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.fetch_sequence("chr1", 0, 4).unwrap(), b"ACGT");
        assert_eq!(reader.fetch_sequence("chr1", 4, 8).unwrap(), b"ACGT");
        assert_eq!(reader.fetch_sequence("chrM", 0, 4).unwrap(), b"NNNN");
        assert_eq!(reader.fetch_sequence("chrM", 4, 8).unwrap(), b"ACGT");
    }

    #[test]
    fn fetch_sequence_uppercases_lowercase_input() {
        // Builder uppercases — verify lowercase input survives.
        let contigs: &[(&str, &[u8])] = &[("1", b"acgtacgt")];
        let (_tmp, reader) = write_test_genome(contigs);
        let seq = reader.fetch_sequence("chr1", 0, 8).unwrap();
        assert_eq!(seq, b"ACGTACGT");
    }

    #[test]
    fn fetch_sequence_raw_returns_same_as_fetch_sequence() {
        // With flat binary, raw and uppercased are identical.
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        let raw = reader.fetch_sequence_raw("chr1", 0, 8).unwrap();
        let upper = reader.fetch_sequence("chr1", 0, 8).unwrap();
        assert_eq!(raw, upper);
    }

    // -- verify_ref --------------------------------------------------------

    #[test]
    fn verify_ref_match_and_mismatch() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert!(reader.verify_ref("chr1", 0, b"ACGT").unwrap());
        assert!(!reader.verify_ref("chr1", 0, b"TTTT").unwrap());

        // Case-insensitive: caller passes lowercase, binary is uppercase.
        assert!(reader.verify_ref("chr1", 0, b"acgt").unwrap());

        // Empty allele trivially verifies.
        assert!(reader.verify_ref("chr1", 0, b"").unwrap());
    }

    #[test]
    fn verify_ref_out_of_bounds_returns_err() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        // pos + ref_allele.len() > length
        let err = reader.verify_ref("chr1", 2, b"ACGT").unwrap_err();
        assert!(matches!(err, VarEffectError::CoordinateOutOfRange { .. }));
    }

    // -- chrom_length ------------------------------------------------------

    #[test]
    fn chrom_length_reports_correct_values() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT"), ("MT", b"NNNNACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);

        assert_eq!(reader.chrom_length("chr1"), Some(16));
        assert_eq!(reader.chrom_length("chrM"), Some(12));
        assert_eq!(reader.chrom_length("chrZZ"), None);
    }

    // -- try_clone ---------------------------------------------------------

    #[test]
    fn try_clone_produces_independent_reader_sharing_mmap() {
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGTACGTACGTACGT")];
        let (_tmp, reader) = write_test_genome(contigs);
        let cloned = reader.try_clone().expect("try_clone must succeed");

        assert_eq!(reader.fetch_base("chr1", 0).unwrap(), b'A');
        assert_eq!(cloned.fetch_base("chr1", 0).unwrap(), b'A');

        // Shared Arc -> same allocation.
        assert!(
            Arc::ptr_eq(&reader.mmap, &cloned.mmap),
            "mmap Arc must be shared across clones"
        );
        assert!(
            Arc::ptr_eq(&reader.contigs, &cloned.contigs),
            "contigs Arc must be shared across clones"
        );
    }

    // -- Index validation --------------------------------------------------

    #[test]
    fn missing_idx_reports_index_not_found() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("no_index.bin");
        std::fs::write(&bin_path, b"ACGT").unwrap();

        let err = FastaReader::open_with_assembly(&bin_path, Assembly::GRCh38).unwrap_err();
        assert!(matches!(err, VarEffectError::IndexNotFound { .. }));
    }

    #[test]
    fn truncated_bin_reports_malformed() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("truncated.bin");
        let idx_path = tmp.path().join("truncated.bin.idx");

        // Write a valid index claiming 100 bytes...
        let index = GenomeBinIndex {
            version: GENOME_BIN_INDEX_VERSION,
            build: "test".into(),
            expected_size: 100,
            contigs: vec![ContigEntry {
                name: "chr1".into(),
                offset: 0,
                length: 100,
            }],
        };
        let idx_bytes = rmp_serde::to_vec(&index).unwrap();
        std::fs::write(&idx_path, &idx_bytes).unwrap();
        // ...but only write 10 bytes to the bin file.
        std::fs::write(&bin_path, [b'A'; 10]).unwrap();

        let err = FastaReader::open_with_assembly(&bin_path, Assembly::GRCh38).unwrap_err();
        assert!(matches!(err, VarEffectError::Malformed(_)));
    }

    #[test]
    fn wrong_index_version_reports_malformed() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("badver.bin");
        let idx_path = tmp.path().join("badver.bin.idx");

        std::fs::write(&bin_path, b"ACGT").unwrap();
        let index = GenomeBinIndex {
            version: 99,
            build: "test".into(),
            expected_size: 4,
            contigs: vec![ContigEntry {
                name: "chr1".into(),
                offset: 0,
                length: 4,
            }],
        };
        let idx_bytes = rmp_serde::to_vec(&index).unwrap();
        std::fs::write(&idx_path, &idx_bytes).unwrap();

        let err = FastaReader::open_with_assembly(&bin_path, Assembly::GRCh38).unwrap_err();
        assert!(matches!(err, VarEffectError::Malformed(_)));
    }

    // -- Multi-contig boundary ---------------------------------------------

    #[test]
    fn multi_contig_boundary_no_cross_contamination() {
        // chr1 ends with TTTT, chr2 starts with AAAA. Fetching the last
        // byte of chr1 and the first byte of chr2 must not cross over.
        let contigs: &[(&str, &[u8])] = &[("chr1", b"GGGGTTTT"), ("chr2", b"AAAACCCC")];
        let (_tmp, reader) = write_test_genome(contigs);

        // Last byte of chr1
        assert_eq!(reader.fetch_base("chr1", 7).unwrap(), b'T');
        // First byte of chr2
        assert_eq!(reader.fetch_base("chr2", 0).unwrap(), b'A');

        // Range fetch at the boundary
        assert_eq!(reader.fetch_sequence("chr1", 4, 8).unwrap(), b"TTTT");
        assert_eq!(reader.fetch_sequence("chr2", 0, 4).unwrap(), b"AAAA");
    }

    // -- Patch aliases -----------------------------------------------------

    #[test]
    fn ncbi_fasta_rejects_ucsc_patch_without_alias_csv() {
        let contigs: &[(&str, &[u8])] = &[
            ("NC_000001.11", b"ACGTACGTACGTACGT"),
            ("NW_009646194.1", b"GGGGCCCC"),
        ];
        let (_tmp, reader) = write_test_genome(contigs);

        let err = reader.fetch_base("chr1_KN196472v1_fix", 0).unwrap_err();
        assert!(
            matches!(err, VarEffectError::ChromNotFound { .. }),
            "expected ChromNotFound, got {err:?}",
        );
    }

    #[test]
    fn ncbi_fasta_resolves_ucsc_patch_via_alias_csv() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("ncbi.bin");
        let idx_path = tmp.path().join("ncbi.bin.idx");
        let contigs: &[(&str, &[u8])] = &[
            ("NC_000001.11", b"ACGTACGTACGTACGT"),
            ("NW_009646194.1", b"GGGGCCCC"),
        ];
        write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();

        let csv_path = write_patch_alias_csv(&tmp);
        let reader = FastaReader::open_with_patch_aliases_and_assembly(
            &bin_path,
            Some(&csv_path),
            Assembly::GRCh38,
        )
        .unwrap();

        assert_eq!(reader.naming, ContigNaming::NcbiRefSeq);
        assert!(reader.patch_aliases.is_some());

        // `chr1_KN196472v1_fix` -> `NW_009646194.1` via alias.
        assert_eq!(reader.fetch_base("chr1_KN196472v1_fix", 0).unwrap(), b'G');
        assert_eq!(
            reader.fetch_sequence("chr1_KN196472v1_fix", 0, 4).unwrap(),
            b"GGGG",
        );
    }

    #[test]
    fn ncbi_fasta_ignores_alias_csv_for_non_ncbi_binary() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("ens.bin");
        let idx_path = tmp.path().join("ens.bin.idx");
        let contigs: &[(&str, &[u8])] = &[("1", b"ACGT"), ("MT", b"NNNN")];
        write_genome_binary(contigs, "test", &bin_path, &idx_path).unwrap();

        let csv_path = write_patch_alias_csv(&tmp);
        let reader = FastaReader::open_with_patch_aliases_and_assembly(
            &bin_path,
            Some(&csv_path),
            Assembly::GRCh38,
        )
        .unwrap();

        assert_eq!(reader.naming, ContigNaming::EnsemblBare);
        assert!(
            reader.patch_aliases.is_none(),
            "patch aliases must not load against a non-NCBI binary",
        );
        assert_eq!(reader.fetch_base("chr1", 0).unwrap(), b'A');
    }

    // -- Builder validation ------------------------------------------------

    #[test]
    fn builder_rejects_non_iupac_bytes() {
        let tmp = TempDir::new().unwrap();
        let bin_path = tmp.path().join("bad.bin");
        let idx_path = tmp.path().join("bad.bin.idx");

        // 'X' is not a valid IUPAC nucleotide code.
        let err =
            write_genome_binary(&[("chr1", b"ACGTXN")], "test", &bin_path, &idx_path).unwrap_err();
        assert!(matches!(err, VarEffectError::Malformed(_)));
    }

    #[test]
    fn builder_accepts_iupac_ambiguity_codes() {
        // The NCBI GRCh38.p14 assembly uses ambiguity codes (M, R, Y, etc.)
        // in some patch-scaffold regions. The builder must accept all 15
        // standard IUPAC nucleotide codes.
        let seq = b"ACGTNRYSWKMBDHV";
        let contigs: &[(&str, &[u8])] = &[("chr1", seq.as_slice())];
        let (_tmp, reader) = write_test_genome(contigs);
        let fetched = reader.fetch_sequence("chr1", 0, seq.len() as u64).unwrap();
        assert_eq!(fetched, seq.as_slice());
    }

    #[test]
    fn builder_round_trips_all_valid_bases() {
        let seq = b"ACGTNNNTTTAAACCCGGG";
        let contigs: &[(&str, &[u8])] = &[("chr1", seq.as_slice())];
        let (_tmp, reader) = write_test_genome(contigs);
        let fetched = reader.fetch_sequence("chr1", 0, seq.len() as u64).unwrap();
        assert_eq!(fetched, seq.as_slice());
    }

    // -- append_idx_extension ----------------------------------------------

    #[test]
    fn append_idx_extension_preserves_existing_ext() {
        let out = append_idx_extension(Path::new("/tmp/GRCh38.bin"));
        assert_eq!(out, PathBuf::from("/tmp/GRCh38.bin.idx"));
    }
}
