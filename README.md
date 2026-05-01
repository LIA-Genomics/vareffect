# vareffect


[![CI](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml)
[![vareffect](https://img.shields.io/crates/v/vareffect.svg?label=vareffect)](https://crates.io/crates/vareffect)
[![vareffect-cli](https://img.shields.io/crates/v/vareffect-cli.svg?label=vareffect-cli)](https://crates.io/crates/vareffect-cli)
[![docs.rs](https://docs.rs/vareffect/badge.svg)](https://docs.rs/vareffect)
[![Apache-2.0](https://img.shields.io/crates/l/vareffect)](#license)


Variant consequence prediction and HGVS notation in Rust, concordant with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/).

vareffect takes a variant (chromosome, position, ref, alt) and returns Sequence Ontology consequence terms, HGVS c./p. notation, and NMD prediction for every overlapping transcript -- the same questions VEP answers, at 100--400x the speed, with zero runtime dependencies. The engine is thread-safe (`Send + Sync`), memory-mapped, and designed to embed directly into Rust pipelines.

This repository contains two crates:

| Crate | Description |
|-------|-------------|
| [`vareffect`](vareffect/) | Core library -- consequence prediction, HGVS, transcript/genome stores |
| [`vareffect-cli`](vareffect-cli/) | CLI tool -- config scaffolding (`init`), validation (`check`), data provisioning (`setup`), and VCF annotation (`annotate`) |

## Quick start

### 1. Install and set up runtime data

```bash
cargo install vareffect-cli

# Scaffold a config file (one-time):
vareffect init

# Downloads GRCh38 + MANE, builds genome binary and transcript store.
# One-time, ~10 minutes, ~3 GB disk. Idempotent.
vareffect setup                        # default: --assembly grch38
vareffect setup --assembly grch37      # GRCh37 + RefSeq Select (~6 GB more)
vareffect setup --assembly all         # both, side-by-side under data/vareffect/

# Validate everything is in place:
vareffect check
```

### 2. Annotate a VCF with the CLI

```bash
vareffect annotate \
  --assembly grch38 \
  --input sample.vcf.gz \
  --output annotated.vcf.gz \
  --fasta data/vareffect/GRCh38.bin \
  --transcripts data/vareffect/transcript_models_grch38.bin \
  --threads 8
```

Output is a standard VCF with a `CSQ` INFO field whose pipe-delimited layout matches VEP `--vcf` output. Downstream tools that parse VEP CSQ fields (SnpSift, bcftools, GEMINI) work unchanged.

### 3. Use as a library

```toml
# Cargo.toml
[dependencies]
vareffect = "0.1"
```

```rust
use std::path::Path;
use vareffect::{Assembly, Consequence, VarEffect};

let ve = VarEffect::builder()
    .with_grch38(
        Path::new("data/vareffect/transcript_models_grch38.bin"),
        Path::new("data/vareffect/GRCh38.bin"),
    )?
    .build()?;

// TP53 c.742C>T (p.Arg248Trp) — chr17, 0-based position
let result = ve.annotate(Assembly::GRCh38, "chr17", 7_674_219, b"C", b"T")?;

for r in &result.consequences {
    println!(
        "{} {}: {} ({})",
        r.transcript,
        r.hgvs_p.as_deref().unwrap_or("-"),
        r.consequences.iter().map(Consequence::as_str).collect::<Vec<_>>().join(","),
        r.impact,
    );
}
// => NM_000546.6 p.Arg248Trp: missense_variant (MODERATE)
```

The builder accepts any combination of GRCh38 / GRCh37 slots in a single
`VarEffect`; pass the matching `Assembly` to each `annotate()` call:

```rust
let ve = VarEffect::builder()
    .with_grch38(grch38_transcripts, grch38_genome)?
    .with_grch37(grch37_transcripts, grch37_genome)?
    .build()?;

let g38 = ve.annotate(Assembly::GRCh38, "chr17", 7_674_219, b"C", b"T")?;
let g37 = ve.annotate(Assembly::GRCh37, "chr17", 7_577_120, b"C", b"T")?;
```

For multi-threaded use, wrap in `Arc`:

```rust
use std::sync::Arc;
let ve = Arc::new(ve);
// Clone the Arc into each worker — zero contention, no interior mutability.
```

## Features

- **24 SO consequence terms** covering SNVs, insertions, deletions, MNVs, and complex indels
- **HGVS c./n. and p. notation** -- substitutions, deletions, insertions, duplications, delins, intronic/UTR offsets, frameshift extension walks
- **HGVS reverse resolution** -- parse `NM_000546.6:c.742C>T` back to genomic coordinates
- **NMD prediction** via the 50-nucleotide rule on truncating variants
- **IMPACT ranking** -- HIGH / MODERATE / LOW / MODIFIER, matching VEP's severity scale
- **Multi-transcript annotation** -- one result per overlapping transcript with MANE Select / MANE Plus Clinical / RefSeq Select tier metadata
- **Memory-mapped genome** -- ~5 ns base fetches via mmap'd flat binary
- **Thread-safe** -- `VarEffect`, `TranscriptStore`, and `FastaReader` are `Send + Sync` with a compile-time assertion
- **Standard + mitochondrial genetic codes** -- `chrM` uses NCBI table 2 automatically
- **Reference allele verification** -- mismatched REF fails fast before annotation

## VarEffect API

All coordinates are 0-based half-open (BED convention). Chromosome names are UCSC-style (`chr17`, `chrM`). Assembly is always passed explicitly — there is no implicit default, because misrouting under the wrong chrom table silently produces off-by-millions annotations.

### Construction

`VarEffect` is built via the typestate `VarEffectBuilder`. Each `with_*` method loads one assembly slot; `build()` finalises. Both slots are independent, so a single `VarEffect` can serve GRCh38 and GRCh37 traffic simultaneously.

| Builder method | Description |
|--------|-------------|
| `VarEffect::builder()` | Start a new builder (returns `VarEffectBuilder`) |
| `.with_grch38(transcripts, genome)` | Load GRCh38 from disk into the GRCh38 slot |
| `.with_grch38_and_patch_aliases(transcripts, genome, aliases)` | GRCh38 with NCBI patch-contig alias support |
| `.with_grch37(transcripts, genome)` | Load GRCh37 from disk into the GRCh37 slot |
| `.with_grch37_and_patch_aliases(transcripts, genome, aliases)` | GRCh37 with patch-contig alias support |
| `.with_handles(assembly, transcripts, fasta)` | Insert a pre-built `(TranscriptStore, FastaReader)` pair (assembly-generic) |
| `.build()` | Finalise into `VarEffect` |

### Annotation

| Method | Description |
|--------|-------------|
| `annotate(assembly, chrom, pos, ref, alt)` | Returns `AnnotationResult { consequences, warnings }` — one `ConsequenceResult` per overlapping transcript |
| `annotate_to_vep_json(assembly, chrom, pos, ref, alt)` | Same, serialised as VEP REST-compatible JSON |
| `resolve_hgvs_c(assembly, hgvs)` | Parse HGVS c. string to `GenomicVariant` (genomic coordinates + alleles) |

### Reference genome

All take `assembly` as the first argument; route to the matching slot. For raw access to the inner `FastaReader`, use `ve.fasta(assembly)`.

| Method | Description |
|--------|-------------|
| `fetch_base(assembly, chrom, pos)` | Single base lookup (~5 ns) |
| `fetch_sequence(assembly, chrom, start, end)` | Sequence slice as `Vec<u8>` |
| `verify_ref(assembly, chrom, pos, ref_allele)` | Check REF matches genome (zero-copy) |
| `chrom_length(assembly, chrom)` | Chromosome length or `None` |

### Indel normalization

| Method | Description |
|--------|-------------|
| `anchor_prepend_indel(assembly, chrom, pos, ref, alt)` | Convert HGVS `"-"` placeholders to VCF anchor-prepended form |
| `left_align_indel(assembly, chrom, pos, ref, alt)` | Left-align to leftmost position (Tan et al. 2015) |

### Transcript queries

All take `assembly` as the first argument. For raw access to the inner `TranscriptStore`, use `ve.transcripts(assembly)`.

| Method | Description |
|--------|-------------|
| `get_by_accession(assembly, accession)` | O(1) lookup by versioned accession (e.g. `NM_000546.6`) |
| `query_overlap(assembly, chrom, start, end)` | O(log n + k) interval-tree overlap query |
| `has_assembly(assembly)` | True if the slot was loaded by the builder |

## Performance

Single-threaded throughput on a modern x86_64 laptop (excludes startup cost):

| Tool | Language | Variants/sec (1 thread) |
|------|----------|------------------------|
| VEP (`--cache`) | Perl | ~200 -- 500 |
| SnpEff | Java | ~2,000 -- 10,000 |
| Nirvana | C# | ~5,000 -- 15,000 |
| bcftools csq | C | ~10,000 -- 50,000 |
| **vareffect** | **Rust** | **~50,000 -- 200,000** |

The gap is almost entirely I/O: VEP reads BGZF-compressed FASTA through block decompression per base; vareffect memory-maps a flat binary and reads bytes directly. The CLI annotator achieves ~270k variants/sec with near-linear multi-thread scaling.

## VEP divergences

vareffect targets concordance with VEP releases 115/116. Validated against 136 hand-curated variants across 6 test suites (6/6 pass as of 2026-04-11). Full details in [`VEP_DIVERGENCES.md`](vareffect/VEP_DIVERGENCES.md).

### Large-scale concordance

At-scale run over 50,000 variants (45,644 compared after skipping rows without a RefSeq transcript): 99.38% consequence concordance, 97.52% HGVS c., 97.72% HGVS p., 99.78% IMPACT, 98.10% protein start. Zero vareffect errors or panics.

```
============================================================
LARGE-SCALE VEP CONCORDANCE REPORT
============================================================
Total TSV rows:                50000
Skipped (no RefSeq tc):           97
Skipped (vareffect error):         0
Skipped (vareffect PANIC):         0  <- hard bugs, see mismatches.log
Transcript not found:           4259
Transcript version mismatch:       0
Compared:                      45644

Consequence concordance:   45361/45644 (99.38%)
  Exact match:             40180       (88.03%)
  Normalised (splice/NMD):  4504       (9.87%)
  Real mismatch:             283       (0.62%)

HGVS c. concordance:       44511/45641 (97.52%)
HGVS p. concordance:       33008/33777 (97.72%)
Impact concordance:        45543/45644 (99.78%)
Protein start concordance: 33447/34096 (98.10%)

TOP MISMATCH PATTERNS (up to 20):
    728x  hgvsc
    365x  protein_start
    338x  hgvsp
    316x  hgvsc+hgvsp
    184x  hgvsp+protein_start
    144x  csq
     66x  hgvsc+protein_start
     53x  csq+impact
     24x  csq+hgvsp
     18x  csq+hgvsp+protein_start
     13x  csq+hgvsp+impact
     13x  impact
      7x  csq+protein_start
      5x  csq+hgvsc+hgvsp+impact
      5x  csq+hgvsc+protein_start
      3x  csq+hgvsc+hgvsp
      3x  csq+hgvsc+impact
      3x  hgvsp+impact+protein_start
      2x  csq+hgvsc
      2x  csq+hgvsc+impact+protein_start

Throughput: 249684 variants/second
Elapsed:    0.2s
============================================================
```

### Intentional divergences

| Area | vareffect behaviour | VEP behaviour |
|------|---------------------|---------------|
| Splice sub-terms | All emitted as `splice_region_variant` | Optionally emits `splice_polypyrimidine_tract_variant`, `splice_donor_region_variant`, `splice_donor_5th_base_variant` |
| Non-coding introns | Emits `intron_variant` | May emit `non_coding_transcript_variant` |
| NMD | `predicts_nmd` boolean via 50-nt rule | Separate `NMD_transcript_variant` SO term |
| HGVS 3' shift | Always on (no toggle) | Configurable via `--shift_hgvs` |

### Not yet implemented

| Feature | Reason |
|---------|--------|
| `mature_miRNA_variant` | Requires miRNA locus track |
| Regulatory terms (`TFBS_ablation`, `TF_binding_site_variant`, `regulatory_region_variant`) | Requires regulatory feature store |
| SV terms (`feature_elongation`, `feature_truncation`, `transcript_amplification`) | Requires segment-level SV input |
| Alternate builds (CHM13) | Build pipeline work, not engine limitation. GRCh38 and GRCh37 are both supported via `Assembly` selector. |
| Multi-allele VCF splitting | Caller's responsibility |
| Canonical transcript selection | Caller filters on tier metadata |

### Out of scope (by design)

No plugin system. No co-located variant lookup (gnomAD, ClinVar, dbSNP). No allele frequency annotation. No built-in VCF I/O in the library.

## License

Licensed under [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)