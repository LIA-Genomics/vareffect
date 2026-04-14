# vareffect

Rust variant consequence prediction and HGVS notation, concordant with
[Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/).
 
[![crates.io](https://img.shields.io/crates/v/vareffect.svg)](https://crates.io/crates/vareffect)
[![docs.rs](https://img.shields.io/docsrs/vareffect)](https://docs.rs/vareffect)
[![CI status](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-genomics/vareffect/actions/workflows/ci.yml)

`vareffect` takes a variant (chromosome, position, reference allele,
alternate allele) and tells you what it does to every transcript that
overlaps it: which protein residues change, whether it disrupts a splice
site, whether it introduces a premature termination codon and is likely to
trigger nonsense-mediated decay, and what the Sequence Ontology consequence
terms are — the same questions VEP answers, returned as strongly typed Rust
values instead of TSV.

The crate is deliberately small. It does one thing — assign SO consequences
and HGVS notation — and leaves regulatory layers, frequency lookups, and
plugin ecosystems to downstream code. If you need a fast, embeddable
consequence predictor inside a Rust pipeline, this is for you. If you need
the full VEP experience including gnomAD enrichment and custom plugins, use
VEP.

## Quick start

Add the library to your project and install the companion CLI that
provisions the reference data:

```toml
# Cargo.toml
[dependencies]
vareffect = "0.1.2"
```

```bash
# One-time data setup: downloads GRCh38, builds the transcript store,
# and writes everything to `data/vareffect/`. Takes ~10 minutes and
# ~3 GB of disk, then you never have to touch it again.
cargo install vareffect-cli
vareffect init
vareffect setup
```

Then annotate a variant against every overlapping transcript:

```rust,no_run
use std::path::Path;
use vareffect::{Consequence, VarEffect};

# fn main() -> Result<(), vareffect::VarEffectError> {
let ve = VarEffect::open(
    Path::new("data/vareffect/transcript_models.bin"),
    Path::new("data/vareffect/GRCh38.bin"),
)?;

// TP53 c.742C>T (p.Arg248Trp) — a well-known hotspot missense variant.
// chr17:7674220 is 0-based (BED / UCSC style).
let results = ve.annotate("chr17", 7_674_220, b"G", b"A")?;

for result in &results {
    println!(
        "{} {}: {} ({})",
        result.transcript,
        result.hgvs_p.as_deref().unwrap_or("-"),
        result
            .consequences
            .iter()
            .map(Consequence::as_str)
            .collect::<Vec<_>>()
            .join(","),
        result.impact,
    );
}
# Ok(())
# }
```

Expected output (one line per overlapping transcript):

```text
NM_000546.6 p.Arg248Trp: missense_variant (MODERATE)
```

For sharing across threads, wrap the loaded `VarEffect` in a
`std::sync::Arc` and clone the `Arc` into each worker — the underlying
stores are `Send + Sync` with zero interior mutability.

## Performance

Single-threaded throughput, measured on a modern x86_64 laptop. All numbers
exclude startup cost (genome load, transcript store parse).

| Tool                          | Language | Variants / sec (1 thread) |
|-------------------------------|----------|---------------------------|
| VEP (`--cache`)               | Perl     |                ~200 – 500 |
| SnpEff                        | Java     |           ~2,000 – 10,000 |
| Nirvana                       | C#       |           ~5,000 – 15,000 |
| bcftools csq                  | C        |          ~10,000 – 50,000 |
| **vareffect**                 | Rust     |         ~50,000 – 200,000 |

The gap is almost entirely I/O: VEP reads BGZF-compressed FASTA through
block decompression per base; `vareffect` memory-maps a flat uppercase
binary and reads bytes directly.

## Features

### Variant consequence prediction
- **24 Sequence Ontology terms** — `missense_variant`, `synonymous_variant`,
  `stop_gained`, `stop_lost`, `start_lost`, `start_retained_variant`,
  `stop_retained_variant`, `frameshift_variant`, `inframe_insertion`,
  `inframe_deletion`, `splice_donor_variant`, `splice_acceptor_variant`,
  `splice_region_variant`, `incomplete_terminal_codon_variant`,
  `5_prime_UTR_variant`, `3_prime_UTR_variant`, `intron_variant`,
  `non_coding_transcript_exon_variant`, `upstream_gene_variant`,
  `downstream_gene_variant`, `intergenic_variant`, `coding_sequence_variant`,
  `transcript_ablation`, `protein_altering_variant`.
- **IMPACT ranking** — `HIGH` / `MODERATE` / `LOW` / `MODIFIER`, matching
  VEP's severity scale. `Consequence` derives `Ord` for sorting.
- **Multi-transcript annotation** — one `ConsequenceResult` per overlapping
  transcript; the caller picks a canonical isoform if desired (every result
  carries MANE Select / MANE Plus Clinical / RefSeq Select tier metadata).
- **NMD prediction** — the 50-nucleotide rule applied to truncating variants
  (`stop_gained`, `frameshift_variant`).
- **Standard and mitochondrial genetic codes** — `chrM` variants translate
  with NCBI table 2 automatically.

### HGVS nomenclature
- **Forward c. / n.** — substitutions, deletions, insertions, duplications,
  delins, intronic offsets (`c.672+1`), 5'/3' UTR offsets (`c.-15`, `c.*42`),
  combined forms (`c.-15+1`).
- **Forward p.** — missense, synonymous, stop gain / loss, start loss
  (`p.Met1?`), frameshift with extension walk (`p.Glu23ValfsTer17`),
  incomplete terminal codon, stop extension (`p.Ter394CysextTer9`). Matches
  VEP's default capitalisation and three-letter code convention.
- **3' normalization** — shift indels to the most 3' equivalent position on
  the coding strand for HGVS notation (matches VEP's `--shift_hgvs`).
- **Reverse c.** — parse an HGVS c. string back to plus-strand 0-based
  genomic coordinates and alleles, round-tripped through the transcript
  store and verified against the genome.

### Variant localization
- Classify variants as CDS exon, intron, 5' / 3' UTR, splice donor /
  acceptor, splice region, upstream / downstream, intergenic.
- Multi-exon indel handling with exon-boundary-spanning logic.
- Reference-allele verification against the loaded genome before
  annotation, so VCFs built against the wrong build fail fast instead of
  silently producing wrong calls.

### Transcript models
- MANE Select, MANE Plus Clinical, and RefSeq Select tiers — your build
  pipeline decides which to ingest.
- Full exon and CDS segment layout with GFF3 phase preserved, so codon
  walks across exon boundaries and frameshift detection are O(1) lookups.
- Interval-tree indexed per chromosome for O(log n + k) overlap queries.
- Strand-aware: correct reverse-complement handling on minus-strand genes.

### Runtime characteristics
- **Memory-mapped reference genome** — base fetches are a pointer
  dereference (~5 ns), not a BGZF block decompression.
- **Thread-safe** — `VarEffect`, `TranscriptStore`, and `FastaReader` are
  all `Send + Sync` with zero interior mutability. A compile-time assertion
  in `lib.rs` guarantees this for every release.
- **Zero external runtime** — no network, no database, no background
  workers. Load two files at startup, share an `Arc<VarEffect>` across your
  worker pool, call `.annotate(...)`.

## What vareffect does not do

- No regulatory / TFBS / motif annotation layer.
- No co-located variant lookup (no gnomAD, ClinVar, dbSNP).
- No allele frequency or population annotation.
- No canonical transcript *selection* — every overlapping transcript is
  returned with its tier metadata; callers decide which to keep.
- No multi-allele VCF splitting — the caller must split comma-separated
  `ALT`s before invoking `annotate`.
- No plugin system.
- No alternate genome builds out of the box — GRCh37 or CHM13 require
  regenerating the transcript and genome binaries with your own build.

See [`VEP_DIVERGENCES.md`](./VEP_DIVERGENCES.md) for the complete list of
intentional divergences from VEP and features that are not yet implemented.

## Setting up the data files

`vareffect` does not ship reference data. You provide two files at runtime:

1. **`transcript_models.bin`** — a MessagePack-serialised
   `Vec<TranscriptModel>` built from a MANE / RefSeq Select GFF3.
2. **`GRCh38.bin`** (or whatever build you use) — a flat uppercase-IUPAC
   binary plus a `.bin.idx` MessagePack sidecar produced by the builder.

### Recommended: use `vareffect-cli`

The companion crate `vareffect-cli` ships a `vareffect` binary that handles
the entire provisioning flow in one command. It downloads the GRCh38
reference FASTA from NCBI, writes the flat-binary genome + index, fetches
the latest MANE release GFF3 + summary, builds the transcript store, and
writes NCBI patch-contig aliases — all under `data/vareffect/`.

```bash
cargo install vareffect-cli

# Scaffold a config file (one-time):
vareffect init

# Full provisioning (GRCh38 genome + MANE transcript models).
vareffect setup

# Validate everything is in place:
vareffect check

# Write runtime files to a custom directory instead of the config default:
vareffect setup --output /data/genomes/vareffect

# Only the reference genome:
vareffect setup --fasta-only

# Only rebuild the transcript model store (reuses an existing genome):
vareffect setup --models-only
```

`setup` is idempotent — source archives are cached in `data/raw/`, the
genome binary is skipped if it already exists, and transcript models are
rebuilt on every run so a new MANE release picks up automatically.

After `vareffect setup` finishes, you have the layout the library expects:

```text
data/vareffect/
  GRCh38.bin               # flat-binary reference genome
  GRCh38.bin.idx           # MessagePack contig index
  transcript_models.bin    # serialised Vec<TranscriptModel>
  patch_chrom_aliases.csv  # UCSC <-> RefSeq patch-contig map
```

### Alternative: roll your own

If you're building a custom store (a different transcript source, a
non-human genome, a subset of the human transcriptome), the
[`fasta::write_genome_binary`] function is public so you can generate the
flat binary yourself:

```rust,no_run
use std::path::Path;
use vareffect::fasta::write_genome_binary;

# fn main() -> Result<(), vareffect::VarEffectError> {
// Uppercase ASCII bytes, one contig per tuple.
let chr_toy: &[u8] = b"ACGTACGTNNN";
let contigs: &[(&str, &[u8])] = &[("chrToy", chr_toy)];

write_genome_binary(
    contigs,
    "toy",                               // build label, stored in the index
    Path::new("out/toy.bin"),            // flat binary
    Path::new("out/toy.bin.idx"),        // MessagePack index sidecar
)?;
# Ok(())
# }
```

Building the transcript store by hand means producing a
`Vec<vareffect::TranscriptModel>` and MessagePack-serialising it with
`rmp-serde`; see `crates/vareffect-cli/src/builders/` in the source tree
for a worked example.

## HGVS reverse resolution

Take an HGVS c. string, resolve it to plus-strand genomic coordinates, and
feed the result straight back into `annotate` — useful when your input is
transcript-relative rather than coordinate-based:

```rust,no_run
# use vareffect::VarEffect;
# fn example(ve: &VarEffect) -> Result<(), vareffect::VarEffectError> {
let gv = ve.resolve_hgvs_c("NM_000546.6:c.742C>T")?;

let results = ve.annotate(&gv.chrom, gv.pos, &gv.ref_allele, &gv.alt_allele)?;
# let _ = results;
# Ok(())
# }
```

`resolve_hgvs_c` supports substitutions, deletions, insertions,
duplications, and delins across CDS, 5' UTR, 3' UTR, and intronic
positions.

## Threading model

`VarEffect`, `TranscriptStore`, and `FastaReader` are all `Send + Sync` —
construct one `VarEffect` at startup, wrap it in `std::sync::Arc`, and
share it across every worker thread or async task. There is no interior
mutability, no contention, and no hidden cloning cost; every thread reads
from the same memory-mapped bytes.

## Coordinate conventions

- **0-based, half-open** (BED / UCSC) for every coordinate in the public
  API. GFF3's 1-based-fully-closed input is converted at build time.
- **UCSC chromosome names** (`chr17`, `chrM`). The `FastaReader` has an
  alias table so it can transparently accept `"17"`, `"NC_000017.11"`, etc.
  Patch-contig lookups additionally need `open_with_patch_aliases` and a
  `refseq,ucsc` alias CSV (written automatically by `vareffect setup`).
- **Uppercase ASCII** allele bytes — no case coercion at call time, so
  lowercase input is a bug in your code, not something we silently fix.

## VEP feature parity

`vareffect` targets near-100% concordance with VEP release 115 / 116 on the
core coding and transcript consequence layer. Every intentional divergence,
not-yet-implemented feature, and out-of-scope design decision is catalogued
in [`VEP_DIVERGENCES.md`](./VEP_DIVERGENCES.md). The validation suite lives
in `tests/vep_concordance_*.rs`.

## Testing

The unit and fast integration tests run with no external data:

```bash
cargo test -p vareffect
```

The VEP concordance tests are `#[ignore]`-gated because they need the
transcript store and reference genome on disk (run `vareffect setup` first
if you haven't already):

```bash
FASTA_PATH="$(pwd)/data/vareffect/GRCh38.bin" CONCORDANCE_THREADS=1 \
    cargo test -p vareffect --release -- --ignored vep_concordance
```

## License

Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or
  <https://www.apache.org/licenses/LICENSE-2.0>)


## Contributing

Contributions are welcome. The most valuable areas for outside help are:

- Adding SO terms currently listed as not-yet-implemented in
  `VEP_DIVERGENCES.md` (splice polypyrimidine-tract, donor region, donor
  5th base).
- Broadening the VEP concordance corpus in `tests/vep_concordance_*.rs`
  with variants that stress a part of the pipeline that isn't already
  covered.
- Alternate genome build support (GRCh37, CHM13) in the `vareffect-cli`
  provisioning flow.

Open an issue before starting on anything larger than a bug fix so we can
agree on scope.
