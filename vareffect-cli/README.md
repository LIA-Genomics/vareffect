# vareffect-cli

[![crates.io](https://img.shields.io/crates/v/vareffect-cli.svg)](https://crates.io/crates/vareffect-cli)
[![CI status](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-genomics/vareffect/actions/workflows/ci.yml)

CLI for the [`vareffect`](https://crates.io/crates/vareffect) variant
consequence prediction engine. Installs as a single binary named
`vareffect` with four capabilities:

1. **`init`** -- scaffold a `vareffect_build.toml` config file at a
   discoverable location
2. **`check`** -- validate that config and reference data files are in
   place and readable
3. **`setup`** -- download GRCh38 (MANE) and/or GRCh37 (RefSeq Select)
   and build the runtime data files that `vareffect` needs. Selectable
   per assembly via `--assembly grch38` (default), `--assembly grch37`,
   or `--assembly all`.
4. **`annotate`** -- annotate a VCF with VEP-compatible CSQ fields
   (~270k variants/sec single-threaded, near-linear multi-thread scaling)

## Install

```bash
cargo install vareffect-cli
```

The binary is named `vareffect`.

## Initialize configuration

Scaffold a `vareffect_build.toml` configuration file:

```bash
vareffect init
```

Writes the default config to `~/.config/vareffect/vareffect_build.toml`
(Linux/macOS). Auto-discovery means `setup`, `check`, and other
subcommands will find it without `--config`.

### Custom locations

```bash
# Write to a specific path:
vareffect init --config /path/to/vareffect_build.toml

# Override the data directory in the generated config:
vareffect init --data-dir /data/genomes/vareffect

# Overwrite an existing config without prompting:
vareffect init --force
```

## Validate setup

Check that config and all reference data files are in place:

```bash
vareffect check
```

Reports pass/fail status for the config file, data directory, genome
binary, genome index, and transcript model store. Exit code 0 if all
pass, 1 if any fail.

## Setup runtime data

Run this first -- it builds the genome and transcript files that
`annotate` requires.

```bash
vareffect setup                       # default: --assembly grch38 (NCBI MANE)
vareffect setup --assembly grch37     # NCBI RefSeq Select on GRCh37.p13
vareffect setup --assembly all        # both, side-by-side under data/vareffect/
```

`--assembly grch38` downloads NCBI's GRCh38.p14 reference FASTA and the
MANE v1.5 GFF3, builds the flat-binary genome and transcript store, and
writes everything to `data/vareffect/` (the default from
`vareffect_build.toml`). `--assembly grch37` does the same for GRCh37.p13
+ RefSeq Select. `--assembly all` builds every assembly the config
defines. Takes ~10 minutes per assembly on broadband. Idempotent --
source archives are cached in `data/raw/`, so repeated runs skip
downloads.

### Custom output directory

Use `--output` to override the output directory from the config file:

```bash
vareffect setup --output /data/genomes/vareffect
```

Runtime files (genome binary, transcript models, patch aliases) are
written to the given directory instead of `[vareffect].output_dir`.

After setup, `data/vareffect/` contains the following per-assembly artifacts
(suffixed with `_grch37` / `_grch38` whenever both can coexist; FASTA binaries
already use the canonical NCBI casing):

| File | Purpose |
|------|---------|
| `GRCh38.bin` + `.idx` (or `GRCh37.bin`) | Flat-binary genome + MessagePack contig index |
| `transcript_models_grch{37,38}.bin` | MANE / RefSeq Select transcript models (MessagePack) |
| `transcript_models_grch{37,38}.manifest.json` | Sibling manifest (records assembly + record count) |
| `patch_chrom_aliases_grch{37,38}.csv` | UCSC <-> RefSeq patch-contig mapping |

### Partial runs

```bash
vareffect setup --fasta-only     # skip transcript model build
vareffect setup --models-only    # skip reference FASTA download
```

Both flags can be combined with `--output`:

```bash
vareffect setup --fasta-only --output /data/genomes/vareffect
```

## Annotate a VCF

```bash
vareffect annotate \
  --assembly grch38 \
  --input sample.vcf.gz \
  --output annotated.vcf.gz \
  --fasta data/vareffect/GRCh38.bin \
  --transcripts data/vareffect/transcript_models_grch38.bin \
  --threads 8
```

Reads a VCF (`.vcf` or `.vcf.gz`), runs consequence prediction on every
variant, and writes an annotated VCF with a `CSQ` INFO field whose
pipe-delimited layout matches VEP's `--vcf` output:

```
Allele|Consequence|IMPACT|SYMBOL|Feature|Feature_type|BIOTYPE|EXON|INTRON|
HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|
STRAND|MANE_SELECT|MANE_PLUS_CLINICAL
```

Downstream tools that parse VEP CSQ fields (SnpSift, bcftools, GEMINI)
work unchanged.

A 100-variant test VCF is included at `tests/data/clinvar_100.vcf`
(sampled from ClinVar across 10 chromosomes):

```bash
vareffect annotate \
  --assembly grch38 \
  --input tests/data/clinvar_100.vcf \
  --output /tmp/annotated.vcf \
  --fasta data/vareffect/GRCh38.bin \
  --transcripts data/vareffect/transcript_models_grch38.bin
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--assembly` | required | `grch38` or `grch37` (must match the data files) |
| `--input` | required | Input VCF (`.vcf` or `.vcf.gz`) |
| `--output` | required | Output VCF (`.vcf` or `.vcf.gz`) |
| `--fasta` | required | `GRCh38.bin` / `GRCh37.bin` flat-binary genome |
| `--transcripts` | required | `transcript_models_grch{37,38}.bin` |
| `--threads` | `1` | Rayon worker threads |
| `--patch-aliases` | none | `patch_chrom_aliases_grch{37,38}.csv` for patch contigs |

### Behaviour

- Chromosome names are normalized automatically (`1` -> `chr1`,
  `MT` -> `chrM`). Both UCSC-prefixed and Ensembl/ClinVar bare names
  work.
- Multi-allelic ALT alleles are annotated independently; CSQ entries are
  comma-separated with the `Allele` field identifying each ALT.
- Intergenic variants pass through without a CSQ field (matching VEP).
- Malformed lines, annotation errors, and panics are logged at WARN and
  the original line is passed through unchanged. The pipeline never
  drops a line.
- `RUST_LOG=debug` enables per-variant reason logging for unannotated
  variants.
- Output `.vcf.gz` is standard gzip (not BGZF). Use `tabix` post-hoc
  if indexing is needed.

## Configuration

URLs and output paths live in `vareffect_build.toml`. Run `vareffect init`
to scaffold one, or it is auto-discovered from:

1. `--config <PATH>`
2. `VAREFFECT_BUILD_CONFIG` env var
3. `crates/vareffect-cli/vareffect_build.toml` (repo checkout)
4. `/etc/vareffect/vareffect_build.toml`
5. `~/.config/vareffect/vareffect_build.toml`

Update `transcript_models_version` and the MANE URLs to target a newer
MANE release.

## `build-transcripts`

Standalone transcript model builder for iterative development:

```bash
vareffect build-transcripts \
  --assembly grch38 \
  --input data/raw/mane_grch38.gff.gz \
  --summary-input data/raw/mane_transcript_models_v1.5.summary.tsv.gz \
  --patch-chrom-aliases data/vareffect/patch_chrom_aliases_grch38.csv \
  --output data/vareffect \
  --version 1.5
```

## Disk requirements

| When | Disk |
|------|------|
| During first `setup` | ~7 GB peak |
| After `setup` (runtime only) | ~3 GB |

Delete `data/raw/` after setup if disk is tight.

## License

Apache License, Version 2.0
([LICENSE-APACHE](LICENSE-APACHE) or
<https://www.apache.org/licenses/LICENSE-2.0>)
