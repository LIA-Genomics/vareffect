# vareffect-cli

[![crates.io](https://img.shields.io/crates/v/vareffect-cli.svg)](https://crates.io/crates/vareffect-cli)
[![CI status](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-genomics/vareffect/actions/workflows/ci.yml)

The command-line tool for the [`vareffect`](https://crates.io/crates/vareffect)
consequence predictor. It builds the runtime data the library needs and
annotates VCFs.

```bash
cargo install vareffect-cli    # installs a binary named `vareffect`
```

| Command | What it does |
|---|---|
| `vareffect init` | Write a `vareffect_build.toml` config file |
| `vareffect setup` | Download the genome and transcripts and build the stores |
| `vareffect check` | Validate the config and the built files |
| `vareffect annotate` | Annotate a VCF with VEP-style `CSQ` fields |
| `vareffect build-transcripts` | Build a transcript store from a GFF3 you already have |

## First run

```bash
vareffect init     # writes vareffect_build.toml to your user config directory
vareffect setup    # ~10 minutes, ~7 GB peak disk, ~3 GB afterwards
vareffect check    # exit code 0 when everything is in place
```

`setup` runs three steps: the GRCh38 genome (NCBI GRCh38.p14, converted to
a memory-mappable flat binary), the patch-contig alias table, and the
transcript store for the configured source. Downloads are cached in
`data/raw/`, so a rerun skips them; the transcript store is rebuilt every
time.

```bash
vareffect setup --output /data/vareffect   # write somewhere other than the config's output_dir
vareffect setup --fasta-only               # genome only
vareffect setup --models-only              # transcript store only (reuses the genome)
```

After setup, `data/vareffect/` holds:

| File | Contents |
|---|---|
| `GRCh38.bin`, `GRCh38.bin.idx` | Flat-binary genome and its contig index |
| `transcript_models*.bin` | Transcript store (name depends on the source, below) |
| `*.manifest.json`, `*.sha256`, `*.excluded.tsv` | Build provenance: source, release, checksums, counts, transcripts not admitted and why |
| `patch_chrom_aliases.csv` | UCSC to RefSeq names for patch contigs |

`vareffect init` options: `--config <PATH>` to write elsewhere,
`--data-dir <DIR>` to set the output directory, `--force` to overwrite.

## Choosing a transcript source

Set the source in the config before running `setup`:

```toml
[vareffect.transcripts]
source = "mane"    # or "refseq" or "gencode"
```

| Source | Transcripts | Store file | Store size |
|---|---|---|---:|
| `mane` (default) | MANE Select and Plus Clinical, v1.5 (19,437) | `transcript_models.bin` | 30 MB |
| `refseq` | Curated `NM_`/`NR_` on the primary assembly, RS_2025_08 (93,432) | `transcript_models_refseq.bin` | 135 MB |
| `gencode` | GENCODE comprehensive, v50 (644,292) | `transcript_models_gencode.bin` | 714 MB |

The RefSeq and GENCODE files are pinned in the config's `[vareffect.refseq]`
and `[vareffect.gencode]` tables (`url`, `release`, `sha256`). `setup` refuses
a download whose SHA-256, release or genome build does not match the pin, so an
upstream change cannot silently alter your store. `check` fails when the
config, the store and its manifest disagree on the source or release.

GENCODE is large: annotating ClinVar with it peaks at about 2.4 GB of memory,
against about 1 GB for MANE.

## Annotating a VCF

```bash
vareffect annotate \
  --input sample.vcf.gz \
  --output annotated.vcf.gz \
  --fasta data/vareffect/GRCh38.bin \
  --transcripts data/vareffect/transcript_models.bin
```

| Flag | Default | Meaning |
|---|---|---|
| `--input` | required | Input VCF (`.vcf` or `.vcf.gz`) |
| `--output` | required | Output VCF; `.vcf.gz` is written as BGZF |
| `--fasta` | required | `GRCh38.bin` |
| `--transcripts` | required | The store to annotate against |
| `--threads` | `0` | Worker threads; `0` uses every core |
| `--patch-aliases` | none | `patch_chrom_aliases.csv`, for patch contigs |

The output adds two header lines and a `CSQ` INFO field per record:

```text
##vareffect_transcripts=mane 1.5
##INFO=<ID=CSQ,...,Format: Allele|Consequence|IMPACT|SYMBOL|Feature|Feature_type|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|MANE_SELECT|MANE_PLUS_CLINICAL">
```

- Each `CSQ` entry is one transcript, with the same fields as VEP's `--vcf`
  output, so bcftools, SnpSift and GEMINI read it unchanged.
- Transcripts within 5 kb of a variant that they do not overlap get an
  `upstream_gene_variant` or `downstream_gene_variant` entry. A record with no
  transcript in range is written without `CSQ`.
- Chromosome names are normalized (`1` to `chr1`, `MT` to `chrM`).
- Multi-allelic records are annotated per ALT; the `Allele` field says which.
- No record is dropped. A malformed line or an annotation error is logged as a
  warning and the line is written unchanged (so, without `CSQ`).
  `RUST_LOG=debug` logs the reason for every unannotated variant.
- `.vcf.gz` output is BGZF, compressed in parallel and indexable with
  `tabix -p vcf`.

## Building a store from your own GFF3

`build-transcripts` builds a store without downloading anything. `--source`
selects the admission rules:

```bash
vareffect build-transcripts --source mane \
  --input MANE.GRCh38.v1.5.refseq_genomic.gff.gz \
  --summary-input MANE.GRCh38.v1.5.summary.txt.gz \
  --patch-chrom-aliases data/vareffect/patch_chrom_aliases.csv \
  --output data/vareffect --version 1.5

vareffect build-transcripts --source refseq \
  --input GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
  --output data/vareffect --version RS_2025_08

vareffect build-transcripts --source gencode \
  --input gencode.v50.annotation.gff3.gz \
  --output data/vareffect --version 50
```

For MANE, `--summary-input` cross-checks every transcript against the MANE
summary file and fails the build on a mismatch. Every build writes the
manifest and `excluded.tsv` listed above.

## Configuration file

`vareffect_build.toml` holds the download URLs, the transcript source and the
output directory. Commands look for it in this order:

1. `--config <PATH>`
2. the `VAREFFECT_BUILD_CONFIG` environment variable
3. `vareffect-cli/vareffect_build.toml` (in a repository checkout)
4. `/etc/vareffect/vareffect_build.toml`
5. `vareffect/vareffect_build.toml` in the user config directory
   (`~/.config` on Linux, `~/Library/Application Support` on macOS)

The file written by `vareffect init` documents every key. Unknown keys inside
`[vareffect]` are an error, so a misspelled table cannot fall back to a
default.

## License

Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
<https://www.apache.org/licenses/LICENSE-2.0>).
