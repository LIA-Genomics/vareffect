# vareffect

[![CI](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml)
[![vareffect](https://img.shields.io/crates/v/vareffect.svg?label=vareffect)](https://crates.io/crates/vareffect)
[![vareffect-cli](https://img.shields.io/crates/v/vareffect-cli.svg?label=vareffect-cli)](https://crates.io/crates/vareffect-cli)
[![docs.rs](https://docs.rs/vareffect/badge.svg)](https://docs.rs/vareffect)
[![Apache-2.0](https://img.shields.io/crates/l/vareffect)](#license)

Variant consequence prediction and HGVS notation in Rust, matching
[Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/) release 115/116.

Give vareffect a variant (chromosome, position, ref, alt) and it returns, for
every transcript that overlaps it or lies within 5 kb: Sequence Ontology
consequence terms, IMPACT, HGVS c./n./p., protein and cDNA positions, and an
NMD prediction. It runs in-process, memory-maps the genome and is safe to share
across threads.

| Crate | What it is |
|---|---|
| [`vareffect`](vareffect/) | The library: consequence prediction, HGVS, transcript and genome stores |
| [`vareffect-cli`](vareffect-cli/) | The `vareffect` binary: builds the runtime data and annotates VCFs |

## Quick start

```bash
cargo install vareffect-cli
vareffect init     # write a config file (once)
vareffect setup    # download GRCh38 + transcripts, build the stores (~10 min, ~3 GB)
vareffect check    # confirm everything is in place

vareffect annotate \
  --input sample.vcf.gz --output annotated.vcf.gz \
  --fasta data/vareffect/GRCh38.bin \
  --transcripts data/vareffect/transcript_models.bin
```

The output VCF carries a VEP-style `CSQ` INFO field, so tools that read VEP
output (bcftools, SnpSift, GEMINI) work unchanged.

As a library:

```toml
[dependencies]
vareffect = "0.6"
```

```rust
use std::path::Path;
use vareffect::VarEffect;

let ve = VarEffect::open(
    Path::new("data/vareffect/transcript_models.bin"),
    Path::new("data/vareffect/GRCh38.bin"),
)?;

// TP53 c.742C>T, p.Arg248Trp. Positions are 0-based.
for r in ve.annotate("chr17", 7_674_220, b"G", b"A")? {
    println!("{} {:?} {:?}", r.transcript, r.consequences, r.hgvs_p);
}
```

## Transcript sources

`vareffect setup` builds one of three transcript sets, chosen in the config
(`[vareffect.transcripts] source`):

| Source | Transcripts | Store file |
|---|---|---|
| `mane` (default) | MANE Select and MANE Plus Clinical (v1.5, 19,437) | `transcript_models.bin` |
| `refseq` | Curated RefSeq `NM_`/`NR_` on the primary assembly (RS_2025_08, 93,432) | `transcript_models_refseq.bin` |
| `gencode` | GENCODE comprehensive (v50, 644,292) | `transcript_models_gencode.bin` |

RefSeq and GENCODE downloads are pinned by SHA-256. Each store records its
source and release, which `vareffect check` and the annotated VCF header
report. See the [CLI README](vareffect-cli/README.md) for details.

## Accuracy

Measured against VEP 116 (REST) on 11,499 variants: ClinVar records, ClinVar
indels and HG002 benchmark calls. Rates are per transcript pair, after the
documented differences in [`VEP_DIVERGENCES.md`](vareffect/VEP_DIVERGENCES.md)
are folded.

| Store | Consequence terms | HGVS c. | HGVS p. | IMPACT |
|---|---:|---:|---:|---:|
| MANE | 99.54% | 99.98% | 99.99% | 100.00% |
| RefSeq | 99.45% | 99.95% | 99.98% | 99.99% |
| GENCODE | 99.33% | 99.97% | 100.00% | 100.00% |

Every known difference is listed in
[`VEP_DIVERGENCES.md`](vareffect/VEP_DIVERGENCES.md), either as intentional or
as a defect to fix.

## Speed

Annotating all 4.43 million GRCh38 ClinVar records against MANE, including
reading the VCF and writing compressed output, takes about 10 s on one thread
(~430,000 records/s) and about 3 s on 8 threads (Apple M5 Max). The
genome is memory-mapped, so a base lookup is a pointer read rather than a
compressed-block decode.

## Documentation

- [Library README](vareffect/README.md): API, result fields, data files, threading
- [CLI README](vareffect-cli/README.md): commands, configuration, VCF output
- [`VEP_DIVERGENCES.md`](vareffect/VEP_DIVERGENCES.md): every difference from VEP
- [`CHANGELOG.md`](CHANGELOG.md): release notes
- [docs.rs/vareffect](https://docs.rs/vareffect): full API reference

## License

Licensed under [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).
