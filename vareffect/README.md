# vareffect

Variant consequence prediction and HGVS notation in Rust, matching
[Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/) release 115/116.

[![crates.io](https://img.shields.io/crates/v/vareffect.svg)](https://crates.io/crates/vareffect)
[![docs.rs](https://img.shields.io/docsrs/vareffect)](https://docs.rs/vareffect)
[![CI status](https://github.com/LIA-Genomics/vareffect/actions/workflows/ci.yml/badge.svg)](https://github.com/LIA-genomics/vareffect/actions/workflows/ci.yml)

`vareffect` tells you what a variant does to each nearby transcript: the
Sequence Ontology consequence terms and their IMPACT, the HGVS c./n./p.
names, the affected cDNA, CDS and protein positions, and whether a truncating
variant is likely to trigger nonsense-mediated decay. Results are typed Rust
values.

It does one job. Regulatory annotation, allele frequencies and plugins are
left to downstream code; if you need those, use VEP.

## Quick start

```toml
[dependencies]
vareffect = "0.6"
```

The library ships no reference data. Build it once with the companion CLI
(about 10 minutes and 3 GB of disk):

```bash
cargo install vareffect-cli
vareffect init
vareffect setup
```

Then annotate:

```rust,no_run
use std::path::Path;
use vareffect::{Consequence, VarEffect};

# fn main() -> Result<(), vareffect::VarEffectError> {
let ve = VarEffect::open(
    Path::new("data/vareffect/transcript_models.bin"),
    Path::new("data/vareffect/GRCh38.bin"),
)?;

// TP53 c.742C>T (p.Arg248Trp). Positions are 0-based.
let results = ve.annotate("chr17", 7_674_220, b"G", b"A")?;

for result in &results {
    let terms: Vec<_> = result.consequences.iter().map(Consequence::as_str).collect();
    println!(
        "{} {}: {} ({})",
        result.transcript,
        result.hgvs_p.as_deref().unwrap_or("-"),
        terms.join(","),
        result.impact,
    );
}
# Ok(())
# }
```

```text
NM_000546.6 p.Arg248Trp: missense_variant (MODERATE)
```

## What a result contains

`annotate` returns one `ConsequenceResult` per transcript that overlaps the
variant or lies within 5,000 bases of it. A transcript the variant does not
overlap gets `upstream_gene_variant` or `downstream_gene_variant` and a
`distance`, as VEP's default `--distance 5000` does. With no transcript in
range, the result is a single `intergenic_variant` entry with an empty
`transcript`.

| Field | Meaning |
|---|---|
| `consequences`, `impact` | SO terms, most severe first, and the highest IMPACT (`HIGH` to `MODIFIER`) |
| `hgvs_c`, `hgvs_p` | HGVS names, 3'-shifted as VEP does (`NM_000546.6:c.742C>T`, `p.Arg248Trp`) |
| `cdna_position`, `cds_position`, `protein_start` (and `_end`) | 1-based positions, VEP style |
| `codons`, `amino_acids`, `exon`, `intron` | As in VEP's output (`Cgg/Tgg`, `R/W`, `7/11`) |
| `distance` | Bases to the transcript, for up/downstream results |
| `is_mane_select`, `is_mane_plus_clinical`, `is_refseq_select` | Transcript tier, for choosing a canonical result |
| `predicts_nmd` | 50-nucleotide rule at the variant site (VEP `NMD.pm` parity) |
| `ptc` | The new termination codon and its own NMD verdict: use this for ACMG PVS1 |

`predicts_nmd` and `ptc` differ for frameshifts, whose new stop lies
downstream of the variant. `ptc` is `Indeterminate` when the stop cannot be
located, never a guessed negative.

## Transcript stores

A store holds one transcript set, built by `vareffect setup` or
`vareffect build-transcripts`:

- **MANE** (default): MANE Select and MANE Plus Clinical.
- **RefSeq**: curated `NM_`/`NR_` transcripts on the primary assembly.
- **GENCODE**: GENCODE comprehensive.

The store records its source and release
(`ve.transcripts().source()`, `ve.transcripts().source_release()`). The file
starts with a versioned header; a store from another format version, including
any written by vareffect 0.5 or earlier, is refused at load time. Rebuild it
with `vareffect setup`.

RefSeq and GENCODE include transcripts MANE does not have, and vareffect
handles them as VEP does:

- **Incomplete CDSs** (GENCODE `cds_start_NF`/`cds_end_NF`, partial RefSeq
  CDSs): numbered in VEP's padded frame; a partial codon gets `p.?`.
- **Short exon gaps** (13 bases or fewer) are frameshift introns: no splice or
  intron terms inside them.
- **RefSeq transcripts that differ from the genome**: numbers the genome
  alignment would shift are withheld rather than reported wrong.
- **Selenocysteine codons**: `TGA` there reads as `Sec`.
- **Pseudoautosomal transcripts**: one copy on chrX and one on chrY.

The details, with VEP examples, are in
[`VEP_DIVERGENCES.md`](./VEP_DIVERGENCES.md).

## Other entry points

- **VEP JSON**: `annotate_to_vep_json` returns the same results shaped like
  VEP REST output.
- **Structural variants**: `annotate_interval` (DEL, DUP, INV),
  `annotate_breakend` (`<BND>`) and `annotate_sv_insertion` (`<INS>`) return
  VEP's region-level terms (`transcript_ablation`, `feature_truncation`, ...)
  with exon-overlap geometry for CNV scoring.
- **HGVS back to coordinates**: `resolve_hgvs_c` maps a transcript-relative
  name to genomic coordinates; `resolve_hgvs_g` takes a `g.` or `m.` name and
  returns the left-aligned, VCF-canonical variant.

```rust,no_run
# use vareffect::VarEffect;
# fn example(ve: &VarEffect) -> Result<(), vareffect::VarEffectError> {
let gv = ve.resolve_hgvs_c("NM_000546.6:c.742C>T")?;
let results = ve.annotate(&gv.chrom, gv.pos, &gv.ref_allele, &gv.alt_allele)?;
# let _ = results;
# Ok(())
# }
```

## Conventions

- **Coordinates** are 0-based and half-open (BED style) everywhere in the API.
- **Chromosomes** use UCSC names (`chr17`, `chrM`); `17` and `NC_000017.11`
  are accepted too. Patch contigs need `VarEffect::open_with_patch_aliases`
  and the alias CSV that `vareffect setup` writes.
- **Alleles** are uppercase plus-strand bytes. The reference allele is checked
  against the genome before annotation; a mismatch is an error.
- **`chrM`** translates with the vertebrate mitochondrial code (NCBI table 2).
- **Threads**: `VarEffect` is `Send + Sync` with no interior mutability. Load
  it once, wrap it in `Arc` and share it.

## Building the data yourself

`vareffect setup` is the supported path: it downloads the GRCh38 FASTA and
the configured transcript source and writes everything to `data/vareffect/`.
To build your own genome or store (another source, another species, a subset),
use the same public functions the CLI uses:

```rust,no_run
use std::path::Path;
use vareffect::fasta::write_genome_binary;
use vareffect::{TranscriptModel, TranscriptSource, encode_transcript_store};

# fn build_models() -> Vec<TranscriptModel> { unimplemented!() }
# fn main() -> Result<(), Box<dyn std::error::Error>> {
// Genome: uppercase ASCII, one contig per tuple.
let contigs: &[(&str, &[u8])] = &[("chrToy", &b"ACGTACGTNNN"[..])];
write_genome_binary(contigs, "toy", Path::new("toy.bin"), Path::new("toy.bin.idx"))?;

// Transcripts: build `TranscriptModel`s with `TranscriptModel::new` and its
// `with_*` setters (an empty store is refused), then encode them.
let models = build_models();
let bytes = encode_transcript_store(Some(TranscriptSource::Mane), Some("1.5"), &models)?;
std::fs::write("transcript_models.bin", &bytes)?;
# Ok(())
# }
```

The GFF3 builder in `vareffect-cli/src/builders/transcript_models/` is a
complete example.

## Limits

No regulatory annotation, no co-located variant or frequency lookup, no
plugin system and no canonical-transcript selection (every result carries its
tier; you choose). Split multi-allelic records before calling `annotate`.
GRCh37 and CHM13 need their own genome and transcript builds. Everything that
differs from VEP is in [`VEP_DIVERGENCES.md`](./VEP_DIVERGENCES.md).

## Testing

```bash
cargo test -p vareffect
```

The VEP concordance suites need the built data and are `#[ignore]`-gated:

```bash
FASTA_PATH="$PWD/data/vareffect/GRCh38.bin" \
    cargo test -p vareffect --release --no-fail-fast --test '*' -- --ignored
```

## License

Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
<https://www.apache.org/licenses/LICENSE-2.0>).

## Contributing

Contributions are welcome, especially VEP concordance cases that exercise an
untested path and alternate genome builds in the CLI. Open an issue before
starting anything larger than a bug fix.
