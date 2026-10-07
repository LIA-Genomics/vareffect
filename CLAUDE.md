# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```bash
make all          # full local gate: harness-check fmt lint test doc-check build
make fmt          # cargo fmt --all
make fmt-check    # check formatting without modifying
make lint         # cargo clippy --workspace --all-targets -- -D warnings (as CI)
make test         # cargo test --workspace
make doc-check    # rustdoc with warnings denied (as CI)
make harness-check # agent harness structure and its tests
make build        # debug build
make release      # release build
make check        # cargo check --workspace
make install-cli  # install vareffect binary locally
```

Run a single test:
```bash
cargo test -p vareffect test_name
```

Integration tests are `#[ignore]`-gated and require runtime data files (`data/vareffect/GRCh38.bin`, `data/vareffect/transcript_models.bin`). To run them:
```bash
FASTA_PATH="$PWD/data/vareffect/GRCh38.bin" cargo test -p vareffect --release --test '*' -- --ignored --nocapture
```
A test that prints `SKIP:` did not run. Details and how to read results: [vareffect-verification](.claude/skills/vareffect-verification/SKILL.md).

Generate data files with `vareffect setup` (one-time, ~10 min, ~3 GB disk).

## Agent workflow

The session agent is the main engineer and the only editor; the `.claude/agents/` reviewers and researcher are read-only. Overview, model policy and controls: [`.claude/README.md`](.claude/README.md). [`.claude/INVARIANTS.md`](.claude/INVARIANTS.md) is the maintainer-owned charter: code shows current behavior, the charter required behavior.

Work through [vareffect-ground](.claude/skills/vareffect-ground/SKILL.md) (scope, tier T0–T3, plan review from T2), [vareffect-review](.claude/skills/vareffect-review/SKILL.md) (snapshot, fresh parallel review, bounded repair) and [vareffect-sequence](.claude/skills/vareffect-sequence/SKILL.md) (several authorized tasks). Spawning reviewers is authorized and expected when the workflow calls for it; only T0 skips independent review. Any change that can move annotation output is measured before and after ([concordance-measurement](.claude/skills/concordance-measurement/SKILL.md)). Accepting a VEP divergence, a transcript-source policy, a breaking change or a release is the maintainer's decision; commits and releases each need their own authorization.

Commit only when authorized: short imperative subject, at most two sentences, only `Co-Authored-By:` trailers (no `Claude-Session:`).

## Architecture

**Workspace layout:** Two crates — `vareffect` (core library) and `vareffect-cli` (data provisioning + VCF annotation CLI).

### Core library (`vareffect/`)

The annotation pipeline flows through three main components:

1. **`TranscriptStore`** (`transcript.rs`) — In-memory store of MANE/RefSeq Select transcript models loaded from MessagePack binary. Indexed by COITree per chromosome for O(log n + k) overlap queries and by accession HashMap for O(1) lookup. Arc-backed for cheap thread sharing.

2. **`FastaReader`** (`fasta.rs`) — Memory-mapped flat binary reference genome (~3.1 GB for GRCh38). Zero-copy random access (~5 ns per base). Maps UCSC/NCBI/Ensembl chromosome naming conventions automatically.

3. **`VarEffect`** (`var_effect.rs`) — Stateful entrypoint bundling TranscriptStore + FastaReader. Main API: `annotate(chrom, pos, ref_allele, alt_allele) → Vec<ConsequenceResult>`. Wrap in `Arc<VarEffect>` for multi-threaded use.

**Annotation pipeline** for a single variant:
- `VarEffect::annotate` queries overlapping transcripts via interval tree
- For each transcript: `locate_variant`/`locate_indel` (`locate/`) classifies position (CDS, intron, UTR, splice site)
- Consequence assignment (`consequence/snv.rs`, `indel.rs`, `complex.rs`) translates ref/alt codons and assigns SO terms
- HGVS notation generated (`hgvs_c.rs` for c./n., `hgvs_p.rs` for p.)
- NMD prediction (`consequence/nmd.rs`) via 50-nucleotide rule on truncating variants

### CLI (`vareffect-cli/`)

- **`setup`** — Downloads GRCh38 FASTA + MANE GFF3, builds flat binary genome and MessagePack transcript store. Idempotent.
- **`annotate`** — Parallel VCF annotation (rayon). Writes CSQ INFO field in VEP `--vcf` format.
- **`models`** — Standalone transcript model builder from GFF3.
- Config in `vareffect_build.toml` (download URLs, output paths).

## Critical Conventions

**All coordinates are 0-based, half-open** (BED/UCSC style). GFF3 input (1-based, fully-closed) is converted at build time. Getting this wrong silently produces off-by-one annotation errors.

**Chromosome names** are UCSC-style (`chr17`, `chrM`). The `chrom` module handles UCSC ↔ RefSeq accession mapping.

**`chrM` uses NCBI genetic code table 2** (vertebrate mitochondrial). Standard chromosomes use table 1. The `codon` module handles this automatically.

## Design Principles

- **No embedded data** — Library ships no reference genomes or transcripts; all data built offline via CLI
- **Thread-safe by default** — `VarEffect`, `TranscriptStore`, `FastaReader` are `Send + Sync` (compile-time assertion in `lib.rs`)
- **VEP concordance** — Targets Ensembl VEP release 115/116 output; intentional divergences documented in `vareffect/VEP_DIVERGENCES.md`
- **Biotype forward compatibility** — `Biotype` enum uses `Other(String)` for unknown labels so new upstream biotypes don't break deserialization

## Code rules

Hard rules for every edit. The full conventions (naming, doc-comment shape, type and function design, ownership, concurrency, imports, testing) are in [idiomatic-rust](.claude/skills/idiomatic-rust/SKILL.md): load it before writing Rust. Reviewers have it preloaded.

- `make all` passes with no errors or warnings before a task is handed off (fmt, clippy `-D warnings` on all targets, tests, rustdoc `-D warnings`).
- Every public function, struct, enum and method has a doc comment covering parameters, return value and errors, with an example for complex APIs.
- No `.unwrap()` in library or production code; `.expect("why this holds")` only for genuine invariants. Fallible operations return `Result`: `thiserror` types in the library, `anyhow` with `.context()` in the CLI. Never swallow an error.
- No `unsafe` without a `// SAFETY:` comment stating the invariant.
- Exhaustive `match`; avoid the catch-all `_`.
- Unit tests for every new function and type; no commented-out tests.
- No wildcard imports except preludes and `use super::*` in tests.
- No emoji or emoji-like unicode, except in tests of multibyte input.
