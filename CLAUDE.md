# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```bash
make fmt          # cargo fmt --all
make fmt-check    # check formatting without modifying
make lint         # cargo clippy --workspace -- -D warnings
make test         # cargo test --workspace
make build        # debug build
make release      # release build
make check        # cargo check --workspace
make install-cli  # install vareffect binary locally
```

Run a single test:
```bash
cargo test -p vareffect test_name
```

Integration tests are `#[ignore]`-gated and require runtime data files (`data/vareffect/GRCh38.bin`, `data/vareffect/transcript_models_grch38.bin`, and the matching `data/vareffect/patch_chrom_aliases_grch38.csv`). The GRCh38 helpers also read the `GRCH38_FASTA` env var (and optionally `GRCH38_TRANSCRIPTS`); GRCh37-suffixed tests use `GRCH37_FASTA` / `GRCH37_TRANSCRIPTS`. To run them:
```bash
GRCH38_FASTA=data/vareffect/GRCh38.bin cargo test -p vareffect -- --ignored
```

Generate data files with `vareffect setup --assembly grch38` (one-time, ~10 min, ~3 GB disk). Add `--assembly grch37` for the GRCh37 build, or `--assembly all` to build both.

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

## Code Style and Formatting

- **MUST** use meaningful, descriptive variable and function names
- **MUST** follow Rust API Guidelines and idiomatic Rust conventions
- **MUST** use 4 spaces for indentation (never tabs)
- **NEVER** use emoji, or unicode that emulates emoji (e.g. ✓, ✗). The only exception is when writing tests and testing the impact of multibyte characters.
- Use snake_case for functions/variables/modules, PascalCase for types/traits, SCREAMING_SNAKE_CASE for constants
- Limit line length to 100 characters (rustfmt default)

## Documentation

- **MUST** include doc comments for all public functions, structs, enums, and methods
- **MUST** document function parameters, return values, and errors
- Keep comments up-to-date with code changes
- Include examples in doc comments for complex functions

Example doc comment:

````rust
/// Calculate the total cost of items including tax.
///
/// # Arguments
///
/// * `items` - Slice of item structs with price fields
/// * `tax_rate` - Tax rate as decimal (e.g., 0.08 for 8%)
///
/// # Returns
///
/// Total cost including tax
///
/// # Errors
///
/// Returns `CalculationError::EmptyItems` if items is empty
/// Returns `CalculationError::InvalidTaxRate` if tax_rate is negative
///
/// # Examples
///
/// ```
/// let items = vec![Item { price: 10.0 }, Item { price: 20.0 }];
/// let total = calculate_total(&items, 0.08)?;
/// assert_eq!(total, 32.40);
/// ```
pub fn calculate_total(items: &[Item], tax_rate: f64) -> Result<f64, CalculationError> {
````

## Type System

- **MUST** leverage Rust's type system to prevent bugs at compile time
- **NEVER** use `.unwrap()` in library code; use `.expect()` only for invariant violations with a descriptive message
- **MUST** use meaningful custom error types with `thiserror`
- Use newtypes to distinguish semantically different values of the same underlying type
- Prefer `Option<T>` over sentinel values

## Error Handling

- **NEVER** use `.unwrap()` in production code paths
- **MUST** use `Result<T, E>` for fallible operations
- **MUST** use `thiserror` for defining error types and `anyhow` for application-level errors
- **MUST** propagate errors with `?` operator where appropriate
- Provide meaningful error messages with context using `.context()` from `anyhow`

## Function Design

- **MUST** keep functions focused on a single responsibility
- **MUST** prefer borrowing (`&T`, `&mut T`) over ownership when possible
- Limit function parameters to 5 or fewer; use a config struct for more
- Return early to reduce nesting
- Use iterators and combinators over explicit loops where clearer

## Struct and Enum Design

- **MUST** keep types focused on a single responsibility
- **MUST** derive common traits: `Debug`, `Clone`, `PartialEq` where appropriate
- Use `#[derive(Default)]` when a sensible default exists
- Prefer composition over inheritance-like patterns
- Use builder pattern for complex struct construction
- Make fields private by default; provide accessor methods when needed

## Testing

- **MUST** write unit tests for all new functions and types
- **MUST** mock external dependencies (APIs, databases, file systems)
- **MUST** use the built-in `#[test]` attribute and `cargo test`
- Follow the Arrange-Act-Assert pattern
- Do not commit commented-out tests
- Use `#[cfg(test)]` modules for test code

## Imports and Dependencies

- **MUST** avoid wildcard imports (`use module::*`) except for preludes, test modules (`use super::*`), and prelude re-exports
- **MUST** document dependencies in `Cargo.toml` with version constraints
- Use `cargo` for dependency management
- Organize imports: standard library, external crates, local modules
- Use `rustfmt` to automate import formatting

## Rust Best Practices

- **NEVER** use `unsafe` unless absolutely necessary; document safety invariants when used
- **MUST** call `.clone()` explicitly on non-`Copy` types; avoid hidden clones in closures and iterators
- **MUST** use pattern matching exhaustively; avoid catch-all `_` patterns when possible
- **MUST** use `format!` macro for string formatting
- Use iterators and iterator adapters over manual loops
- Use `enumerate()` instead of manual counter variables
- Prefer `if let` and `while let` for single-pattern matching

## Memory and Performance

- **MUST** avoid unnecessary allocations; prefer `&str` over `String` when possible
- **MUST** use `Cow<'_, str>` when ownership is conditionally needed
- Use `Vec::with_capacity()` when the size is known
- Prefer stack allocation over heap when appropriate
- Use `Arc` and `Rc` judiciously; prefer borrowing

## Concurrency

- **MUST** use `Send` and `Sync` bounds appropriately
- **MUST** prefer `tokio` for async runtime in async applications
- **MUST** use `rayon` for CPU-bound parallelism
- Avoid `Mutex` when `RwLock` or lock-free alternatives are appropriate
- Use channels (`mpsc`, `crossbeam`) for message passing

## Security

- **NEVER** store secrets, API keys, or passwords in code. Only store them in `.env`.
    - Ensure `.env` is declared in `.gitignore`.
- **MUST** use environment variables for sensitive configuration via `dotenvy` or `std::env`
- **NEVER** log sensitive information (passwords, tokens, PII)
- Use `secrecy` crate for sensitive data types

## Tools

- **MUST** use `rustfmt` for code formatting
- **MUST** use `clippy` for linting and follow its suggestions
- **MUST** ensure code compiles with no warnings (use `-D warnings` flag in CI, not `#![deny(warnings)]` in source)
- Use `cargo` for building, testing, and dependency management
- Use `cargo test` for running tests
- Use `cargo doc` for generating documentation
- Use `make all` after every task and make sure there are no errors and warning.

---

**Remember:** Prioritize clarity and maintainability over cleverness.