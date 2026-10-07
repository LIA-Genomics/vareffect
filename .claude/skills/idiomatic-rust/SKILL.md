---
name: idiomatic-rust
description: >-
  vareffect's Rust conventions — naming and style, rustdoc, types and errors, function and
  type design, ownership and performance, concurrency, imports, testing, comments. Load
  before writing or reviewing Rust; the hard rules are repeated in CLAUDE.md.
---

# Rust conventions

The full rule set for vareffect code. `CLAUDE.md` repeats only the hard rules that apply on every edit; everything else lives here. Follow the Rust API Guidelines; where a rule below conflicts with the surrounding code, match the surrounding code and say so. Prioritize clarity and maintainability over cleverness.

## Style and naming

- Meaningful, descriptive names. `snake_case` for functions, variables and modules, `PascalCase` for types and traits, `SCREAMING_SNAKE_CASE` for constants.
- `rustfmt` formatting: 4-space indentation, 100-column lines.
- No emoji or emoji-like unicode (e.g. check marks) anywhere, except in tests that exercise multibyte input.
- Use `format!` for string formatting; `if let` / `while let` for single-pattern matching; `enumerate()` instead of manual counters.

## Documentation and comments

- **Every public function, struct, enum and method has a doc comment** (published library; rustdoc builds with `-D warnings`). Document parameters, return value and errors; add an example for complex APIs. Keep docs in step with code. The shape (illustrative function, not a real API):

  ````rust
  /// Translate a codon under the given NCBI genetic code table.
  ///
  /// # Arguments
  ///
  /// * `codon` - Three uppercase bases on the coding strand
  /// * `table` - Genetic code (table 1 nuclear, table 2 for `chrM`)
  ///
  /// # Returns
  ///
  /// The one-letter amino acid, `*` for a stop codon.
  ///
  /// # Errors
  ///
  /// Returns `VarEffectError::Malformed` if `codon` is not three valid bases.
  ````
- **Inline comments** explain **why**, not **what**: non-obvious code, or a citation (a spec clause, VEP behavior, a paper), at most 2-3 sentences. Match the density of the surrounding code; delete comments that restate it.

## Types and API design

- Let the type system prevent bugs: enums over stringly-typed fields, newtypes for semantically different values of one underlying type, `Option<T>` over sentinel values.
- **Exhaustive `match`; avoid the catch-all `_`** so new variants surface as compile errors.
- Derive `Debug`, `Clone`, `PartialEq` where appropriate, `Default` when a sensible default exists.
- Types and functions have a single responsibility. Private fields by default, with accessors; builders for complex construction; composition over inheritance-like patterns.
- At most 5 function parameters (a config struct beyond that); return early to reduce nesting.
- Public types that may grow fields are `#[non_exhaustive]` with a constructor (see `ConsequenceResult`); changing a public shape or meaning is a breaking change (`.claude/INVARIANTS.md` V-6).

## Errors

- `Result<T, E>` for every fallible operation, propagated with `?`.
- `thiserror` error types in the library (`vareffect/src/error.rs`); `anyhow` with `.context("what we were doing")` in the CLI.
- **No `.unwrap()` in library code or production paths.** `.expect("...")` only for a genuine invariant, with a message saying why it holds.
- **Never silently swallow an error**: a discarded `Err` (`let _ = fallible();`) is a defect unless deliberate and commented. An error must not become an annotation (V-3).

## Ownership and performance

The per-variant path is hot (memory-mapped genome, interval-tree query, codon translation): allocations there are measurable.

- Prefer borrowing (`&T`, `&mut T`, `&str`) over ownership; `Cow<'_, str>` when ownership is conditional.
- `Vec::with_capacity` when the size is known; prefer stack over heap where it is natural.
- Iterators and combinators over manual index loops where clearer.
- `.clone()` on non-`Copy` types is explicit — no hidden clones in closures and iterators.
- `Arc` only where sharing is real; never double-wrap a type that already shares (`TranscriptStore` is Arc-backed).

## Concurrency

- `rayon` for CPU-bound parallelism (the CLI's VCF annotation). The library has no async runtime; do not add one.
- `VarEffect`, `TranscriptStore` and `FastaReader` stay `Send + Sync` (compile-time assertion in `vareffect/src/lib.rs`).
- Prefer `RwLock` or lock-free structures over `Mutex` when reads dominate; channels for message passing. Parallel output keeps input order (V-2).

## Unsafe

No `unsafe` unless necessary. Every `unsafe` block states the safety invariant it relies on in a `// SAFETY:` comment, and the change is reviewed at T2.

## Imports and dependencies

- No wildcard imports, except preludes and `use super::*` in test modules.
- Order: standard library, external crates, local modules (let `rustfmt` keep it tidy).
- Dependencies declared with version constraints in `Cargo.toml`, workspace-level where shared. A dependency change can change annotation output or store compatibility (routes to the specialist).

## Testing

- Unit tests for every new function and type, in `#[cfg(test)]` modules, with the built-in `#[test]`; Arrange-Act-Assert.
- Mock external dependencies (network, file system); use small synthetic fixtures rather than the real genome where possible.
- Tests needing the built data (`data/vareffect/`) are `#[ignore]`-gated and locate it relative to the workspace root; name the opt-in command (see `vareffect-verification`).
- Expected annotation values come from recorded VEP output, never from vareffect's own output or memory.
- No commented-out tests.

## Tools

`rustfmt`, `clippy` with `-D warnings` (in CI and `make lint`, never `#![deny(warnings)]` in source), and `make all` before handing off a task with no errors or warnings.
