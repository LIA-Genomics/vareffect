---
name: vareffect-verification
description: Choose and run the right vareffect checks — fast targeted gates while iterating, CI-equivalent gates at handoff, and the data-gated VEP concordance suites for anything that can move annotation output. Use before review and after every repair.
---

# Verification

The main engineer runs all builds and checks; reviewers inspect evidence and request targeted checks. An `#[ignore]`d data-gated test that was not run is not a pass: report it as a limit.

## While iterating (fast)

- `cargo fmt --all`, then `cargo clippy -p <crate> --all-targets -- -D warnings` and `cargo test -p <crate>` for each changed crate. A library change also needs `vareffect-cli` (it depends on `vareffect`).
- One test: `cargo test -p vareffect <test_name>`.

## Before review (handoff gates)

- `make all` — format, harness check, clippy (`--all-targets`, as CI), tests, rustdoc with warnings denied (as CI), build. It is the local equivalent of CI's `check` job.
- Not covered locally, report as a limit unless run: CI's MSRV job (`cargo +1.94.1 test --workspace`, if that toolchain is installed), the Windows and macOS matrix, and `cargo deny check advisories bans licenses` (if `cargo-deny` is installed). Run them when the change touches `Cargo.toml`, `Cargo.lock`, platform-specific code (paths, mmap, line endings) or uses a newly stabilized language or std feature.
- Harness-only changes: `make harness-check`; do not rebuild unrelated Rust. It does not prove the client loaded the agents; see `.claude/README.md` for the smoke procedure.

## Annotation output (T2/T3)

Any change the route sends to `vareffect-bioinformatics` runs the data-gated suites against the built data in `data/vareffect/` (`vareffect setup` builds it):

```bash
FASTA_PATH="$PWD/data/vareffect/GRCh38.bin" \
    cargo test -p vareffect --release --test '*' -- --ignored --nocapture
```

`--test '*'` runs the integration suites only; `--tests` would also run the library's `#[ignore]`d unit tests. The suites load `data/vareffect/transcript_models.bin` and the genome named by `FASTA_PATH`. Record both files' SHA-256 (`shasum -a 256`); a contrast is valid only when they are identical across the before and after runs, and a store is stale when `types.rs`, `chrom.rs` or the builders changed after `transcript_models.manifest.json`'s `built_at` (`git log -1 --format=%cI -- <paths>`); rebuild with `vareffect setup`, or report it as a limit. The genome has no manifest; its hash is the only provenance. Know the baseline: run the same command at the comparison base, since a failure that already exists there is pre-existing, not a pass and not introduced.

Any test that prints `SKIP:` returned early and reports `ok` — `vep_large_concordance` does so whenever its ground-truth TSV is absent, which is always on `main`. Treat it as not run (a limit), never a pass; read the output, not only the result line. That harness passes on a rate threshold, so even a real run shows movement only through its mismatch log.

That covers the curated `vep_concordance_*` suites (per-variant equality against recorded VEP REST output) and the location, FASTA and normalization tests. Then run the before/after contrast in `concordance-measurement`. A builder change also rebuilds the store from a pinned source and states the record-count and exclusion deltas (V-5).

Adding a variant to a curated suite: record its expected output from the Ensembl VEP REST API (`refseq=1&hgvs=1`, plus `numbers=1` for exon numbers) and keep the request and the recorded response in the fixture comment (`vareffect/VEP_DIVERGENCES.md`, Validation methodology). Never write an expected value from memory or from vareffect's own output.

## Evidence

Record commands, exit results, skipped checks, environment limits, the data files' SHA-256 and the snapshot tree id. Create `.vareffect-tasks/<task>/` for concordance outputs; `snapshot --diff-out` creates it too. After a fix, re-run the affected checks — not unrelated builds.
