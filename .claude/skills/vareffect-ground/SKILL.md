---
name: vareffect-ground
description: Scope a vareffect task, assign its risk tier, ground it in current code, and get independent plan review where the tier requires it. Use at the start of any nontrivial vareffect change, before editing.
---

# Ground a task

Read `.claude/INVARIANTS.md`, `CLAUDE.md` and the relevant code. Code establishes current behavior; the charter establishes required behavior, so a disagreement may be an implementation defect.

## 1. Scope and route

Fix an explicit Git comparison and run `python3 scripts/vareffect_harness.py route --base <ref>` (or `<left>...<right>` for a merge-base comparison). Its reviewer list is a path-based **minimum**; add `vareffect-bioinformatics` for semantic effects the paths do not show — a shared helper the annotator calls, a serializer default, a CLI flag that changes which transcripts load.

## 2. Assign the tier (record it with a one-line reason; when unsure, go one tier up)

| Tier | When | Plan review | Change review |
| --- | --- | --- | --- |
| **T0 trivial** | Spelling or formatting in documentation or comments only, changing no stated behavior | none | self-review, recorded |
| **T1 routine** | The route lists no specialist and nothing can change annotation output, a public contract or a data format: CLI ergonomics, build tooling that cannot change the compiled annotator, internal refactor of non-annotation code, tests of non-annotation code | none | one fresh `vareffect-reviewer`, default model |
| **T2 protected** | The route names `vareffect-bioinformatics` but output is meant to stay identical (refactor, performance work on the hot path); or a public API or store-format change; `unsafe`, mmap or concurrency; release or CI workflow; dependency upgrades; harness changes; cross-crate work | fresh `vareffect-reviewer` with `model: opus` plus `vareffect-bioinformatics`, in parallel | same, plus proof that output did not move (concordance-measurement, "no movement expected") |
| **T3 annotation semantics** | Intends to change what vareffect emits for some variant — consequence terms, HGVS strings, NMD/PTC, exon/intron numbering, normalization, transcript selection — or what data it builds from: transcript source, release, admission rules, genome build | as T2 | as T2, plus the measured before/after contrast (concordance-measurement), `VEP_DIVERGENCES.md` and `CHANGELOG.md` updated, and every maintainer decision listed |

A task never moves down a tier during review; any `route_to` a reviewer raises moves it to at least T2. A performance change that unexpectedly moves output is T3.

Engineering reviewer: `vareffect-reviewer`. Specialist: `vareffect-bioinformatics`, pinned to Opus — never pass a `model` for it (a hook blocks downgrades).

## 3. Ground

Inspect the implementation yourself. Load `vareffect-verification` for checks. For a concrete unknown — a format convention, VEP's behavior on a case, an HGVS rule, a crate API — spawn `vareffect-researcher`; for pure code location the built-in `Explore` agent is enough. Do not commission a general research panel. Never plan from remembered format or VEP behavior: verify it.

Audits, concordance runs and investigations **document** findings; fixing them is separate authorized work.

## 4. Plan and review (T2, T3)

Draft the approach: verified existing symbols vs proposed new ones, error and availability handling (V-3), public-contract and store-format compatibility (V-6), the concordance measurement and which variants are expected to move, tests and observable acceptance criteria. A plan the user supplied still gets this review.

Spawn the plan reviewers **in one message** so they run in parallel. Give each the task, requirements, plan and evidence by path (see vareffect-review §3 for the packet); each verifies consequential claims itself. The main engineer cannot substitute self-review.

Integrate findings in at most **two plan revision rounds**; stop on no progress. For a revision, continue the same plan reviewers with `SendMessage` rather than spawning new ones. Unavailable reviewers, unverified dependencies and maintainer decisions stay open — never a clean pass. Once required review is complete and nothing blocking remains, start implementing without asking for routine approval.

Then implement, verify (`vareffect-verification`) and run `vareffect-review` with **fresh** reviewers — plan reviewers never review the change.
