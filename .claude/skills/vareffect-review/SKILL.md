---
name: vareffect-review
description: Independently review a vareffect change — route it, freeze the reviewed bytes, dispatch fresh parallel reviewers, integrate permitted fixes within bounded repair loops, and hand unresolved decisions to the user. Use after implementing and verifying any T1+ vareffect change.
---

# Review changes

Reviewers are read-only; the main engineer owns fixes, checks and task state. Obligations and the report format are in `vareffect-findings`.

## 1. Route and tier

Re-run `python3 scripts/vareffect_harness.py route --base <comparison>` and confirm the tier from vareffect-ground (scope may have grown). Add `vareffect-bioinformatics` for semantic effects the paths miss. Reconcile any task record with current Git state.

## 2. Verify and freeze

Run the checks from `vareffect-verification` and record exact commands, results and limits. Then freeze:

```
python3 scripts/vareffect_harness.py snapshot --base <comparison> --diff-out .vareffect-tasks/<task>/review.diff
```

It prints `HEAD` and one tree id covering the whole worktree (tracked, staged, unstaged and untracked non-ignored files) without touching the real index or files, and writes the exact binary diff from the comparison start to that tree as the durable record of what was reviewed. It refuses files that look like secrets or sequencing data (`.env` variants, keys and credentials, VCF/BAM/CRAM/FASTQ). For a confirmed-synthetic fixture pass `--allow <path>` for that one path and say so in the packet — never a blanket override. Record the tree id; re-running `snapshot` later tells you whether the reviewed bytes changed, and `git diff --stat <old-tree> <new-tree>` shows which files did. Any change invalidates affected review.

## 3. Build one review packet

Write `.vareffect-tasks/<task>/packet.md` (gitignored) with: the original requirement and authorization, tier and reason, comparison and snapshot tree id, the `review.diff` path, changed-file list, verification commands and results, the concordance contrast where required, and known limits. Reviewer prompts then contain only the mode, the role focus for this change, and the packet path — never paste the diff into several prompts.

## 4. Dispatch

Spawn every required reviewer **in a single message** so they run in parallel: a **fresh** `vareffect-reviewer` (`model: opus` for T2/T3) and, when required, `vareffect-bioinformatics` (no `model` parameter). Never reuse a plan-review instance. A missing reviewer, partial review, or an unavailable required check is **incomplete review**, never a clean pass.

## 5. Integrate findings

Check each finding against the evidence yourself. Unverified claims become open questions. Re-request a malformed finding once, then report it as incomplete review. Protected findings (vareffect-findings) cannot be voted away, downgraded to avoid a maintainer decision, or dismissed by the implementer. A `route_to` raises the task to at least T2.

**Closing** a protected finding — repaired, refuted or discarded — needs a **fresh `vareffect-bioinformatics` instance** to confirm it against the current snapshot. Record each confirmation (agent id, tree id) in the task record. An unavailable specialist leaves the task incomplete.

## 6. Repair (bounded)

- At most **two attempts per finding** and **two repair rounds per task**. A second attempt needs new evidence, recorded before editing. Stop on no progress or identical findings in consecutive rounds.
- After each repair, re-run affected checks (including the concordance run when output paths changed), `route` and `snapshot`.
- Protected findings are re-checked by a fresh specialist instance, so the check is not anchored to the reviewer's own proposed fix. Other findings: **continue the reviewer that raised them** with `SendMessage` (mode `re-review`, repair diff and new tree id). Spawn fresh instead when the repair grew beyond what that reviewer saw or the instance is unavailable.

## 7. Unresolved findings

- introduced by this task and exhausted → the task is **incomplete**. Nothing is committed around an exhausted blocker without the user's explicit decision.
- pre-existing and outside scope → list under **"Pre-existing issues"** in the report, with enough detail to file an issue. Opening a GitHub issue is outward-facing: only at the user's direction.
- introduced `nit` → fix it in the first pass or drop it; never a repair round.
- introduced `minor`, deliberately not worth a round → state it in the report; the task stays complete.

A protected finding introduced by this task is repaired, or the task is incomplete; it is never dropped as minor.

**Maintainer decisions go to the user.** List each under **"Decisions for you"**: the question, realistic options, a recommendation, what is held, and the evidence. Every unresolved blocker, charter conflict, breaking change and accepted output movement appears there too. A blanket "do what you recommend" covers engineering choices, not accepting a VEP divergence or a breaking release.

## 8. Report

State what changed, checks run with results, the concordance contrast (or why none was needed), review coverage (who reviewed which snapshot, specialist confirmations), **Decisions for you**, **Pre-existing issues**, repair counts and every limit. Commits and releases need their own authorization; agent review and passing checks grant neither.
