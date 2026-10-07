---
name: vareffect-sequence
description: Carry an authorized sequence of vareffect tasks through grounding, implementation and independent review with resumable per-task records. Use when the user authorizes several vareffect tasks or milestones in one go.
---

# Sequence tasks

Follow the user-authorized order and scope. Explain dependency changes; do not invent extra work or commit boundaries. For each task run [vareffect-ground](../vareffect-ground/SKILL.md), implement and verify, then [vareffect-review](../vareffect-review/SKILL.md). Proceed automatically through authorized tasks and report each checkpoint.

Stop dependent work on an open maintainer decision, unavailable required review, repeated repair failure or no progress; unaffected authorized work may continue. Never call an incomplete task done. Commits and releases need their own authorization — a sequence request grants neither.

## Task records

Keep one gitignored directory per task, `.vareffect-tasks/<task-id>/`, holding `record.md`, the review `packet.md`, `review.diff` and any concordance outputs. Records are required for sequences, for any task with more than one review round, and whenever a finding or a maintainer decision is left open; a single T0/T1 task may rely on its final report.

Record: scope and authorization; worktree; comparison and resolved commits; tier and reason; acceptance criteria; plan-review rounds; decisions with source evidence; reviewers (role, model, mode, agent id for `SendMessage` continuation, snapshot tree id); each finding's stable signature, status, confirmations and repair count (with the new evidence for attempt two); concordance runs (data hashes, commits, result paths); open decisions for the user; verification commands, results and limits; next action.

Keep a separate record and comparison per task so earlier uncommitted work stays visible without being mistaken for new work.

## Resume

Reconcile the record with Git: `git rev-parse HEAD`, `git status --short`, and `python3 scripts/vareffect_harness.py snapshot` against the recorded tree id. Changed bytes invalidate affected review. Carry every finding and repair count forward; never reset limits by renaming tasks or findings, and never silently adopt a different base. If state is missing or inconsistent, reconstruct from Git and prior evidence, mark unknown counts and decisions unresolved, and treat a finding whose repair count is unknown as exhausted.
