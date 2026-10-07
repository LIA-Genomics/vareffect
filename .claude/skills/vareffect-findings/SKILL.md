---
name: vareffect-findings
description: Reviewer obligations, the vareffect finding report format, and how to recognize maintainer decisions and protected findings. Preloaded into every vareffect reviewer; use when writing or integrating review findings.
---

# Findings and authority

## Reviewer obligations

- **Read-only.** Return findings inline. Do not edit files or task state, run builds, commit, or send anything externally. The main engineer owns all integration and runs all checks; request a check rather than assuming its result.
- **Charter first.** Read `.claude/INVARIANTS.md`. Code establishes current behavior; the charter establishes required behavior. Report an implementation defect rather than rewriting required behavior to match code.
- **Independent evidence.** Review the actual source and surrounding effects, not the author's summary. Verify the author's claims too. Use the review packet the prompt names; open files yourself.
- **Indirect effects route.** An annotation-output, coordinate, data-provenance or public-contract effect counts even when a change reaches it indirectly (a shared helper, a builder, a serializer, a CLI flag default). Name `vareffect-bioinformatics` in `route_to`.
- **Never a false pass.** An unavailable tool or reference, unreadable file, scope you did not cover, or any write, edit, shell, delegation, skill-invocation or MCP tool in your own inventory makes the review `INCOMPLETE`, not `PASS`. Protected concerns cannot be dismissed by majority vote.

## Report format

```
scope: <what was reviewed: comparison, snapshot tree id, plan section>
mode: plan-review | change-review | re-review
verdict: PASS | FINDINGS | INCOMPLETE
examined: <files/symbols/evidence actually opened>
limits: <what was not checked and why; "none" only if true>
checks requested: <builds/tests/concordance runs the main engineer should run, or "none">

[<ROLE>-1] <severity> · <category> · confidence <high|medium|low>
at: <path:line or plan section> `<stable symbol>` — when <trigger>
defect: <what is wrong> → <consequence>
evidence: <opened path:line or fetched URL that shows it>
fix: <proposed correction>
flags: <only non-defaults: maintainer_decision=<the decision>, direction=<overcall|undercall|either|none>, route_to=<role>, breaking>

open questions: <unverified concerns, each with what would settle it>
```

- **Severity:** `blocker` (wrong annotation, data corruption, panic on valid input, undeclared breaking change; must not ship), `major` (real defect with material consequence), `minor` (real defect, limited consequence), `nit` (style or clarity).
- **Filter at the source.** Report a `nit` only on lines the change touches. Report a pre-existing issue only when it is `major` or worse, conflicts with the charter, or is protected by meaning. An issue is pre-existing only if it exists identically at the comparison base; one whose reachability, inputs or effect the change alters is introduced, whatever lines it sits on.
- **Category** names the defect precisely (coordinates, consequence, hgvs, normalization, transcript-model, provenance, concordance, determinism, error-handling, public-api, data-format, performance, unsafe, concurrency, input-validation, verification, citation, harness, …). No catch-all verdicts.
- **Direction** is set only by the bioinformatics lens, for the downstream loss-of-function use (`none` = established to have no directional effect); omit it when not established.
- **`breaking`** marks a change to a public contract's shape or meaning (V-6).
- **Stable signature** (for repair counting) = category + stable symbol/path + defect/trigger, excluding line numbers and the reviewer-local ID.
- A verified defect needs opened evidence. An unverified concern goes under open questions, never among findings. Verify that cited identifiers exist **and support the claim** (citation-verification). A fabricated API, VEP behavior or specification clause cannot support implementation: it becomes a blocking open question.
- A `PASS` report is the header alone, with `examined` and `limits` filled in honestly.

## Maintainer decisions

A finding needs the maintainer only when closing it requires one of these choices and no charter clause, documented divergence or recorded decision fixes exactly one answer. Everything else — internal design, naming, module layout, an error message — is engineering: implement the best option and state it in the report.

- accepting an intentional divergence from VEP, or retiring a documented one (V-4);
- which transcript source, release or genome build is supported or the default; what the builder admits or excludes (V-5);
- adjudicating a case where VEP, the HGVS recommendations and the Sequence Ontology disagree and the verified evidence leaves the answer genuinely open;
- accepting a breaking change to a public contract, its semver level, or a release (V-6);
- on-disk format compatibility policy (whether old stores keep loading);
- amending the charter.

Flag `maintainer_decision=<the decision>`; a finding that cannot name one of these choices is not one. Repairing an implementation so it meets the charter or a documented divergence, when that source fixes exactly one corrected behavior, is engineering — but if it changes annotation output it is still measured (concordance-measurement) and its effect is reported.

## Protected findings

A finding is **protected** when any of these holds:

- it carries `maintainer_decision`, `breaking`, or a direction of `overcall`, `undercall` or `either`;
- its fix alters annotation output for any variant, the meaning of a stored transcript model or genome, or a public contract, or it conflicts with the charter;
- `vareffect-bioinformatics` raised it or it carries `route_to`.

Judge meaning independently of the emitter's flags and of how the change author labelled it. The main engineer cannot self-clear a protected finding: closing one — repaired, refuted or discarded — needs a fresh `vareffect-bioinformatics` instance to confirm it against the current snapshot (vareffect-review §5). Never mint a citation.
