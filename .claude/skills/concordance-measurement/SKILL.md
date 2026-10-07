---
name: concordance-measurement
description: Measure what a change does to annotation output before shipping — pinned data, baseline worktree, per-variant before/after joins, movement buckets checked against VEP, and what the numbers do not establish. Use for any change that can move consequences, HGVS, NMD/PTC or the transcript set.
---

# Measuring an annotation change

A change to annotation logic carries two questions:

1. **Is it defensible?** Do the HGVS recommendations, the Sequence Ontology definition, VEP's documented or observed behavior, or the source format's specification support it? Answered by reading.
2. **What does it do?** Which variants change output, from what to what, and are the new outputs right? **Answered only by running it.**

A reviewer who confirms that a fix matches its spec sentence has established that it is defensible. A fix for one case can still move hundreds of others; only a run over a corpus shows that. **If vareffect can be run over the corpus, an output change is not reviewed until it has been.** A change meant to leave output unchanged (refactor, performance) is measured too: the expected movement is zero, and anything else is a finding.

## Corpora available

- **Curated suites** — `vareffect/tests/vep_concordance_*.rs`, about 140 variants with recorded VEP REST output. Exact per-variant equality; they exercise named edge cases, not the distribution.
- **Large ClinVar corpus** — `vareffect/tests/vep_large_concordance.rs` compares about 50k ClinVar variants against a ground-truth TSV (vep_ground_truth.tsv in the tests' data directory) and writes a mismatch log. It holds one RefSeq transcript per variant and passes on a concordance *rate* threshold after folding VEP's granular splice terms, so its exit status is not per-variant concordance: movement comes from its mismatch log or from the join below. On `main` the TSV is absent — it was committed only on the unmerged `grch37-support` branch — and the test then prints `SKIP` and reports `ok`. Until it lands, the large-corpus contrast is unavailable; the curated suites alone are a limit, not a pass.
- **Ad-hoc corpus** — for a targeted change, a deterministic set of variants that exercises it (for example every variant in the affected transcripts, or a ClinVar subset filtered by consequence), with VEP output fetched from the REST API or a pinned local VEP. Save the corpus and the VEP output under `.vareffect-tasks/<task>/`.

**Record the reference.** Every VEP fetch, curated or ad hoc, records beside its output the VEP release (`GET https://rest.ensembl.org/info/software` — the public REST server moves to each new release), the assembly, the full request (endpoint and parameters, documented parameters only) and the date. A release outside the 115/116 target, or one that differs between baseline and change, is a limit: a VEP-side change would otherwise read as a vareffect regression or fix.

## The procedure

1. **Fix one corpus** and reuse the same file for both runs.
2. **Run both builds end-to-end**: the baseline from a worktree pinned at the comparison base (`git worktree add`), not by reverting in place; the change from the snapshot.
3. **Pin the inputs.** `GRCh38.bin`, `transcript_models.bin` and the ground truth must be byte-identical across runs: compare SHA-256, not version strings. Hashing pins the bytes, not the reference version: the recorded VEP release (above) pins that. A builder change is the exception: then the store is the variable — build both from the same pinned source file and record both hashes.
4. **Join per variant** on the input exactly as submitted to both tools — VCF `CHROM`, 1-based `POS`, `REF`, `ALT` — plus the versioned transcript accession. Never join on VEP's echoed trimmed `start`/`allele_string`: indels would fall into the transcript gained/lost bucket. Report:
   - agreement with VEP before and after;
   - the **movement**, bucketed by field and from→to (consequence set, `hgvs_c`, `hgvs_p`, `predicts_nmd`/PTC), each bucket annotated with whether the new value agrees with VEP;
   - new panics, errors and indeterminate results, and transcripts gained or lost.
5. **State the trade in one line**: what was fixed, what moved that was not meant to, and how many of those now disagree with VEP.
6. **Name the variants.** Give a handful of concrete examples per bucket; a count is not reviewable, a variant is.

## Reading the result honestly

- **VEP is the reference, not the truth.** Where VEP contradicts the HGVS recommendations or SO, agreeing with it may be wrong. That is a candidate intentional divergence — a maintainer decision, recorded in `VEP_DIVERGENCES.md` (V-4), never a silent choice.
- **Transcript and version mismatch are not concordance.** Variants whose VEP transcript is missing from the store, or present at another version, are a separate bucket; do not count them as agreements or as defects.
- **ClinVar is not a random sample.** It over-represents clinically reported genes and variant classes; say which classes the corpus does not exercise.
- **The contrast is two builds on one corpus**, never a claim about accuracy in general.
- **Check that the changed branch has a test that exercises it.** A rule changed with no test through it means a green suite said nothing about the change.
- **Downstream direction.** For movements into or out of loss-of-function terms (`stop_gained`, `frameshift_variant`, splice donor/acceptor, `start_lost`) or NMD/PTC, state the direction for a PVS1 consumer: overcall, undercall or either.

## What the measurement does and does not authorize

It records a cost and a benefit; it does not decide. Accepting output movement that diverges from VEP, or that changes a public contract's meaning, is a maintainer decision made with the numbers in hand. A good number does not rescue a change with no basis in the specification; overfitting to one corpus is its own failure.
