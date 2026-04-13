# VEP divergences

This document catalogues every known point where `vareffect`'s output
differs from Ensembl VEP, every VEP feature that is not yet implemented,
and every feature that is intentionally out of scope for the core crate.

## Overview

- **Target VEP version:** releases 115 and 116.
- **Genome build:** GRCh38 (GRCh37 and CHM13 are not validated; the
  transcript and genome binaries would need to be rebuilt from
  build-specific sources).
- **Transcript source:** MANE Select, MANE Plus Clinical, and RefSeq
  Select. Variants on accessions outside the loaded store are reported as
  `intergenic_variant` even if they would overlap a non-loaded transcript.
- **Validation date:** last full concordance run 2026-04-11 — 6 / 6 test
  files pass (see [Validation methodology](#validation-methodology) for the
  per-file variant counts).

## Intentional divergences

These are cases where `vareffect` deliberately produces different output
from VEP. Each is documented as an intentional choice rather than a bug.

### Splice region sub-terms not emitted

VEP optionally emits three fine-grained splice region sub-terms via plugin
extensions or newer SO vocabularies:

- `splice_polypyrimidine_tract_variant` (acceptor intronic positions −3 to
  −17)
- `splice_donor_region_variant` (donor intronic positions +3, +4, +6)
- `splice_donor_5th_base_variant` (donor intronic position +5)

`vareffect` classifies all of these as the generic `splice_region_variant`
plus the surrounding context term (`intron_variant` for deep intronic
positions). These sub-terms are not part of VEP's default core SO output —
they are produced by plugin extensions or newer releases that have not
been adopted widely in clinical pipelines.

If you need sub-term granularity, post-process the `ConsequenceResult` and
the `IndelLocation` splice offsets from `locate_indel` to reclassify.

### Non-coding transcript classification

For an intronic variant on a non-coding transcript (`NR_*` accession),
VEP may emit `non_coding_transcript_variant`. `vareffect` emits the more
specific `intron_variant` and leaves the non-coding context implicit in
the transcript accession prefix.

This is a lossless divergence — the caller can recover the non-coding
distinction from `TranscriptModel::biotype`, but if your downstream
pipeline keys on the SO term string it will need a one-line adapter.

### NMD prediction via 50-nucleotide rule only

`vareffect` applies the 50-nucleotide NMD rule to `stop_gained` and
`frameshift_variant` consequences and exposes the result as
`ConsequenceResult::predicts_nmd`. VEP emits a separate
`NMD_transcript_variant` SO term when the *transcript itself* is biotyped
as NMD. `vareffect` does not emit this biotype-level term; if your store
contains explicitly NMD-flagged transcripts, filter them in the caller.

The rule checks the distance from the premature termination codon to the
last exon-exon junction in the CDS. Variants in single-exon genes
(`PRNP`, the mitochondrial genome) are correctly reported as
`predicts_nmd = false` because there is no junction.

### HGVS 3' normalization is always on

VEP has both `--shift_hgvs` (3' shift, default since release 109) and an
undocumented legacy mode. `vareffect` always applies 3' normalization on
the coding strand, matching VEP's current default. There is no toggle to
disable it.

## Not yet implemented

These are features `vareffect` intends to cover but has not yet. Pull
requests welcome.

- **`mature_miRNA_variant`** — variant overlaps a miRNA hairpin. Requires
  a miRNA locus track that is not currently loaded into `TranscriptStore`.
- **`TFBS_ablation`, `TF_binding_site_variant`, `regulatory_region_variant`** —
  the entire regulatory / motif annotation layer. Requires a
  RegulatoryFeatureStore with transcription factor binding sites,
  enhancers, and promoters. See the [Out of scope](#out-of-scope-by-design)
  section below for why this is not on the core roadmap.
- **`feature_elongation`, `feature_truncation`, `transcript_amplification`** —
  structural-variant consequence terms. `transcript_ablation` is
  implemented for variants that delete an entire transcript; elongation
  and truncation require a different input surface (segment-level SV
  calls, not VCF rows).
- **`--shift_3prime` equivalent** — a toggle to switch off 3' normalization
  for callers who need VEP's pre-release-109 behaviour.
- **Alternate genome builds** — GRCh37, CHM13, and non-human assemblies.
  The crate is build-agnostic at runtime; the work is all upstream, in the
  transcript-model and genome-binary build pipelines.
- **Multi-allele VCF splitting** — VCF `ALT` columns can carry multiple
  comma-separated alternate alleles. `vareffect` expects one ref / one alt
  per `annotate` call; callers must split beforehand. A convenience
  wrapper that fans out multi-allele rows would be a welcome contribution.
- **Canonical transcript *selection*** — every `ConsequenceResult` carries
  MANE Select / MANE Plus Clinical / RefSeq Select metadata, but the
  caller is responsible for filtering to the chosen tier. A small helper
  method on `VarEffect` to return only the canonical result would be
  backward-compatible.

## Out of scope (by design)

These are features `vareffect` will not implement as part of the core
crate, because they belong to a different layer of the pipeline.

- **Plugin system.** VEP has a Perl plugin architecture with ~60 plugins
  in the wild (CADD, REVEL, SpliceAI, LoFtool, …). `vareffect` is a
  library, not an extension platform — downstream code can wrap results
  and apply plugins as a separate step.
- **Co-located variant lookup.** No gnomAD, dbSNP, ClinVar integration,
  no allele-frequency annotation, no known-variant matching. These are
  enrichment steps that happen after consequence prediction, and coupling
  them into the predictor would force every consumer to carry network
  dependencies.
- **Custom annotation files.** No BED / GFF / VCF overlay interface. Use a
  dedicated tool (bcftools, bedtools) and merge on the output side.
- **VCF I/O.** `vareffect` takes variants as four arguments (chrom, pos,
  ref, alt). There is no built-in VCF parser — use `rust-bio`, `noodles`,
  or `bcftools pipe`.

## SO term coverage matrix

VEP produces roughly 35 distinct SO consequence terms. This table shows
which are covered, which are not, and the reasoning.

| SO term                                   | Status | Note |
|-------------------------------------------|--------|------|
| `missense_variant`                        | yes    | Covered by the SNV annotator. |
| `synonymous_variant`                      | yes    | |
| `stop_gained`                             | yes    | Drives NMD prediction. |
| `stop_lost`                               | yes    | Extension distance computed via 3'UTR stop-scan. |
| `start_lost`                              | yes    | Covers SNV and partial-deletion start codon loss. |
| `start_retained_variant`                  | yes    | Synonymous change at the start codon. |
| `stop_retained_variant`                   | yes    | |
| `frameshift_variant`                      | yes    | Drives NMD prediction. |
| `inframe_insertion`                       | yes    | |
| `inframe_deletion`                        | yes    | |
| `protein_altering_variant`                | yes    | Catch-all for length-preserving complex changes. |
| `coding_sequence_variant`                 | yes    | Ambiguous-codon fallback. |
| `incomplete_terminal_codon_variant`       | yes    | CDS length not divisible by 3. |
| `splice_donor_variant`                    | yes    | Intronic ±1, ±2 on the donor side. |
| `splice_acceptor_variant`                 | yes    | Intronic ±1, ±2 on the acceptor side. |
| `splice_region_variant`                   | yes    | Intronic ±3..±8 and exonic 1..3. |
| `splice_polypyrimidine_tract_variant`     | no     | Intentional divergence — emitted as `splice_region_variant`. |
| `splice_donor_region_variant`             | no     | Intentional divergence — plugin-extension term. |
| `splice_donor_5th_base_variant`           | no     | Intentional divergence — plugin-extension term. |
| `5_prime_UTR_variant`                     | yes    | |
| `3_prime_UTR_variant`                     | yes    | |
| `intron_variant`                          | yes    | |
| `non_coding_transcript_exon_variant`      | yes    | |
| `non_coding_transcript_variant`           | no     | Intentional divergence — non-coding context is preserved via `TranscriptModel::biotype`. |
| `upstream_gene_variant`                   | yes    | |
| `downstream_gene_variant`                 | yes    | |
| `intergenic_variant`                      | yes    | No overlapping transcript in the loaded store. |
| `transcript_ablation`                     | yes    | Variant deletes an entire transcript. |
| `transcript_amplification`                | no     | Not yet implemented — SV-shaped term. |
| `feature_elongation`                      | no     | Not yet implemented — SV-shaped term. |
| `feature_truncation`                      | no     | Not yet implemented — SV-shaped term. |
| `NMD_transcript_variant`                  | no     | Intentional divergence — biotype-based, exposed via `predicts_nmd` instead. |
| `mature_miRNA_variant`                    | no     | Not yet implemented — requires a miRNA locus track. |
| `TF_binding_site_variant`                 | no     | Out of scope — regulatory layer. |
| `TFBS_ablation`                           | no     | Out of scope — regulatory layer. |
| `regulatory_region_variant`               | no     | Out of scope — regulatory layer. |

## Validation methodology

The `tests/vep_concordance_*.rs` suite holds hand-curated variants whose
expected outputs were recorded from the Ensembl VEP REST API on real
clinical transcripts (TP53, BRCA1, BRCA2, APC, BRAF, EGFR, CFTR, ERBB2,
PRNP, and others). Each test file opens the transcript store and FASTA
reader once, iterates its `VARIANTS` fixture, and asserts per-variant
equality on `hgvs_c`, `hgvs_p`, the consequence subset, and (where
relevant) the `predicts_nmd` flag.

| Test file                              | Variants | Focus |
|----------------------------------------|---------:|-------|
| `vep_concordance_snv.rs`               |       20 | SNV consequences across missense, synonymous, stop gain / loss, start loss / retained, splice donor / acceptor, splice region |
| `vep_concordance_indel.rs`             |       28 | Frameshift, inframe insertions / deletions, boundary-spanning indels, splice-overlap indels |
| `vep_concordance_hgvs.rs`              |       20 | HGVS c. forward formatting (substitution, del, ins, dup, delins, intronic offsets, UTR offsets) |
| `vep_concordance_hgvs_p.rs`            |       30 | HGVS p. forward formatting for every consequence type |
| `vep_concordance_hgvs_reverse.rs`      |       20 | HGVS c. → genomic coordinate round-trip for every position type |
| `vep_concordance_normalization.rs`     |       18 | HGVS 3' normalization, intergenic classification, NMD 50-nt rule |
| **Total**                              |  **136** | |

As of 2026-04-11 all six test files pass (6 / 6 test functions, 136 / 136
variants). The suite is `#[ignore]`-gated because it requires the
transcript store and genome binary on disk; run with:

```bash
FASTA_PATH=/absolute/path/to/GRCh38.bin \
    cargo test -p vareffect --release -- --ignored vep_concordance
```

When adding a new variant to a fixture, record its expected output by
querying the Ensembl VEP REST API with `refseq=1&hgvs=1` (or `&numbers=1`
for exon number assertions) and paste the response into the fixture
comment so future reviewers can cross-check against the recorded ground
truth.

## Known edge cases

- **NMD with final-exon 3'UTR tails.** The 50-nucleotide rule is applied
  to the distance from the PTC to the last exon-exon CDS junction. For
  transcripts where the final exon's CDS portion is under 50 nt but the
  3'UTR extends well beyond, the rule is conservative — a variant just
  inside the last coding exon may be flagged as NMD-predicted when VEP
  would let it escape. This matches the ClinGen-recommended
  interpretation of the rule but may produce false positives in rare
  transcript architectures.
- **Patch-contig lookups.** The reference genome binary can be built with
  or without NCBI patch contigs. When you need patch lookups, use
  `VarEffect::open_with_patch_aliases` and supply a `refseq,ucsc` alias
  CSV — otherwise variants on patch contigs return `ChromNotFound`.
- **IUPAC ambiguity codes.** The NCBI GRCh38.p14 assembly contains
  ambiguity bases (`M`, `R`, `Y`, etc.) in a few patch regions. The
  genome reader preserves them, but codon translation treats any
  non-ACGTN base as an unknown amino acid (`Xaa`). Consensus-based
  patch contigs should not contain these bases.
- **Mitochondrial variants.** `chrM` uses NCBI translation table 2. The
  `codon` module switches tables automatically based on the chromosome
  name; variants on `chrMT` or `MT` aliases rely on the reader's
  chromosome-name alias table to resolve to `chrM` first.
