# VEP divergences

This document catalogues every known point where `vareffect`'s output
differs from Ensembl VEP, every VEP feature that is not yet implemented,
and every feature that is intentionally out of scope for the core crate.

## Overview

- **Target VEP version:** releases 115 and 116.
- **Genome builds:** GRCh38 and GRCh37, selected per call via the
  `Assembly` enum. CHM13 is not yet supported. UCSC `hg19` / `hg38`
  aliases are explicitly rejected at parse time — UCSC `hg19` chrM
  (`NC_001807`) differs from GRCh37 chrMT (`NC_012920.1`, the rCRS) by
  ~10 bases and would silently mis-annotate every chrM variant.
- **Transcript source:** **GRCh38** uses MANE Select + MANE Plus
  Clinical from the NCBI MANE v1.5 release. **GRCh37** uses NCBI RefSeq
  Select from `GCF_000001405.25_GRCh37.p13`; MANE is GRCh38-only and
  has no GRCh37 equivalent. Variants on accessions outside the loaded
  store are reported as `intergenic_variant` even if they would overlap
  a non-loaded transcript.
- **Transcript-vs-genome divergence (GRCh37 only):** ~5 % of NCBI
  RefSeq Select transcripts on GRCh37 carry sequence that differs from
  the reference assembly at one or more positions. NCBI flags these in
  the GFF3 `Note=` attribute on the affected mRNA / CDS rows; vareffect
  detects the flag at build time, persists it on
  [`TranscriptModel::genome_transcript_divergent`], and surfaces a
  structured [`Warning::DivergentTranscript`] on every `annotate(...)`
  call whose chosen transcript carries it. HGVS positions emitted
  against a divergent transcript may not map back to the same genomic
  position they would against the reference; clinical callers should
  consider falling back to a non-divergent transcript before reporting.
  GRCh38 MANE transcripts are curated to exclude divergence by
  construction, so this warning is silent there.
- **Translational exceptions (selenocysteine / pyrrolysine
  readthrough):** parsed from the GFF3 `transl_except=` attribute and
  recorded on [`TranscriptModel::translational_exception`]. Distinct
  from divergence: a `transl_except` codon is a known biological
  mechanism, not a sequence disagreement, and HGVS positions remain
  reliable — vareffect therefore does NOT raise a clinical warning.
- **Cross-validation:** GRCh38 builds are second-source-checked against
  the MANE summary TSV at build time, with mismatches failing the
  build. GRCh37 builds are checked against UCSC's hg19 `ncbiRefSeq.txt`
  + `ncbiRefSeqSelect.txt` pair. UCSC re-derives NCBI's annotation
  release, so the GRCh37 cross-check catches GFF3 attribute-parser
  regressions and coordinate-conversion drift but does not catch
  divergence between two independent biological databases.

  **chrM is excluded from the GRCh37 UCSC cross-check.** UCSC `hg19`
  chrM is `NC_001807` (the original 1981 Anderson reference), while
  NCBI GRCh37 chrMT is `NC_012920.1` (the rCRS). The two references
  differ by ~10 bp plus indels, so naive coordinate comparison would
  systematically false-positive on every chrM transcript. The UCSC
  parser skips chrM rows at parse time with a build-log warning; the
  ~37 affected mitochondrial transcripts (MT-RNR1, MT-CO1, …) are
  validated separately via downstream ClinVar concordance.
- **Validation date:** last full concordance run 2026-04-30.
  - **GRCh38:** 6 / 6 spot-check files pass (136 / 136 hand-curated
    variants — see [Validation methodology](#validation-methodology) for
    the per-file breakdown). Large-scale concordance against VEP REST on
    ~50,000 ClinVar variants (`vep_large_concordance_grch38.rs`) is in
    progress with the threshold currently set at ≥95 % consequence
    concordance; tightening to ≥99 % once divergences are triaged.
  - **GRCh37 ClinVar self-concordance** (vareffect-GRCh37 vs
    vareffect-GRCh38 on 5,000 dual-coordinate ClinVar pairs,
    `grch37_clinvar_concordance.rs`): **100.00 %** consequence,
    **99.98 %** HGVS.c, **100.00 %** HGVS.p across 4,698
    shared-transcript comparisons, with **902 / 19,262 (4.68 %)**
    divergent transcripts excluded by construction (within tolerance
    of NCBI's published ~5 % figure for RefSeq Select on GRCh37).
  - **GRCh37 VEP REST concordance** on 10,000 stratified ClinVar
    variants captured from `https://grch37.rest.ensembl.org`
    (`vep_large_concordance_grch37.rs`): **99.33 %** consequence
    concordance across 5,220 in-store comparisons (89.64 % strict-equal
    + 8.31 % normalised splice/NMD folding + 0.67 % real mismatch).
    HGVS.c **96.99 %**, HGVS.p **96.53 %**, IMPACT **99.67 %**, protein
    start **97.68 %**. Zero panics, zero annotation errors.
  - **GRCh37 spot-check tier**: `vep_concordance_grch37_snv.rs` ships
    with **3 hand-curated variants** validated against VEP REST. Five
    additional spot-check files (`indel`, `hgvs`, `hgvs_p`,
    `hgvs_reverse`, `normalization`) — mirroring the GRCh38 layout —
    will land as separate files once their fixtures are curated. Use
    the same `https://grch37.rest.ensembl.org/vep/human/region/...`
    REST endpoint as the GRCh38 spot-checks, with `assembly=GRCh37`.

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
- **CHM13 and non-human assemblies.** The crate accepts an `Assembly`
  selector; adding new variants requires a hardcoded NC_* accession
  table per assembly plus a transcript-source pipeline. GRCh38 and
  GRCh37 are both supported.
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

### GRCh38 hand-curated spot-checks

| Test file                                     | Variants | Focus |
|-----------------------------------------------|---------:|-------|
| `vep_concordance_grch38_snv.rs`               |       20 | SNV consequences across missense, synonymous, stop gain / loss, start loss / retained, splice donor / acceptor, splice region |
| `vep_concordance_grch38_indel.rs`             |       28 | Frameshift, inframe insertions / deletions, boundary-spanning indels, splice-overlap indels |
| `vep_concordance_grch38_hgvs.rs`              |       20 | HGVS c. forward formatting (substitution, del, ins, dup, delins, intronic offsets, UTR offsets) |
| `vep_concordance_grch38_hgvs_p.rs`            |       30 | HGVS p. forward formatting for every consequence type |
| `vep_concordance_grch38_hgvs_reverse.rs`      |       20 | HGVS c. → genomic coordinate round-trip for every position type |
| `vep_concordance_grch38_normalization.rs`     |       18 | HGVS 3' normalization, intergenic classification, NMD 50-nt rule |
| **GRCh38 spot-check total**                   |  **136** | |

### GRCh37 hand-curated spot-checks

| Test file                                     | Variants | Focus |
|-----------------------------------------------|---------:|-------|
| `vep_concordance_grch37_snv.rs`               |        3 | SNV consequences (in progress; target ~50) |
| **GRCh37 spot-check total**                   |    **3** | additional files (`indel`, `hgvs`, `hgvs_p`, `hgvs_reverse`, `normalization`) will mirror the GRCh38 layout once curated |

### Statistical / large-scale concordance

| Test file                                | Rows | Focus |
|------------------------------------------|-----:|-------|
| `vep_large_concordance_grch38.rs`        | ~50,000 | Statistical consequence concordance vs VEP REST on stratified ClinVar variants. Threshold ≥95 %. |
| `vep_large_concordance_grch37.rs`        |  9,986 | Same shape on GRCh37; 5,220 in-store comparisons → 99.33 % consequence concordance. |
| `grch37_clinvar_concordance.rs`          |  5,000 | GRCh37 self-concordance vs vareffect-GRCh38 on dual-coord ClinVar pairs. Threshold ≥99 % per metric. |

The suite is `#[ignore]`-gated because it requires the transcript stores
and genome binaries on disk; run with:

```bash
GRCH38_FASTA=/abs/path/to/GRCh38.bin \
GRCH37_FASTA=/abs/path/to/GRCh37.bin \
GRCH37_TRANSCRIPTS=/abs/path/to/transcript_models_grch37.bin \
  cargo test -p vareffect --release -- --ignored
```

or via `make test-ignored` for the env-var-laden form. When adding a new
variant to a hand-curated fixture, record its expected output from the
relevant Ensembl VEP REST endpoint (`rest.ensembl.org` for GRCh38,
`grch37.rest.ensembl.org` for GRCh37) with
`refseq=1&hgvs=1&shift_hgvs=1&numbers=1` and paste the response into
the fixture comment so future reviewers can cross-check against the
recorded ground truth.

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
  or without NCBI patch contigs. When you need patch lookups, build the
  `VarEffect` via `VarEffect::builder().with_grch38_and_patch_aliases(...)`
  (or the GRCh37 variant) and supply the `refseq,ucsc` alias CSV that
  `vareffect setup` writes as `patch_chrom_aliases_grch{37,38}.csv` —
  otherwise variants on patch contigs return `ChromNotFound`.
- **IUPAC ambiguity codes.** The NCBI GRCh38.p14 assembly contains
  ambiguity bases (`M`, `R`, `Y`, etc.) in a few patch regions. The
  genome reader preserves them, but codon translation treats any
  non-ACGTN base as an unknown amino acid (`Xaa`). Consensus-based
  patch contigs should not contain these bases.
- **Mitochondrial variants.** `chrM` uses NCBI translation table 2. The
  `codon` module switches tables automatically based on the chromosome
  name; variants on `chrMT` or `MT` aliases rely on the reader's
  chromosome-name alias table to resolve to `chrM` first.
