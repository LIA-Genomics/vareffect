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

`predicts_nmd` measures the distance from the **variant's own CDS
position** to the last exon-exon junction, not from the termination codon.
The two coincide for `stop_gained`, but for `frameshift_variant` the
termination codon lies downstream of the variant and may sit in a later
exon — so `predicts_nmd` can report NMD for a frameshift whose actual
termination codon escapes it.

Measuring at the variant site follows `NMD.pm`, but the agreement is
partial. Matching: the variant site as the anchor, the last-junction
distance rule, and the intronless case. Diverging:

- `NMD.pm` rule 3 escapes a variant in the **first ~100 coding bases**
  (`$variant_coding_region <= 101`). This crate does not implement it, so a
  5'-proximal truncating variant is reported NMD-predicted where the plugin
  reports escape.
- `NMD.pm` walks **all transcript exons**; this crate uses CDS segments only.
- `NMD.pm`'s penultimate-exon window is the 3'-most **52 genomic bases**
  (`variant_exon_check`: `$diff_end = $coding_region_end - 51`, tested
  inclusively at both ends), against the 50 coding bases used here. The 50
  matches the ClinGen SVI wording ("within the 3'-most 50 nucleotides"), so
  `ptc` is deliberately ClinGen-aligned rather than `NMD.pm`-aligned.

`ConsequenceResult::ptc` carries that codon's own residue position and an
NMD verdict computed there. Consumers applying the ClinGen SVI PVS1
decision tree (Abou Tayoun et al. 2018, PMID 30192042,
doi:10.1002/humu.23626) must read `ptc`, not `predicts_nmd`.

**Which escape rules `ptc` applies.** Only the last-junction rule, because
that is the only rule the SVI decision tree uses; `AutoPVS1`, the published
reference implementation of that flowchart (Xiang et al. 2020, PMID
32442321, doi:10.1002/humu.24051), also implements no other. Two further
rules that are well established in the NMD literature — start-proximal
escape within 150 nt of the coding start site, and escape from a
termination codon inside an exon longer than 407 bp, both after Lindeboom
et al. (PMID 27618451, doi:10.1038/ng.3664; PMID 31659324,
doi:10.1038/s41588-019-0517-5) and both implemented by `aenmd` — are
deliberately **not** applied. Either would withdraw PVS1's full strength
from a large class of established loss-of-function alleles, which is a
guideline decision rather than an annotation one. Neither `predicts_nmd`
nor `ptc` models translation re-initiation; `ptc.protein_position` is
everything a consumer needs to layer these rules itself.

**Which base of the termination codon `ptc` measures from.** The 3'-most
one. The SVI rule is written as "the premature termination codon occurs in
the 3' most exon or within the 3'-most 50 nucleotides of the penultimate
exon"; a codon spans three bases, so "occurs within" is read as *overlaps*,
and the codon's last base is what decides it. Anchoring on the first base
instead would shrink the escape window to 48 nucleotides and predict NMD
for two codon positions the guideline places outside it — an overcall,
since predicting NMD is what earns PVS1 its full strength. Across the
GRCh38 ClinVar release this moves 202 of 4,666,135 annotation rows, every
one of them from NMD-predicted to NMD-escaping.

Variants in single-exon genes (`PRNP`, the mitochondrial genome) are
correctly reported as `predicts_nmd = false`, and `ptc` is NMD-negative,
because there is no junction.

### HGVS 3' normalization is always on

VEP has both `--shift_hgvs` (3' shift, default since release 109) and an
undocumented legacy mode. `vareffect` always applies 3' normalization on
the coding strand, matching VEP's current default. There is no toggle to
disable it.

### Structural variants are annotated at the region level

`annotate_interval` (DEL / DUP / INV), `annotate_breakend` (`<BND>`), and
`annotate_sv_insertion` (`<INS>`) emit the per-transcript SO term **set** VEP
assigns to a structural variant: the overlapped sub-regions
(`coding_sequence_variant`, `5_prime_UTR_variant`, `3_prime_UTR_variant`,
`intron_variant`, `non_coding_transcript_exon_variant`) plus a copy-number
headline where applicable — `transcript_ablation` / `transcript_amplification`
for a whole-transcript deletion / duplication, `feature_truncation` for a
partial deletion overlapping an exon or any breakend breakpoint, and
`feature_elongation` for an intragenic, exonic duplication or insertion.
Inversions are copy-number-neutral and carry no headline term (matching
[ensembl-vep#79](https://github.com/Ensembl/ensembl-vep/issues/79)). The term
sets are validated against the VEP REST API in `tests/vep_concordance_sv.rs`.

Two deliberate divergences, both because the SV path is pure interval geometry
(no FASTA, no codon translation):

- **No precise codon consequences at breakpoints.** Where an SV breakpoint
  removes a stop/start codon or hits a splice site, VEP may add `stop_lost`,
  `start_lost`, `stop_gained`, `frameshift_variant`, or a splice term;
  `vareffect` reports only the region-level set.
- **Symbolic insertions are annotated at the insertion point**, not an
  `SVLEN`-projected reference footprint, so the touched sub-regions (and
  `feature_elongation`) can differ when the point and footprint straddle a
  region boundary.

Breakend mate pairing (`MATEID` / `EVENT`) is not performed — each breakpoint is
annotated independently, matching VEP. `Breakend::parse` decodes the VCF breakend
ALT grammar so callers can extract and annotate the mate locus.

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
| `transcript_ablation`                     | yes    | Deletion spanning a whole transcript (`annotate_interval`, DEL + full span). |
| `transcript_amplification`                | yes    | Duplication spanning a whole transcript (`annotate_interval`, DUP + full span). |
| `feature_elongation`                      | yes    | Intragenic, exonic duplication or insertion (`annotate_interval` DUP / `annotate_sv_insertion`). |
| `feature_truncation`                      | yes    | Partial deletion overlapping an exon, or a breakend breakpoint (`annotate_interval` DEL / `annotate_breakend`). |
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

- **A deletion spanning a CDS/UTR boundary is not a frameshift.** VEP's
  `frameshift` predicate opens with `return 0 unless defined $bvfo->cds_start
  && defined $bvfo->cds_end`, so a deletion that begins in the 5'UTR or runs
  into the 3'UTR never draws `frameshift_variant` — however many CDS bases it
  removes. This crate replicates that on both ends: `{5_prime_UTR_variant,
  start_lost}` at the 5' boundary and `{stop_lost, 3_prime_UTR_variant}` at the
  3'. Neither result carries an NMD prediction or a termination codon, since
  the reading frame is undefined. A deletion beginning at CDS offset 1 or 2
  leaves `cds_start` defined and does draw `frameshift_variant` alongside
  `start_lost`. Still divergent on these variants: `hgvs_p` (VEP emits
  `p.Met1_?5`-style range notation, this crate emits nothing) and
  `protein_start` (VEP `None`, this crate `Some(1)`).
- **NMD is measured at the variant site, not the termination codon.**
  `predicts_nmd` applies the 50-nucleotide rule to the variant's own CDS
  position, for VEP `NMD.pm` parity. For a frameshift the termination
  codon is downstream and can fall in a later exon, so `predicts_nmd` may
  report NMD where the real termination codon escapes it.
  `ConsequenceResult::ptc` reports that codon's position and an NMD
  verdict computed there; it is the field to use for PVS1. Worked
  example: BRCA1 `NM_007294.4` c.5266dup (`p.Gln1756ProfsTer74`,
  `predicts_nmd = true`) terminates at residue 1829 (c.5485), inside the
  last coding exon c.5468-5592 — an NMD escape.
- **NMD with final-exon 3'UTR tails.** The junction set is built from CDS
  segments only, so when a transcript's 3'UTR spans exons beyond the
  CDS-end exon, the junction used here sits **upstream** of the mRNA's true
  last junction. Distances are therefore understated and real NMD can be
  missed — the residual is under-reporting NMD, not over-reporting it. A
  variant at or past the last coding exon always escapes
  (`nmd.rs`, `cds_position >= last_junction_cds_pos`), so the opposite
  failure cannot occur.
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
