# VEP divergences

Where `vareffect`'s output differs from Ensembl VEP, and which VEP features
it does not cover.

- **Target:** VEP releases 115 and 116, GRCh38. GRCh37 and CHM13 are not
  validated (the transcript and genome binaries would need rebuilding).
- **Transcripts:** MANE Select, MANE Plus Clinical and RefSeq Select. A
  variant overlapping only transcripts outside the loaded store is
  `intergenic_variant`.

Contents:

1. [Summary](#1-summary)
2. [Consequence terms](#2-consequence-terms)
3. [Incomplete transcripts](#3-incomplete-transcripts)
4. [NMD](#4-nmd)
5. [Structural variants](#5-structural-variants)
6. [Coverage and scope](#6-coverage-and-scope)
7. [Validation](#7-validation)
8. [Known edge cases](#8-known-edge-cases)

## 1. Summary

| Area | VEP | `vareffect` | Section |
|---|---|---|---|
| Splice sub-terms | `splice_polypyrimidine_tract_variant`, `splice_donor_region_variant`, `splice_donor_5th_base_variant` | `splice_region_variant` | [2.1](#21-splice-region-sub-terms) |
| Intronic variant, non-coding transcript | `non_coding_transcript_variant` | `intron_variant` | [2.2](#22-intronic-variants-of-non-coding-transcripts) |
| Deletion across a CDS/UTR boundary | `p.Met1_?5`-style HGVS p., no `protein_start` | no HGVS p., `protein_start` 1 | [2.3](#23-deletions-across-a-cdsutr-boundary) |
| Partial codon, HGVS p. | `p.Ter…` | `p.?` | [3.1](#31-partial-first-and-last-codons) |
| Deletion removing a CDS end with no stop | `stop_lost` / `stop_retained_variant` | no stop term | [3.1](#31-partial-first-and-last-codons) |
| Splice-region windows near a frameshift intron | spill up to 7 exon bases (SNVs), stop at 8 (indels) | 3 exon bases (SNVs), up to 16 (indel acceptor side) | [3.2](#32-frameshift-introns) |
| NMD | `NMD.pm` rules, biotype term `NMD_transcript_variant` | last-junction rule only; `ptc` field for PVS1 | [4](#4-nmd) |
| Structural variants | codon and splice terms at breakpoints | region-level terms only | [5](#5-structural-variants) |
| HGVS 3' shift | `--shift_hgvs` toggle | always on | [6.2](#62-not-yet-implemented) |

## 2. Consequence terms

### 2.1 Splice region sub-terms

VEP 116 emits three finer splice-region terms:
`splice_polypyrimidine_tract_variant` (acceptor intron −3 to −17),
`splice_donor_region_variant` (donor intron +3, +4, +6) and
`splice_donor_5th_base_variant` (donor intron +5). `vareffect` reports all
three as `splice_region_variant` (plus `intron_variant` where applicable).
To recover the finer term, reclassify from the splice offsets that
`locate_indel` exposes.

For an insertion on an exon short enough that two introns' windows reach
it, VEP sets the splice-region flag intron by intron, so its result depends on
the order it visits them; `vareffect` takes either intron's region (accepted, maintainer decision 2026-10-09: LOW term, rare).

### 2.2 Intronic variants of non-coding transcripts

VEP may give `non_coding_transcript_variant`; `vareffect` gives
`intron_variant`, and the non-coding context is in
`TranscriptModel::biotype`. Inside a frameshift intron `vareffect` gives
`non_coding_transcript_variant` as VEP does.

### 2.3 Deletions across a CDS/UTR boundary

VEP's `frameshift` predicate needs both `cds_start` and `cds_end`, so a
deletion starting in the 5'UTR or running into the 3'UTR never draws
`frameshift_variant`. `vareffect` matches: `{5_prime_UTR_variant,
start_lost}` and `{stop_lost, 3_prime_UTR_variant}`, with no NMD prediction
or termination codon. A deletion starting at CDS offset 1 or 2 does draw
`frameshift_variant` with `start_lost`.

**Difference:** VEP writes `p.Met1_?5`-style HGVS p. and no
`protein_start`; `vareffect` writes no HGVS p. and `protein_start` 1.

## 3. Incomplete transcripts

These cases do not occur in MANE except TMEM247 NM_001424184.1; they matter
for GENCODE and full RefSeq transcripts.

### 3.1 Partial first and last codons

A CDS can start mid-codon (GENCODE `cds_start_NF`, GFF3 phase 1 or 2 on the
first CDS row) or end mid-codon (`cds_end_NF`, or a stop completed only by
polyadenylation, as in MT-ND4 and TMEM247). `vareffect` follows VEP 116:

| Aspect | Behavior (as VEP 116) |
|---|---|
| HGVS `c.` | The missing bases of a partial first codon are counted only for the exonic position of a substitution (`c.3` is the first genomic base after a phase-1 start); everywhere else genomic bases are counted (`c.191+6T>A`, `c.4del` for the base whose substitution is `c.5`). `resolve_hgvs_c` reads both. |
| cDNA | Transcript bases only. |
| Codons, amino acids | `nnT/nnA` for a padded first codon, `Ta/Aa` for a partial last codon; amino acid `X`. A left-over partial codon translates to `X`, unless the complete codons are a lone stop. |
| Padded first codon | `coding_sequence_variant` (an MNV too, never `synonymous_variant`). No `start_lost` or `start_retained_variant` without a complete first codon. |
| Partial last codon | `incomplete_terminal_codon_variant` + `coding_sequence_variant`; `protein_altering_variant` when a deletion removes exactly the partial codon; `stop_gained` when an insertion creates a stop (`X/*`, `X/W*X`). |
| Insertion before the partial codon | `inframe_insertion` + `incomplete_terminal_codon_variant`, inserted bases translated on their own: `GGG` `-/G`; `G`, `GG` `-/X` + `coding_sequence_variant`; `TGA`, `TGAG` `-/*`, `TGAGGG` `-/*G`, `GGGTGAG` `-/G*X`, each with `stop_gained`. |
| Minus-strand insertion after the last CDS base | `3_prime_UTR_variant` with a 3'UTR (ENST00000654733.1 `c.*1dup`), `downstream_gene_variant` without one (ENST00000518351.2). |
| Deletion ending inside the partial codon | Frameshift or in-frame deletion (`c.658_659del` after 3' shifting). |
| Stop-codon terms | None on a CDS that does not end in a stop codon on the genome. |

**Differences:**

- **HGVS p. is `p.?`.** VEP prints `p.Ter1=`, `p.Ter309del`, `p.TerLeu1=`,
  naming a stop codon that is not there; for a stop an insertion creates at
  the partial codon it writes `p.Ter220=`, `p.Ter220delinsTer` or nothing.
  The residue is unknown, so `vareffect` writes `p.?`.
- **No stop term on a deletion that removes the CDS end** (maintainer
  decision, 2026-10-07). For a deletion starting before the partial codon
  VEP reports `frameshift_variant,stop_lost` (`c.656_659del`) or
  `stop_retained_variant,3_prime_UTR_variant` (`c.656_*2del`). There is no
  stop to lose or retain: `vareffect` reports
  `incomplete_terminal_codon_variant` + `coding_sequence_variant` (+
  `3_prime_UTR_variant` when the deletion reaches the 3'UTR), with no HGVS
  p.
- **HGVS on deletions at the CDS end.** VEP omits or clips it for some
  (`c.656_*2del` written as `c.656_659del`; nothing for a deletion of the
  partial codon into the 3'UTR, or of ENST00000518351.2's whole partial
  codon). `vareffect` writes the deletion it annotates, with `p.?` when only
  the partial codon is removed and no HGVS p. when the deletion reaches the
  3'UTR.
- **No `FLAGS` field.** VEP lists `cds_start_NF` / `cds_end_NF` there.

### 3.2 Frameshift introns

An exon gap of 13 bases or fewer is a frameshift intron, a correction of the
genome against the transcript (VEP: `abs(intron_end - intron_start) <= 12`
on 1-based bounds). `vareffect` follows VEP 116:

| Variant | Result (as VEP 116) |
|---|---|
| An end inside the gap, or an insertion between two gap bases | `coding_sequence_variant` (`non_coding_transcript_variant` on a non-coding transcript), intron number, no splice term from the gap. A real intron's canonical site in the same footprint keeps its term. |
| Deletion spanning the gap | Normal frameshift or in-frame deletion of its coding bases, no splice term from the gap. |
| Coding base next to a one-base gap | VEP's two-base canonical windows spill onto it: `splice_acceptor_variant` before the gap, `splice_donor_variant` after it (HIGH), for SNVs and indels (ENST00000400999.7 `c.228C>A`, `c.229del`). |
| Insertion at a gap edge | Frameshift plus that splice term. |

**Differences:**

- **Splice-region windows** (maintainer decision, 2026-10-07). VEP
  computes its splice-region windows from the gap bounds, so they spill
  onto exon bases. `vareffect` does the same for indels on the donor side
  (`c.233_235del`: `splice_donor_region_variant`), but:
  - for an SNV it gives exonic `splice_region_variant` only to the three
    bases next to the gap, so VEP's `c.232T>A`
    `splice_donor_5th_base_variant` is not reproduced;
  - for an indel its 17-base acceptor window reaches up to 16 bases into
    the upstream exon where VEP stops at 8 (`c.218_220del`,
    `c.212_214del`, ENST00000436261.7 `c.1465_1467del`: VEP
    `inframe_deletion` only, `vareffect` adds `splice_region_variant`).

  Accepted because `splice_region_variant` is LOW and not a loss-of-function
  term, and MANE 1.5 has no frameshift intron.
- **Codon case.** For a gap-spanning deletion VEP uppercases as many codon
  bases as the genomic footprint is long (`tCCGAg/tag`); `vareffect`
  uppercases the coding bases removed (`tCCGag/tag`).

## 4. NMD

`vareffect` reports NMD in two fields:

| Field | Measured at | Use |
|---|---|---|
| `predicts_nmd` | the variant's own CDS position, as VEP's `NMD.pm` | VEP parity |
| `ptc` | the termination codon (residue and NMD verdict) | ClinGen SVI PVS1 (Abou Tayoun et al. 2018, PMID 30192042) |

For a frameshift the termination codon lies downstream and may sit in a
later exon, so `predicts_nmd` can report NMD where the real termination
codon escapes it (BRCA1 NM_007294.4 c.5266dup terminates at c.5485, in the
last coding exon). Single-exon genes (`PRNP`, the mitochondrial genome)
have no junction: `predicts_nmd` is false and `ptc` is NMD-negative.

**Differences from `NMD.pm`:**

- No escape for a variant in the first ~100 coding bases (`NMD.pm`:
  `$variant_coding_region <= 101`).
- Junctions come from CDS segments, not all exons. When a 3'UTR spans
  further exons the junction used sits upstream of the true last one: NMD
  can be missed, never overcalled.
- A frameshift intron is an exon boundary for `predicts_nmd` (as in
  `NMD.pm`) but not for `ptc`: no exon-junction complex sits at a genome
  correction.
- `ptc`'s penultimate-exon window is 50 coding bases (ClinGen SVI: "within
  the 3'-most 50 nucleotides of the penultimate exon"); `NMD.pm` uses 52
  genomic bases.
- No `NMD_transcript_variant`: it is a transcript biotype term. Filter
  NMD-biotype transcripts in the caller.

**How `ptc` decides:**

- **Only the last-junction rule**, the only one the SVI tree uses (as does
  `AutoPVS1`, Xiang et al. 2020, PMID 32442321). Start-proximal escape
  (150 nt) and the long-exon rule (> 407 bp) of Lindeboom et al. (PMIDs
  27618451, 31659324) are deliberately not applied: they would withdraw
  PVS1 from many established loss-of-function alleles, a guideline decision
  rather than an annotation one. Re-initiation is not modeled;
  `ptc.protein_position` lets a consumer add these rules.
- **The codon's 3'-most base decides.** "Occurs within the 3'-most 50
  nucleotides" is read as "overlaps". Using the first base would shrink the
  window to 48 nucleotides and overcall NMD; across GRCh38 ClinVar this
  moves 202 of 4,666,135 rows, all from NMD to escape.

## 5. Structural variants

`annotate_interval` (DEL, DUP, INV), `annotate_breakend` (`<BND>`) and
`annotate_sv_insertion` (`<INS>`) give VEP's per-transcript term set: the
overlapped regions plus a copy-number term (`transcript_ablation`,
`transcript_amplification`, `feature_truncation`, `feature_elongation`).
Inversions carry no copy-number term
([ensembl-vep#79](https://github.com/Ensembl/ensembl-vep/issues/79)).
Breakend mates are annotated independently, as in VEP; `Breakend::parse`
decodes the ALT grammar to find the mate locus. Validated against VEP REST
in `tests/vep_concordance_sv.rs`.

**Differences** (the SV path is pure interval geometry, with no FASTA and
no translation):

- No codon or splice consequences at breakpoints (VEP may add `stop_lost`,
  `start_lost`, `frameshift_variant` or a splice term).
- A symbolic insertion is annotated at its insertion point, not an
  `SVLEN`-projected footprint.

## 6. Coverage and scope

### 6.1 SO terms

| SO term | Status | Note |
|---|---|---|
| `missense_variant` | yes | |
| `synonymous_variant` | yes | |
| `stop_gained` | yes | Drives NMD prediction. |
| `stop_lost` | yes | Extension found by a 3'UTR stop scan. |
| `start_lost` | yes | SNVs and partial start-codon deletions. |
| `start_retained_variant` | yes | |
| `stop_retained_variant` | yes | |
| `frameshift_variant` | yes | Drives NMD prediction. |
| `inframe_insertion` | yes | |
| `inframe_deletion` | yes | |
| `protein_altering_variant` | yes | Length-preserving complex changes. |
| `coding_sequence_variant` | yes | Ambiguous codons, partial codons, frameshift introns of coding transcripts. |
| `incomplete_terminal_codon_variant` | yes | CDS length not divisible by 3. |
| `splice_donor_variant` | yes | Intron +1, +2; the exon base after a one-base frameshift intron. |
| `splice_acceptor_variant` | yes | Intron −1, −2; the exon base before a one-base frameshift intron. |
| `splice_region_variant` | yes | Intron ±3..±8 (acceptor to −17), exon 1..3. |
| `splice_polypyrimidine_tract_variant` | no | Reported as `splice_region_variant` (2.1). |
| `splice_donor_region_variant` | no | Reported as `splice_region_variant` (2.1). |
| `splice_donor_5th_base_variant` | no | Reported as `splice_region_variant` (2.1). |
| `5_prime_UTR_variant` | yes | |
| `3_prime_UTR_variant` | yes | |
| `intron_variant` | yes | |
| `non_coding_transcript_exon_variant` | yes | |
| `non_coding_transcript_variant` | partial | Frameshift introns of non-coding transcripts only (2.2). |
| `upstream_gene_variant` | yes | |
| `downstream_gene_variant` | yes | |
| `intergenic_variant` | yes | No overlapping transcript in the store. |
| `transcript_ablation` | yes | DEL spanning a whole transcript. |
| `transcript_amplification` | yes | DUP spanning a whole transcript. |
| `feature_elongation` | yes | Intragenic exonic DUP or `<INS>`. |
| `feature_truncation` | yes | Partial DEL over an exon, or a breakend. |
| `NMD_transcript_variant` | no | Biotype term; see `predicts_nmd` and `ptc` (4). |
| `mature_miRNA_variant` | no | Not yet implemented. |
| `TF_binding_site_variant` | no | Out of scope. |
| `TFBS_ablation` | no | Out of scope. |
| `regulatory_region_variant` | no | Out of scope. |

### 6.2 Not yet implemented

- **`mature_miRNA_variant`**: needs a miRNA locus track.
- **A toggle to switch off HGVS 3' normalization.** `vareffect` always
  shifts 3' on the coding strand, as VEP's default `--shift_hgvs` (since
  release 109).
- **Other genome builds** (GRCh37, CHM13, non-human): the work is in the
  build pipelines; the runtime is build-agnostic.
- **Multi-allelic VCF rows**: `annotate` takes one ALT; split beforehand.
- **Canonical transcript selection**: results carry the MANE / RefSeq
  Select flags; filtering is the caller's job.

### 6.3 Out of scope

- **Regulatory annotation** (`TFBS_ablation`, `TF_binding_site_variant`,
  `regulatory_region_variant`): a separate layer with its own feature data.
- **Plugins** (CADD, REVEL, SpliceAI, ...): wrap results downstream.
- **Co-located variant lookup** (gnomAD, dbSNP, ClinVar, allele
  frequencies): an enrichment step after consequence prediction.
- **Custom annotation overlays** (BED, GFF, VCF): use bcftools or bedtools.
- **VCF parsing in the library**: `vareffect` takes chrom, pos, ref, alt
  (the CLI reads VCF).

## 7. Validation

Expected values come from the Ensembl VEP REST API on real transcripts. The
suites are `#[ignore]`-gated because they need the store and genome:

```bash
FASTA_PATH=/absolute/path/to/GRCh38.bin \
    cargo test -p vareffect --release -- --ignored
```

| Test file | Variants | Focus |
|---|---:|---|
| `vep_concordance_snv.rs` | 20 | SNV consequences, splice sites and regions |
| `vep_concordance_indel.rs` | 28 | Frameshift, in-frame, boundary and splice indels |
| `vep_concordance_hgvs.rs` | 20 | HGVS c. formatting |
| `vep_concordance_hgvs_p.rs` | 30 | HGVS p. formatting |
| `vep_concordance_hgvs_reverse.rs` | 20 | HGVS c. to genomic coordinates |
| `vep_concordance_normalization.rs` | 18 | 3' normalization, intergenic, NMD rule |
| `incomplete_cds_concordance.rs` | 224 | Sections 3.1 and 3.2, on GENCODE, RefSeq and MANE models in `tests/data/` |

The first six (TP53, BRCA1, BRCA2, APC, BRAF, EGFR, CFTR, ERBB2, PRNP and
others) all passed on 2026-04-11, 136 / 136 variants. The incomplete-CDS
suite passed on 2026-10-07 with the normalizations listed in its header.
When adding a variant, record VEP's output (`refseq=1&hgvs=1`, `&numbers=1`
for exon numbers) next to it.

## 8. Known edge cases

- **Patch contigs.** Lookups on patch contigs need
  `VarEffect::open_with_patch_aliases` with a `refseq,ucsc` alias CSV;
  otherwise they return `ChromNotFound`.
- **IUPAC codes.** A few GRCh38.p14 patch regions contain ambiguity bases
  (`M`, `R`, `Y`); they are kept, and a codon containing one translates to
  `Xaa`.
- **Mitochondrion.** `chrM` uses NCBI translation table 2; `chrMT` and `MT`
  resolve to `chrM` through the reader's alias table.
