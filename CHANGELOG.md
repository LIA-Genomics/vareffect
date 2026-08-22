# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For releases before `0.5.0`, see the git history.

## [0.5.0]

### Added

- `ConsequenceResult::ptc`, a `PtcStatus` reporting the premature termination
  codon's residue position and an NMD verdict computed **at that codon**. This
  is the input the ClinGen SVI PVS1 decision tree asks for (Abou Tayoun et al.
  2018, PMID 30192042, doi:10.1002/humu.23626). Populated for `stop_gained`
  and `frameshift_variant`; `PtcStatus::NotApplicable` otherwise.
- `PtcStatus`, a four-state enum. `NoStopCodon` (HGVS `fsTer?` — the alternate
  frame reaches the end of the 3'UTR with no stop) is kept distinct from
  `Indeterminate` (the codon could not be located) because the two carry
  opposite clinical meanings; `PtcStatus::nmd_at_ptc` returns `None` rather
  than `Some(false)` for absences of an answer.
- `ConsequenceResult::new`, the constructor for building a minimal result
  outside this crate.

### Changed

- **Breaking.** `ConsequenceResult` is now `#[non_exhaustive]` and has gained a
  public field. Two distinct things stop compiling outside this crate:
  - **Construction** — struct-literal *and* functional-update syntax
    (`ConsequenceResult { .., ..other }`). Build via `ConsequenceResult::new`
    and assign the fields you need.
  - **Destructuring** — any field pattern must now end in `..`, *including*
    one that already names every field. `let ConsequenceResult { transcript,
    .. } = r;` is fine; the same pattern without `..` is not.

  Marking the struct `#[non_exhaustive]` means future field additions will no
  longer be breaking.
- **Breaking.** `ConsequenceResult` derives `PartialEq`, so two results that
  compared equal on 0.4.0 may differ on 0.5.0 once `ptc` is populated.
- The termination-codon NMD verdict is measured from the **3'-most base** of
  that codon. ClinGen SVI states the rule as "the premature termination codon
  occurs in the 3' most exon or within the 3'-most 50 nucleotides of the
  penultimate exon"; a codon spans three bases, so "occurs within" is read as
  *overlaps*, and its last base is what decides it. Anchoring on the first base
  would shrink the escape window to 48 nucleotides and predict NMD for two
  codon positions the guideline places outside it — an overcall, since
  predicting NMD is what earns PVS1 its full strength. Only the last-junction
  rule is applied; the start-proximal (150 nt) and long-exon (407 bp) escape
  rules of Lindeboom et al. (PMID 27618451, PMID 31659324) are deliberately
  omitted, because the SVI tree does not use them and applying them would
  withdraw PVS1's full strength from a large class of established
  loss-of-function alleles. Across the GRCh38 ClinVar release this moves 202 of
  4,666,135 annotation rows, all from NMD-predicted to NMD-escaping.
- `predicts_nmd` is documented honestly. Its value is unchanged except on the
  start-codon path noted under **Fixed** below, where it stopped asserting NMD
  against a destroyed reading frame. It measures the 50-nucleotide rule at
  the **variant site**, matching Ensembl VEP's `NMD.pm` convention. For
  `stop_gained` the variant site is the termination codon, so it remains
  correct there. For `frameshift_variant` the termination codon lies
  downstream, often in a later exon, and this field ignores that displacement.

### Fixed

- An `inframe_deletion` overlapping the reference terminator no longer reports
  `stop_gained`. Deleting a codon next to the stop can fuse the flanking bases
  into a fresh terminator, but that is the normal stop re-formed one codon
  earlier, not a premature one — and the HIGH impact it carried was a
  truncation claim about a protein that terminates where it always did. VEP
  guards its own predicate the same way (`VariationEffect.pm::stop_gained`:
  `($alt_pep =~ /\*/) and ($ref_pep !~ /\*/)`); four of the five sibling call
  sites in this crate already did. Two rows in the GRCh38 ClinVar release:
  `HRAS` `NM_005343.4` p.Ser189del and `USH1C` `NM_005709.4` p.Phe552del, both
  deletions of the protein's final residue.
- A boundary-spanning deletion that removes part of the initiator codon now
  reports `start_lost` alongside `frameshift_variant`, emits `p.Met1?`, and
  withholds both the NMD prediction and the termination codon. It previously
  named a downstream termination codon and asserted NMD against a reading frame
  that no longer exists. VEP fires `start_lost` and `frameshift` independently
  here (`VariationEffect.pm::start_lost` via `_ins_del_start_altered`), which
  the pure-CDS path in this crate already matched. 41 rows in the GRCh38
  ClinVar release.

### Migration

If you use `predicts_nmd` for ACMG PVS1, switch to `ptc`. The two disagree
whenever a frameshift's termination codon lands in a different exon than the
variant. Worked example: BRCA1 `NM_007294.4` c.5266dup reports
`predicts_nmd = true`, but terminates at residue 1829 (c.5485), inside the last
coding exon (c.5468-5592) — `ptc` correctly reports the NMD escape.

`annotate_to_vep_json` output and the CLI's CSQ columns are unchanged; neither
carries `predicts_nmd` or `ptc`, since neither has a VEP counterpart.
