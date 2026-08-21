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
  public field. External struct-literal *and* functional-update syntax
  (`..other`) no longer compile; build via `ConsequenceResult::new` and assign
  the fields you need. Pattern matches using `..` are unaffected. Marking the
  struct `#[non_exhaustive]` now means future field additions will not be
  breaking.
- **Breaking.** `ConsequenceResult` derives `PartialEq`, so two results that
  compared equal on 0.4.0 may differ on 0.5.0 once `ptc` is populated.
- `predicts_nmd` is documented honestly. **Its value is unchanged for every
  input** — only the doc comment moved. It measures the 50-nucleotide rule at
  the **variant site**, matching Ensembl VEP's `NMD.pm` convention. For
  `stop_gained` the variant site is the termination codon, so it remains
  correct there. For `frameshift_variant` the termination codon lies
  downstream, often in a later exon, and this field ignores that displacement.

### Migration

If you use `predicts_nmd` for ACMG PVS1, switch to `ptc`. The two disagree
whenever a frameshift's termination codon lands in a different exon than the
variant. Worked example: BRCA1 `NM_007294.4` c.5266dup reports
`predicts_nmd = true`, but terminates at residue 1829 (c.5485), inside the last
coding exon (c.5468-5592) — `ptc` correctly reports the NMD escape.

`annotate_to_vep_json` output and the CLI's CSQ columns are unchanged; neither
carries `predicts_nmd` or `ptc`, since neither has a VEP counterpart.
