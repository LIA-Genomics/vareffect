# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For releases before `0.5.0`, see the git history.

## [0.6.0]

### Added

- A versioned transcript store format. `transcript_models.bin` now opens
  with a 12-byte header (magic `VEFXTXM\0` plus a `u32` little-endian format
  version) followed by the MessagePack `Vec<TranscriptModel>`. The version is
  checked before the payload is decoded.
- `TranscriptStore::from_bytes`, for callers that read the file themselves
  (for example to hash exactly the bytes they load).
  `TranscriptStore::load_from_path` now reads the file and calls it.
- `encode_transcript_store`, the one supported way to write a store, and
  `StoreFormat` (both re-exported at the crate root), plus
  `transcript::{store_format, TRANSCRIPT_STORE_FORMAT_VERSION,
  TRANSCRIPT_STORE_HEADER_LEN}`.
- `VarEffectError::UnsupportedStoreFormat`, `VarEffectError::Serialize` and
  `VarEffectError::EmptyStore`.
- Store manifests record `format_version`, `payload_sha256` (SHA-256 of the
  bytes after the header) and `source_file_sha256` (SHA-256 of the GFF3 the
  store was built from).
- `vareffect check` loads the transcript store the way `annotate` does and
  reports its format version and record count, so an unversioned,
  other-version, truncated or empty store fails the check.
- Incomplete CDS support in the engine, following VEP 116. GFF3 phase on
  the first CDS row is read: `LocateIndex::cds_start_pad` gives the bases
  missing from a 5'-incomplete first codon. CDS positions, codon numbers and
  protein positions count them; cDNA positions do not; HGVS `c.` counts them
  only for the exonic position of a substitution, as VEP does. Variants in a
  partial first or last codon (SNVs, MNVs, deletions, insertions) get VEP's
  terms, codons and amino acid `X`, with HGVS p. `p.?`. A deletion that
  removes the end of a CDS without a stop codon gets no stop term. It has no
  HGVS p. when it starts before the partial codon (`c.656_659del`) or runs
  from the partial codon into the 3'UTR (`c.658_*2del`); deleting exactly
  the partial codon is `protein_altering_variant` with `p.?`. An insertion
  inside the partial codon translates every codon it completes
  (`stop_gained` for `X/W*X`). An
  insertion between the last complete codon and the partial codon is VEP's
  `inframe_insertion` + `incomplete_terminal_codon_variant` with the
  inserted bases translated (`stop_gained` when they contain a stop,
  `coding_sequence_variant` when the translation has an `X`). A
  minus-strand insertion after the last CDS base is `3_prime_UTR_variant`,
  or `downstream_gene_variant` when the transcript has no 3'UTR (see
  `VEP_DIVERGENCES.md`).
- Frameshift introns (exon gaps of `locate::FRAMESHIFT_INTRON_MAX_GAP` = 13
  bases or fewer), following VEP 116: a variant inside one is
  `coding_sequence_variant` or `non_coding_transcript_variant`, a deletion
  across one is a normal coding change, and on a one-base gap the splice
  windows spill onto the flanking coding bases.
- `Consequence::NonCodingTranscriptVariant`,
  `VariantLocation::FrameshiftIntron` and `locate::FRAMESHIFT_INTRON_MAX_GAP`.

### Changed

- **Breaking.** Stores written by vareffect 0.5 and earlier are refused with
  `UnsupportedStoreFormat { found: StoreFormat::Unversioned, .. }`. Rebuild
  them with `vareffect setup` (the transcript store rebuilds in seconds from
  cached sources). A 0.5 reader given a 0.6 store fails at byte 0, so the CLI
  and every library consumer upgrade together.
- **Breaking.** The new `VarEffectError` variants break exhaustive matches on
  the error type. The enum stays exhaustive on purpose, so downstream code is
  told about every new failure mode at compile time.
- **Breaking.** `VariantLocation::CdsExon::cds_offset` (and `codon_number`,
  `codon_position`, `IndelRegion::Cds` offsets) count in VEP's padded frame:
  on a 5'-incomplete CDS the first genomic coding base is offset
  `cds_start_pad()`, not 0. `resolve_hgvs_c` maps `c.N` on such a transcript
  as VEP numbers it, which differs by the pad between substitutions and
  other variant types, and refuses a substitution inside the pad. Unchanged
  for every complete CDS, which includes all of MANE 1.5.
- **Breaking.** `locate_variant` returns `VariantLocation::FrameshiftIntron`
  instead of a splice or intron location inside an exon gap of 13 bases or
  fewer, and indels touching such a gap change terms as described under
  Added. MANE 1.5 has no such gap.
- **Breaking.** `VariantLocation`, `IndelLocation` and `IndelRegion` are now
  `#[non_exhaustive]`. The new `Consequence` variant breaks exhaustive
  matches on `Consequence`, and `Consequence::severity_rank` values from
  `UpstreamGeneVariant` on move up by one.
- **Breaking.** `start_lost`, `start_retained_variant` and `p.Met1?` are no
  longer emitted on a CDS without a complete first codon, and the
  stop-codon deletion pattern (`stop_lost` + `3_prime_UTR_variant`) needs a
  CDS that ends in a stop codon on the genome. A variant in a 3' partial
  codon now reports VEP's terms, codons, amino acid `X`, protein position
  and `p.?`; before, SNVs, MNVs and delins got only
  `incomplete_terminal_codon_variant`, and deletions and insertions a
  frameshift call. In MANE 1.5 this changes only NM_001424184.1 (TMEM247),
  whose stop is completed by polyadenylation; no ClinVar GRCh38 record
  moves.
- `ConsequenceResult::ptc` skips frameshift introns when it locates the last
  exon-exon junction: no exon-junction complex sits at a sequence
  correction. `predicts_nmd` keeps VEP's NMD plugin behavior and counts
  them. Unchanged on MANE 1.5, which has none.
- **Breaking.** `LocateIndex::build` returns `Malformed` for a transcript
  with more than one CDS segment in one exon (a ribosomal-slippage CDS), a
  CDS phase above 2 or a CDS segment pointing at a missing exon, and
  `TranscriptStore` construction fails with it, instead of keeping only the
  last segment, computing with the bad phase or panicking. MANE 1.5 has
  none of these.
- **Breaking.** A store with no transcripts, or with bytes after the payload,
  is an error at load time. The builder refuses a GFF3 that admits no
  transcripts instead of writing an empty store, which would annotate every
  variant as overlapping nothing.
- The whole-file SHA-256 of a rebuilt store changes because of the header.
  The models do not: for the MANE 1.5 GRCh38 store, `payload_sha256` equals
  the 0.5 file's SHA-256 (`edff7c1416b2989a882186360c1782fefb33274edb6efe0441f5236a0ec38b84`).
  To link provenance records across this upgrade, compute
  `sha256(bytes[TRANSCRIPT_STORE_HEADER_LEN..])` over the bytes you loaded and
  compare it with the 0.5 file SHA. That link holds only for this
  header-only transition; across later format versions, key on
  (`format_version`, `payload_sha256`).

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
- A deletion spanning the 5'UTR/CDS boundary into the initiator codon now
  reports `{5_prime_UTR_variant, start_lost}` regardless of how many CDS bases
  it removes, and no longer names a termination codon or predicts NMD against a
  reading frame that no longer exists. Two things were wrong: a 1-2 base
  deletion took the frameshift path and asserted a downstream PTC (41 rows),
  and the `5_prime_UTR_variant` term was missing throughout (190 rows). VEP
  withholds `frameshift_variant` on this geometry however long the deletion —
  its `frameshift` predicate short-circuits on the undefined `cds_start`
  (`VariationEffect.pm::frameshift`), the same short-circuit this crate already
  replicated at the 3' end for `stop_lost`. Verified against recorded VEP
  output for seven such variants, e.g. `ALPL` `NM_000478.6`
  chr1:21554077 `GCACCATGATTT>G`, which removes 7 CDS bases and is still not a
  frameshift to VEP. A deletion beginning at CDS offset 1 or 2 leaves
  `cds_start` defined and does draw both terms; that case keeps the
  pure-CDS-path treatment.

### Migration

If you use `predicts_nmd` for ACMG PVS1, switch to `ptc`. The two disagree
whenever a frameshift's termination codon lands in a different exon than the
variant. Worked example: BRCA1 `NM_007294.4` c.5266dup reports
`predicts_nmd = true`, but terminates at residue 1829 (c.5485), inside the last
coding exon (c.5468-5592) — `ptc` correctly reports the NMD escape.

`annotate_to_vep_json` output and the CLI's CSQ columns are unchanged; neither
carries `predicts_nmd` or `ptc`, since neither has a VEP counterpart.
