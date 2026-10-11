# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For releases before `0.5.0`, see the git history.

## [0.6.0]

RefSeq and GENCODE transcript stores, upstream/downstream gene variants, and
a broad pass of VEP 116 concordance fixes. Stores from 0.5 must be rebuilt
with `vareffect setup`.

### Highlights

- **Three transcript sources.** `vareffect setup` builds MANE (default), the
  full curated RefSeq set (RS_2025_08, 93,432 transcripts) or GENCODE
  comprehensive (v50, 644,292), chosen with
  `[vareffect.transcripts] source = "mane" | "refseq" | "gencode"`. RefSeq and
  GENCODE downloads are pinned by SHA-256, release and genome build. Every
  store records its source and release, which `vareffect check` validates and
  `vareffect annotate` writes as `##vareffect_transcripts=<source> <release>`.
- **Upstream and downstream gene variants.** `annotate` also returns
  transcripts within 5,000 bases of a variant, as VEP's default
  `--distance 5000` does, with the new `ConsequenceResult::distance`.
- **Transcripts MANE never had.** Incomplete CDSs, short frameshift introns,
  RefSeq transcripts that differ from the genome, selenocysteine codons and
  pseudoautosomal copies are annotated as VEP 116 does (see
  `VEP_DIVERGENCES.md` 3).

### Breaking changes

- **Store format.** `transcript_models.bin` now starts with a 12-byte header
  (magic `VEFXTXM\0`, `u32` format version 3) before a MessagePack payload
  holding the source, release and models. Stores from 0.5 and earlier are
  refused with `UnsupportedStoreFormat`; a 0.5 reader fails on a 0.6 store, so
  the CLI and library upgrade together. Empty stores and stores with trailing
  bytes are refused. Identify a store by (`format_version`, `payload_sha256`)
  from its manifest; the file SHA-256 changes.
- **`annotate` returns more results.** Transcripts within 5 kb that the
  variant does not overlap now appear with `upstream_gene_variant` or
  `downstream_gene_variant` (MODIFIER), and `intergenic_variant` means no
  transcript in range. A consumer that picks one transcript per variant, or
  treats any result as "in a gene", should drop results whose only terms are
  these two. In `vareffect annotate` output the CSQ layout is unchanged, but
  records near a transcript gain entries. The structural-variant entry points
  are unchanged.
- **`TranscriptModel`** is `#[non_exhaustive]` (build it with
  `TranscriptModel::new` and the `with_*` setters), `tier` is
  `Option<TranscriptTier>`, and `encode_transcript_store` takes the source and
  release: `encode_transcript_store(source, source_release, transcripts)`.
- **`TranscriptStore`** allows one accession on chrX and chrY
  (pseudoautosomal copies); `get_by_accession` returns the chrX copy. Building
  a store fails on an `exon_count` that disagrees with `exons`, a CDS phase
  above 2, two CDS segments in one exon (ribosomal slippage), or a non-zero
  first phase on a CDS marked start-complete.
- **Enums.** `VariantLocation`, `IndelLocation` and `IndelRegion` are
  `#[non_exhaustive]`. New variants break exhaustive matches:
  `Consequence::NonCodingTranscriptVariant` (`severity_rank` from
  `UpstreamGeneVariant` on moves up by one),
  `VariantLocation::FrameshiftIntron`, and `VarEffectError::UnsupportedStoreFormat`,
  `Serialize`, `EmptyStore` and `CdsNumberingUnavailable`.
- **Coordinates on incomplete CDSs.** On a 5'-incomplete CDS, CDS offsets,
  codon numbers and `resolve_hgvs_c` use VEP's padded frame (the first coding
  base is offset `LocateIndex::cds_start_pad()`). Unchanged for complete
  CDSs, which includes all of MANE.
- **Insertions** are located at their insertion point: on the minus strand
  `IndelRegion::Cds` offsets move one base later, an insertion at an
  exon/intron junction has no exon or intron number, and
  `crosses_exon_boundary` is always `false` for an insertion.
- **Config.** An unknown key inside `[vareffect]` or its tables is an error,
  so a misspelled table cannot silently fall back to MANE.

### Added

- `TranscriptSource` (`Mane`, `RefSeq`, `Gencode`),
  `TranscriptStore::{source, source_release, try_from_data, from_bytes,
  get_all_by_accession}`, `encode_transcript_store`, `StoreFormat` and
  `transcript::{store_format, TRANSCRIPT_STORE_FORMAT_VERSION,
  TRANSCRIPT_STORE_HEADER_LEN}`.
- `TranscriptModel` fields `cds_start_complete`, `cds_end_complete`,
  `genome_mismatch`, `cds_alignment_gap`, `utr_numbering_breaks` and
  `selenocysteine_codons`; `ConsequenceResult::distance` (also in the VEP JSON
  as `distance`); `locate::FRAMESHIFT_INTRON_MAX_GAP`.
- `vareffect build-transcripts --source mane|refseq|gencode`. Every build
  writes `<store>.excluded.tsv` (transcripts not admitted, and why) and a
  manifest with the source, counts, skipped rows, `format_version`,
  `payload_sha256` and `source_file_sha256`.
- `vareffect check` loads the store as `annotate` does and fails when the
  config, the store and its manifest disagree on source or release.

### Changed: annotation output, following VEP 116

- **Incomplete CDSs** (GENCODE `cds_start_NF`/`cds_end_NF`, partial RefSeq
  CDSs): a partial first or last codon gets VEP's terms, amino acid `X` and
  `p.?`; no start terms without a complete start codon; no stop terms without
  a stop codon on the genome.
- **Frameshift introns** (exon gaps of 13 bases or fewer): a variant inside one
  is `coding_sequence_variant` or `non_coding_transcript_variant`, with no
  splice or intron term. `ptc` skips them when it locates the last junction.
- **RefSeq transcripts that differ from the genome**: on a model flagged
  `cds_alignment_gap`, HGVS, positions, codons and amino acids are withheld
  and `ptc` is `Indeterminate`; past a `utr_numbering_breaks` bound, cDNA
  positions and HGVS naming them are withheld. Terms, impact, exon numbers and
  `predicts_nmd` are kept.
- **Selenocysteine**: `TGA` at a listed codon reads as `U` (`Sec`) in the
  reference peptide (MANE reads them from `transl_except`). A Sec codon SNV is
  `missense_variant` or `stop_gained` instead of `stop_lost`, and
  `codon::aa_three_letter` maps `U` to `"Sec"`.
- **Insertions** are classified at their insertion point (fixing a one-base
  error on the minus strand, TP53 `c.523_524insT` `p.Arg175LeufsTer6`).
  HGVS is 3'-shifted on the genome, across exon/intron junctions, with the
  inserted bases rotated. In-frame insertions follow VEP's codon window,
  terms and HGVS p. (`p.Gln52_Trp53insTer`, `p.Pro626_Gly628dup`).
- **Deletions**: HGVS is 3'-shifted on the genome, across exon/intron
  junctions (`c.1497+1del`); a coding deletion shifted out of the CDS has no
  HGVS p.
- **Frameshifts** show VEP's codon window (`gat/gaTt`, `D/DX`), add
  `stop_gained` when the window reads a new stop, and follow VEP's
  start-codon rules (`start_retained_variant` when `ATG` is kept). Through the
  stop codon they are extensions (`p.Ter426CysextTer4`, `p.Ter362AsnextTer?`);
  one that keeps the stop is `p.Ter298=`, and `stop_retained_variant` (LOW,
  no PTC) when the submitted variant reaches into the stop codon.
- **Alleles with `N`**: no HGVS or amino acids, and peptide terms become
  `coding_sequence_variant`, as in VEP.

### Fixed

- Variants near a 3' partial codon no longer fail with `Malformed`, which
  dropped every transcript's annotation for the variant.
- HGVS p. of a frameshift deletion just before the stop codon is no longer the
  malformed `p.Ter362?fsTer?` when the transcript has a 3'UTR to read.

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
