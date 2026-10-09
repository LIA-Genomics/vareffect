# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For releases before `0.5.0`, see the git history.

## [0.6.0]

### Added

- A versioned transcript store format. `transcript_models.bin` now opens
  with a 12-byte header (magic `VEFXTXM\0` plus a `u32` little-endian format
  version, 2 in this release) followed by a MessagePack payload: the
  transcript source, its release and the `Vec<TranscriptModel>`. The version
  is checked before the payload is decoded.
- `TranscriptSource` (`Mane`, `RefSeq`, `Gencode`; `#[non_exhaustive]`), and
  `TranscriptStore::source`, `TranscriptStore::source_release` and
  `TranscriptStore::try_from_data`. `vareffect setup` records MANE and its
  release in the store.
- `TranscriptModel::new` and `with_*` setters. New model fields:
  `cds_start_complete` / `cds_end_complete` (GENCODE `cds_start_NF` /
  `cds_end_NF`, RefSeq partial CDS ranges), `genome_mismatch`,
  `cds_alignment_gap` and `selenocysteine_codons`. MANE models are complete,
  genome-matching and carry no selenocysteine codons yet.
- Multi-locus accessions: one chrX and one chrY copy of a pseudoautosomal
  transcript. `TranscriptStore::get_all_by_accession` returns every copy in
  canonical chromosome order (`chrom::canonical_order_key`, then
  `tx_start`). The builder admits a repeated accession only as one chrX and
  one chrY copy, each overlapping a pseudoautosomal region, and refuses any
  other.
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
- `vareffect build-transcripts --source mane|refseq|gencode` builds a store
  from the MANE GFF3 (default, unchanged output), NCBI RefSeq's full GRCh38
  annotation (curated `NM_` / `NR_` on the primary assembly, written as
  `transcript_models_refseq.bin`) or GENCODE comprehensive
  (`transcript_models_gencode.bin`). The CDS rows are taken as they are:
  `stop_codon` rows are never added to them (GENCODE's CDS already contains
  the stop). RefSeq and GENCODE models record CDS completeness (a RefSeq range marker
  counts only on the CDS's first or last row; an interior one is a coverage
  break), RefSeq's genome-mismatch flags (in-row `Gap` operations and
  transcript stretches aligned nowhere between `cDNA_match` blocks) and both
  sources' selenocysteine codons. Every
  build writes `<stem>.excluded.tsv`, the candidate transcripts not
  admitted and why (RefSeq alternate loci and patches, ribosomal-slippage
  CDSs with two segments in one exon, ...), and the manifest gains
  `source`, `counts`, `skipped_rows`, `anomalies`, `excluded_file` and
  `excluded_sha256`. The builder loads the encoded store back before
  writing it.

### Changed

- **Breaking.** `TranscriptModel` is `#[non_exhaustive]` (build it with
  `TranscriptModel::new`), `TranscriptModel::tier` is an
  `Option<TranscriptTier>` (`None` for a transcript without a MANE or RefSeq
  Select tag), and `encode_transcript_store` takes the source and release:
  `encode_transcript_store(source, source_release, transcripts)`.
- **Breaking.** `TranscriptStore` accepts an accession on two chromosomes;
  the uniqueness rule is now per (accession, chromosome).
  `get_by_accession`, and with it HGVS c. resolution, returns the first copy
  in canonical chromosome order (chrX for a PAR transcript). Construction
  also fails when `exon_count` disagrees with `exons`, or when a CDS with a
  non-zero first phase is marked `cds_start_complete`.
- **Breaking.** On a CDS marked `cds_start_complete = false`, codon 1 is not a
  start codon even when it reads ATG: no `start_lost`,
  `start_retained_variant` or `p.Met1?` (VEP 116 ENST00000616125.5 `c.1A>G`
  is `missense_variant` `p.Met1Val`). On `cds_end_complete = false` the CDS
  has no stop codon on the genome. A model flagged `cds_alignment_gap` gets
  no HGVS c. or p., cDNA / CDS / protein positions, codons or amino acids,
  its termination codon is `PtcStatus::Indeterminate` and `predicts_nmd` is
  false (consequence terms, impact and exon / intron numbers are kept), and
  `resolve_hgvs_c` refuses it with the new
  `VarEffectError::CdsNumberingUnavailable` (see `VEP_DIVERGENCES.md` 3.3).
  No MANE 1.5 output changes.
- **Breaking.** Stores written by vareffect 0.5 and earlier are refused with
  `UnsupportedStoreFormat { found: StoreFormat::Unversioned, .. }`, and
  format 1 stores from development builds of 0.6 with
  `found: StoreFormat::Version(1)`. Rebuild
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
- **Breaking.** `locate_indel` classifies an insertion at its insertion
  point, as VEP 116 does: on the minus strand `IndelRegion::Cds` offsets move
  one base later (the insertion goes before `cds_offset_start` on both
  strands); an insertion at an exon/intron junction is the mRNA insertion
  between the two exons, with no exon or intron number; one next to a UTR
  base is a UTR insertion; the splice flags follow VEP's two-flank rule; and
  `crosses_exon_boundary` is always `false` for an insertion.
- **Breaking.** A store with no transcripts, or with bytes after the payload,
  is an error at load time. The builder refuses a GFF3 that admits no
  transcripts instead of writing an empty store, which would annotate every
  variant as overlapping nothing.
- The SHA-256 of a rebuilt store changes: the payload now holds the source,
  the release and the new model fields. Identify a store by
  (`format_version`, `payload_sha256`) from its manifest, never by
  `payload_sha256` alone; the 0.5 file SHA does not link to any 0.6 value.
  The MANE 1.5 GRCh38 models themselves are unchanged: every 0.5 field
  holds the same value, with `tier` wrapped in `Some`.

### Fixed

- Frameshifts follow VEP 116's codon window (`TranscriptVariationAllele`
  `codon` / `peptide`): an insertion inside a codon shows that codon and its
  peptide (`gat/gaTt`, `D/DX`), one between codons shows no reference codon
  (`-/T`, `-/X`), and a peptide that is a lone stop gets no `X` (`CE/*`).
  A frameshift whose window reads a new stop is also `stop_gained` (2,713
  ClinVar x MANE 1.5 entries), except a stop inserted right before the
  stop codon. At the start codon the frameshift path now
  applies VEP's `_ins_del_start_altered` and `_inv_start_altered`: an indel
  that keeps or rebuilds `ATG`, or only shifts it into the 5'UTR, is
  `start_retained_variant` (with `start_lost` when `ATG` no longer opens the
  CDS), and only an insertion inside codon 1 touches it (84 entries).
- An alternate allele with an `N` (or another base outside A, C, G, T) is
  annotated as VEP 116 does: the terms read from the peptide (`missense`,
  `stop_gained`, `protein_altering`, ...) become `coding_sequence_variant`,
  `frameshift_variant` stays, and there is no HGVS c., HGVS p. or amino acid
  change. Codons and positions are kept, for SNVs too. Before, such alleles could
  be called `stop_gained`, `missense_variant` or `protein_altering_variant`
  (6 of 668 ClinVar N alleles) and got an HGVS c. A frameshift's
  termination codon is `PtcStatus::Indeterminate`.
- Insertions are annotated at their insertion point (VEP 116). Before, an
  insertion was classified by a single genomic base, which on the minus
  strand lies on the wrong side of it: minus-strand coding insertions were one
  base early (TP53 `c.523_524insT` was `p.Arg175SerfsTer6`, VEP
  `p.Arg175LeufsTer6`), and each exon/intron or UTR/CDS edge was wrong on one
  strand (plus: `splice_donor_variant` at the donor junction,
  `frameshift_variant` + `start_lost` before the first coding base; minus:
  `splice_acceptor_variant` at the acceptor junction, a frameshift after the
  last coding base). Canonical-site insertions no longer carry
  `coding_sequence_variant` and `intron_variant`, and splice terms follow
  VEP's insertion predicates. In-frame insertion HGVS p. now sees the residue
  before a split codon (NM_006231.4 `c.2465_2467dup`:
  `p.Lys822_Gly823insGlu`). Over GRCh38 ClinVar on MANE 1.5, 15,549 of
  4,667,273 transcript rows move, all insertions; no HGVS c. changes. In a
  stratified sample (up to 20 variants per term-change bucket; 484 variants,
  517 transcript rows), terms and HGVS p. both match VEP 116 in 390 rows,
  against none before.
- Insertion HGVS follows VEP 116's 3' shift. The shift is on the genome, so
  it crosses exon/intron junctions (NM_000548.5 `c.336+1_336+2dup`, was
  `c.336_336+1insGT`), reads at most 1000 bases and moves at most 1001 minus
  the insertion length, as VEP's. The inserted bases rotate with the shift,
  which also finds duplications the unrotated bases hid (NM_000051.4
  `c.-32_-31dup`, was `c.-31_-31+1insGT`), and a duplication of the bases 3'
  of the insertion is checked before one of the 5' bases, as VEP does. HGVS
  p. is written at the shifted point: none when a coding insertion shifts
  out of the CDS, into the intron or past the stop codon (VEP; NM_003036.4
  `c.2152_*3dup`), and an in-frame insertion translates the rotated bases
  (NM_001378454.1 `p.His2513_Ala2514insAsn`, was `insIle`). Consequences,
  NMD and PTC stay on the unshifted insertion. Over GRCh38 ClinVar on MANE
  1.5, 12,634 of 4,667,273 transcript rows move, all insertions, in HGVS c.
  (12,605) and HGVS p. (4,437) only; the 63 rows shifted past the stop codon
  all match VEP 116. In a stratified sample (up to 20 per change bucket, 385
  rows), HGVS c. matches VEP 116 in
  all 385 (none before) and HGVS p. in 371 (232 before); the other 14 are
  in-frame protein duplications and stop insertions VEP writes differently
  (a separate fix) and one frameshift at the stop codon, none of which
  matched before either.
- In-frame insertions follow VEP 116's codon window, terms and HGVS p.
  Between two codons only the inserted codon is shown (NM_000059.4
  `c.44_45insATT`: `-/TTA`, `-/L`, was `attttt/attTTAttt`, `IF/ILF`) and the
  protein position spans both flanking residues. An insertion inside a codon
  whose peptide neither starts nor ends with the reference residue is
  `protein_altering_variant` (NM_000546.6 `c.523_524insTTT`, was
  `inframe_insertion`), and a stop-carrying one may be `stop_gained` alone;
  a stop inserted just before the native stop is `stop_retained_variant`
  (NM_003036.4 `c.2152_*3dup`, was `stop_gained`), and an insertion before
  the stop codon whose bases carry a stop is now `stop_gained` with the
  termination codon past the native one (NM_000546.6
  `c.1179_1180insCTGTAA`). Inside codon 1, VEP's start-codon rules apply
  (`c.2_3insGAT`: `start_retained_variant`, was `inframe_insertion`;
  `c.1_2insCGA`: `start_lost,start_retained_variant`). HGVS p. shifts the
  inserted peptide C-terminal across the whole protein before checking for a
  duplication (NM_001854.4 `p.Pro626_Gly628dup`, was `p.Gly625_Pro627dup`),
  writes a stop insertion as an insertion (NM_000546.6
  `p.Gln52_Trp53insTer`, was `p.Gln53Ter`), gives a stop-codon insertion
  VEP's `delins` (`p.Ter394delinsTrpGly`, was none) and no notation for a
  non-duplicating insertion just before the stop codon. An insertion with
  an `N` is `coding_sequence_variant` with no amino acids or HGVS p., as in
  VEP (was translated, e.g. a false `stop_gained`). Over GRCh38 ClinVar
  insertions on MANE 1.5, 4,113 transcript rows move, all in-frame CDS
  insertions: codons and protein position (4,013), amino acids (4,113),
  HGVS p. (2,062), terms (562). In a stratified VEP 116 sample (up to 20 per
  change bucket, 460 rows) codons, amino acids and protein position match
  in all 460 (330, 285 and 330 before), terms in 458 (177 before) and HGVS
  p. in 453 (157 before); the others are `start_lost` HGVS p. and two
  `N`-containing minus-strand insertions VEP calls `stop_retained_variant`
  (see `VEP_DIVERGENCES.md`).

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
