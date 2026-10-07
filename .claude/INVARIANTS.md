# vareffect invariants — the binding charter

Human-owned required behavior. Agents quote these obligations and flag conflicts; they do not amend them without an explicitly authorized human revision. Code establishes current behavior, not permission to violate this charter: a disagreement may be an implementation defect.

Status: initial draft written by an agent when the harness was adapted from LIA (2026-10-07). It binds once the maintainer ratifies it; until then treat it as the intended policy and report conflicts the same way. IDs are stable and never reused.

## V-1 · Coordinates

Every genomic coordinate inside the library is 0-based, half-open (BED/UCSC). Conversion from 1-based genomic inputs (VCF `POS`, GFF3) happens once, at the boundary that reads them, and the inverse conversion once at the boundary that writes them. HGVS `c.`/`n.`/`p.` positions are transcript- or protein-relative, not genomic: they are mapped through the transcript model, never by a genomic offset. Exons are stored in 5'→3' transcript order. Failure produces silent off-by-one annotation that every downstream tool inherits.

References: `CLAUDE.md` (Critical Conventions), [transcript model types](../vareffect/src/types.rs), [GFF3 builder](../vareffect-cli/src/builders/transcript_models/mod.rs).

## V-2 · Determinism

The same variant, transcript store and genome produce byte-identical results regardless of thread count, scheduling, hash-map iteration order or platform. Output order is defined, not incidental. Failure makes concordance runs, diffs and downstream records irreproducible.

References: [annotation entry point](../vareffect/src/var_effect.rs), [parallel VCF annotation](../vareffect-cli/src/annotate.rs).

## V-3 · Unknown is not negative

A result that could not be computed — unreadable reference, transcript or codon not located, reference mismatch, unsupported variant shape — stays distinguishable from a computed negative answer through the public API and every output format. It is an error, `None`, or an explicit indeterminate state, never `false`, an empty consequence set or `intergenic_variant` by default. Failure turns an availability problem into an assertion a clinical consumer acts on (for example a missed or invented loss-of-function call).

References: [`PtcStatus`](../vareffect/src/consequence/mod.rs), [errors](../vareffect/src/error.rs).

## V-4 · VEP concordance

Output targets Ensembl VEP release 115/116 on the same transcripts. Every known difference is either a defect to fix or an intentional divergence recorded, with its reason, in [`VEP_DIVERGENCES.md`](../vareffect/VEP_DIVERGENCES.md). Adopting an intentional divergence is a maintainer decision. Failure leaves users unable to tell a bug from a design choice.

## V-5 · Transcript and reference provenance

The transcript set, its source release and the genome build an annotation was produced with are knowable from the built data (manifest and checksum), and the builder never silently drops, rewrites or re-coordinates source records: every exclusion is counted and reported. A store or genome built by an incompatible version fails to load loudly rather than deserializing into a different meaning. Failure attributes annotations to a transcript model the source never contained.

References: [store build sidecars](../vareffect-cli/src/common.rs), [transcript store](../vareffect/src/transcript.rs).

## V-6 · Public contract

The public Rust API, the VEP-compatible JSON, the VCF `CSQ` field layout and the on-disk store formats are contracts. Changing their shape or meaning is a breaking change: semver-signalled and recorded in [`CHANGELOG.md`](../CHANGELOG.md). A change of meaning under an unchanged shape is still breaking. `VarEffect`, `TranscriptStore` and `FastaReader` stay `Send + Sync`. Failure breaks consumers silently.

## V-7 · No embedded data

The library ships no reference genome, transcript set or derived table that changes annotation results; all such data is built offline by the CLI from declared, checksummed sources. Specification constants are exempt — genetic code tables, SO terms, and assembly accession maps — provided each names its source and release in code. Failure hides provenance and pins users to stale data.

## V-8 · Representation and normalization

API inputs are plus-strand genomic alleles on the loaded assembly. HGVS output is 3'-shifted on the transcript strand (always on), while consequences are computed at the submitted position, as VEP does by default. The same variant given in different equivalent representations (left-aligned, anchor-padded, trimmed) yields the same HGVS; any consequence difference between representations is a documented divergence or a defect. Failure lets a reviewer "fix" representation dependence in a way that silently diverges from VEP.

References: [normalization](../vareffect/src/normalize.rs), [left alignment](../vareffect/src/left_align.rs), [`VEP_DIVERGENCES.md`](../vareffect/VEP_DIVERGENCES.md) ("HGVS 3' normalization is always on").

## V-9 · Chromosome and build

Chromosome names resolve through the declared alias rules (UCSC, RefSeq accession, Ensembl, `chrM`/`MT`, patch aliases) or fail; an unrecognized chromosome is an error, never an empty result. The REF allele is checked against the loaded genome before annotation, and `chrM` uses NCBI genetic code table 2. Failure turns a build or naming mismatch into silent absence of annotation.

References: [chromosome mapping](../vareffect/src/chrom.rs), [genetic codes](../vareffect/src/codon.rs).

## Known conflicts with current code

These clauses state required behavior that current code does not yet meet. Each is a defect to fix as separate, measured work (or a reason to reword the clause on ratification), not something a review of unrelated work should block on. Remove an entry only when its fix lands.

- **V-3 (VCF output).** An annotation error, REF mismatch or caught panic passes the VCF line through without `CSQ`, byte-identical to an intergenic variant; only stderr counters differ (`vareffect-cli/src/annotate.rs`, module docs and the `catch_unwind` branch).
- **V-3 / V-9 (structural variants).** `annotate_interval`, `annotate_breakend` and `annotate_sv_insertion` do not check that the chromosome exists, and the overlap query returns an empty result for an unknown chromosome, so an unrecognized name yields the documented "trustworthy" empty result. Fixing it changes contract meaning (V-6).
- **V-5 (load and exclusions).** The transcript store is a bare MessagePack vector with no format version, and the library never reads its checksum or manifest; only the genome index carries a version. The genome build writes no manifest or checksum. The GFF3 builder skips malformed rows (short rows, unparsable coordinates, unknown strand, genes without an ID) without counting them.
- **V-2 (platform), open.** Per-transcript result order follows the interval tree's callback order, whose callback type differs by SIMD backend; whether order is identical across platforms is unverified.

## The standing human-authority limit

Releasing, accepting an intentional VEP divergence, choosing the authoritative transcript source or release, and accepting a breaking change are maintainer decisions. Agent review and passing checks replace none of them. vareffect output feeds clinical interpretation downstream (for example PVS1 through LIA); an agent's judgment about that effect is a review aid, never clinical sign-off.
