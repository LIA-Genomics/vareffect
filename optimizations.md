# vareffect performance optimizations

Concrete, ROI-ordered optimization opportunities surfaced while profiling the
GRCh37/GRCh38 large-concordance harness. Each entry references the actual code
location and proposes a specific change.

## Context

Measured baseline (single-threaded release build, warm OS page cache):

| Configuration | GRCh38 (45 644 vars) | GRCh37 (5 220 vars) |
|---|---|---|
| VRS on (`annotate_with_options(.., emit_vrs_ids: true)`) | 49 570 v/s | 1 975 v/s |
| VRS off (default `annotate(..)`) | 110 268 v/s | 22 314 v/s |
| Pre-VRS hardcoded baseline | 134 219 v/s | 27 737 v/s |

The GRCh37 set has the same SNV/ins/del/delins ratio (50/10/35/5) as GRCh38 but
~11x larger average deletion size (13.1 bp vs 1.2 bp). Per-base cost in the
indel paths is what drives the gap — these notes are about cutting that cost.

The single largest residual VRS overhead is the per-chromosome SHA-512 cache
fill (~5.7 s for all 25 GRCh38 primaries, paid once per `FastaReader`).
Amortizes invisibly on long jobs but dominates short ones — see the VRS section
at the end.

---

## Indel hot-path optimizations

### 1. `left_align_indel`: O(n^2) -> O(n) via cursor

**File:** `vareffect/src/var_effect.rs:480-507`

The hot loops use `Vec::<u8>::remove(0)` to strip shared prefix bytes:

```rust
while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
    r.remove(0);   // O(n) memmove of every remaining byte
    a.remove(0);   // same
    pos += 1;
}
```

For a 1 kb deletion with 10 bp of shared prefix, that is 10 x 1000 = 10 000
byte shifts to do work that should be 10 cursor increments. The trailing-suffix
loop (uses `pop()`) is fine; the prefix loop is the killer.

**Fix.** Cursor-based slicing — borrow the input bytes once and walk start/end
indices, allocating `String`s only at the end from the surviving slices.

```rust
let r_bytes = ref_allele.as_bytes();
let a_bytes = alt_allele.as_bytes();
let (mut r_lo, mut r_hi) = (0usize, r_bytes.len());
let (mut a_lo, mut a_hi) = (0usize, a_bytes.len());
// suffix-equality: decrement r_hi/a_hi
// prefix-equality: increment r_lo/a_lo
```

| Effort | ~30 min |
|---|---|
| Risk | Low (pure refactor of an internal hot loop) |
| Expected GRCh37 throughput delta | +20-40 % on the indel half of the test set |

---

### 2. `compute_3prime_shift_deletion`: bounded probe instead of whole-region fetch

**File:** `vareffect/src/normalize.rs:51-73` (plus strand) and `:80-104`
(minus strand)

The shift loop itself is correctly O(shift). The **fetch** is over-eager:

```rust
let fetch_end = if let Some((_, exon_end)) = containing_exon_range(...) {
    exon_end                     // entire exon, can be 100s of kb
} else if let Some((_, intron_end)) = containing_intron_range(...) {
    intron_end                   // entire intron, often >50 kb
} ...;
let seq = fasta.fetch_sequence(chrom, del_start, fetch_end)?;
```

A 1 bp deletion mid-intron forces a 50+ kb fetch to scan ~10 bp of repeat
shift. mmap is cheap per byte but this still pays the `Vec::with_capacity` +
memcpy of 50 kb plus cache pollution.

**Fix.** Doubling probe — fetch 64 bp past `del_end`; if the scan reaches the
fetch end without breaking AND there is still room, double and retry. Bound
wasted fetch to ~2x actual shift distance:

```rust
let mut probe = 64u64;
loop {
    let bound = (del_end + probe).min(exon_end);
    let seq = fasta.fetch_sequence(chrom, del_start, bound)?;
    let shift = scan_shift(&seq, del_len);
    let exhausted = (del_start + del_len + shift as u64) >= bound;
    if !exhausted || bound == exon_end {
        return Ok(shift);
    }
    probe *= 2;
}
```

Mirror change on the minus-strand branch. Most variants converge on the first
64 bp probe.

| Effort | ~1 hr |
|---|---|
| Risk | Low (output is identical; only the read window changes) |
| Expected throughput delta | +10-20 % across both builds |

---

### 3. `fetch_cds_sequence` minus-strand: drop one allocation

**File:** `vareffect/src/consequence/helpers.rs:189-198`

For minus-strand transcripts the inner loop is three full passes over the
data:

```rust
let chunk = fasta.fetch_sequence(chrom, gstart, gend)?;   // alloc 1
let rc = crate::codon::reverse_complement(&chunk);        // alloc 2
seq.extend_from_slice(&rc);                               // copy 3
```

`rc` is gratuitous — complement-while-reversing into `seq` directly:

```rust
let chunk = fasta.fetch_sequence(chrom, gstart, gend)?;
seq.reserve(chunk.len());
for &b in chunk.iter().rev() {
    seq.push(crate::codon::complement(b));
}
```

One pass, zero extra allocation. For a 1 kb minus-strand deletion that is
~1 kb less heap traffic per call, multiplied by the number of overlapping
transcripts (3-5 typically).

| Effort | ~15 min |
|---|---|
| Risk | Low |
| Expected throughput delta | +5-10 % on minus-strand indels |

---

### 4. Borrowed `FastaReader::fetch_sequence`

**File:** `vareffect/src/fasta.rs:540`

Today `fetch_sequence` always allocates an owned `Vec<u8>` even though the
underlying storage is a memory-mapped contiguous slice the kernel has already
page-cached. The `verify_ref` method already does this borrow-and-compare
inline.

**Fix.** Add a sibling `fetch_sequence_borrowed(&self, ..) -> Result<&[u8], _>`
returning a slice into the mmap. Move the read-only callers in the indel path
to it:

- `compute_3prime_shift_deletion` (`normalize.rs`) — only scans, never owns. Pure win.
- `is_duplication` (`hgvs_c.rs`) — same. Pure win.
- `fetch_cds_sequence` plus-strand path — borrow, then `extend_from_slice` into the contiguous output. Still one copy into the codon buffer, but no intermediate `Vec`.

Minus-strand stays on the owned form (the in-place complement needs a write
target).

| Effort | ~2-3 hr (broader API change) |
|---|---|
| Risk | Medium (lifetime plumbing across multiple call sites) |
| Expected throughput delta | +5-15 % spread across the indel paths |

---

### 5. Cross-transcript redundancy in the dispatcher

**File:** `vareffect/src/consequence/mod.rs::annotate` (per-transcript loop)

For a locus with N overlapping transcripts on one strand, the FASTA fetch for
the deletion + 3' shift is essentially repeated N times — once per transcript
— even though they read the same plus-strand bytes. Per-transcript work that
**must** stay per-transcript: exon/intron boundary detection, codon-frame
computation, HGVS notation. Per-transcript work that is actually
**locus-level**: the genomic-base fetches.

**Fix.** Hoist strand-invariant FASTA reads (deletion bases, 3' shift probe
bytes) one level up into the per-variant dispatcher. Pass the resulting slice
down to each transcript pass. Keep the strand-aware codon stitch in
`fetch_cds_sequence`.

Bigger refactor — only worth it after #1-4 land. Delivers the most on dense
multi-transcript regions (HLA, immunoglobulin loci, alt-spliced disease genes).

| Effort | ~1 day |
|---|---|
| Risk | Medium-high (touches the dispatcher contract) |
| Expected throughput delta | +10-30 % on multi-transcript loci |

---

### 6. Cap display strings for very large indels

**Files:** `vareffect/src/codon.rs::format_codons_indel`,
`format_amino_acids_indel`

These build display strings whose size scales with deletion length. VEP
truncates with `...` past a threshold; vareffect today does not. A 1500 bp
inframe deletion produces a ~1500-char codon string per CSQ row — allocated,
formatted, and threaded through the rest of the pipeline.

**Fix.** Cap output to match VEP behavior (typically 100 codons / 33 amino
acids), eliding the middle with `...`.

| Effort | ~30 min |
|---|---|
| Risk | Low (improves VEP concordance on large-indel rows) |
| Expected throughput delta | Marginal CPU; meaningful memory saving on outliers |

---

## VRS-side optimizations (orthogonal — only relevant when VRS is enabled)

The VRS pipeline currently costs ~13 us per variant on warm cache, which is
why GRCh38 throughput drops from 134k -> 49k v/s when VRS is enabled. The
opt-in `AnnotateOptions::emit_vrs_v1` / `emit_vrs_v2` flags eliminate this
for clients that do not need VRS, but the cost is still worth attacking
when it is on.

### V1. Eager parallel SQ-digest fill — DONE

Shipped as `VarEffect::warm_vrs_cache(assembly)`. Fans the 25 primary-contig
SHA-512 hashes across rayon's pool out of the hot path. Lazy fill remains
the default; eager warm is opt-in for callers about to enter a VRS-emitting
batch.

### V2. Per-schema VRS opt-in — DONE

Shipped as independent `emit_vrs_v1` / `emit_vrs_v2` flags on
`AnnotateOptions` (replacing the single master `emit_vrs_ids`). Each
schema's serialize step is gated on its own flag; the shared upstream
(VOCA + SQ digest) runs once even when both are on.

### V3. Direct-to-bytes serialization

**File:** `vareffect/src/vrs/serialize.rs`

Today both schemas build a `serde_json::Value` tree, serialize it to bytes,
then SHA-512 the bytes. The Value tree allocates `BTreeMap`s, `String`s, and
`Vec`s for fields whose canonical form is statically known.

**Fix.** Write the canonical form directly into a thread-local `Vec<u8>` with
`write!` (or just `extend_from_slice` of literals + `itoa` for the integers).
Pre-size, reuse the buffer. Cuts VRS cost by another ~30-40 %. The
`vrs-python` golden-vector tests already pin the exact expected bytes, so
regressions here are caught by unit tests.

---

## Suggested attack order

For the indel-heavy GRCh37 path the user is investigating, the order is:

1. **#1** (`left_align_indel` cursor) — highest single ROI, smallest change.
2. Measure GRCh37 large-concordance throughput.
3. **#2** (bounded 3' shift probe) if step 2 is not enough.
4. **#3** (drop minus-strand rev-complement Vec) — bundle with #2.
5. **#4** (borrowed `fetch_sequence`) only if profiling still shows
   FASTA-fetch alloc as a top hit.
6. **#5** (cross-transcript hoist) only after #1-4 — risk/reward shifts
   only when the per-transcript cost is already minimized.
7. **#6** any time; correctness benefit independent of perf.

For VRS-on workflows, V1 and V2 have shipped. **V3** (direct-to-bytes
serialization) is the remaining polish layer.

## Validation

Every change here has the same validation contract:

```bash
make fmt && make lint && cargo test -p vareffect

# GRCh38 perf
GRCH38_FASTA=data/vareffect/GRCh38.bin cargo test --release -p vareffect \
  --test vep_large_concordance_grch38 -- --ignored --nocapture vep_large_concordance_grch38

# GRCh37 perf (the dataset most affected by indel-path changes)
GRCH37_FASTA=data/vareffect/GRCh37.bin cargo test --release -p vareffect \
  --test vep_large_concordance_grch37 -- --ignored --nocapture
```

Concordance rates (consequence / HGVS c. / HGVS p.) are the correctness
guardrail — they must not move. Throughput is the performance signal.
