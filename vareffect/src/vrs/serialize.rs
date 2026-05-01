//! Recursive `ga4gh_serialize` implementation for VRS 1.3 and VRS 2.0.
//!
//! VRS computed identifiers use a content-addressed digest (`sha512t24u`)
//! over a canonicalized JSON serialization of the typed object. The
//! canonicalization is recursive: nested **identifiable** objects (those
//! with their own GA4GH digest) are hashed first and substituted as
//! their **bare 32-character digest string** (no `ga4gh:` prefix, no
//! type prefix) in the parent object before the parent is hashed.
//!
//! Per the VRS spec and `vrs-python` reference implementation
//! (`ga4gh/core/_internal/identifiers.py`):
//!
//! - JSON canonicalization: keys sorted in unicode order, no
//!   insignificant whitespace, RFC 7159 escape semantics.
//! - Pre-existing `ga4gh:NS.xxx` CURIEs (e.g. `sequence_id` /
//!   `refgetAccession`) are reduced to the bare digest portion before
//!   serialization.
//! - `_id` and any underscore-prefixed keys are stripped.
//! - `null` values are stripped.
//! - Unordered arrays are sorted lexicographically by element.
//!
//! For our VRS allele use case the inputs are always ASCII (canonical
//! base symbols and base64url digests), so the escape edge cases do
//! not apply.
//!
//! The canonical bytes are produced by writing literal segments and
//! `itoa`-formatted integers directly into a thread-local `Vec<u8>`,
//! bypassing any `serde_json::Value` tree construction. This is sound
//! only because every input string is ASCII alphanumeric (or `-` / `_`
//! for base64url, or empty for pure deletions) — no character ever
//! requires JSON escaping. Any future field added to the canonical form
//! that takes user-controlled string content **must** be JSON-escaped,
//! or the digests will diverge from `vrs-python` and break downstream
//! interop. The pinned-byte tests below catch this drift on first run.

use std::cell::RefCell;

use super::digest::{ga4gh_identify, sha512t24u};
use super::sq_digest::bare_digest_from_curie;

// Per-thread scratch buffer used to build canonical VRS blobs without
// per-call heap churn.
//
// SAFETY-INVARIANT: no callee transitively reached while a
// `borrow_mut()` of this cell is held may itself attempt to compute a
// VRS ID (i.e. re-enter `vrs1_allele_id` / `vrs2_allele_id`), or the
// nested `borrow_mut` will panic at runtime. Audited callees today —
// `sha512t24u`, `ga4gh_identify`, `bare_digest_from_curie` — are all
// VRS-free.
thread_local! {
    static VRS_BUF: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

/// Append the minimal-decimal ASCII form of `n` to `buf`.
///
/// Uses `itoa::Buffer` (a 40-byte stack array) so this function performs
/// zero heap allocation. Output is byte-equivalent to `serde_json`'s
/// integer formatter, which itself routes through `itoa`.
#[inline]
fn write_u64(buf: &mut Vec<u8>, n: u64) {
    let mut itoa_buf = itoa::Buffer::new();
    buf.extend_from_slice(itoa_buf.format(n).as_bytes());
}

/// Build the VRS 1.3 SequenceLocation canonical serialization blob.
///
/// `sq_curie` is the full `ga4gh:SQ.<digest>` form; the canonical blob
/// uses only the bare 32-char digest portion in the `sequence_id`
/// field (per the CURIE-reduction rule).
fn write_vrs1_location(buf: &mut Vec<u8>, sq_curie: &str, start: u64, end: u64) {
    debug_assert!(
        sq_curie.is_ascii(),
        "VRS location writer requires ASCII sq_curie"
    );
    let bare_sq = bare_digest_from_curie(sq_curie);
    buf.extend_from_slice(br#"{"interval":{"end":{"type":"Number","value":"#);
    write_u64(buf, end);
    buf.extend_from_slice(br#"},"start":{"type":"Number","value":"#);
    write_u64(buf, start);
    buf.extend_from_slice(br#"},"type":"SequenceInterval"},"sequence_id":""#);
    buf.extend_from_slice(bare_sq.as_bytes());
    buf.extend_from_slice(br#"","type":"SequenceLocation"}"#);
}

/// Build the VRS 1.3 Allele canonical serialization blob.
///
/// `location_digest` is the bare 32-char `sha512t24u` digest of the
/// SequenceLocation, computed by the caller via the recursive enref
/// step. The Allele's `location` field substitutes that bare digest
/// string in place of the inlined object.
fn write_vrs1_allele(buf: &mut Vec<u8>, location_digest: &str, alt: &str) {
    debug_assert!(
        location_digest.is_ascii(),
        "VRS allele writer requires ASCII location_digest"
    );
    debug_assert!(
        alt.is_ascii(),
        "VRS allele writer requires ASCII alt; got non-ASCII bytes"
    );
    buf.extend_from_slice(br#"{"location":""#);
    buf.extend_from_slice(location_digest.as_bytes());
    buf.extend_from_slice(br#"","state":{"sequence":""#);
    buf.extend_from_slice(alt.as_bytes());
    buf.extend_from_slice(br#"","type":"LiteralSequenceExpression"},"type":"Allele"}"#);
}

/// Compute the full VRS 1.3 Allele identifier (`ga4gh:VA.<digest>`)
/// for a fully-justified normalized variant.
///
/// Inputs:
/// - `sq_curie`: `ga4gh:SQ.<digest>` for the chromosome.
/// - `start`, `end`: 0-based interbase coordinates after VOCA
///   normalization. `start == end` for pure insertions.
/// - `alt`: ASCII alt sequence after VOCA normalization. Empty string
///   for pure deletions.
pub(super) fn vrs1_allele_id(sq_curie: &str, start: u64, end: u64, alt: &str) -> String {
    VRS_BUF.with(|cell| {
        let mut buf = cell.borrow_mut();
        buf.clear();
        write_vrs1_location(&mut buf, sq_curie, start, end);
        let location_digest = sha512t24u(&buf);
        buf.clear();
        write_vrs1_allele(&mut buf, &location_digest, alt);
        ga4gh_identify("VA", &buf)
    })
}

/// Build the VRS 2.0 SequenceLocation canonical serialization blob.
///
/// Differs from 1.3 in: no `SequenceInterval` wrapper, `start`/`end`
/// are direct integers on the location, and `sequence_id` is replaced
/// by a `sequenceReference` object carrying `refgetAccession` (which
/// uses bare `SQ.<digest>` form per refget v2, no `ga4gh:` namespace).
fn write_vrs2_location(buf: &mut Vec<u8>, sq_curie: &str, start: u64, end: u64) {
    debug_assert!(
        sq_curie.is_ascii(),
        "VRS location writer requires ASCII sq_curie"
    );
    let bare_sq = bare_digest_from_curie(sq_curie);
    buf.extend_from_slice(br#"{"end":"#);
    write_u64(buf, end);
    buf.extend_from_slice(br#","sequenceReference":{"refgetAccession":"SQ."#);
    buf.extend_from_slice(bare_sq.as_bytes());
    buf.extend_from_slice(br#"","type":"SequenceReference"},"start":"#);
    write_u64(buf, start);
    buf.extend_from_slice(br#","type":"SequenceLocation"}"#);
}

/// Build the VRS 2.0 Allele canonical serialization blob.
///
/// Same Allele shape as 1.3 (location -> bare digest string, state,
/// type) — kept as a separate function so future schema drift in the
/// 2.0 Allele lands locally. The location's pre-digest blob differs
/// (see [`write_vrs2_location`]).
fn write_vrs2_allele(buf: &mut Vec<u8>, location_digest: &str, alt: &str) {
    debug_assert!(
        location_digest.is_ascii(),
        "VRS allele writer requires ASCII location_digest"
    );
    debug_assert!(
        alt.is_ascii(),
        "VRS allele writer requires ASCII alt; got non-ASCII bytes"
    );
    buf.extend_from_slice(br#"{"location":""#);
    buf.extend_from_slice(location_digest.as_bytes());
    buf.extend_from_slice(br#"","state":{"sequence":""#);
    buf.extend_from_slice(alt.as_bytes());
    buf.extend_from_slice(br#"","type":"LiteralSequenceExpression"},"type":"Allele"}"#);
}

/// Compute the full VRS 2.0 Allele identifier (`ga4gh:VA.<digest>`).
pub(super) fn vrs2_allele_id(sq_curie: &str, start: u64, end: u64, alt: &str) -> String {
    VRS_BUF.with(|cell| {
        let mut buf = cell.borrow_mut();
        buf.clear();
        write_vrs2_location(&mut buf, sq_curie, start, end);
        let location_digest = sha512t24u(&buf);
        buf.clear();
        write_vrs2_allele(&mut buf, &location_digest, alt);
        ga4gh_identify("VA", &buf)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vrs1_location_blob_is_canonical_sorted_compact_json() {
        // Manually-constructed expected blob with keys in sorted order
        // and no whitespace. Any drift in field order or spacing
        // breaks the digest contract.
        let mut buf = Vec::new();
        write_vrs1_location(
            &mut buf,
            "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            44908821,
            44908822,
        );
        let expected = r#"{"interval":{"end":{"type":"Number","value":44908822},"start":{"type":"Number","value":44908821},"type":"SequenceInterval"},"sequence_id":"IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceLocation"}"#;
        assert_eq!(std::str::from_utf8(&buf).unwrap(), expected);
    }

    #[test]
    fn vrs1_allele_blob_substitutes_location_as_bare_digest_string() {
        // The recursive enref step: `location` must be a bare digest
        // STRING in the parent's blob, not an inlined object.
        let mut buf = Vec::new();
        write_vrs1_allele(&mut buf, "esDSArZQC-Sx-96ZZzHnzAVNOc439oE5", "T");
        let expected = r#"{"location":"esDSArZQC-Sx-96ZZzHnzAVNOc439oE5","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#;
        assert_eq!(std::str::from_utf8(&buf).unwrap(), expected);
    }

    #[test]
    fn vrs2_location_blob_uses_sequence_reference_object() {
        // 2.0 differs from 1.3: no SequenceInterval wrapper, and
        // `sequenceReference` is an object not a `sequence_id` string.
        let mut buf = Vec::new();
        write_vrs2_location(
            &mut buf,
            "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            44908821,
            44908822,
        );
        let expected = r#"{"end":44908822,"sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceReference"},"start":44908821,"type":"SequenceLocation"}"#;
        assert_eq!(std::str::from_utf8(&buf).unwrap(), expected);
    }

    #[test]
    fn allele_ids_have_correct_curie_shape() {
        let id1 = vrs1_allele_id("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 100, 101, "A");
        let id2 = vrs2_allele_id("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 100, 101, "A");
        assert!(id1.starts_with("ga4gh:VA."), "id1 = {id1}");
        assert!(id2.starts_with("ga4gh:VA."), "id2 = {id2}");
        assert_eq!(id1.len(), "ga4gh:VA.".len() + 32);
        assert_eq!(id2.len(), "ga4gh:VA.".len() + 32);
        // Same input under different schemas must produce different IDs.
        assert_ne!(
            id1, id2,
            "VRS 1.3 and 2.0 IDs must differ for the same allele"
        );
    }

    #[test]
    fn empty_alt_is_valid_for_pure_deletion() {
        // Pure deletion: alt sequence is the empty string. Must not
        // panic and must produce a stable digest.
        let id1 = vrs1_allele_id("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 100, 105, "");
        let id2 = vrs2_allele_id("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 100, 105, "");
        assert!(id1.starts_with("ga4gh:VA."));
        assert!(id2.starts_with("ga4gh:VA."));
    }

    #[test]
    fn pure_insertion_uses_zero_width_interval() {
        // Pure insertion: start == end (zero-width interbase interval),
        // alt holds the inserted bases.
        let id = vrs1_allele_id("ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 100, 100, "AT");
        assert!(id.starts_with("ga4gh:VA."));
    }

    //
    // Inputs and expected bytes/digests below are copied verbatim from
    // vrs-python's `tests/test_vrs.py::test_vr` (commit on `main` as of
    // 2026-04-30):
    //
    //   start=55181319, end=55181320, refget=SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul
    //   alt="T"  ->  ga4gh:VA.Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE
    //
    // These tests pin our canonical-form output to the upstream
    // implementation byte-for-byte; any drift in field order, JSON
    // spacing, integer formatting, or refgetAccession framing breaks
    // them. They are the single strongest guarantee that `vareffect`'s
    // VRS IDs interoperate with anyvar / ClinGen / ClinVar / MAVEDB.

    const VRS_PY_SQ_CURIE: &str = "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul";
    const VRS_PY_START: u64 = 55181319;
    const VRS_PY_END: u64 = 55181320;
    const VRS_PY_ALT: &str = "T";
    const VRS_PY_LOC_DIGEST_2_0: &str = "_G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd";
    const VRS_PY_VA_DIGEST_2_0: &str = "Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE";

    #[test]
    fn vrs2_location_blob_matches_vrs_python_bytes() {
        let mut buf = Vec::new();
        write_vrs2_location(&mut buf, VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END);
        let expected = br#"{"end":55181320,"sequenceReference":{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"},"start":55181319,"type":"SequenceLocation"}"#;
        assert_eq!(&buf, expected);
    }

    #[test]
    fn vrs2_location_digest_matches_vrs_python() {
        let mut buf = Vec::new();
        write_vrs2_location(&mut buf, VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END);
        assert_eq!(sha512t24u(&buf), VRS_PY_LOC_DIGEST_2_0);
    }

    #[test]
    fn vrs2_allele_blob_matches_vrs_python_bytes() {
        let mut buf = Vec::new();
        write_vrs2_allele(&mut buf, VRS_PY_LOC_DIGEST_2_0, VRS_PY_ALT);
        let expected = br#"{"location":"_G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#;
        assert_eq!(&buf, expected);
    }

    #[test]
    fn vrs2_allele_id_matches_vrs_python() {
        let id = vrs2_allele_id(VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END, VRS_PY_ALT);
        assert_eq!(id, format!("ga4gh:VA.{VRS_PY_VA_DIGEST_2_0}"));
    }

    // Second 2.0 vector — the docstring example in
    // `vrs-python/src/ga4gh/core/identifiers.py::ga4gh_identify`:
    //   start=44908821, end=44908822, same refget
    //   ga4gh:SL.4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT
    #[test]
    fn vrs2_location_digest_matches_second_vrs_python_vector() {
        let mut buf = Vec::new();
        write_vrs2_location(&mut buf, VRS_PY_SQ_CURIE, 44908821, 44908822);
        assert_eq!(sha512t24u(&buf), "4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT");
    }

    /// Microbenchmark for the canonical-byte writer + digest path.
    ///
    /// Times `vrs2_allele_id` over a fixed iteration count using
    /// realistic inputs. Run manually with:
    ///
    /// ```sh
    /// cargo test --release -p vareffect --lib \
    ///   vrs::serialize::tests::microbench_vrs2_allele_id_throughput \
    ///   -- --ignored --nocapture
    /// ```
    ///
    /// Not part of normal CI — gated `#[ignore]` so it never runs by
    /// default. The gating tests above are what actually pin
    /// correctness; this is purely a perf observability tool.
    #[test]
    #[ignore]
    fn microbench_vrs2_allele_id_throughput() {
        const ITERS: u64 = 1_000_000;
        let sq = VRS_PY_SQ_CURIE;
        let alt = VRS_PY_ALT;
        // Warmup — primes the thread-local buffer.
        for _ in 0..1_000 {
            let _ = vrs2_allele_id(sq, VRS_PY_START, VRS_PY_END, alt);
        }
        let t0 = std::time::Instant::now();
        let mut sink = 0usize;
        for i in 0..ITERS {
            // Vary inputs so the optimizer can't fold the call.
            let id = vrs2_allele_id(sq, VRS_PY_START + i, VRS_PY_END + i, alt);
            sink = sink.wrapping_add(id.len());
        }
        let elapsed = t0.elapsed();
        let per_call_ns = elapsed.as_nanos() as f64 / ITERS as f64;
        let v_per_sec = ITERS as f64 / elapsed.as_secs_f64();
        std::hint::black_box(sink);
        println!(
            "vrs2_allele_id microbench: {ITERS} calls in {elapsed:?} \
             ({per_call_ns:.1} ns/call, {v_per_sec:.0} calls/sec)"
        );
    }

    #[test]
    fn buffer_reuse_across_calls_produces_identical_ids() {
        // Same input twice; the second call reuses the dirty thread-local
        // buffer (post-clear). Catches "did someone forget a buf.clear()".
        let id1 = vrs2_allele_id(VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END, VRS_PY_ALT);
        let id2 = vrs2_allele_id(VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END, VRS_PY_ALT);
        assert_eq!(id1, id2);
        assert_eq!(id1, format!("ga4gh:VA.{VRS_PY_VA_DIGEST_2_0}"));

        // Cross-schema reuse: VRS 1.3 then 2.0 with different inputs,
        // ensuring the buffer is correctly cleared between schemas.
        let id_v1 = vrs1_allele_id(VRS_PY_SQ_CURIE, 100, 101, "A");
        let id_v2 = vrs2_allele_id(VRS_PY_SQ_CURIE, 100, 101, "A");
        assert_ne!(id_v1, id_v2);
        assert!(id_v1.starts_with("ga4gh:VA."));
        assert!(id_v2.starts_with("ga4gh:VA."));
    }
}
