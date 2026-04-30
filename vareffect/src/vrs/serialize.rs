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
//! not apply. `serde_json::Map` in default configuration is backed by
//! `BTreeMap`, which iterates in sorted key order, giving us canonical
//! key ordering for free.

use serde_json::{Value, json};

use super::digest::{ga4gh_identify, sha512t24u};
use super::sq_digest::bare_digest_from_curie;

/// Build the VRS 1.3 SequenceLocation canonical serialization blob.
///
/// `sq_curie` is the full `ga4gh:SQ.<digest>` form; the canonical blob
/// uses only the bare 32-char digest portion in the `sequence_id`
/// field (per the CURIE-reduction rule).
fn vrs1_location_blob(sq_curie: &str, start: u64, end: u64) -> Vec<u8> {
    let bare_sq = bare_digest_from_curie(sq_curie);
    let val: Value = json!({
        "interval": {
            "end":   {"type": "Number", "value": end},
            "start": {"type": "Number", "value": start},
            "type":  "SequenceInterval",
        },
        "sequence_id": bare_sq,
        "type": "SequenceLocation",
    });
    serde_json::to_vec(&val).expect("VRS 1.3 location JSON serialization")
}

/// Build the VRS 1.3 Allele canonical serialization blob.
///
/// `location_digest` is the bare 32-char `sha512t24u` digest of the
/// SequenceLocation, computed by the caller via the recursive enref
/// step. The Allele's `location` field substitutes that bare digest
/// string in place of the inlined object.
fn vrs1_allele_blob(location_digest: &str, alt: &str) -> Vec<u8> {
    let val: Value = json!({
        "location": location_digest,
        "state": {
            "sequence": alt,
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    });
    serde_json::to_vec(&val).expect("VRS 1.3 allele JSON serialization")
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
    let location_blob = vrs1_location_blob(sq_curie, start, end);
    let location_digest = sha512t24u(&location_blob);
    let allele_blob = vrs1_allele_blob(&location_digest, alt);
    ga4gh_identify("VA", &allele_blob)
}

/// Build the VRS 2.0 SequenceLocation canonical serialization blob.
///
/// Differs from 1.3 in: no `SequenceInterval` wrapper, `start`/`end`
/// are direct integers on the location, and `sequence_id` is replaced
/// by a `sequenceReference` object carrying `refgetAccession` (which
/// uses bare `SQ.<digest>` form per refget v2, no `ga4gh:` namespace).
fn vrs2_location_blob(sq_curie: &str, start: u64, end: u64) -> Vec<u8> {
    let bare_sq = bare_digest_from_curie(sq_curie);
    // refget convention: `refgetAccession` is `SQ.<digest>` (no `ga4gh:` prefix).
    let refget_acc = format!("SQ.{bare_sq}");
    let val: Value = json!({
        "end": end,
        "sequenceReference": {
            "refgetAccession": refget_acc,
            "type": "SequenceReference",
        },
        "start": start,
        "type": "SequenceLocation",
    });
    serde_json::to_vec(&val).expect("VRS 2.0 location JSON serialization")
}

/// Build the VRS 2.0 Allele canonical serialization blob.
///
/// Same Allele shape as 1.3 (location → bare digest string, state,
/// type), but the location's pre-digest blob differs (see
/// [`vrs2_location_blob`]).
fn vrs2_allele_blob(location_digest: &str, alt: &str) -> Vec<u8> {
    let val: Value = json!({
        "location": location_digest,
        "state": {
            "sequence": alt,
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    });
    serde_json::to_vec(&val).expect("VRS 2.0 allele JSON serialization")
}

/// Compute the full VRS 2.0 Allele identifier (`ga4gh:VA.<digest>`).
pub(super) fn vrs2_allele_id(sq_curie: &str, start: u64, end: u64, alt: &str) -> String {
    let location_blob = vrs2_location_blob(sq_curie, start, end);
    let location_digest = sha512t24u(&location_blob);
    let allele_blob = vrs2_allele_blob(&location_digest, alt);
    ga4gh_identify("VA", &allele_blob)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vrs1_location_blob_is_canonical_sorted_compact_json() {
        // Manually-constructed expected blob with keys in sorted order
        // and no whitespace. Any drift in field order or spacing
        // breaks the digest contract.
        let blob = vrs1_location_blob(
            "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            44908821,
            44908822,
        );
        let expected = r#"{"interval":{"end":{"type":"Number","value":44908822},"start":{"type":"Number","value":44908821},"type":"SequenceInterval"},"sequence_id":"IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceLocation"}"#;
        assert_eq!(std::str::from_utf8(&blob).unwrap(), expected);
    }

    #[test]
    fn vrs1_allele_blob_substitutes_location_as_bare_digest_string() {
        // The recursive enref step: `location` must be a bare digest
        // STRING in the parent's blob, not an inlined object.
        let blob = vrs1_allele_blob("esDSArZQC-Sx-96ZZzHnzAVNOc439oE5", "T");
        let expected = r#"{"location":"esDSArZQC-Sx-96ZZzHnzAVNOc439oE5","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#;
        assert_eq!(std::str::from_utf8(&blob).unwrap(), expected);
    }

    #[test]
    fn vrs2_location_blob_uses_sequence_reference_object() {
        // 2.0 differs from 1.3: no SequenceInterval wrapper, and
        // `sequenceReference` is an object not a `sequence_id` string.
        let blob = vrs2_location_blob(
            "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
            44908821,
            44908822,
        );
        let expected = r#"{"end":44908822,"sequenceReference":{"refgetAccession":"SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl","type":"SequenceReference"},"start":44908821,"type":"SequenceLocation"}"#;
        assert_eq!(std::str::from_utf8(&blob).unwrap(), expected);
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

    // -----------------------------------------------------------------
    // Cross-implementation correctness (vrs-python golden vectors)
    // -----------------------------------------------------------------
    //
    // Inputs and expected bytes/digests below are copied verbatim from
    // vrs-python's `tests/test_vrs.py::test_vr` (commit on `main` as of
    // 2026-04-30):
    //
    //   start=55181319, end=55181320, refget=SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul
    //   alt="T"  →  ga4gh:VA.Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE
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
        let blob = vrs2_location_blob(VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END);
        let expected = br#"{"end":55181320,"sequenceReference":{"refgetAccession":"SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul","type":"SequenceReference"},"start":55181319,"type":"SequenceLocation"}"#;
        assert_eq!(&blob, expected);
    }

    #[test]
    fn vrs2_location_digest_matches_vrs_python() {
        let blob = vrs2_location_blob(VRS_PY_SQ_CURIE, VRS_PY_START, VRS_PY_END);
        assert_eq!(sha512t24u(&blob), VRS_PY_LOC_DIGEST_2_0);
    }

    #[test]
    fn vrs2_allele_blob_matches_vrs_python_bytes() {
        let blob = vrs2_allele_blob(VRS_PY_LOC_DIGEST_2_0, VRS_PY_ALT);
        let expected = br#"{"location":"_G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd","state":{"sequence":"T","type":"LiteralSequenceExpression"},"type":"Allele"}"#;
        assert_eq!(&blob, expected);
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
        let blob = vrs2_location_blob(VRS_PY_SQ_CURIE, 44908821, 44908822);
        assert_eq!(sha512t24u(&blob), "4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT");
    }
}
