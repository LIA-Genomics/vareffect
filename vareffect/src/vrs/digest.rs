//! GA4GH `sha512t24u` digest construction.
//!
//! Truncated SHA-512 digest with URL-safe base64 encoding, used as the
//! identifier suffix in every `ga4gh:NS.<digest>` CURIE. Defined by the
//! GA4GH digest specification:
//!
//! 1. Compute the SHA-512 of the input bytes.
//! 2. Truncate to the first 24 bytes (192 bits).
//! 3. Encode with URL-safe base64 alphabet (RFC 4648 §5), no padding.
//!
//! 24 bytes is exactly divisible by 6, so the encoding is always 32
//! ASCII characters with no `=` padding.

use base64::Engine;
use base64::engine::general_purpose::URL_SAFE_NO_PAD;
use sha2::{Digest, Sha512};

/// Compute the GA4GH `sha512t24u` digest for an arbitrary byte string.
///
/// Returns the 32-character base64url-encoded suffix (no namespace, no
/// type prefix). The returned string is always exactly 32 ASCII bytes.
///
/// # Examples
///
/// ```ignore
/// use crate::vrs::digest::sha512t24u;
/// assert_eq!(sha512t24u(b""), "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
/// ```
pub(super) fn sha512t24u(blob: &[u8]) -> String {
    let full = Sha512::digest(blob);
    URL_SAFE_NO_PAD.encode(&full[..24])
}

/// Construct a full GA4GH CURIE identifier for the given type prefix
/// and pre-serialized blob.
///
/// `type_prefix` is the GA4GH digest namespace (`"VA"` for Allele,
/// `"VSL"` for VRS 1.3 SequenceLocation, `"SL"` for VRS 2.0
/// SequenceLocation, `"SQ"` for Sequence, etc.).
pub(super) fn ga4gh_identify(type_prefix: &str, blob: &[u8]) -> String {
    format!("ga4gh:{type_prefix}.{}", sha512t24u(blob))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_matches_vrs_python_reference() {
        // Cross-checked against `vrs-python` 0.8 `sha512t24u(b"")`.
        // Empty digest is the canonical regression vector for the
        // sha512t24u algorithm — any drift in either truncation length
        // or alphabet choice breaks this.
        assert_eq!(sha512t24u(b""), "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
    }

    #[test]
    fn output_is_always_32_bytes_no_padding() {
        for input in [b"".as_ref(), b"a", b"hello world", b"\x00\xff\x80\x01"] {
            let out = sha512t24u(input);
            assert_eq!(out.len(), 32, "wrong length for input {input:?}");
            assert!(!out.contains('='), "padding present for input {input:?}");
            assert!(
                out.bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_'),
                "non-url-safe character in output {out}"
            );
        }
    }

    #[test]
    fn ga4gh_identify_assembles_curie_with_namespace() {
        let id = ga4gh_identify("VA", b"");
        assert_eq!(id, "ga4gh:VA.z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
    }
}
