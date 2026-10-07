//! On-disk format of the transcript model store (`transcript_models.bin`).
//!
//! A store file is a fixed 12-byte header followed by a MessagePack payload:
//!
//! ```text
//! [magic: b"VEFXTXM\0"][format version: u32 little-endian][Vec<TranscriptModel>]
//! ```
//!
//! The header is checked before the ~27 MB payload is decoded, so a store
//! written by an incompatible version fails loudly instead of deserializing
//! into a different meaning. Stores written before the header existed are a
//! bare MessagePack array; their first byte is an array marker, never `b'V'`,
//! so they are recognized exactly as [`StoreFormat::Unversioned`].
//!
//! The public items are re-exported from [`crate::transcript`].

use std::fmt;

use serde::Deserialize;

use crate::error::VarEffectError;
use crate::types::TranscriptModel;

/// Current transcript store format version, written into every store header
/// and required at load time.
///
/// Bump it on any change to:
///
/// - the header or the payload container;
/// - the fields, variants or encoding of any serialized type
///   ([`TranscriptModel`] and everything it contains), **including additive
///   `Option` fields and new enum variants**: serde ignores unknown fields and
///   defaults a missing `Option` to `None`, so an additive change decodes
///   silently in an older reader with the new information dropped;
/// - the meaning of any field or variant (accession namespace, coordinate or
///   chromosome convention, tier semantics);
/// - any store-level invariant readers rely on (accession uniqueness,
///   ordering, admission policy).
///
/// The schema-fingerprint test in this module catches shape changes; changes
/// of meaning alone are not caught mechanically and must follow this list.
/// The version says whether a reader can decode a store. It is not
/// provenance: which transcript source or release a store holds is recorded
/// in its manifest.
pub const TRANSCRIPT_STORE_FORMAT_VERSION: u32 = 1;

/// Length in bytes of the store header that precedes the MessagePack payload.
pub const TRANSCRIPT_STORE_HEADER_LEN: usize = 12;

/// Magic bytes opening every versioned store.
const MAGIC: [u8; 8] = *b"VEFXTXM\0";

const _: () = assert!(MAGIC.len() + size_of::<u32>() == TRANSCRIPT_STORE_HEADER_LEN);

/// Bytes reserved per record when sizing the encode buffer (the MANE 1.5
/// store averages about 1,380 bytes per transcript).
const ENCODED_BYTES_PER_RECORD: usize = 1_500;

/// What the opening bytes of a transcript store file say about its format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StoreFormat {
    /// A versioned store with this format version.
    Version(u32),
    /// A store written before the header existed: a bare MessagePack array.
    Unversioned,
    /// The file opens with (part of) the magic but ends before the header does.
    TruncatedHeader,
    /// Empty, or neither the magic nor a MessagePack array.
    Unrecognized,
}

impl fmt::Display for StoreFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Version(version) => write!(f, "format version {version}"),
            Self::Unversioned => f.write_str("an unversioned store from before format versioning"),
            Self::TruncatedHeader => f.write_str("a truncated store header"),
            Self::Unrecognized => f.write_str("a file that is not a transcript store"),
        }
    }
}

/// Classify the opening bytes of a transcript store file.
///
/// Only the first [`TRANSCRIPT_STORE_HEADER_LEN`] bytes are inspected, so a
/// caller can pass just the head of a large file.
///
/// # Arguments
///
/// * `bytes` - The file contents, or at least its first 12 bytes
///
/// # Returns
///
/// The detected [`StoreFormat`]; never fails.
///
/// # Examples
///
/// ```
/// use vareffect::transcript::{StoreFormat, store_format};
///
/// assert_eq!(store_format(b"VEFXTXM\0\x01\x00\x00\x00"), StoreFormat::Version(1));
/// assert_eq!(store_format(&[0x91]), StoreFormat::Unversioned);
/// assert_eq!(store_format(b""), StoreFormat::Unrecognized);
/// ```
pub fn store_format(bytes: &[u8]) -> StoreFormat {
    if bytes.len() >= TRANSCRIPT_STORE_HEADER_LEN && bytes[..MAGIC.len()] == MAGIC {
        let mut version = [0u8; 4];
        version.copy_from_slice(&bytes[MAGIC.len()..TRANSCRIPT_STORE_HEADER_LEN]);
        return StoreFormat::Version(u32::from_le_bytes(version));
    }
    let opening = &bytes[..bytes.len().min(MAGIC.len())];
    if !bytes.is_empty() && MAGIC.starts_with(opening) {
        return StoreFormat::TruncatedHeader;
    }
    match bytes.first() {
        // MessagePack fixarray, array16 and array32 markers.
        Some(0x90..=0x9f | 0xdc | 0xdd) => StoreFormat::Unversioned,
        Some(_) | None => StoreFormat::Unrecognized,
    }
}

/// Encode transcript models as a versioned store file.
///
/// Writes the header for [`TRANSCRIPT_STORE_FORMAT_VERSION`] followed by the
/// models as MessagePack with named fields, the form
/// [`crate::TranscriptStore::from_bytes`] reads. Third-party builders use this
/// to write stores from models they built; it is not a way to upgrade an
/// existing store of another version.
///
/// # Arguments
///
/// * `transcripts` - The models to store, in the order they should be loaded
///
/// # Returns
///
/// The complete file contents.
///
/// # Errors
///
/// * [`VarEffectError::EmptyStore`] if `transcripts` is empty.
/// * [`VarEffectError::Serialize`] if MessagePack encoding fails.
pub fn encode_transcript_store(transcripts: &[TranscriptModel]) -> Result<Vec<u8>, VarEffectError> {
    if transcripts.is_empty() {
        return Err(VarEffectError::EmptyStore);
    }
    let mut buffer = Vec::with_capacity(
        TRANSCRIPT_STORE_HEADER_LEN + transcripts.len() * ENCODED_BYTES_PER_RECORD,
    );
    buffer.extend_from_slice(&MAGIC);
    buffer.extend_from_slice(&TRANSCRIPT_STORE_FORMAT_VERSION.to_le_bytes());
    rmp_serde::encode::write_named(&mut buffer, transcripts)?;
    Ok(buffer)
}

/// Decode a versioned store file into its transcript models.
///
/// Requires the current format version, a payload that is consumed exactly
/// (trailing bytes are an error) and at least one model.
///
/// # Errors
///
/// * [`VarEffectError::UnsupportedStoreFormat`] if the header is missing,
///   truncated or carries another version.
/// * [`VarEffectError::Deserialize`] if the payload is not a valid model list.
/// * [`VarEffectError::Malformed`] if bytes follow the payload.
/// * [`VarEffectError::EmptyStore`] if the payload holds no models.
pub(crate) fn decode_transcript_store(
    bytes: &[u8],
) -> Result<Vec<TranscriptModel>, VarEffectError> {
    let found = store_format(bytes);
    if found != StoreFormat::Version(TRANSCRIPT_STORE_FORMAT_VERSION) {
        return Err(VarEffectError::UnsupportedStoreFormat {
            found,
            expected: TRANSCRIPT_STORE_FORMAT_VERSION,
        });
    }
    // A `&mut &[u8]` reader advances past exactly what the decoder consumed,
    // so whatever remains afterwards is trailing data.
    let mut remaining = &bytes[TRANSCRIPT_STORE_HEADER_LEN..];
    let transcripts = {
        let mut deserializer = rmp_serde::Deserializer::new(&mut remaining);
        Vec::<TranscriptModel>::deserialize(&mut deserializer)?
    };
    if !remaining.is_empty() {
        return Err(VarEffectError::Malformed(format!(
            "{} trailing bytes after the transcript store payload",
            remaining.len()
        )));
    }
    if transcripts.is_empty() {
        return Err(VarEffectError::EmptyStore);
    }
    Ok(transcripts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::plus_strand_coding;
    use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptTier};
    use sha2::{Digest, Sha256};

    /// SHA-256 of [`fingerprint_fixture`] encoded at the current version.
    const SCHEMA_FINGERPRINT: &str =
        "51aa8fa35034c03a8898ce01e47c6beba53930064f9a7538ca12eb0fcda58e72";

    fn one_store() -> Vec<TranscriptModel> {
        vec![plus_strand_coding()]
    }

    // Each enum's variants come from a successor chain written as an
    // exhaustive match (no `_` arm): a new variant fails to compile until it
    // gets an arm, which forces a decision about its place in the chain. An
    // arm no other arm points to compiles but stays out of the fingerprint.

    fn next_tier(tier: Option<TranscriptTier>) -> Option<TranscriptTier> {
        match tier {
            None => Some(TranscriptTier::ManeSelect),
            Some(TranscriptTier::ManeSelect) => Some(TranscriptTier::ManePlusClinical),
            Some(TranscriptTier::ManePlusClinical) => Some(TranscriptTier::RefSeqSelect),
            Some(TranscriptTier::RefSeqSelect) => None,
        }
    }

    fn next_strand(strand: Option<Strand>) -> Option<Strand> {
        match strand {
            None => Some(Strand::Plus),
            Some(Strand::Plus) => Some(Strand::Minus),
            Some(Strand::Minus) => None,
        }
    }

    fn next_biotype(biotype: Option<&Biotype>) -> Option<Biotype> {
        match biotype {
            None => Some(Biotype::ProteinCoding),
            Some(Biotype::ProteinCoding) => Some(Biotype::NonCodingRna),
            Some(Biotype::NonCodingRna) => Some(Biotype::LncRna),
            Some(Biotype::LncRna) => Some(Biotype::AntisenseRna),
            Some(Biotype::AntisenseRna) => Some(Biotype::SnoRna),
            Some(Biotype::SnoRna) => Some(Biotype::SnRna),
            Some(Biotype::SnRna) => Some(Biotype::RnaseMrpRna),
            Some(Biotype::RnaseMrpRna) => Some(Biotype::TelomeraseRna),
            Some(Biotype::TelomeraseRna) => Some(Biotype::VaultRna),
            Some(Biotype::VaultRna) => Some(Biotype::Other("misc_RNA".into())),
            Some(Biotype::Other(_)) => Some(Biotype::Unknown),
            Some(Biotype::Unknown) => None,
        }
    }

    fn all_tiers() -> Vec<TranscriptTier> {
        std::iter::successors(next_tier(None), |tier| next_tier(Some(*tier))).collect()
    }

    fn all_strands() -> Vec<Strand> {
        std::iter::successors(next_strand(None), |strand| next_strand(Some(*strand))).collect()
    }

    fn all_biotypes() -> Vec<Biotype> {
        std::iter::successors(next_biotype(None), |biotype| next_biotype(Some(biotype))).collect()
    }

    /// A record owned by this test alone, as a full struct literal (no
    /// `..base`), so a new field in any serialized type fails to compile here
    /// rather than slipping past the fingerprint.
    fn fingerprint_record(
        index: usize,
        tier: TranscriptTier,
        strand: Strand,
        biotype: Biotype,
    ) -> TranscriptModel {
        let offset = index as u64 * 10_000;
        TranscriptModel {
            accession: format!("NM_{:06}.1", index + 1),
            protein_accession: Some(format!("NP_{:06}.1", index + 1)),
            gene_symbol: format!("GENE{index}"),
            hgnc_id: Some(format!("HGNC:{}", index + 1)),
            ensembl_accession: Some(format!("ENST{:011}.1", index + 1)),
            chrom: "chr1".into(),
            strand,
            tx_start: offset + 1_000,
            tx_end: offset + 5_000,
            cds_genomic_start: Some(offset + 1_500),
            cds_genomic_end: Some(offset + 4_500),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: offset + 1_000,
                    genomic_end: offset + 2_000,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: offset + 4_000,
                    genomic_end: offset + 5_000,
                },
            ],
            cds_segments: vec![
                CdsSegment {
                    exon_index: 0,
                    genomic_start: offset + 1_500,
                    genomic_end: offset + 2_000,
                    phase: 0,
                },
                CdsSegment {
                    exon_index: 1,
                    genomic_start: offset + 4_000,
                    genomic_end: offset + 4_500,
                    phase: 1,
                },
            ],
            tier,
            biotype,
            exon_count: 2,
        }
    }

    /// One record per enum variant, cycling the shorter lists.
    fn fingerprint_fixture() -> Vec<TranscriptModel> {
        let tiers = all_tiers();
        let strands = all_strands();
        let biotypes = all_biotypes();
        let count = tiers.len().max(strands.len()).max(biotypes.len());
        (0..count)
            .map(|index| {
                fingerprint_record(
                    index,
                    tiers[index % tiers.len()],
                    strands[index % strands.len()],
                    biotypes[index % biotypes.len()].clone(),
                )
            })
            .collect()
    }

    #[test]
    fn variant_chains_cover_every_variant() {
        assert_eq!(all_tiers().len(), 3);
        assert_eq!(all_strands().len(), 2);
        assert_eq!(all_biotypes().len(), 11);
    }

    fn reencode(bytes: &[u8]) -> Vec<u8> {
        encode_transcript_store(&decode_transcript_store(bytes).expect("decode")).expect("encode")
    }

    #[test]
    fn header_bytes_are_pinned() {
        let encoded = encode_transcript_store(&one_store()).expect("encode");
        assert_eq!(
            &encoded[..TRANSCRIPT_STORE_HEADER_LEN],
            b"VEFXTXM\0\x01\x00\x00\x00"
        );
    }

    #[test]
    fn round_trip_reencodes_identically() {
        let encoded = encode_transcript_store(&fingerprint_fixture()).expect("encode");
        assert_eq!(reencode(&encoded), encoded);
    }

    #[test]
    fn schema_fingerprint_is_unchanged() {
        let encoded = encode_transcript_store(&fingerprint_fixture()).expect("encode");
        let fingerprint = hex::encode(Sha256::digest(&encoded));
        assert_eq!(
            fingerprint, SCHEMA_FINGERPRINT,
            "the encoded store bytes changed. If a serialized type changed (a field, a \
             variant or its encoding), bump TRANSCRIPT_STORE_FORMAT_VERSION and update \
             SCHEMA_FINGERPRINT. If only this fixture or the rmp-serde version changed, \
             update SCHEMA_FINGERPRINT alone after checking the payload still decodes in \
             the previous release"
        );
    }

    #[test]
    fn legacy_named_and_compact_arrays_are_unversioned() {
        let named = rmp_serde::to_vec_named(&one_store()).expect("named");
        let compact = rmp_serde::to_vec(&one_store()).expect("compact");
        for legacy in [named, compact] {
            assert_eq!(store_format(&legacy), StoreFormat::Unversioned);
            assert!(matches!(
                decode_transcript_store(&legacy),
                Err(VarEffectError::UnsupportedStoreFormat {
                    found: StoreFormat::Unversioned,
                    expected: TRANSCRIPT_STORE_FORMAT_VERSION,
                })
            ));
        }
    }

    #[test]
    fn every_legacy_array_marker_is_unversioned() {
        // A real legacy store opens with array16 (MANE holds 19,437 records).
        assert_eq!(store_format(&[0xdc, 0x4b, 0xed]), StoreFormat::Unversioned);
        assert_eq!(
            store_format(&[0xdd, 0, 0, 0x4b, 0xed]),
            StoreFormat::Unversioned
        );
        let sixteen: Vec<TranscriptModel> = (0..16).map(|_| plus_strand_coding()).collect();
        let legacy = rmp_serde::to_vec_named(&sixteen).expect("encode");
        assert_eq!(legacy[0], 0xdc);
        assert_eq!(store_format(&legacy), StoreFormat::Unversioned);
    }

    #[test]
    fn header_without_payload_is_a_deserialize_error() {
        assert!(matches!(
            decode_transcript_store(b"VEFXTXM\0\x01\x00\x00\x00"),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    #[test]
    fn other_versions_are_rejected() {
        let mut encoded = encode_transcript_store(&one_store()).expect("encode");
        encoded[8..12].copy_from_slice(&2u32.to_le_bytes());
        assert!(matches!(
            decode_transcript_store(&encoded),
            Err(VarEffectError::UnsupportedStoreFormat {
                found: StoreFormat::Version(2),
                ..
            })
        ));
    }

    #[test]
    fn truncated_header_is_named() {
        for length in [1, 5, 8, 11] {
            let header = &b"VEFXTXM\0\x01\x00\x00\x00"[..length];
            assert_eq!(
                store_format(header),
                StoreFormat::TruncatedHeader,
                "{length} bytes"
            );
        }
    }

    #[test]
    fn empty_and_foreign_bytes_are_unrecognized() {
        assert_eq!(store_format(b""), StoreFormat::Unrecognized);
        assert_eq!(store_format(b"genome"), StoreFormat::Unrecognized);
        assert_eq!(store_format(&[0x00, 0x01, 0x02]), StoreFormat::Unrecognized);
        assert!(matches!(
            decode_transcript_store(b""),
            Err(VarEffectError::UnsupportedStoreFormat {
                found: StoreFormat::Unrecognized,
                ..
            })
        ));
    }

    #[test]
    fn garbage_payload_is_a_deserialize_error() {
        let mut bytes = b"VEFXTXM\0\x01\x00\x00\x00".to_vec();
        bytes.extend_from_slice(b"\xc1not messagepack");
        assert!(matches!(
            decode_transcript_store(&bytes),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let mut encoded = encode_transcript_store(&one_store()).expect("encode");
        encoded.extend_from_slice(b"extra");
        assert!(matches!(
            decode_transcript_store(&encoded),
            Err(VarEffectError::Malformed(message)) if message.contains("5 trailing bytes")
        ));
    }

    #[test]
    fn empty_stores_are_rejected_on_both_sides() {
        assert!(matches!(
            encode_transcript_store(&[]),
            Err(VarEffectError::EmptyStore)
        ));
        let mut empty = b"VEFXTXM\0\x01\x00\x00\x00".to_vec();
        empty.push(0x90);
        assert!(matches!(
            decode_transcript_store(&empty),
            Err(VarEffectError::EmptyStore)
        ));
    }

    #[test]
    fn compact_payload_behind_a_valid_header_decodes() {
        let mut bytes = b"VEFXTXM\0\x01\x00\x00\x00".to_vec();
        bytes.extend(rmp_serde::to_vec(&one_store()).expect("compact"));
        let decoded = decode_transcript_store(&bytes).expect("compact payload");
        assert_eq!(decoded[0].accession, one_store()[0].accession);
    }

    #[test]
    fn pre_header_readers_fail_on_versioned_stores() {
        let encoded = encode_transcript_store(&one_store()).expect("encode");
        assert!(rmp_serde::from_slice::<Vec<TranscriptModel>>(&encoded).is_err());
    }

    /// Pins the serde behavior the bump rule depends on: an added field is
    /// ignored and a missing `Option` field decodes as `None`, both silently.
    #[test]
    fn additive_changes_decode_silently() {
        let mut record = serde_json::to_value(plus_strand_coding()).expect("to json");
        let fields = record.as_object_mut().expect("struct is a map");
        fields.remove("protein_accession");
        fields.insert("field_from_a_newer_version".into(), serde_json::json!(1));
        let mut bytes = b"VEFXTXM\0\x01\x00\x00\x00".to_vec();
        bytes.extend(rmp_serde::to_vec_named(&vec![record]).expect("encode map"));
        let decoded = decode_transcript_store(&bytes).expect("silently accepted");
        assert_eq!(decoded[0].protein_accession, None);
    }
}
