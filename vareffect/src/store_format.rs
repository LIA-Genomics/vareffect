//! On-disk format of the transcript model store (`transcript_models.bin`).
//!
//! A store file is a fixed 12-byte header followed by a MessagePack payload:
//!
//! ```text
//! [magic: b"VEFXTXM\0"][format version: u32 little-endian][payload]
//! ```
//!
//! Since format version 2 the payload is a map with named fields:
//! `source` (a [`TranscriptSource`] label or nil), `source_release` (a
//! string or nil) and `transcripts` (the `Vec<TranscriptModel>`). Version 1
//! stores held the bare model array and are refused like any other version.
//!
//! The header is checked before the ~27 MB payload is decoded, so a store
//! written by an incompatible version fails loudly instead of deserializing
//! into a different meaning. Stores written before the header existed are a
//! bare MessagePack array; their first byte is an array marker, never `b'V'`,
//! so they are recognized exactly as [`StoreFormat::Unversioned`].
//!
//! The public items are re-exported from [`crate::transcript`].

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::error::VarEffectError;
use crate::types::{TranscriptModel, TranscriptSource};

/// Current transcript store format version, written into every store header
/// and required at load time.
///
/// Bump it on any change to:
///
/// - the header or the payload container (its fields, including `source`
///   and `source_release`);
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
/// provenance: the payload names the transcript source and release, and the
/// manifest records the input file. Version 2 added that container, the
/// optional tier, the CDS completeness and genome-mismatch flags, the
/// selenocysteine codons and multi-locus accessions (one chrX and one chrY
/// PAR copy).
pub const TRANSCRIPT_STORE_FORMAT_VERSION: u32 = 2;

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

/// The decoded payload of a store file.
#[derive(Debug, Deserialize)]
pub(crate) struct StoreData {
    /// Annotation source the store was built from, if recorded. The key is
    /// required (nil allowed): the encoder always writes it.
    #[serde(deserialize_with = "Option::deserialize")]
    pub(crate) source: Option<TranscriptSource>,
    /// Release of that source (e.g. `"1.5"`, `"50"`), if recorded. Required
    /// key, nil allowed.
    #[serde(deserialize_with = "Option::deserialize")]
    pub(crate) source_release: Option<String>,
    /// The transcript models, in load order.
    pub(crate) transcripts: Vec<TranscriptModel>,
}

/// Borrowed view of [`StoreData`], so encoding does not copy the models.
/// Field names and order must match [`StoreData`].
#[derive(Serialize)]
struct StoreDataRef<'a> {
    source: Option<TranscriptSource>,
    source_release: Option<&'a str>,
    transcripts: &'a [TranscriptModel],
}

/// Encode transcript models as a versioned store file.
///
/// Writes the header for [`TRANSCRIPT_STORE_FORMAT_VERSION`] followed by the
/// payload as MessagePack with named fields, the form
/// [`crate::TranscriptStore::from_bytes`] reads. Third-party builders use this
/// to write stores from models they built; it is not a way to upgrade an
/// existing store of another version.
///
/// # Arguments
///
/// * `source` - Annotation source of the models, `None` if not one of the
///   [`TranscriptSource`] variants
/// * `source_release` - Release of that source (e.g. `"1.5"`), if known
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
///
/// # Examples
///
/// ```
/// use vareffect::{Exon, Strand, TranscriptModel, TranscriptSource, TranscriptStore};
///
/// let exon = Exon { exon_number: 1, genomic_start: 100, genomic_end: 400 };
/// let model = TranscriptModel::new("NR_000001.1", "GENE1", "chr1", Strand::Plus, vec![exon]);
/// let bytes = vareffect::encode_transcript_store(
///     Some(TranscriptSource::RefSeq),
///     Some("RS_2025_08"),
///     &[model],
/// )?;
/// let store = TranscriptStore::from_bytes(&bytes)?;
/// assert_eq!(store.source(), Some(TranscriptSource::RefSeq));
/// assert_eq!(store.source_release(), Some("RS_2025_08"));
/// # Ok::<(), vareffect::VarEffectError>(())
/// ```
pub fn encode_transcript_store(
    source: Option<TranscriptSource>,
    source_release: Option<&str>,
    transcripts: &[TranscriptModel],
) -> Result<Vec<u8>, VarEffectError> {
    if transcripts.is_empty() {
        return Err(VarEffectError::EmptyStore);
    }
    let mut buffer = Vec::with_capacity(
        TRANSCRIPT_STORE_HEADER_LEN + transcripts.len() * ENCODED_BYTES_PER_RECORD,
    );
    buffer.extend_from_slice(&MAGIC);
    buffer.extend_from_slice(&TRANSCRIPT_STORE_FORMAT_VERSION.to_le_bytes());
    let payload = StoreDataRef {
        source,
        source_release,
        transcripts,
    };
    rmp_serde::encode::write_named(&mut buffer, &payload)?;
    Ok(buffer)
}

/// Decode a versioned store file into its payload.
///
/// Requires the current format version, a payload that is consumed exactly
/// (trailing bytes are an error) and at least one model.
///
/// # Errors
///
/// * [`VarEffectError::UnsupportedStoreFormat`] if the header is missing,
///   truncated or carries another version.
/// * [`VarEffectError::Deserialize`] if the payload is not a valid store
///   payload.
/// * [`VarEffectError::Malformed`] if bytes follow the payload.
/// * [`VarEffectError::EmptyStore`] if the payload holds no models.
pub(crate) fn decode_transcript_store(bytes: &[u8]) -> Result<StoreData, VarEffectError> {
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
    let data = {
        let mut deserializer = rmp_serde::Deserializer::new(&mut remaining);
        StoreData::deserialize(&mut deserializer)?
    };
    if !remaining.is_empty() {
        return Err(VarEffectError::Malformed(format!(
            "{} trailing bytes after the transcript store payload",
            remaining.len()
        )));
    }
    if data.transcripts.is_empty() {
        return Err(VarEffectError::EmptyStore);
    }
    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::plus_strand_coding;
    use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptTier};
    use sha2::{Digest, Sha256};

    /// SHA-256 over [`fingerprint_fixture`] encoded once per source chain
    /// entry at the current version.
    const SCHEMA_FINGERPRINT: &str =
        "fcd9c8c2272541e357f087c8acea5a4121c49aeded2f4e1cfee9a2480a31f058";

    /// The header of a current-version store.
    const HEADER: &[u8; TRANSCRIPT_STORE_HEADER_LEN] = b"VEFXTXM\0\x02\x00\x00\x00";

    fn one_store() -> Vec<TranscriptModel> {
        vec![plus_strand_coding()]
    }

    fn encode(transcripts: &[TranscriptModel]) -> Result<Vec<u8>, VarEffectError> {
        encode_transcript_store(None, None, transcripts)
    }

    /// `HEADER` followed by `payload` encoded as a named MessagePack map.
    fn with_header(payload: &serde_json::Value) -> Vec<u8> {
        let mut bytes = HEADER.to_vec();
        bytes.extend(rmp_serde::to_vec_named(payload).expect("encode payload"));
        bytes
    }

    // Each enum's variants come from a successor chain written as an
    // exhaustive match (no `_` arm): a new variant fails to compile until it
    // gets an arm, which forces a decision about its place in the chain. An
    // arm no other arm points to compiles but stays out of the fingerprint.

    /// Tiers including `None`; the chain ends when it returns to `None`.
    fn next_tier(tier: Option<TranscriptTier>) -> Option<Option<TranscriptTier>> {
        match tier {
            None => Some(Some(TranscriptTier::ManeSelect)),
            Some(TranscriptTier::ManeSelect) => Some(Some(TranscriptTier::ManePlusClinical)),
            Some(TranscriptTier::ManePlusClinical) => Some(Some(TranscriptTier::RefSeqSelect)),
            Some(TranscriptTier::RefSeqSelect) => None,
        }
    }

    fn next_source(source: Option<TranscriptSource>) -> Option<Option<TranscriptSource>> {
        match source {
            None => Some(Some(TranscriptSource::Mane)),
            Some(TranscriptSource::Mane) => Some(Some(TranscriptSource::RefSeq)),
            Some(TranscriptSource::RefSeq) => Some(Some(TranscriptSource::Gencode)),
            Some(TranscriptSource::Gencode) => None,
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

    fn all_tiers() -> Vec<Option<TranscriptTier>> {
        std::iter::successors(Some(None), |tier| next_tier(*tier)).collect()
    }

    fn all_sources() -> Vec<Option<TranscriptSource>> {
        std::iter::successors(Some(None), |source| next_source(*source)).collect()
    }

    fn all_strands() -> Vec<Strand> {
        std::iter::successors(next_strand(None), |strand| next_strand(Some(*strand))).collect()
    }

    fn all_biotypes() -> Vec<Biotype> {
        std::iter::successors(next_biotype(None), |biotype| next_biotype(Some(biotype))).collect()
    }

    /// A record owned by this test alone, as a full struct literal (no
    /// `..base`), so a new field in any serialized type fails to compile here
    /// rather than slipping past the fingerprint. The boolean flags and the
    /// selenocysteine list alternate with `index`.
    fn fingerprint_record(
        index: usize,
        tier: Option<TranscriptTier>,
        strand: Strand,
        biotype: Biotype,
    ) -> TranscriptModel {
        let offset = index as u64 * 10_000;
        let odd = index % 2 == 1;
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
            cds_start_complete: !odd,
            cds_end_complete: odd,
            genome_mismatch: odd,
            cds_alignment_gap: !odd,
            selenocysteine_codons: if odd {
                vec![offset + 1_600, offset + 4_100]
            } else {
                Vec::new()
            },
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

    /// The fixture encoded once per source, alternating a release and none.
    fn fingerprint_stores() -> Vec<Vec<u8>> {
        let fixture = fingerprint_fixture();
        all_sources()
            .into_iter()
            .enumerate()
            .map(|(index, source)| {
                let release = (index % 2 == 0).then_some("1.5");
                encode_transcript_store(source, release, &fixture).expect("encode")
            })
            .collect()
    }

    #[test]
    fn variant_chains_cover_every_variant() {
        assert_eq!(all_tiers().len(), 4);
        assert_eq!(all_sources().len(), 4);
        assert_eq!(all_strands().len(), 2);
        assert_eq!(all_biotypes().len(), 11);
    }

    fn reencode(bytes: &[u8]) -> Vec<u8> {
        let data = decode_transcript_store(bytes).expect("decode");
        encode_transcript_store(
            data.source,
            data.source_release.as_deref(),
            &data.transcripts,
        )
        .expect("encode")
    }

    #[test]
    fn header_bytes_are_pinned() {
        let encoded = encode(&one_store()).expect("encode");
        assert_eq!(&encoded[..TRANSCRIPT_STORE_HEADER_LEN], HEADER);
    }

    #[test]
    fn round_trip_reencodes_identically() {
        for encoded in fingerprint_stores() {
            assert_eq!(reencode(&encoded), encoded);
        }
    }

    #[test]
    fn source_and_release_round_trip() {
        let encoded =
            encode_transcript_store(Some(TranscriptSource::Gencode), Some("50"), &one_store())
                .expect("encode");
        let data = decode_transcript_store(&encoded).expect("decode");
        assert_eq!(data.source, Some(TranscriptSource::Gencode));
        assert_eq!(data.source_release.as_deref(), Some("50"));
        assert_eq!(data.transcripts.len(), 1);
        let bare = decode_transcript_store(&encode(&one_store()).expect("encode")).expect("decode");
        assert_eq!((bare.source, bare.source_release), (None, None));
    }

    #[test]
    fn schema_fingerprint_is_unchanged() {
        let mut hasher = Sha256::new();
        for encoded in fingerprint_stores() {
            hasher.update(&encoded);
        }
        let fingerprint = hex::encode(hasher.finalize());
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
            decode_transcript_store(HEADER),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    #[test]
    fn other_versions_are_rejected() {
        for version in [1u32, 3] {
            let mut encoded = encode(&one_store()).expect("encode");
            encoded[8..12].copy_from_slice(&version.to_le_bytes());
            assert!(matches!(
                decode_transcript_store(&encoded),
                Err(VarEffectError::UnsupportedStoreFormat {
                    found: StoreFormat::Version(found),
                    ..
                }) if found == version
            ));
        }
    }

    #[test]
    fn version_one_payload_behind_the_current_header_is_a_deserialize_error() {
        let mut bytes = HEADER.to_vec();
        bytes.extend(rmp_serde::to_vec_named(&one_store()).expect("encode"));
        assert!(matches!(
            decode_transcript_store(&bytes),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    #[test]
    fn truncated_header_is_named() {
        for length in [1, 5, 8, 11] {
            assert_eq!(
                store_format(&HEADER[..length]),
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
        let mut bytes = HEADER.to_vec();
        bytes.extend_from_slice(b"\xc1not messagepack");
        assert!(matches!(
            decode_transcript_store(&bytes),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let mut encoded = encode(&one_store()).expect("encode");
        encoded.extend_from_slice(b"extra");
        assert!(matches!(
            decode_transcript_store(&encoded),
            Err(VarEffectError::Malformed(message)) if message.contains("5 trailing bytes")
        ));
    }

    #[test]
    fn empty_stores_are_rejected_on_both_sides() {
        assert!(matches!(encode(&[]), Err(VarEffectError::EmptyStore)));
        let empty = with_header(&serde_json::json!({
            "source": null,
            "source_release": null,
            "transcripts": [],
        }));
        assert!(matches!(
            decode_transcript_store(&empty),
            Err(VarEffectError::EmptyStore)
        ));
    }

    #[test]
    fn compact_payload_behind_a_valid_header_decodes() {
        let mut bytes = HEADER.to_vec();
        let payload = StoreDataRef {
            source: Some(TranscriptSource::Mane),
            source_release: Some("1.5"),
            transcripts: &one_store(),
        };
        bytes.extend(rmp_serde::to_vec(&payload).expect("compact"));
        let decoded = decode_transcript_store(&bytes).expect("compact payload");
        assert_eq!(decoded.source, Some(TranscriptSource::Mane));
        assert_eq!(decoded.transcripts[0].accession, one_store()[0].accession);
    }

    #[test]
    fn pre_header_readers_fail_on_versioned_stores() {
        let encoded = encode(&one_store()).expect("encode");
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
        let bytes = with_header(&serde_json::json!({
            "source": "mane",
            "source_release": null,
            "transcripts": [record],
            "container_field_from_a_newer_version": true,
        }));
        let decoded = decode_transcript_store(&bytes).expect("silently accepted");
        assert_eq!(decoded.transcripts[0].protein_accession, None);
    }

    /// The container's `source`, `source_release` and the model's `tier`
    /// keys are required even though their values may be nil, so a
    /// truncated or foreign payload cannot load as "source unknown" or
    /// "no tier".
    #[test]
    fn missing_container_or_tier_keys_are_deserialize_errors() {
        let record = serde_json::to_value(plus_strand_coding()).expect("to json");
        let full = serde_json::json!({
            "source": null,
            "source_release": null,
            "transcripts": [record.clone()],
        });
        assert!(decode_transcript_store(&with_header(&full)).is_ok());
        for key in ["source", "source_release"] {
            let mut payload = full.clone();
            payload.as_object_mut().expect("map").remove(key);
            assert!(
                matches!(
                    decode_transcript_store(&with_header(&payload)),
                    Err(VarEffectError::Deserialize(_))
                ),
                "{key}"
            );
        }
        let mut no_tier = record;
        no_tier.as_object_mut().expect("map").remove("tier");
        let payload = serde_json::json!({
            "source": null,
            "source_release": null,
            "transcripts": [no_tier],
        });
        assert!(matches!(
            decode_transcript_store(&with_header(&payload)),
            Err(VarEffectError::Deserialize(_))
        ));
    }

    /// A missing required field (here a v2 flag) is an error, not a default.
    #[test]
    fn missing_required_field_is_a_deserialize_error() {
        let mut record = serde_json::to_value(plus_strand_coding()).expect("to json");
        let fields = record.as_object_mut().expect("struct is a map");
        fields.remove("cds_start_complete");
        let bytes = with_header(&serde_json::json!({
            "source": null,
            "source_release": null,
            "transcripts": [record],
        }));
        assert!(matches!(
            decode_transcript_store(&bytes),
            Err(VarEffectError::Deserialize(_))
        ));
    }
}
