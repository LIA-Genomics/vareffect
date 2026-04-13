//! In-memory transcript store with per-chromosome interval tree indexing.
//!
//! The store is built from a MessagePack-serialized `Vec<TranscriptModel>`
//! produced offline by `vareffect-cli`. At load time it constructs:
//!
//! 1. A shared `Arc<[TranscriptModel]>` slice for O(1) index lookup and cheap
//!    cross-thread sharing.
//! 2. A `HashMap<chrom, COITree>` — one interval tree per chromosome, whose
//!    payload is the index of the corresponding record in the slice.
//! 3. A `HashMap<accession, index>` for O(1) direct lookup by accession.
//!
//! Query coordinates are **0-based half-open** `[start, end)`. Internally,
//! coitrees uses end-inclusive `[first, last]`, so we convert both at
//! construction and at query time.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use coitrees::{COITree, IntervalNode, IntervalTree};

use crate::error::VarEffectError;
use crate::locate::LocateIndex;
use crate::types::TranscriptModel;

/// In-memory transcript model store indexed by genomic interval and by
/// accession.
///
/// Construct via [`TranscriptStore::load_from_path`] (production) or
/// [`TranscriptStore::from_transcripts`] (tests and build-time callers).
///
/// `Debug` is implemented manually (the `COITree` field does not derive
/// `Debug`); it prints a compact summary (`len`, chromosome count, first few
/// accessions) rather than dumping the full interval tree.
pub struct TranscriptStore {
    /// Backing slice. `Arc` lets multiple threads share the data without
    /// copying; `[T]` (not `Vec<T>`) because we never mutate after
    /// construction.
    transcripts: Arc<[TranscriptModel]>,
    /// Precomputed locate indices, one per transcript, same indexing as
    /// `transcripts`. Built at construction time for O(log n) exon lookup
    /// and O(1) CDS offset computation at query time.
    locate_indices: Arc<[LocateIndex]>,
    /// chrom → coitrees interval tree. The node metadata (`usize`) is the
    /// index into `transcripts`. The second type parameter (`u32`) is the
    /// internal tree index type; `u32` is smaller than the default `usize`
    /// and is fine because no chromosome has > 4B transcripts.
    trees: HashMap<String, COITree<usize, u32>>,
    /// Full accession (with version, e.g. `"NM_006772.2"`) → index into
    /// `transcripts`. Version-less lookup is intentionally *not* provided
    /// here; resolving an unversioned HGVS input to "the latest version"
    /// belongs with the caller, not the raw store.
    by_accession: HashMap<String, usize>,
}

impl TranscriptStore {
    /// Load and index a MessagePack-serialized `Vec<TranscriptModel>` from
    /// disk.
    ///
    /// # Errors
    ///
    /// * [`VarEffectError::Io`] if the file cannot be read.
    /// * [`VarEffectError::Deserialize`] if the MessagePack payload is
    ///   malformed or the wrong type.
    /// * [`VarEffectError::Malformed`] if any record violates a coordinate
    ///   invariant (`tx_end <= tx_start` or `tx_end > i32::MAX`).
    pub fn load_from_path(path: &Path) -> Result<Self, VarEffectError> {
        let bytes = std::fs::read(path).map_err(|source| VarEffectError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let transcripts: Vec<TranscriptModel> = rmp_serde::from_slice(&bytes)?;
        Self::try_from_transcripts(transcripts)
    }

    /// Build a store from an owned `Vec<TranscriptModel>`.
    ///
    /// Panics on invariant violations — use [`TranscriptStore::try_from_transcripts`]
    /// in contexts where malformed input is possible.
    pub fn from_transcripts(transcripts: Vec<TranscriptModel>) -> Self {
        Self::try_from_transcripts(transcripts)
            .expect("TranscriptStore::from_transcripts: invariant violation in input")
    }

    /// Fallible version of [`TranscriptStore::from_transcripts`].
    ///
    /// Used by the MessagePack loader so corrupt on-disk data surfaces as an
    /// error rather than a panic.
    pub fn try_from_transcripts(transcripts: Vec<TranscriptModel>) -> Result<Self, VarEffectError> {
        // Per-chromosome scratch buckets so we can `COITree::new(&nodes)` once
        // per chromosome without re-sorting the global vec.
        let mut scratch: HashMap<String, Vec<IntervalNode<usize, u32>>> = HashMap::new();
        let mut by_accession: HashMap<String, usize> = HashMap::with_capacity(transcripts.len());

        for (idx, tx) in transcripts.iter().enumerate() {
            // Invariant: tx_start < tx_end. Rejecting `==` catches zero-length
            // transcripts, which would also confuse the `end - 1` conversion.
            if tx.tx_end <= tx.tx_start {
                return Err(VarEffectError::Malformed(format!(
                    "transcript {} has tx_start={} >= tx_end={}",
                    tx.accession, tx.tx_start, tx.tx_end
                )));
            }

            // coitrees stores coordinates as i32. Human chromosome 1 maxes out
            // at ~249M, well under i32::MAX (~2.1B), but we assert symmetry on
            // both endpoints so a pathological GFF3 input (e.g., a patch
            // sequence with unusual coordinates) fails loudly instead of
            // silently wrapping during the `as i32` cast. `vareffect-cli` already
            // filters obvious garbage; this is a defense-in-depth check.
            if tx.tx_start > i32::MAX as u64 {
                return Err(VarEffectError::Malformed(format!(
                    "transcript {} has tx_start={} exceeding i32::MAX",
                    tx.accession, tx.tx_start
                )));
            }
            if tx.tx_end > i32::MAX as u64 {
                return Err(VarEffectError::Malformed(format!(
                    "transcript {} has tx_end={} exceeding i32::MAX",
                    tx.accession, tx.tx_end
                )));
            }

            // Insert a duplicate-accession guard. MANE is supposed to have
            // unique (accession, version) tuples, but defending against dup
            // keys means the builder cannot silently drop records.
            if by_accession.insert(tx.accession.clone(), idx).is_some() {
                return Err(VarEffectError::Malformed(format!(
                    "duplicate accession {} in transcript store",
                    tx.accession
                )));
            }

            // Convert half-open [start, end) to end-inclusive [start, end-1]
            // for coitrees.
            let start_i32 = tx.tx_start as i32;
            let end_i32 = (tx.tx_end - 1) as i32;
            let node = IntervalNode::new(start_i32, end_i32, idx);
            scratch.entry(tx.chrom.clone()).or_default().push(node);
        }

        let mut trees: HashMap<String, COITree<usize, u32>> = HashMap::with_capacity(scratch.len());
        for (chrom, nodes) in scratch {
            let tree = COITree::new(&nodes);
            trees.insert(chrom, tree);
        }

        let locate_indices: Vec<LocateIndex> = transcripts
            .iter()
            .map(LocateIndex::build)
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            transcripts: Arc::from(transcripts.into_boxed_slice()),
            locate_indices: Arc::from(locate_indices.into_boxed_slice()),
            trees,
            by_accession,
        })
    }

    /// Return every transcript whose `tx_start..tx_end` overlaps the
    /// half-open genomic interval `[start, end)` on `chrom`.
    ///
    /// Unknown chromosomes return an empty vec (not an error). Callers that
    /// need to distinguish "no chromosome" from "no overlap" should check
    /// [`TranscriptStore::has_chrom`] first.
    ///
    /// # Panics
    ///
    /// Panics if `end > i32::MAX as u64` or `end == 0` — both indicate a
    /// caller bug rather than bad data. Human genomic coordinates never
    /// approach these bounds.
    pub fn query_overlap(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Vec<(&TranscriptModel, &LocateIndex)> {
        assert!(
            end > 0,
            "query_overlap: end must be > 0 (half-open interval)"
        );
        // `start` is validated symmetrically with `end` to guarantee the
        // `start as i32` cast below cannot silently wrap. Human genomic
        // coordinates never approach i32::MAX, so a violation indicates a
        // caller bug.
        assert!(
            start <= i32::MAX as u64,
            "query_overlap: start={start} exceeds i32::MAX"
        );
        assert!(
            end <= i32::MAX as u64,
            "query_overlap: end={end} exceeds i32::MAX"
        );
        assert!(
            start < end,
            "query_overlap: start={start} must be < end={end}"
        );

        let Some(tree) = self.trees.get(chrom) else {
            return Vec::new();
        };

        let mut results = Vec::new();
        tree.query(start as i32, (end - 1) as i32, |node| {
            // coitrees uses different callback types per SIMD backend:
            // neon (ARM): &Interval<&usize>, nosimd (x86_64): &IntervalNode<usize, _>.
            // Explicit &usize annotation + deref coercion handles both portably.
            #[allow(clippy::needless_borrow)]
            let metadata: &usize = &node.metadata;
            let idx = *metadata;
            if let Some(tx) = self.transcripts.get(idx) {
                results.push((tx, &self.locate_indices[idx]));
            }
        });
        results
    }

    /// Direct lookup by full accession with version, e.g. `"NM_006772.2"`.
    /// Returns both the transcript model and its precomputed [`LocateIndex`].
    pub fn get_by_accession(&self, accession: &str) -> Option<(&TranscriptModel, &LocateIndex)> {
        self.by_accession
            .get(accession)
            .map(|&idx| (&self.transcripts[idx], &self.locate_indices[idx]))
    }

    /// Total number of transcripts in the store.
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// `true` if the store is empty.
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Read-only access to the full backing slice. Useful for
    /// aggregate statistics (e.g., counting tiers, protein-coding
    /// transcripts) without going through the interval tree.
    pub fn transcripts(&self) -> &[TranscriptModel] {
        &self.transcripts
    }

    /// `true` if the store has at least one transcript on `chrom`.
    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.trees.contains_key(chrom)
    }
}

impl std::fmt::Debug for TranscriptStore {
    /// Manually implemented because `COITree` does not derive `Debug`. Prints
    /// a compact summary — `len`, chromosome count, and the first few
    /// accessions — instead of dumping the full interval tree.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Snapshot of the first few accessions. `BTreeMap`-style sort is not
        // worth doing here; the vec is already sorted by accession at build
        // time, so `take(3)` gives a stable preview.
        let accessions: Vec<&str> = self
            .transcripts
            .iter()
            .take(3)
            .map(|t| t.accession.as_str())
            .collect();
        f.debug_struct("TranscriptStore")
            .field("len", &self.transcripts.len())
            .field("chromosomes", &self.trees.len())
            .field("first_accessions", &accessions)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Biotype, CdsSegment, Exon, Strand, TranscriptTier};

    fn plus_strand_3_exon() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_000001.1".into(),
            protein_accession: Some("NP_000001.1".into()),
            gene_symbol: "FAKE1".into(),
            hgnc_id: Some("HGNC:1".into()),
            ensembl_accession: Some("ENST00000000001.1".into()),
            chrom: "chr1".into(),
            strand: Strand::Plus,
            tx_start: 1_000,
            tx_end: 5_000,
            cds_genomic_start: Some(1_500),
            cds_genomic_end: Some(4_500),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: 1_000,
                    genomic_end: 2_000,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: 3_000,
                    genomic_end: 3_500,
                },
                Exon {
                    exon_number: 3,
                    genomic_start: 4_000,
                    genomic_end: 5_000,
                },
            ],
            cds_segments: vec![
                CdsSegment {
                    exon_index: 0,
                    genomic_start: 1_500,
                    genomic_end: 2_000,
                    phase: 0,
                },
                CdsSegment {
                    exon_index: 1,
                    genomic_start: 3_000,
                    genomic_end: 3_500,
                    phase: 2,
                },
                CdsSegment {
                    exon_index: 2,
                    genomic_start: 4_000,
                    genomic_end: 4_500,
                    phase: 0,
                },
            ],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 3,
        }
    }

    fn minus_strand_5_exon() -> TranscriptModel {
        // Minus strand: exon 1 (transcript 5') is at the highest genomic
        // coordinate. Exon vector is in transcript order, so genomic
        // coordinates decrease as the index increases.
        TranscriptModel {
            accession: "NM_000002.1".into(),
            protein_accession: Some("NP_000002.1".into()),
            gene_symbol: "FAKE2".into(),
            hgnc_id: Some("HGNC:2".into()),
            ensembl_accession: None,
            chrom: "chr17".into(),
            strand: Strand::Minus,
            tx_start: 10_000,
            tx_end: 20_000,
            cds_genomic_start: Some(11_000),
            cds_genomic_end: Some(19_000),
            exons: vec![
                Exon {
                    exon_number: 1,
                    genomic_start: 19_000,
                    genomic_end: 20_000,
                },
                Exon {
                    exon_number: 2,
                    genomic_start: 17_000,
                    genomic_end: 18_000,
                },
                Exon {
                    exon_number: 3,
                    genomic_start: 15_000,
                    genomic_end: 16_000,
                },
                Exon {
                    exon_number: 4,
                    genomic_start: 13_000,
                    genomic_end: 14_000,
                },
                Exon {
                    exon_number: 5,
                    genomic_start: 10_000,
                    genomic_end: 11_000,
                },
            ],
            cds_segments: vec![],
            tier: TranscriptTier::ManePlusClinical,
            biotype: Biotype::ProteinCoding,
            exon_count: 5,
        }
    }

    fn noncoding_1_exon() -> TranscriptModel {
        TranscriptModel {
            accession: "NR_000001.1".into(),
            protein_accession: None,
            gene_symbol: "FAKE3".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr2".into(),
            strand: Strand::Plus,
            tx_start: 100,
            tx_end: 500,
            cds_genomic_start: None,
            cds_genomic_end: None,
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: 100,
                genomic_end: 500,
            }],
            cds_segments: vec![],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::NonCodingRna,
            exon_count: 1,
        }
    }

    fn patch_sequence_transcript() -> TranscriptModel {
        TranscriptModel {
            accession: "NM_000003.1".into(),
            protein_accession: Some("NP_000003.1".into()),
            gene_symbol: "FAKE4".into(),
            hgnc_id: Some("HGNC:4".into()),
            ensembl_accession: None,
            chrom: "NW_025791820.1".into(),
            strand: Strand::Plus,
            tx_start: 500,
            tx_end: 1_500,
            cds_genomic_start: Some(600),
            cds_genomic_end: Some(1_400),
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: 500,
                genomic_end: 1_500,
            }],
            cds_segments: vec![CdsSegment {
                exon_index: 0,
                genomic_start: 600,
                genomic_end: 1_400,
                phase: 0,
            }],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 1,
        }
    }

    fn build_sample_store() -> TranscriptStore {
        TranscriptStore::from_transcripts(vec![
            plus_strand_3_exon(),
            minus_strand_5_exon(),
            noncoding_1_exon(),
            patch_sequence_transcript(),
        ])
    }

    #[test]
    fn store_length_and_emptiness() {
        let store = build_sample_store();
        assert_eq!(store.len(), 4);
        assert!(!store.is_empty());
    }

    #[test]
    fn query_overlap_returns_plus_strand_transcript() {
        let store = build_sample_store();
        let hits = store.query_overlap("chr1", 1_500, 1_501);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0.accession, "NM_000001.1");
    }

    #[test]
    fn query_overlap_returns_minus_strand_with_reverse_genomic_order() {
        let store = build_sample_store();
        let hits = store.query_overlap("chr17", 12_000, 12_001);
        assert_eq!(hits.len(), 1);
        let (tx, _) = hits[0];
        assert_eq!(tx.strand, Strand::Minus);
        assert!(tx.exons[0].genomic_start > tx.exons[1].genomic_start);
        assert_eq!(tx.exons[0].exon_number, 1);
        assert_eq!(tx.exons.last().unwrap().exon_number, 5);
    }

    #[test]
    fn query_overlap_unknown_chrom_returns_empty() {
        let store = build_sample_store();
        assert!(store.query_overlap("chr99", 1, 1_000).is_empty());
    }

    #[test]
    fn query_overlap_exclusive_end_does_not_match() {
        let store = build_sample_store();
        let hits = store.query_overlap("chr1", 500, 1_000);
        assert!(
            hits.iter().all(|(tx, _)| tx.accession != "NM_000001.1"),
            "exclusive-end semantics violated: [500, 1000) should not touch tx_start=1000"
        );
    }

    #[test]
    fn query_overlap_boundary_inclusive_start() {
        let store = build_sample_store();
        let hits = store.query_overlap("chr1", 1_000, 1_001);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0.accession, "NM_000001.1");
    }

    #[test]
    fn get_by_accession_round_trips_every_record() {
        let store = build_sample_store();
        assert_eq!(
            store.get_by_accession("NM_000001.1").unwrap().0.gene_symbol,
            "FAKE1"
        );
        assert_eq!(
            store.get_by_accession("NM_000002.1").unwrap().0.gene_symbol,
            "FAKE2"
        );
        assert_eq!(
            store.get_by_accession("NR_000001.1").unwrap().0.gene_symbol,
            "FAKE3"
        );
        assert_eq!(
            store.get_by_accession("NM_000003.1").unwrap().0.gene_symbol,
            "FAKE4"
        );
        assert!(store.get_by_accession("nonexistent").is_none());
    }

    #[test]
    fn patch_sequence_transcripts_remain_queryable_by_accession() {
        let store = build_sample_store();
        let (tx, _) = store.get_by_accession("NM_000003.1").unwrap();
        assert_eq!(tx.chrom, "NW_025791820.1");
        assert!(crate::chrom::is_patch_sequence(&tx.chrom));
    }

    #[test]
    fn messagepack_round_trip_preserves_all_fields() {
        let original = vec![
            plus_strand_3_exon(),
            minus_strand_5_exon(),
            noncoding_1_exon(),
        ];
        let bytes = rmp_serde::to_vec_named(&original).expect("serialize");
        let decoded: Vec<TranscriptModel> = rmp_serde::from_slice(&bytes).expect("deserialize");
        let store = TranscriptStore::from_transcripts(decoded);

        let (tx, _) = store.get_by_accession("NM_000001.1").unwrap();
        assert_eq!(tx.strand, Strand::Plus);
        assert_eq!(tx.exons.len(), 3);
        assert_eq!(tx.exons[0].genomic_start, 1_000);
        assert_eq!(tx.exons[2].genomic_end, 5_000);
        assert_eq!(tx.tier, TranscriptTier::ManeSelect);
        assert_eq!(tx.cds_segments.len(), 3);
        assert_eq!(tx.cds_segments[0].phase, 0);
        assert_eq!(tx.cds_segments[1].phase, 2);
        assert_eq!(tx.cds_segments[0].exon_index, 0);
        assert_eq!(tx.cds_segments[2].exon_index, 2);
        assert_eq!(tx.biotype, Biotype::ProteinCoding);

        let (minus_tx, _) = store.get_by_accession("NM_000002.1").unwrap();
        assert_eq!(minus_tx.tier, TranscriptTier::ManePlusClinical);
        assert_eq!(minus_tx.ensembl_accession, None);

        let (nc, _) = store.get_by_accession("NR_000001.1").unwrap();
        assert_eq!(nc.cds_genomic_start, None);
        assert_eq!(nc.cds_genomic_end, None);
        assert_eq!(nc.protein_accession, None);
        assert_eq!(nc.biotype, Biotype::NonCodingRna);
        assert!(nc.cds_segments.is_empty());
    }

    #[test]
    fn debug_impl_prints_compact_summary() {
        // With the manual Debug impl, TranscriptStore formats to a summary
        // instead of dumping the interval tree. Guard against a future
        // refactor that accidentally derives Debug and blows up on COITree.
        let store = build_sample_store();
        let formatted = format!("{store:?}");
        assert!(formatted.contains("TranscriptStore"));
        assert!(formatted.contains("len"));
        assert!(formatted.contains("first_accessions"));
    }

    #[test]
    fn rejects_zero_length_transcript() {
        let bad = vec![TranscriptModel {
            accession: "BAD.1".into(),
            protein_accession: None,
            gene_symbol: "BAD".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr1".into(),
            strand: Strand::Plus,
            tx_start: 100,
            tx_end: 100,
            cds_genomic_start: None,
            cds_genomic_end: None,
            exons: vec![],
            cds_segments: vec![],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::Unknown,
            exon_count: 0,
        }];
        // Explicit match instead of `unwrap_err` because TranscriptStore
        // owns a COITree that does not implement Debug.
        match TranscriptStore::try_from_transcripts(bad) {
            Ok(_) => panic!("expected Malformed error for zero-length transcript"),
            Err(VarEffectError::Malformed(msg)) => assert!(msg.contains("tx_start=100")),
            Err(other) => panic!("expected Malformed, got {other:?}"),
        }
    }

    #[test]
    fn rejects_duplicate_accession() {
        let dup = vec![plus_strand_3_exon(), plus_strand_3_exon()];
        match TranscriptStore::try_from_transcripts(dup) {
            Ok(_) => panic!("expected Malformed error for duplicate accession"),
            Err(VarEffectError::Malformed(msg)) => assert!(msg.contains("duplicate")),
            Err(other) => panic!("expected Malformed, got {other:?}"),
        }
    }

    #[test]
    fn has_chrom_reports_populated_chromosomes() {
        let store = build_sample_store();
        assert!(store.has_chrom("chr1"));
        assert!(store.has_chrom("chr17"));
        assert!(store.has_chrom("NW_025791820.1"));
        assert!(!store.has_chrom("chrX"));
    }

    // Compile-time sanity check: u16 holds exon counts for every known
    // human transcript. TTN, the longest known, has ~363 exons.
    const _: () = assert!(u16::MAX as usize > 363);
}
