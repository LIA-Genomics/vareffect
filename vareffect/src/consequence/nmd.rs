//! NMD prediction via the 50-nucleotide rule.
//!
//! Predicts whether a premature termination codon triggers nonsense-mediated
//! mRNA decay. Used to populate [`super::ConsequenceResult::predicts_nmd`].
//!
//! Uses CDS-only junctions via [`LocateIndex`]. For MANE Select transcripts
//! where the stop codon is in the last exon (>99% of cases), this is
//! identical to VEP's transcript-level junction check. For the rare case
//! where 3'UTR spans additional exons, the CDS-only approach may miss
//! some NMD predictions (conservative: reports NMD escape, which leads to
//! PVS1 downgrade rather than overcalling pathogenicity).

use crate::locate::LocateIndex;

/// Returns `true` if a PTC at the given CDS position is predicted to
/// trigger NMD (>50 nt upstream of the last exon-exon junction).
///
/// Returns `false` (NMD escape) if:
/// - Transcript has only 1 CDS segment (no exon-exon junction in CDS)
/// - PTC is in the last coding exon
/// - PTC is within 50 nt of the last exon-exon junction
///
/// # Arguments
///
/// * `cds_position` -- 1-based CDS position of the variant site (not the
///   downstream frameshift stop). Matches VEP's convention of using the
///   variant position rather than the predicted termination site.
/// * `index` -- Precomputed locate index for the transcript.
pub(super) fn predicts_nmd(cds_position: u32, index: &LocateIndex) -> bool {
    // Find the last coding exon's CDS segment index.
    let cds_end_exon = match index.cds_end_exon_idx() {
        Some(idx) => idx,
        None => return false, // non-coding transcript
    };

    let last_seg_idx = match index.exon_to_cds_seg().get(cds_end_exon) {
        Some(Some(seg)) => *seg as usize,
        _ => return false,
    };

    // Single CDS segment -> no exon-exon junction within CDS.
    if last_seg_idx == 0 {
        return false;
    }

    // Last junction = boundary between penultimate and last CDS segments.
    // cumulative_cds[i] = total CDS bases in segments 0..i (prefix sum),
    // so cumulative_cds[last_seg_idx] is the CDS offset where the last
    // segment begins (0-based). Convert to 1-based for comparison.
    let last_junction_cds_pos = index.cumulative_cds()[last_seg_idx] + 1;

    // PTC in or past the last coding exon -> NMD escape.
    if cds_position >= last_junction_cds_pos {
        return false;
    }

    // Strictly more than 50 nt upstream -> NMD predicted.
    last_junction_cds_pos - cds_position > 50
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locate::LocateIndex;
    use crate::test_fixtures::{minus_strand_coding, plus_strand_coding};

    /// Helper: build LocateIndex for a transcript fixture.
    fn build_index(tx: &crate::types::TranscriptModel) -> LocateIndex {
        LocateIndex::build(tx).unwrap()
    }

    // plus_strand_coding(): 3 exons, CDS segs [500, 500, 500],
    // cumulative_cds = [0, 500, 1000, 1500].
    // Last CDS seg idx = 2, last junction CDS pos = 1000 + 1 = 1001.

    #[test]
    fn nmd_ptc_far_upstream() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 100, junction at 1001 -> distance 901 > 50
        assert!(predicts_nmd(100, &idx));
    }

    #[test]
    fn nmd_ptc_in_last_exon() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 1100, junction at 1001 -> in last exon
        assert!(!predicts_nmd(1100, &idx));
    }

    #[test]
    fn nmd_ptc_within_50nt() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 960, junction at 1001 -> distance 41 <= 50
        assert!(!predicts_nmd(960, &idx));
    }

    #[test]
    fn nmd_ptc_exactly_at_50nt() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 951, junction at 1001 -> distance 50 = 50, NOT > 50
        assert!(!predicts_nmd(951, &idx));
    }

    #[test]
    fn nmd_ptc_at_51nt() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 950, junction at 1001 -> distance 51 > 50
        assert!(predicts_nmd(950, &idx));
    }

    #[test]
    fn nmd_single_exon() {
        // stop_gained_transcript() has 1 exon -> 1 CDS segment -> false.
        // Inline a minimal single-exon transcript rather than importing
        // the test-only builder from consequence::tests.
        use crate::types::*;
        let tx = TranscriptModel {
            accession: "NM_SINGLE.1".into(),
            protein_accession: Some("NP_SINGLE.1".into()),
            gene_symbol: "SINGLE".into(),
            hgnc_id: None,
            ensembl_accession: None,
            chrom: "chr1".into(),
            strand: Strand::Plus,
            tx_start: 0,
            tx_end: 300,
            cds_genomic_start: Some(100),
            cds_genomic_end: Some(200),
            exons: vec![Exon {
                exon_number: 1,
                genomic_start: 0,
                genomic_end: 300,
            }],
            cds_segments: vec![CdsSegment {
                exon_index: 0,
                genomic_start: 100,
                genomic_end: 200,
                phase: 0,
            }],
            tier: TranscriptTier::ManeSelect,
            biotype: Biotype::ProteinCoding,
            exon_count: 1,
            genome_transcript_divergent: false,
            translational_exception: None,
        };
        let idx = build_index(&tx);
        assert!(!predicts_nmd(10, &idx));
    }

    #[test]
    fn nmd_minus_strand() {
        // minus_strand_coding(): 3 exons, CDS segs [1500, 2000, 1000],
        // cumulative_cds = [0, 1500, 3500, 4500].
        // Last seg idx = 2, last junction CDS pos = 3500 + 1 = 3501.
        let tx = minus_strand_coding();
        let idx = build_index(&tx);
        // CDS pos 100, junction at 3501 -> distance 3401 > 50
        assert!(predicts_nmd(100, &idx));
        // CDS pos 3500, junction at 3501 -> distance 1 <= 50
        assert!(!predicts_nmd(3500, &idx));
    }
}
