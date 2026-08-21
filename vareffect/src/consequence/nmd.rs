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

use crate::hgvs_p::{CdsEdit, PtcSite};
use crate::locate::LocateIndex;

use super::PtcStatus;

/// A termination codon must sit strictly more than this many nucleotides
/// upstream of the last exon-exon junction to trigger NMD.
const NMD_MIN_UPSTREAM_DISTANCE_NT: u32 = 50;

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
    let Some(last_junction_cds_pos) = last_junction_cds_pos(index) else {
        return false;
    };

    // PTC in or past the last coding exon -> NMD escape.
    if cds_position >= last_junction_cds_pos {
        return false;
    }

    // Strictly more than 50 nt upstream -> NMD predicted.
    last_junction_cds_pos - cds_position > NMD_MIN_UPSTREAM_DISTANCE_NT
}

/// 1-based reference-CDS position at which the last coding exon begins.
///
/// `None` when the transcript has no CDS exon-exon junction at all: a
/// non-coding transcript, a missing segment mapping, or a single CDS segment
/// (an intronless transcript deposits no exon-junction complex).
fn last_junction_cds_pos(index: &LocateIndex) -> Option<u32> {
    let cds_end_exon = index.cds_end_exon_idx()?;

    let last_seg_idx = match index.exon_to_cds_seg().get(cds_end_exon) {
        Some(Some(seg)) => *seg as usize,
        _ => return None,
    };

    // Single CDS segment -> no exon-exon junction within CDS.
    if last_seg_idx == 0 {
        return None;
    }

    // cumulative_cds[i] = total CDS bases in segments 0..i (prefix sum), so
    // cumulative_cds[last_seg_idx] is the 0-based CDS offset where the last
    // segment begins. Convert to 1-based for comparison.
    Some(index.cumulative_cds()[last_seg_idx] + 1)
}

/// Attach an NMD verdict to a located termination codon.
///
/// Unlike [`predicts_nmd`], which measures at the variant site for Ensembl VEP
/// `NMD.pm` parity, this evaluates the 50-nucleotide rule **at the termination
/// codon** -- the input the ClinGen SVI PVS1 decision tree asks for (Abou
/// Tayoun et al. 2018, PMID 30192042, doi:10.1002/humu.23626).
///
/// # Coordinate spaces
///
/// [`PtcSite::Residue`] is numbered in the **alternate** protein, while
/// `index`'s junctions are **reference**-CDS coordinates. Rather than mapping
/// the codon back into reference space -- undefined when the codon lies inside
/// inserted sequence -- this projects the junction *forward* into alternate
/// space, where the comparison is always well defined:
///
/// - A junction downstream of the edit shifts by [`CdsEdit::delta`].
/// - A junction at or before the edit does not move. The termination codon is
///   always at or after the edit, so it necessarily sits in the last exon and
///   escapes.
///
/// Shifting unconditionally would push the junction past a termination codon
/// already inside the last exon and manufacture a spurious NMD prediction --
/// an overcall.
///
/// # Arguments
///
/// * `site` -- the located termination codon, from the frameshift formatter.
/// * `edit` -- the 3'-normalized CDS edit `site` is relative to.
/// * `index` -- precomputed locate index for the transcript.
pub(super) fn resolve_ptc(site: PtcSite, edit: &CdsEdit, index: &LocateIndex) -> PtcStatus {
    let residue = match site {
        PtcSite::Indeterminate => return PtcStatus::Indeterminate,
        PtcSite::NoStopCodon => return PtcStatus::NoStopCodon,
        PtcSite::Residue(r) => r,
    };

    PtcStatus::At {
        protein_position: residue,
        nmd_at_ptc: nmd_at_residue(residue, edit, index),
    }
}

/// Locate a premature termination codon in a translated alternate window.
///
/// `alt_aas` is the alternate translation beginning at 0-based codon
/// `first_codon`; `ref_aas` is the reference translation of the same window.
///
/// Returns [`PtcStatus::Indeterminate`] when the reference window already
/// contained a stop, because the `*` found in `alt_aas` may then be the
/// transcript's normal termination codon rather than a premature one. Most
/// callers already guard `StopGained` on `!ref_aas.contains(&b'*')`, but
/// `annotate_cds_inframe_deletion` does not -- an in-frame deletion
/// overlapping the normal stop reaches here with the reference stop still in
/// view, and naming it a PTC would be a false truncation claim.
///
/// # Arguments
///
/// * `alt_aas`, `ref_aas` -- translated alternate and reference windows.
/// * `first_codon` -- 0-based codon index the windows begin at.
/// * `edit` -- the CDS edit the windows describe.
/// * `index` -- precomputed locate index for the transcript.
pub(super) fn ptc_from_alt_window(
    alt_aas: &[u8],
    ref_aas: &[u8],
    first_codon: u32,
    edit: &CdsEdit,
    index: &LocateIndex,
) -> PtcStatus {
    if ref_aas.contains(&b'*') {
        return PtcStatus::Indeterminate;
    }
    let Some(offset) = alt_aas.iter().position(|&aa| aa == b'*') else {
        return PtcStatus::Indeterminate;
    };
    // alt_aas[i] is the residue at 0-based codon `first_codon + i`.
    let residue = first_codon + offset as u32 + 1;
    resolve_ptc(PtcSite::Residue(residue), edit, index)
}

/// Evaluate the 50-nucleotide rule at a 1-based alternate-protein residue.
fn nmd_at_residue(residue: u32, edit: &CdsEdit, index: &LocateIndex) -> bool {
    let Some(junction_ref) = last_junction_cds_pos(index) else {
        return false;
    };

    // First base of the terminating codon, 1-based in the alternate CDS.
    let ptc_pos = i64::from(residue) * 3 - 2;

    // Project the junction into alternate-CDS space. Only a junction strictly
    // downstream of the replaced span moves.
    let junction_alt = if junction_ref > edit.last_replaced_pos() {
        i64::from(junction_ref) + edit.delta()
    } else {
        i64::from(junction_ref)
    };

    ptc_pos < junction_alt && junction_alt - ptc_pos > i64::from(NMD_MIN_UPSTREAM_DISTANCE_NT)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locate::LocateIndex;
    use crate::test_fixtures::{minus_strand_coding, plus_strand_coding, single_exon_coding};

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
        };
        let idx = build_index(&tx);
        assert!(!predicts_nmd(10, &idx));
    }

    // ---- resolve_ptc: the 50-nt rule measured at the termination codon ----
    //
    // plus_strand_coding() again: last junction at reference CDS pos 1001.

    /// Edit with no length change, anchored well upstream of the junction.
    fn neutral_edit(start: u32) -> CdsEdit {
        CdsEdit {
            start,
            del_len: 1,
            ins_len: 1,
        }
    }

    #[test]
    fn resolve_ptc_passes_through_non_positions() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        assert_eq!(
            resolve_ptc(PtcSite::Indeterminate, &neutral_edit(10), &idx),
            PtcStatus::Indeterminate
        );
        assert_eq!(
            resolve_ptc(PtcSite::NoStopCodon, &neutral_edit(10), &idx),
            PtcStatus::NoStopCodon
        );
    }

    #[test]
    fn resolve_ptc_upstream_residue_predicts_nmd() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // Residue 100 -> CDS 298; junction 1001 -> distance 703 > 50.
        assert_eq!(
            resolve_ptc(PtcSite::Residue(100), &neutral_edit(297), &idx),
            PtcStatus::At {
                protein_position: 100,
                nmd_at_ptc: true,
            }
        );
    }

    #[test]
    fn resolve_ptc_last_exon_residue_escapes() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // Residue 400 -> CDS 1198, past the junction at 1001.
        assert_eq!(
            resolve_ptc(PtcSite::Residue(400), &neutral_edit(1197), &idx),
            PtcStatus::At {
                protein_position: 400,
                nmd_at_ptc: false,
            }
        );
    }

    #[test]
    fn resolve_ptc_shifts_junction_for_upstream_deletion() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // A 100 bp deletion at CDS offset 880 pulls the junction back to 901
        // in the alternate transcript. The termination codon at alt CDS 880
        // is then only 21 nt upstream of it -> escape.
        //
        // Comparing the alternate-space codon against the *reference* junction
        // would give 1001 - 880 = 121 and wrongly predict NMD.
        let edit = CdsEdit {
            start: 880,
            del_len: 100,
            ins_len: 0,
        };
        assert_eq!(
            resolve_ptc(PtcSite::Residue(294), &edit, &idx),
            PtcStatus::At {
                protein_position: 294,
                nmd_at_ptc: false,
            }
        );
    }

    #[test]
    fn resolve_ptc_does_not_shift_junction_for_downstream_edit() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // The edit sits inside the last exon, so the junction does not move.
        // The termination codon is downstream of it either way -> escape.
        //
        // Applying the +99 delta unconditionally would place the codon at
        // 913, i.e. 88 nt "upstream" of the junction, and manufacture an NMD
        // prediction for a variant that plainly escapes. This test fails if
        // the junction-projection guard is removed.
        let edit = CdsEdit {
            start: 1010,
            del_len: 1,
            ins_len: 100,
        };
        assert_eq!(
            resolve_ptc(PtcSite::Residue(338), &edit, &idx),
            PtcStatus::At {
                protein_position: 338,
                nmd_at_ptc: false,
            }
        );
    }

    #[test]
    fn resolve_ptc_boundary_is_strictly_greater_than_50() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // No length change, so the junction stays at 1001.
        // Residue 318 -> CDS 952; distance 49 -> escape.
        assert!(matches!(
            resolve_ptc(PtcSite::Residue(318), &neutral_edit(951), &idx),
            PtcStatus::At {
                nmd_at_ptc: false,
                ..
            }
        ));
        // Residue 317 -> CDS 949; distance 52 -> NMD.
        assert!(matches!(
            resolve_ptc(PtcSite::Residue(317), &neutral_edit(948), &idx),
            PtcStatus::At {
                nmd_at_ptc: true,
                ..
            }
        ));
    }

    #[test]
    fn resolve_ptc_single_exon_never_predicts_nmd() {
        let tx = single_exon_coding();
        let idx = build_index(&tx);
        assert_eq!(
            resolve_ptc(PtcSite::Residue(5), &neutral_edit(12), &idx),
            PtcStatus::At {
                protein_position: 5,
                nmd_at_ptc: false,
            }
        );
    }

    #[test]
    fn ptc_from_alt_window_locates_stop_and_rejects_reference_stop() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        let edit = neutral_edit(297);
        // Window starts at 0-based codon 99, so alt_aas[1] is residue 101.
        assert_eq!(
            ptc_from_alt_window(b"MK*", b"MKR", 99, &edit, &idx),
            PtcStatus::At {
                protein_position: 102,
                nmd_at_ptc: true,
            }
        );
        // Reference already terminated here: the `*` found may be the normal
        // stop, so refuse to call it premature.
        assert_eq!(
            ptc_from_alt_window(b"MK*", b"M*R", 99, &edit, &idx),
            PtcStatus::Indeterminate
        );
        // No stop in the alternate window at all.
        assert_eq!(
            ptc_from_alt_window(b"MKR", b"MKQ", 99, &edit, &idx),
            PtcStatus::Indeterminate
        );
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
