//! NMD prediction via the 50-nucleotide rule.
//!
//! Predicts whether a premature termination codon triggers nonsense-mediated
//! mRNA decay. Populates [`super::ConsequenceResult::predicts_nmd`] (measured
//! at the variant site, for Ensembl VEP `NMD.pm` parity) and
//! [`super::ConsequenceResult::ptc`] (measured at the termination codon, for
//! ACMG PVS1).
//!
//! Uses CDS-only junctions via [`LocateIndex`]. For MANE Select transcripts
//! where the stop codon is in the last exon (>99% of cases), this is
//! identical to VEP's transcript-level junction check. For the rare case
//! where 3'UTR spans additional exons, the CDS-only approach may miss
//! some NMD predictions (conservative: reports NMD escape, which leads to
//! PVS1 downgrade rather than overcalling pathogenicity).
//!
//! # Which escape rules this implements, and which it does not
//!
//! The termination-codon verdict implements the **last-junction rule only**,
//! because that is the rule the ClinGen SVI PVS1 decision tree is written
//! against (Abou Tayoun et al. 2018, PMID 30192042, doi:10.1002/humu.23626),
//! and PVS1 is what this field exists to feed. `AutoPVS1`, the published
//! reference implementation of that flowchart (Xiang et al. 2020, PMID
//! 32442321, doi:10.1002/humu.24051), likewise implements no other escape
//! rule.
//!
//! Three further escape rules are well described in the NMD literature and are
//! **deliberately not applied here**:
//!
//! - **Start-proximal escape.** A termination codon close to the coding start
//!   site can be bypassed by translational re-initiation. `aenmd` puts the
//!   default cutoff at 150 nt downstream of the coding start site, after
//!   Lindeboom et al. (PMID 27618451, doi:10.1038/ng.3664; PMID 31659324,
//!   doi:10.1038/s41588-019-0517-5). Ensembl VEP's `NMD.pm` plugin carries a
//!   coarser variant of it — escape when the variant's `cds_end` is `<= 101`.
//! - **Long-exon escape.** A termination codon inside an exon longer than
//!   407 bp, same sources.
//! - **Intronless transcripts**, which this module *does* honour, but as a
//!   consequence of there being no CDS junction rather than as a named rule.
//!
//! Applying the first two would withdraw NMD — and with it PVS1's full
//! strength — from every truncating variant in the first 150 coding bases or
//! in a long exon, a large and clinically well-established class of
//! loss-of-function alleles. That is a guideline change, not an annotation
//! detail, so it belongs to the consumer: this module reports the termination
//! codon's position, which is everything a consumer needs to layer either rule
//! itself.

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
    // Both `Absent` and `Unknown` yield `false` here, preserving this
    // function's exact pre-existing behavior. Only `nmd_at_residue`
    // distinguishes them, because only it can report "no answer".
    let LastJunction::At(last_junction_cds_pos) = last_junction_cds_pos(index) else {
        return false;
    };

    // PTC in or past the last coding exon -> NMD escape.
    if cds_position >= last_junction_cds_pos {
        return false;
    }

    // Strictly more than 50 nt upstream -> NMD predicted.
    last_junction_cds_pos - cds_position > NMD_MIN_UPSTREAM_DISTANCE_NT
}

/// Where the last coding exon begins, or why that could not be determined.
///
/// [`LastJunction::Absent`] and [`LastJunction::Unknown`] are deliberately
/// distinct: the first is a real answer (there is no junction, so NMD cannot
/// be triggered), the second is missing data. Collapsing the second into a
/// negative NMD verdict would turn an absent index into a clinical assertion.
enum LastJunction {
    /// 1-based reference-CDS position at which the last coding exon begins.
    At(u32),
    /// The transcript has no CDS exon-exon junction: it is non-coding, or its
    /// CDS lies in a single segment, so no exon-junction complex is deposited.
    Absent,
    /// The exon-to-CDS-segment mapping is missing, so the junction cannot be
    /// located. Never a verdict.
    Unknown,
}

/// Locate the last CDS exon-exon junction.
fn last_junction_cds_pos(index: &LocateIndex) -> LastJunction {
    let Some(cds_end_exon) = index.cds_end_exon_idx() else {
        return LastJunction::Absent; // non-coding transcript
    };

    let last_seg_idx = match index.exon_to_cds_seg().get(cds_end_exon) {
        Some(Some(seg)) => *seg as usize,
        // The CDS-end exon carries no CDS segment, or is out of range: the
        // index disagrees with itself. Missing data, not "no junction".
        _ => return LastJunction::Unknown,
    };

    // Single CDS segment -> no exon-exon junction within CDS.
    if last_seg_idx == 0 {
        return LastJunction::Absent;
    }

    // cumulative_cds[i] = total CDS bases in segments 0..i (prefix sum), so
    // cumulative_cds[last_seg_idx] is the 0-based CDS offset where the last
    // segment begins. Convert to 1-based for comparison.
    LastJunction::At(index.cumulative_cds()[last_seg_idx] + 1)
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

    match nmd_at_residue(residue, edit, index) {
        Some(nmd_at_ptc) => PtcStatus::At {
            protein_position: residue,
            nmd_at_ptc,
        },
        // The codon was located but no NMD verdict could be reached. Reporting
        // the position with a guessed verdict would be worse than withholding.
        None => PtcStatus::Indeterminate,
    }
}

/// Resolve the PTC for a result whose consequence set may or may not be
/// truncating.
///
/// Folds the `StopGained` gate in once so the six coding call sites cannot
/// drift apart: a change to what counts as truncating lands in one place.
///
/// # Arguments
///
/// * `consequences` -- the finalized consequence set for this result.
/// * `alt_aas`, `ref_aas` -- translated alternate and reference windows.
/// * `first_codon` -- 0-based codon index the windows begin at.
/// * `edit` -- the CDS edit the windows describe.
/// * `index` -- precomputed locate index for the transcript.
pub(super) fn ptc_for_window(
    consequences: &[super::Consequence],
    alt_aas: &[u8],
    ref_aas: &[u8],
    first_codon: u32,
    edit: &CdsEdit,
    index: &LocateIndex,
) -> PtcStatus {
    if !consequences.contains(&super::Consequence::StopGained) {
        return PtcStatus::NotApplicable;
    }
    ptc_from_alt_window(alt_aas, ref_aas, first_codon, edit, index)
}

/// Locate a premature termination codon in a translated alternate window.
///
/// `alt_aas` is the alternate translation beginning at 0-based codon
/// `first_codon`; `ref_aas` is the reference translation of the same window.
///
/// Returns [`PtcStatus::Indeterminate`] when the reference window already
/// contained a stop, because the `*` found in `alt_aas` may then be the
/// transcript's normal termination codon rather than a premature one. Every
/// caller now guards `StopGained` on `!ref_aas.contains(&b'*')` as well, so
/// this is a second line of defence rather than the only one: a future caller
/// that forgets the guard gets a withheld verdict instead of a false
/// truncation claim.
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
///
/// # The anchor base
///
/// Measured from the **3'-most base** of the termination codon. ClinGen SVI
/// states the rule as "NMD is not predicted to occur if the premature
/// termination codon occurs in the 3' most exon or within the 3'-most 50
/// nucleotides of the penultimate exon" (Abou Tayoun et al. 2018). A codon
/// spans three bases, so "occurs within" is read as *overlaps*: if any base of
/// the codon lies in that window the codon escapes, and the 3'-most base is
/// the one that decides it. Anchoring on the codon's first base instead would
/// shrink the escape window to 48 nucleotides and predict NMD for two codon
/// positions the guideline places outside it — an overcall, since predicting
/// NMD is what earns PVS1 its full strength.
///
/// # Returns
///
/// `None` when no verdict can be reached — the junction index is
/// self-inconsistent, or the junction falls inside the edited span so its
/// position in the alternate transcript is undefined. Never guess in those
/// cases: an invented verdict here becomes a PVS1 strength downstream.
fn nmd_at_residue(residue: u32, edit: &CdsEdit, index: &LocateIndex) -> Option<bool> {
    let junction_ref = match last_junction_cds_pos(index) {
        LastJunction::At(pos) => pos,
        // No junction at all -> nothing to measure against, and NMD genuinely
        // cannot be triggered. That is an answer.
        LastJunction::Absent => return Some(false),
        LastJunction::Unknown => return None,
    };

    // 3'-most base of the terminating codon, 1-based in the alternate CDS.
    let ptc_last_base = i64::from(residue) * 3;

    // Project the junction into alternate-CDS space.
    let junction_alt = if junction_ref > edit.last_replaced_pos() {
        // Downstream of the replaced span: it shifts by the length change.
        i64::from(junction_ref) + edit.delta()
    } else if junction_ref > edit.start {
        // Inside the replaced span: the bases carrying this junction are gone,
        // so it has no well-defined alternate-transcript position. Currently
        // unreachable -- a contiguous deletion spanning two CDS segments must
        // cover the canonical splice dinucleotides and is diverted to the
        // splice path first (`complex.rs`, `overlaps_splice_canonical`) -- but
        // withhold rather than silently treat it as unmoved.
        return None;
    } else {
        // At or before the edit: unmoved. The termination codon is at or after
        // the edit, so it necessarily lies in the last exon and escapes.
        i64::from(junction_ref)
    };

    // A separate "is it upstream of the junction" test would be redundant:
    // exceeding the 50-nt distance already implies it, and the arithmetic is
    // signed so a codon inside the last exon simply yields a negative
    // distance.
    Some(junction_alt - ptc_last_base > i64::from(NMD_MIN_UPSTREAM_DISTANCE_NT))
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
        // Residue 100 -> last base CDS 300; junction 1001 -> distance 701.
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
        // Residue 400 -> last base CDS 1200, past the junction at 1001.
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
        // in the alternate transcript. Residue 294 ends at alt CDS 882, only
        // 19 nt upstream of it -> escape.
        //
        // Comparing the alternate-space codon against the *reference* junction
        // would give 1001 - 882 = 119 and wrongly predict NMD.
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
        // Applying the +99 delta unconditionally would move the junction to
        // 1100, leaving the codon at 1014 an apparent 86 nt upstream of it,
        // and manufacture an NMD prediction for a variant that plainly
        // escapes. This test fails if the junction-projection guard is
        // removed.
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

    /// The escape window is the 3'-most 50 nucleotides of the penultimate
    /// coding exon -- CDS 951..=1000 here -- and a codon escapes as soon as it
    /// reaches into it. Residue 317 ends exactly on 951, the first base of the
    /// window, so it escapes; residue 316 ends on 948 and does not.
    ///
    /// This is the pair that pins the anchor base. Measuring from the codon's
    /// first base instead would put residue 317 at CDS 949, 52 nt from the
    /// junction, and predict NMD for a codon the guideline places inside the
    /// escape window.
    #[test]
    fn resolve_ptc_boundary_is_the_3prime_most_base_of_the_codon() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // No length change, so the junction stays at 1001.
        // Residue 317 -> last base CDS 951; distance 50, not > 50 -> escape.
        assert!(matches!(
            resolve_ptc(PtcSite::Residue(317), &neutral_edit(948), &idx),
            PtcStatus::At {
                nmd_at_ptc: false,
                ..
            }
        ));
        // Residue 316 -> last base CDS 948; distance 53 -> NMD.
        assert!(matches!(
            resolve_ptc(PtcSite::Residue(316), &neutral_edit(945), &idx),
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
    fn resolve_ptc_minus_strand() {
        // minus_strand_coding(): CDS segments [1500, 2000, 1000],
        // cumulative_cds = [0, 1500, 3500, 4500] -> last junction at 3501.
        // The projection is pure CDS space, so strand does not enter it --
        // this pins that, since every other resolve_ptc test is plus-strand.
        let tx = minus_strand_coding();
        let idx = build_index(&tx);

        // Residue 100 -> last base CDS 300; distance 3201 -> NMD.
        assert_eq!(
            resolve_ptc(PtcSite::Residue(100), &neutral_edit(297), &idx),
            PtcStatus::At {
                protein_position: 100,
                nmd_at_ptc: true,
            }
        );
        // Residue 1200 -> last base CDS 3600, past the junction -> escape.
        assert_eq!(
            resolve_ptc(PtcSite::Residue(1200), &neutral_edit(3597), &idx),
            PtcStatus::At {
                protein_position: 1200,
                nmd_at_ptc: false,
            }
        );
        // An upstream deletion drags the junction back with it: a 300 bp
        // deletion at CDS offset 3000 puts the junction at 3201, leaving the
        // codon ending at alt CDS 3183 only 18 nt upstream -> escape.
        let edit = CdsEdit {
            start: 3000,
            del_len: 300,
            ins_len: 0,
        };
        assert_eq!(
            resolve_ptc(PtcSite::Residue(1061), &edit, &idx),
            PtcStatus::At {
                protein_position: 1061,
                nmd_at_ptc: false,
            }
        );
    }

    #[test]
    fn resolve_ptc_withholds_when_junction_is_inside_the_edit() {
        let tx = plus_strand_coding();
        let idx = build_index(&tx);
        // The deletion spans CDS 951..1050, swallowing the junction at 1001.
        // That junction has no position in the alternate transcript, so no
        // verdict is available and the codon position must not be reported
        // with a guessed one.
        let edit = CdsEdit {
            start: 950,
            del_len: 100,
            ins_len: 0,
        };
        assert_eq!(
            resolve_ptc(PtcSite::Residue(318), &edit, &idx),
            PtcStatus::Indeterminate
        );
    }

    #[test]
    fn resolve_ptc_reports_escape_for_an_intronless_transcript() {
        // Absent (no junction) is a real answer, unlike Unknown.
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
