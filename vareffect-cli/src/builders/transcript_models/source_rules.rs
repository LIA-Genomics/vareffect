//! Per-source admission and attribute rules for the transcript builder:
//! RefSeq alignment gaps and translation exceptions, RefSeq biotypes, and
//! the CDS checks that exclude a transcript.

use anyhow::{Context, Result, bail};
use vareffect::{Biotype, CdsSegment, Exon, Strand};

/// One gap of a RefSeq `cDNA_match` alignment, in 0-based genomic
/// coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum GapOp {
    /// Transcript bases the genome lacks (`I`), inserted between genomic
    /// bases `at - 1` and `at`.
    Insertion { at: u64 },
    /// Genomic bases `start..end` the transcript lacks (`D`).
    Deletion { start: u64, end: u64 },
}

/// Parse the `Gap` attribute of a `cDNA_match` row into its gaps.
///
/// NCBI writes the operations in transcript order: ascending from the row
/// start on the plus strand, descending from the row end on the minus
/// strand (checked on NM_001371720.2, `M163 D1 M108`, whose exon split lies
/// 163 bases below the row end).
///
/// # Arguments
///
/// * `gap` - The decoded `Gap` value, e.g. `"M174 I3 M3"`
/// * `start` - Row start, 0-based inclusive
/// * `end` - Row end, 0-based exclusive
/// * `strand` - Row strand (column 7)
///
/// # Returns
///
/// The `I` and `D` gaps; `M` operations only advance the position.
///
/// # Errors
///
/// Returns an error for an operation other than `M`, `I` or `D`, a
/// non-numeric length, or operations that run past the row.
pub(super) fn parse_gap_ops(gap: &str, start: u64, end: u64, strand: Strand) -> Result<Vec<GapOp>> {
    // `Strand` is `#[non_exhaustive]`; a GFF3 row is `+` or `-`.
    let plus = matches!(strand, Strand::Plus);
    let mut ops = Vec::new();
    let mut pos = if plus { start } else { end };
    for op in gap.split_whitespace() {
        let (kind, length) = op.split_at(1);
        let length: u64 = length
            .parse()
            .with_context(|| format!("Gap operation {op:?} has no numeric length"))?;
        match kind {
            "I" => ops.push(GapOp::Insertion { at: pos }),
            "M" | "D" => {
                let (low, high) = if plus {
                    (pos, pos.saturating_add(length))
                } else {
                    (pos.saturating_sub(length), pos)
                };
                if low < start || high > end {
                    bail!("Gap {gap:?} runs past its row {start}..{end}");
                }
                if kind == "D" {
                    ops.push(GapOp::Deletion {
                        start: low,
                        end: high,
                    });
                }
                pos = if plus { high } else { low };
            }
            other => bail!("Gap operation {other:?} in {gap:?} is not M, I or D"),
        }
    }
    Ok(ops)
}

/// `true` when an alignment gap makes the genome-built CDS numbering differ
/// from the transcript's own.
///
/// An `I` gap strictly inside the CDS span always does. A `D` gap does
/// unless RefSeq split the exon there (the deleted bases lie in an intron
/// of the model, a frameshift intron): a `D` gap inside an exon leaves
/// genomic bases in the model that the transcript lacks.
///
/// # Arguments
///
/// * `ops` - The transcript's alignment gaps on this locus
/// * `cds_start` - Lowest CDS base, 0-based inclusive
/// * `cds_end` - CDS end, 0-based exclusive
/// * `exons` - The model's exons, in any order
pub(super) fn gap_moves_cds_numbering(
    ops: &[GapOp],
    cds_start: u64,
    cds_end: u64,
    exons: &[Exon],
) -> bool {
    ops.iter().any(|op| match *op {
        GapOp::Insertion { at } => cds_start < at && at < cds_end,
        GapOp::Deletion { start, end } => {
            start < cds_end
                && cds_start < end
                && exons
                    .iter()
                    .any(|exon| start < exon.genomic_end && exon.genomic_start < end)
        }
    })
}

/// The 1-based model cDNA position of the first exonic base at or after
/// genomic base `pos` in transcript order, or `None` past the last exon.
///
/// # Arguments
///
/// * `exons` - The model's exons in transcript order
/// * `plus` - Transcript strand
/// * `pos` - 0-based genomic base
pub(super) fn cdna_at_or_after(exons: &[Exon], plus: bool, pos: u64) -> Option<u32> {
    let mut before = 0u64;
    for exon in exons {
        if plus && pos < exon.genomic_end {
            return u32::try_from(before + pos.saturating_sub(exon.genomic_start) + 1).ok();
        }
        if !plus && pos >= exon.genomic_start {
            return u32::try_from(before + (exon.genomic_end - 1).saturating_sub(pos) + 1).ok();
        }
        before += exon.genomic_end - exon.genomic_start;
    }
    None
}

/// One alignment gap outside the CDS span, in model cDNA positions
/// (1-based, transcript order).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct UtrGapBreak {
    /// First position numbered differently from the transcript: the base
    /// after an `I` gap, the first base of a `D` gap.
    pub(super) shifted_from: u32,
    /// Last position whose count back to a later base crosses the gap: the
    /// base before an `I` gap (0 at the transcript start), the last base of a
    /// `D` gap.
    pub(super) crossed_through: u32,
}

/// The alignment gaps outside the CDS span (anywhere on a non-coding
/// transcript) that move model numbering: an `I` gap at or inside the
/// transcript's 5' end and before its 3' end, a `D` gap overlapping an exon.
///
/// # Arguments
///
/// * `ops` - The transcript's alignment gaps on this locus
/// * `cds` - The CDS span, 0-based half-open, or `None` for a non-coding
///   transcript
/// * `exons` - The model's exons in transcript order
/// * `plus` - Transcript strand
pub(super) fn utr_gap_breaks(
    ops: &[GapOp],
    cds: Option<(u64, u64)>,
    exons: &[Exon],
    plus: bool,
) -> Vec<UtrGapBreak> {
    let tx_start = exons.iter().map(|e| e.genomic_start).min().unwrap_or(0);
    let tx_end = exons.iter().map(|e| e.genomic_end).max().unwrap_or(0);
    let outside_cds = |start: u64, end: u64| cds.is_none_or(|(cs, ce)| end <= cs || ce <= start);
    ops.iter()
        .filter_map(|op| match *op {
            GapOp::Insertion { at } => {
                let inside = if plus {
                    tx_start <= at && at < tx_end
                } else {
                    tx_start < at && at <= tx_end
                };
                if !inside || !outside_cds(at, at) {
                    return None;
                }
                let shifted_from = cdna_at_or_after(exons, plus, if plus { at } else { at - 1 })?;
                Some(UtrGapBreak {
                    shifted_from,
                    crossed_through: shifted_from - 1,
                })
            }
            GapOp::Deletion { start, end } => {
                let exonic: u64 = exons
                    .iter()
                    .map(|e| {
                        end.min(e.genomic_end)
                            .saturating_sub(start.max(e.genomic_start))
                    })
                    .sum();
                if exonic == 0 || !outside_cds(start, end) {
                    return None;
                }
                let shifted_from =
                    cdna_at_or_after(exons, plus, if plus { start } else { end - 1 })?;
                Some(UtrGapBreak {
                    shifted_from,
                    crossed_through: shifted_from + u32::try_from(exonic).ok()? - 1,
                })
            }
        })
        .collect()
}

/// Selenocysteine codons and the count of other entries in a RefSeq
/// `transl_except` attribute.
#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct TranslExcept {
    /// 0-based genomic start (lowest base) of each `aa:Sec` codon.
    pub(super) selenocysteine: Vec<u64>,
    /// Entries for any other amino acid (e.g. `aa:Met` alternative starts),
    /// which the store does not record.
    pub(super) other: usize,
}

/// Parse a decoded RefSeq `transl_except` value such as
/// `(pos:25802093..25802095,aa:Sec),(pos:complement(49358132..49358134),aa:Sec)`.
///
/// # Arguments
///
/// * `raw` - The percent-decoded attribute value
///
/// # Returns
///
/// The Sec codon starts (lowest 1-based position minus one, on either
/// strand) and the number of non-Sec entries.
///
/// # Errors
///
/// Returns an error for an entry without `pos:` and `aa:` parts or without a
/// position.
pub(super) fn parse_transl_except(raw: &str) -> Result<TranslExcept> {
    let mut parsed = TranslExcept::default();
    for entry in raw.split("),(") {
        let entry = entry.trim_start_matches('(').trim_end_matches(')');
        let Some((position, amino_acid)) = entry.rsplit_once(",aa:") else {
            bail!("transl_except entry {entry:?} has no `aa:` part");
        };
        let Some(position) = position.strip_prefix("pos:") else {
            bail!("transl_except entry {entry:?} has no `pos:` part");
        };
        if amino_acid != "Sec" {
            parsed.other += 1;
            continue;
        }
        let lowest = position
            .split(|c: char| !c.is_ascii_digit())
            .filter(|digits| !digits.is_empty())
            .map(str::parse::<u64>)
            .collect::<Result<Vec<_>, _>>()
            .with_context(|| format!("transl_except position {position:?}"))?
            .into_iter()
            .min();
        let Some(lowest) = lowest.filter(|&base| base > 0) else {
            bail!("transl_except entry {entry:?} has no 1-based position");
        };
        parsed.selenocysteine.push(lowest - 1);
    }
    Ok(parsed)
}

/// Biotype of a RefSeq transcript.
///
/// A CDS makes it protein-coding, whatever its row type or prefix. Without
/// one, the specific RNA row types map to their biotype; the generic rows
/// (`transcript`, `primary_transcript`, `ncRNA`, a CDS-less `mRNA`) take the
/// gene's `gene_biotype` when it names a non-coding class, else
/// [`Biotype::NonCodingRna`]. Other row types (`rRNA`, `Y_RNA`, ...) are kept
/// verbatim.
///
/// # Arguments
///
/// * `feature_type` - GFF3 column 3 of the transcript row
/// * `has_cds` - Whether the transcript has CDS rows
/// * `gene_biotype` - The parent gene's `gene_biotype`, if any
pub(super) fn refseq_biotype(
    feature_type: &str,
    has_cds: bool,
    gene_biotype: Option<&str>,
) -> Biotype {
    if has_cds {
        return Biotype::ProteinCoding;
    }
    match feature_type {
        "lnc_RNA" => Biotype::LncRna,
        "antisense_RNA" => Biotype::AntisenseRna,
        "snoRNA" => Biotype::SnoRna,
        "snRNA" => Biotype::SnRna,
        "RNase_MRP_RNA" => Biotype::RnaseMrpRna,
        "telomerase_RNA" => Biotype::TelomeraseRna,
        "vault_RNA" => Biotype::VaultRna,
        "transcript" | "primary_transcript" | "ncRNA" | "mRNA" => match gene_biotype {
            Some(label) if !label.is_empty() && label != "protein_coding" => {
                Biotype::from_label(label)
            }
            Some(_) | None => Biotype::NonCodingRna,
        },
        other => Biotype::from_label(other),
    }
}

/// `true` when no two CDS segments share an exon. RefSeq's ribosomal
/// slippage CDSs split one exon into two overlapping segments, which the
/// store cannot index.
///
/// # Arguments
///
/// * `segments` - CDS segments in transcript order
pub(super) fn one_cds_per_exon(segments: &[CdsSegment]) -> bool {
    segments
        .windows(2)
        .all(|pair| pair[0].exon_index != pair[1].exon_index)
}

/// `true` when every interior CDS phase agrees with the CDS length before
/// it, counting the first segment's pad of `(3 - phase) % 3` bases.
///
/// # Arguments
///
/// * `segments` - CDS segments in transcript order
pub(super) fn phases_consistent(segments: &[CdsSegment]) -> bool {
    let Some(first) = segments.first() else {
        return true;
    };
    let mut cumulative = u64::from((3 - first.phase % 3) % 3);
    for (index, segment) in segments.iter().enumerate() {
        if index > 0 && u64::from(segment.phase) != (3 - cumulative % 3) % 3 {
            return false;
        }
        cumulative += segment.genomic_end - segment.genomic_start;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exon(start: u64, end: u64) -> Exon {
        Exon {
            exon_number: 1,
            genomic_start: start,
            genomic_end: end,
        }
    }

    fn segment(exon_index: u16, start: u64, end: u64, phase: u8) -> CdsSegment {
        CdsSegment {
            exon_index,
            genomic_start: start,
            genomic_end: end,
            phase,
        }
    }

    #[test]
    fn gap_ops_walk_up_on_plus_and_down_on_minus() {
        let plus = parse_gap_ops("M10 I2 M5 D3 M4", 100, 122, Strand::Plus).expect("plus");
        assert_eq!(
            plus,
            vec![
                GapOp::Insertion { at: 110 },
                GapOp::Deletion {
                    start: 115,
                    end: 118
                },
            ]
        );
        // NM_001371720.2: row 155190380..155190651 (1-based), M163 D1 M108;
        // the exons are 155190380-155190487 and 155190489-155190651.
        let minus =
            parse_gap_ops("M163 D1 M108", 155_190_379, 155_190_651, Strand::Minus).expect("minus");
        assert_eq!(
            minus,
            vec![GapOp::Deletion {
                start: 155_190_487,
                end: 155_190_488
            }]
        );
        let insertion = parse_gap_ops("M5 I1 M5", 0, 10, Strand::Minus).expect("minus I");
        assert_eq!(insertion, vec![GapOp::Insertion { at: 5 }]);
    }

    #[test]
    fn utr_gap_breaks_bound_both_sides_of_each_gap() {
        // Two exons 100..200 and 300..400, CDS 150..350.
        let exons = [
            Exon {
                exon_number: 1,
                genomic_start: 100,
                genomic_end: 200,
            },
            Exon {
                exon_number: 2,
                genomic_start: 300,
                genomic_end: 400,
            },
        ];
        let cds = Some((150, 350));
        let gap = |shifted_from, crossed_through| UtrGapBreak {
            shifted_from,
            crossed_through,
        };
        let ops = [
            GapOp::Insertion { at: 100 },
            GapOp::Insertion { at: 150 },
            GapOp::Insertion { at: 200 },
            GapOp::Deletion {
                start: 380,
                end: 382,
            },
        ];
        // Before the first base; right at the CDS start (cDNA 51); the CDS
        // insertion is no UTR break; a 2-base 3'UTR deletion at cDNA 181.
        assert_eq!(
            utr_gap_breaks(&ops, cds, &exons, true),
            vec![gap(1, 0), gap(51, 50), gap(181, 182)]
        );
        let minus: Vec<Exon> = exons.iter().rev().cloned().collect();
        assert_eq!(cdna_at_or_after(&minus, false, 399), Some(1));
        assert_eq!(cdna_at_or_after(&minus, false, 250), Some(101));
        // Minus strand, insertion at the CDS end (its 5' start): cDNA 51.
        let minus_cds = [GapOp::Insertion { at: 350 }];
        assert_eq!(
            utr_gap_breaks(&minus_cds, cds, &minus, false),
            vec![gap(51, 50)]
        );
    }

    #[test]
    fn gap_ops_reject_malformed_input() {
        assert!(parse_gap_ops("M5 X1", 0, 10, Strand::Plus).is_err());
        assert!(parse_gap_ops("Mx", 0, 10, Strand::Plus).is_err());
        assert!(parse_gap_ops("M11", 0, 10, Strand::Plus).is_err());
    }

    #[test]
    fn insertions_move_numbering_only_strictly_inside_the_cds() {
        let exons = [exon(100, 200)];
        let inside = [GapOp::Insertion { at: 150 }];
        let at_edge = [GapOp::Insertion { at: 120 }];
        assert!(gap_moves_cds_numbering(&inside, 120, 180, &exons));
        assert!(!gap_moves_cds_numbering(&at_edge, 120, 180, &exons));
    }

    #[test]
    fn deletions_move_numbering_unless_they_are_exon_splits() {
        let split = [exon(100, 150), exon(151, 200)];
        let deletion = [GapOp::Deletion {
            start: 150,
            end: 151,
        }];
        assert!(!gap_moves_cds_numbering(&deletion, 120, 180, &split));
        let whole = [exon(100, 200)];
        assert!(gap_moves_cds_numbering(&deletion, 120, 180, &whole));
        assert!(!gap_moves_cds_numbering(&deletion, 160, 180, &whole));
    }

    #[test]
    fn transl_except_reads_sec_on_both_strands_and_counts_others() {
        let parsed = parse_transl_except(
            "(pos:25802093..25802095,aa:Sec),(pos:complement(49358132..49358134),aa:Sec),\
             (pos:1000..1002,aa:Met)",
        )
        .expect("parse");
        assert_eq!(parsed.selenocysteine, vec![25_802_092, 49_358_131]);
        assert_eq!(parsed.other, 1);
        let split = parse_transl_except("(pos:join(10..11,20..20),aa:Sec)").expect("join");
        assert_eq!(split.selenocysteine, vec![9]);
    }

    #[test]
    fn transl_except_rejects_malformed_entries() {
        assert!(parse_transl_except("(pos:10..12)").is_err());
        assert!(parse_transl_except("(10..12,aa:Sec)").is_err());
        assert!(parse_transl_except("(pos:x,aa:Sec)").is_err());
    }

    #[test]
    fn refseq_biotype_follows_cds_then_row_type_then_gene() {
        assert_eq!(
            refseq_biotype("transcript", true, None),
            Biotype::ProteinCoding
        );
        assert_eq!(refseq_biotype("lnc_RNA", false, None), Biotype::LncRna);
        assert_eq!(
            refseq_biotype("transcript", false, Some("protein_coding")),
            Biotype::NonCodingRna
        );
        assert_eq!(
            refseq_biotype("primary_transcript", false, Some("miRNA")),
            Biotype::Other("miRNA".to_string())
        );
        assert_eq!(refseq_biotype("mRNA", false, None), Biotype::NonCodingRna);
        assert_eq!(
            refseq_biotype("rRNA", false, Some("rRNA")),
            Biotype::Other("rRNA".to_string())
        );
    }

    #[test]
    fn slippage_cds_shares_an_exon() {
        let slippage = [segment(0, 100, 150, 0), segment(0, 149, 200, 2)];
        assert!(!one_cds_per_exon(&slippage));
        let normal = [segment(0, 100, 150, 0), segment(1, 300, 400, 1)];
        assert!(one_cds_per_exon(&normal));
    }

    #[test]
    fn phases_follow_the_padded_cumulative_length() {
        // 50 bases at phase 0 leave 2 bases over: the next phase is 1.
        assert!(phases_consistent(&[
            segment(0, 100, 150, 0),
            segment(1, 300, 400, 1)
        ]));
        assert!(!phases_consistent(&[
            segment(0, 100, 150, 0),
            segment(1, 300, 400, 2)
        ]));
        // Phase 2 pads one base: 1 + 50 = 51, so the next phase is 0.
        assert!(phases_consistent(&[
            segment(0, 100, 150, 2),
            segment(1, 300, 400, 0)
        ]));
        assert!(phases_consistent(&[]));
    }
}
