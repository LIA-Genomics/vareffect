//! VCF-style left-alignment (5' shift-then-trim) of indel alleles.
//!
//! Implements the Tan et al. 2015 normalization algorithm used by `vt normalize`
//! and `bcftools norm`: it produces the unique **leftmost, parsimonious**
//! representation of a variant. This is the VCF-canonical direction — the
//! opposite axis to the HGVS 3' (rightward) shift in [`crate::normalize`]. The
//! concordance target is `bcftools norm` / biocommons `hgvs_to_vcf`
//! (`shuffle_direction=5`) / gnomAD.
//!
//! Transcript- and strand-free: it consults only the reference genome, so it
//! takes a bare [`FastaReader`] rather than the full [`crate::VarEffect`].
//! [`crate::VarEffect::left_align_indel`] is a thin `&self` wrapper over this
//! function.

use crate::error::VarEffectError;
use crate::fasta::FastaReader;

/// Left-align a VCF-style variant to the leftmost equivalent position.
///
/// Implements the Tan et al. 2015 normalization algorithm used by
/// `vt normalize` and `bcftools norm`: shift-then-trim produces the unique
/// left-aligned parsimonious representation.
///
/// SNVs and MNVs pass through unchanged (the shift loop's rightmost-base
/// comparison exits immediately when alleles differ). Complex inputs like
/// `ref=ACT, alt=AT` are handled naturally: the matching rightmost `T`
/// triggers the shift, exposing the underlying deletion.
///
/// # Arguments
///
/// * `fasta` - Reference genome, used to prepend the base immediately 5' of
///   the event during the shift.
/// * `chrom` - Chromosome name accepted by the FASTA reader (UCSC, bare, or
///   RefSeq naming -- the reader's alias table handles translation).
/// * `pos_1based` - 1-based VCF position of the variant.
/// * `ref_allele` - Reference allele (uppercase ACGTN, no `"-"` placeholders).
/// * `alt_allele` - Alternate allele (uppercase ACGTN, no `"-"` placeholders).
///
/// # Returns
///
/// * `Ok(None)` - No normalization needed (SNV, MNV, or already leftmost).
/// * `Ok(Some((pos, ref, alt)))` - Left-aligned representation with new
///   1-based position and parsimonious alleles.
///
/// # Errors
///
/// * [`VarEffectError::ChromNotFound`] - Chromosome not in the FASTA index.
/// * [`VarEffectError::CoordinateOutOfRange`] - Position exceeds chromosome
///   length during left-extension.
/// * [`VarEffectError::InvalidAllele`] - Allele bytes are not valid UTF-8
///   (should not occur with valid genomic input).
pub(crate) fn left_align_indel(
    fasta: &FastaReader,
    chrom: &str,
    pos_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<Option<(u64, String, String)>, VarEffectError> {
    let orig_pos = pos_1based;
    let orig_ref = ref_allele;
    let orig_alt = alt_allele;

    let mut pos = pos_1based;
    let mut r: Vec<u8> = ref_allele.as_bytes().to_vec();
    let mut a: Vec<u8> = alt_allele.as_bytes().to_vec();

    // Step 1: Shift loop (Tan et al. 2015).
    //
    // Compare rightmost bases. While they match, drop them and shift left.
    // When trimming empties an allele, prepend the reference base immediately
    // 5' to BOTH alleles and step the position back by one. Extending both
    // alleles preserves the indel (the ref/alt length difference); extending
    // only the emptied allele would collapse an insertion/deletion into a
    // substitution.
    //
    // Coordinate convention: `fetch_base` takes a 0-based position. After
    // `pos -= 1`, `pos` is the new 1-based position, so the 0-based index for
    // the base AT `pos` is `pos - 1`.
    loop {
        if pos <= 1 || r.is_empty() || a.is_empty() {
            break;
        }
        let trimmed = if r.last() == a.last() {
            r.pop();
            a.pop();
            true
        } else {
            false
        };
        if r.is_empty() || a.is_empty() {
            pos -= 1;
            let prepend = fasta.fetch_base(chrom, pos - 1)?;
            r.insert(0, prepend);
            a.insert(0, prepend);
        } else if !trimmed {
            break;
        }
    }

    // Step 2: Left-trim for parsimony.
    //
    // Skips shared leading bases via an index (no O(n) `remove(0)`), keeping
    // at least one base per allele (the VCF anchor requirement).
    let mut prefix_skip = 0usize;
    while r.len() - prefix_skip > 1 && a.len() - prefix_skip > 1 && r[prefix_skip] == a[prefix_skip]
    {
        prefix_skip += 1;
        pos += 1;
    }

    // Step 3: Return None if nothing changed, Some if normalized.
    //
    // Alleles contain only ASCII bytes (ACGTN from the FASTA reader or the
    // original input validated upstream). `from_utf8` cannot fail on valid
    // ASCII but we propagate rather than panic.
    let new_ref = std::str::from_utf8(&r[prefix_skip..])
        .map_err(|_| VarEffectError::InvalidAllele)?
        .to_string();
    let new_alt = std::str::from_utf8(&a[prefix_skip..])
        .map_err(|_| VarEffectError::InvalidAllele)?
        .to_string();

    if pos == orig_pos && new_ref == orig_ref && new_alt == orig_alt {
        Ok(None)
    } else {
        Ok(Some((pos, new_ref, new_alt)))
    }
}

#[cfg(test)]
mod left_align_unit_tests {
    //! CI-runnable (no 3.1 GB FASTA) exact-output tests for `left_align_indel`,
    //! covering the shapes whose ref/alt share a rightmost base — duplications
    //! and homopolymer/tandem indels — which is where the per-allele-extend bug
    //! collapsed indels into substitutions. Run against a tiny synthetic genome.

    use super::*;
    use crate::fasta::write_genome_binary;
    use tempfile::TempDir;

    /// Synthetic contig "1" (queried as "chr1"):
    ///   1-based pos: 1 2 3 4 5 6 7 8 9 10
    ///   base:        A A T G G G G T A A
    /// The `GGGG` run (pos 4-7) anchored by `T` at pos 3 exercises left-shifting.
    const CONTIG: &[u8] = b"AATGGGGTAA";

    /// Build a `FastaReader` over a synthetic in-memory genome (temp `.bin` +
    /// `.bin.idx`) — `left_align_indel` only touches the FASTA.
    fn fasta_with(contigs: &[(&str, &[u8])]) -> (TempDir, FastaReader) {
        let tmp = TempDir::new().expect("tempdir");
        let bin = tmp.path().join("g.bin");
        let idx = tmp.path().join("g.bin.idx");
        write_genome_binary(contigs, "test", &bin, &idx).expect("write synthetic genome");
        let fasta = FastaReader::open(&bin).expect("open synthetic genome");
        (tmp, fasta)
    }

    /// A right-shifted 1 bp duplication inside the `GGGG` run must left-align to
    /// the `T` anchor at pos 3 and REMAIN an insertion (`T>TG`). Pre-fix the
    /// buggy per-allele extension collapsed this to a substitution.
    #[test]
    fn duplication_left_aligns_preserving_insertion() {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        let result = left_align_indel(&fasta, "chr1", 6, "G", "GG")
            .expect("left_align_indel should not error");
        assert_eq!(result, Some((3, "T".to_string(), "TG".to_string())));
    }

    /// The symmetric 1 bp tandem deletion must left-align to the same anchor and
    /// REMAIN a deletion (`TG>T`).
    #[test]
    fn tandem_deletion_left_aligns_preserving_deletion() {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        let result = left_align_indel(&fasta, "chr1", 6, "GG", "G")
            .expect("left_align_indel should not error");
        assert_eq!(result, Some((3, "TG".to_string(), "T".to_string())));
    }

    /// Feeding the already-left-aligned form back in is a no-op (idempotent).
    #[test]
    fn left_aligned_insertion_is_idempotent() {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        let result = left_align_indel(&fasta, "chr1", 3, "T", "TG")
            .expect("left_align_indel should not error");
        assert_eq!(result, None);
    }

    /// An SNV (differing rightmost bases) never enters the shift loop → `None`.
    #[test]
    fn snv_passes_through() {
        let (_tmp, fasta) = fasta_with(&[("1", CONTIG)]);
        let result = left_align_indel(&fasta, "chr1", 8, "T", "A")
            .expect("left_align_indel should not error");
        assert_eq!(result, None);
    }
}
