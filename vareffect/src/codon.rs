//! Codon translation, DNA complement, and VEP-style display formatting.
//!
//! All functions operate on **uppercase ASCII bytes** (`b'A'`, `b'C'`, `b'G'`,
//! `b'T'`). The codon lookup tables use an internal 6-bit encoding but callers
//! never see it — pass raw ASCII in, get ASCII out.
//!
//! # Genetic codes
//!
//! Two translation tables are provided:
//!
//! - **Standard (NCBI table 1)** — used for all autosomal and sex-chromosome
//!   transcripts. [`translate_codon`] uses this table.
//! - **Vertebrate mitochondrial (NCBI table 2)** — used for chrM transcripts.
//!   Four codons differ: `TGA→W`, `AGA→*`, `AGG→*`, `ATA→M`.
//!
//! [`translate_codon_for_transcript`] dispatches to the correct table based on
//! an `is_mitochondrial` flag (derived from `transcript.chrom == "chrM"`).

// ---------------------------------------------------------------------------
// Internal encoding
// ---------------------------------------------------------------------------

/// Map an ASCII DNA base to a 2-bit index: A=0, C=1, G=2, T=3.
/// Returns `None` for any non-ACGT byte (N, ambiguity codes, lowercase).
fn base_to_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Pack a 3-base codon into a 6-bit index (0–63) for table lookup.
/// Returns `None` if any base is not in {A, C, G, T}.
fn codon_to_index(codon: &[u8; 3]) -> Option<usize> {
    let a = base_to_index(codon[0])?;
    let b = base_to_index(codon[1])?;
    let c = base_to_index(codon[2])?;
    Some(a * 16 + b * 4 + c)
}

// ---------------------------------------------------------------------------
// Standard genetic code (NCBI translation table 1)
// ---------------------------------------------------------------------------

/// Standard genetic code. Index = `base1*16 + base2*4 + base3` where
/// A=0, C=1, G=2, T=3. Each entry is the one-letter amino acid code;
/// `b'*'` marks stop codons.
///
/// Layout (row = first two bases, column = third base):
/// ```text
///        A    C    G    T
/// AA:    K    N    K    N
/// AC:    T    T    T    T
/// AG:    R    S    R    S
/// AT:    I    I    M    I
/// CA:    Q    H    Q    H
/// CC:    P    P    P    P
/// CG:    R    R    R    R
/// CT:    L    L    L    L
/// GA:    E    D    E    D
/// GC:    A    A    A    A
/// GG:    G    G    G    G
/// GT:    V    V    V    V
/// TA:    *    Y    *    Y
/// TC:    S    S    S    S
/// TG:    *    C    W    C
/// TT:    L    F    L    F
/// ```
static STANDARD_TABLE: [u8; 64] = [
    // AA AC AG AT  (first base = A)
    b'K', b'N', b'K', b'N', // AA{A,C,G,T}
    b'T', b'T', b'T', b'T', // AC{A,C,G,T}
    b'R', b'S', b'R', b'S', // AG{A,C,G,T}
    b'I', b'I', b'M', b'I', // AT{A,C,G,T}
    // CA CC CG CT  (first base = C)
    b'Q', b'H', b'Q', b'H', // CA{A,C,G,T}
    b'P', b'P', b'P', b'P', // CC{A,C,G,T}
    b'R', b'R', b'R', b'R', // CG{A,C,G,T}
    b'L', b'L', b'L', b'L', // CT{A,C,G,T}
    // GA GC GG GT  (first base = G)
    b'E', b'D', b'E', b'D', // GA{A,C,G,T}
    b'A', b'A', b'A', b'A', // GC{A,C,G,T}
    b'G', b'G', b'G', b'G', // GG{A,C,G,T}
    b'V', b'V', b'V', b'V', // GT{A,C,G,T}
    // TA TC TG TT  (first base = T)
    b'*', b'Y', b'*', b'Y', // TA{A,C,G,T}
    b'S', b'S', b'S', b'S', // TC{A,C,G,T}
    b'*', b'C', b'W', b'C', // TG{A,C,G,T}
    b'L', b'F', b'L', b'F', // TT{A,C,G,T}
];

// ---------------------------------------------------------------------------
// Vertebrate mitochondrial genetic code (NCBI translation table 2)
// ---------------------------------------------------------------------------

/// Vertebrate mitochondrial code. Differs from the standard table at 4 codons:
/// - `TGA` → `W` (Trp, not stop)
/// - `AGA` → `*` (stop, not Arg)
/// - `AGG` → `*` (stop, not Arg)
/// - `ATA` → `M` (Met, not Ile)
static MITO_TABLE: [u8; 64] = [
    // AA AC AG AT  (first base = A)
    b'K', b'N', b'K', b'N', // AA{A,C,G,T}
    b'T', b'T', b'T', b'T', // AC{A,C,G,T}
    b'*', b'S', b'*', b'S', // AG{A,C,G,T}  ← AGA=*, AGG=*
    b'M', b'I', b'M', b'I', // AT{A,C,G,T}  ← ATA=M
    // CA CC CG CT  (first base = C)
    b'Q', b'H', b'Q', b'H', // CA{A,C,G,T}
    b'P', b'P', b'P', b'P', // CC{A,C,G,T}
    b'R', b'R', b'R', b'R', // CG{A,C,G,T}
    b'L', b'L', b'L', b'L', // CT{A,C,G,T}
    // GA GC GG GT  (first base = G)
    b'E', b'D', b'E', b'D', // GA{A,C,G,T}
    b'A', b'A', b'A', b'A', // GC{A,C,G,T}
    b'G', b'G', b'G', b'G', // GG{A,C,G,T}
    b'V', b'V', b'V', b'V', // GT{A,C,G,T}
    // TA TC TG TT  (first base = T)
    b'*', b'Y', b'*', b'Y', // TA{A,C,G,T}
    b'S', b'S', b'S', b'S', // TC{A,C,G,T}
    b'W', b'C', b'W', b'C', // TG{A,C,G,T}  ← TGA=W (not stop), TGG=W
    b'L', b'F', b'L', b'F', // TT{A,C,G,T}
];

// ---------------------------------------------------------------------------
// Public translation API
// ---------------------------------------------------------------------------

/// Translate a 3-base codon to a single amino acid character using the
/// standard genetic code (NCBI table 1).
///
/// # Arguments
///
/// * `codon` — Three uppercase ASCII DNA bases (e.g., `b"ATG"`)
///
/// # Returns
///
/// One-letter amino acid code as a `u8`:
/// - Standard amino acids: `b'A'`..`b'Y'`
/// - Stop codons (TAA, TAG, TGA): `b'*'`
/// - Ambiguous codons (containing N or other non-ACGT bases): `b'X'`
///
/// # Examples
///
/// ```
/// use vareffect::codon::translate_codon;
/// assert_eq!(translate_codon(b"ATG"), b'M');
/// assert_eq!(translate_codon(b"TAA"), b'*');
/// assert_eq!(translate_codon(b"NNN"), b'X');
/// ```
pub fn translate_codon(codon: &[u8; 3]) -> u8 {
    match codon_to_index(codon) {
        Some(idx) => STANDARD_TABLE[idx],
        None => b'X',
    }
}

/// Translate a 3-base codon using the vertebrate mitochondrial genetic code
/// (NCBI table 2).
///
/// Differs from [`translate_codon`] at four codons:
/// `TGA→W`, `AGA→*`, `AGG→*`, `ATA→M`.
///
/// # Examples
///
/// ```
/// use vareffect::codon::translate_codon_mito;
/// assert_eq!(translate_codon_mito(b"TGA"), b'W');
/// assert_eq!(translate_codon_mito(b"AGA"), b'*');
/// ```
pub fn translate_codon_mito(codon: &[u8; 3]) -> u8 {
    match codon_to_index(codon) {
        Some(idx) => MITO_TABLE[idx],
        None => b'X',
    }
}

/// Translate a codon using the appropriate genetic code for a transcript.
///
/// Dispatches to the vertebrate mitochondrial code (NCBI table 2) when
/// `is_mitochondrial` is true, otherwise uses the standard code (table 1).
///
/// # Arguments
///
/// * `codon` — Three uppercase ASCII DNA bases
/// * `is_mitochondrial` — `true` for chrM transcripts
///
/// # Examples
///
/// ```
/// use vareffect::codon::translate_codon_for_transcript;
/// // Standard code: TGA is stop
/// assert_eq!(translate_codon_for_transcript(b"TGA", false), b'*');
/// // Mitochondrial code: TGA is Trp
/// assert_eq!(translate_codon_for_transcript(b"TGA", true), b'W');
/// ```
pub fn translate_codon_for_transcript(codon: &[u8; 3], is_mitochondrial: bool) -> u8 {
    if is_mitochondrial {
        translate_codon_mito(codon)
    } else {
        translate_codon(codon)
    }
}

// ---------------------------------------------------------------------------
// DNA complement
// ---------------------------------------------------------------------------

/// Complement a single DNA base. `A↔T`, `C↔G`. Non-ACGT bytes (e.g., `N`)
/// pass through unchanged.
///
/// # Examples
///
/// ```
/// use vareffect::codon::complement;
/// assert_eq!(complement(b'A'), b'T');
/// assert_eq!(complement(b'N'), b'N');
/// ```
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Complement a DNA sequence in place. Each base is replaced with its
/// Watson-Crick complement; non-ACGT bytes are left unchanged.
pub fn complement_in_place(seq: &mut [u8]) {
    for base in seq.iter_mut() {
        *base = complement(*base);
    }
}

/// Return the reverse complement of a DNA sequence.
///
/// # Examples
///
/// ```
/// use vareffect::codon::reverse_complement;
/// assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

// ---------------------------------------------------------------------------
// Amino acid display helpers
// ---------------------------------------------------------------------------

/// Convert a one-letter amino acid code to its three-letter abbreviation.
///
/// Returns `"Ter"` for stop (`*`) and the standard three-letter code for
/// the 20 standard amino acids. Any other byte (including the IUPAC `X`
/// placeholder and invalid input) maps to `"Xaa"`, the IUPAC "any / unknown
/// amino acid" symbol, so this function never panics.
///
/// # Examples
///
/// ```
/// use vareffect::codon::aa_three_letter;
/// assert_eq!(aa_three_letter(b'M'), "Met");
/// assert_eq!(aa_three_letter(b'*'), "Ter");
/// assert_eq!(aa_three_letter(b'X'), "Xaa");
/// ```
pub fn aa_three_letter(one_letter: u8) -> &'static str {
    match one_letter {
        b'A' => "Ala",
        b'C' => "Cys",
        b'D' => "Asp",
        b'E' => "Glu",
        b'F' => "Phe",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'K' => "Lys",
        b'L' => "Leu",
        b'M' => "Met",
        b'N' => "Asn",
        b'P' => "Pro",
        b'Q' => "Gln",
        b'R' => "Arg",
        b'S' => "Ser",
        b'T' => "Thr",
        b'V' => "Val",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'*' => "Ter",
        // `X` is IUPAC for "any amino acid"; any other byte (including
        // garbage input) falls through to the same sentinel — a public
        // helper should never panic on caller-supplied data.
        _ => "Xaa",
    }
}

// ---------------------------------------------------------------------------
// VEP-style display formatting
// ---------------------------------------------------------------------------

/// Format ref/alt codons with VEP's capitalization convention.
///
/// The changed base is uppercase; unchanged bases are lowercase. The two
/// codons are separated by `/`.
///
/// # Arguments
///
/// * `ref_codon` — Reference codon (3 uppercase ASCII bases)
/// * `alt_codon` — Alternate codon (3 uppercase ASCII bases)
/// * `changed_pos` — 0-based position of the changed base (0, 1, or 2)
///
/// # Examples
///
/// ```
/// use vareffect::codon::format_codons;
/// assert_eq!(format_codons(b"CGT", b"TGT", 0), "Cgt/Tgt");
/// assert_eq!(format_codons(b"CGT", b"CAT", 1), "cGt/cAt");
/// assert_eq!(format_codons(b"CGT", b"CGA", 2), "cgT/cgA");
/// ```
pub fn format_codons(ref_codon: &[u8; 3], alt_codon: &[u8; 3], changed_pos: u8) -> String {
    debug_assert!(changed_pos < 3, "codon position must be 0, 1, or 2");
    let mut result = String::with_capacity(7); // "xxx/xxx"
    for i in 0..3u8 {
        if i == changed_pos {
            result.push(ref_codon[i as usize] as char);
        } else {
            result.push((ref_codon[i as usize] as char).to_ascii_lowercase());
        }
    }
    result.push('/');
    for i in 0..3u8 {
        if i == changed_pos {
            result.push(alt_codon[i as usize] as char);
        } else {
            result.push((alt_codon[i as usize] as char).to_ascii_lowercase());
        }
    }
    result
}

/// Format amino acid change for VEP display.
///
/// - Synonymous (ref == alt): single letter, e.g. `"R"`
/// - Non-synonymous: `"R/W"` format
/// - Stop gained: `"R/*"`
///
/// # Examples
///
/// ```
/// use vareffect::codon::format_amino_acids;
/// assert_eq!(format_amino_acids(b'R', b'W'), "R/W");
/// assert_eq!(format_amino_acids(b'R', b'R'), "R");
/// assert_eq!(format_amino_acids(b'R', b'*'), "R/*");
/// ```
pub fn format_amino_acids(ref_aa: u8, alt_aa: u8) -> String {
    if ref_aa == alt_aa {
        String::from(ref_aa as char)
    } else {
        format!("{}/{}", ref_aa as char, alt_aa as char)
    }
}

// ---------------------------------------------------------------------------
// Indel helpers
// ---------------------------------------------------------------------------

/// Translate a DNA sequence to amino acids, codon by codon.
///
/// Uses the appropriate genetic code based on `is_mitochondrial`.
///
/// # Arguments
///
/// * `seq` — Uppercase ASCII DNA bytes. Length must be divisible by 3.
/// * `is_mitochondrial` — `true` for chrM transcripts (NCBI table 2)
///
/// # Returns
///
/// `Vec<u8>` of one-letter amino acid codes (same encoding as
/// [`translate_codon`]).
///
/// # Errors
///
/// Returns [`crate::VarEffectError::Malformed`] if `seq.len() % 3 != 0`.
pub fn translate_sequence(
    seq: &[u8],
    is_mitochondrial: bool,
) -> Result<Vec<u8>, crate::VarEffectError> {
    if !seq.len().is_multiple_of(3) {
        return Err(crate::VarEffectError::Malformed(format!(
            "translate_sequence: sequence length {} is not divisible by 3",
            seq.len(),
        )));
    }
    let mut aas = Vec::with_capacity(seq.len() / 3);
    for codon_bytes in seq.chunks_exact(3) {
        let codon: &[u8; 3] = codon_bytes
            .try_into()
            .expect("chunks_exact(3) always yields a 3-byte slice");
        aas.push(translate_codon_for_transcript(codon, is_mitochondrial));
    }
    Ok(aas)
}

/// Format ref/alt codon sequences for an indel with VEP's capitalisation
/// convention.
///
/// Bases in the ref that are deleted (positions `[changed_start, changed_end)`
/// within `ref_seq`) are uppercase; all other ref bases are lowercase. In the
/// alt, bases that were inserted (positions `[changed_start,
/// changed_start + inserted_len)`) are uppercase; flanking bases are lowercase.
///
/// The two sequences are separated by `/`.
///
/// # Arguments
///
/// * `ref_seq` — Codon-aligned reference CDS bases
/// * `alt_seq` — Codon-aligned alternate CDS bases (after deletion/insertion)
/// * `changed_start` — 0-based index in `ref_seq` where the change begins
/// * `changed_end` — 0-based exclusive end in `ref_seq` where the change ends
///   (for deletions: the range of deleted bases; for insertions: typically
///   `changed_start` since nothing is deleted in the ref)
///
/// # Examples
///
/// ```
/// use vareffect::codon::format_codons_indel;
/// // 3bp deletion at positions 3-5: "atgGAC/atg"
/// assert_eq!(
///     format_codons_indel(b"ATGGAC", b"ATG", 3, 6),
///     "atgGAC/atg"
/// );
/// // 3bp insertion at position 1: "agc/aGACgc"
/// assert_eq!(
///     format_codons_indel(b"AGC", b"AGACGC", 1, 1),
///     "agc/aGACgc"
/// );
/// ```
pub fn format_codons_indel(
    ref_seq: &[u8],
    alt_seq: &[u8],
    changed_start: usize,
    changed_end: usize,
) -> String {
    let mut result = String::with_capacity(ref_seq.len() + alt_seq.len() + 2);

    // Ref side: uppercase for deleted bases [changed_start, changed_end)
    if ref_seq.is_empty() {
        // Pure insertion with no ref codon context — VEP uses "-".
        result.push('-');
    } else {
        for (i, &b) in ref_seq.iter().enumerate() {
            if i >= changed_start && i < changed_end {
                result.push(b as char);
            } else {
                result.push((b as char).to_ascii_lowercase());
            }
        }
    }
    result.push('/');

    if alt_seq.is_empty() {
        // Complete deletion — VEP uses "-" for the empty alt side.
        result.push('-');
    } else {
        // Alt side: bases corresponding to the unchanged prefix and suffix are
        // lowercase; inserted bases (those that don't correspond to the ref
        // flanking regions) are uppercase.
        let prefix_len = changed_start;
        let suffix_len = ref_seq.len() - changed_end;
        let alt_suffix_start = alt_seq.len().saturating_sub(suffix_len);

        for (i, &b) in alt_seq.iter().enumerate() {
            if i < prefix_len || i >= alt_suffix_start {
                result.push((b as char).to_ascii_lowercase());
            } else {
                result.push(b as char);
            }
        }
    }
    result
}

/// Format amino acid change for indels.
///
/// Shows the full ref and alt amino acid sequences separated by `/`.
/// Single-letter codes are used. If ref and alt are identical (synonymous
/// indel at the protein level), returns a single copy. An empty side is
/// represented as `-` (VEP convention).
///
/// # Examples
///
/// ```
/// use vareffect::codon::format_amino_acids_indel;
/// assert_eq!(format_amino_acids_indel(b"RR", b"R"), "RR/R");
/// assert_eq!(format_amino_acids_indel(b"R", b"RDR"), "R/RDR");
/// assert_eq!(format_amino_acids_indel(b"M", b""), "M/-");
/// assert_eq!(format_amino_acids_indel(b"", b"X"), "-/X");
/// ```
pub fn format_amino_acids_indel(ref_aas: &[u8], alt_aas: &[u8]) -> String {
    if ref_aas == alt_aas {
        ref_aas.iter().map(|&b| b as char).collect()
    } else {
        let mut result = String::with_capacity(ref_aas.len() + alt_aas.len() + 2);
        if ref_aas.is_empty() {
            result.push('-');
        } else {
            for &b in ref_aas {
                result.push(b as char);
            }
        }
        result.push('/');
        if alt_aas.is_empty() {
            result.push('-');
        } else {
            for &b in alt_aas {
                result.push(b as char);
            }
        }
        result
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify every one of the 64 standard-code codons against the NCBI
    /// translation table 1 reference.
    #[test]
    fn translate_all_64_codons() {
        // (codon_bytes, expected_aa)
        let expected: &[(&[u8; 3], u8)] = &[
            (b"TTT", b'F'),
            (b"TTC", b'F'),
            (b"TTA", b'L'),
            (b"TTG", b'L'),
            (b"TCT", b'S'),
            (b"TCC", b'S'),
            (b"TCA", b'S'),
            (b"TCG", b'S'),
            (b"TAT", b'Y'),
            (b"TAC", b'Y'),
            (b"TAA", b'*'),
            (b"TAG", b'*'),
            (b"TGT", b'C'),
            (b"TGC", b'C'),
            (b"TGA", b'*'),
            (b"TGG", b'W'),
            (b"CTT", b'L'),
            (b"CTC", b'L'),
            (b"CTA", b'L'),
            (b"CTG", b'L'),
            (b"CCT", b'P'),
            (b"CCC", b'P'),
            (b"CCA", b'P'),
            (b"CCG", b'P'),
            (b"CAT", b'H'),
            (b"CAC", b'H'),
            (b"CAA", b'Q'),
            (b"CAG", b'Q'),
            (b"CGT", b'R'),
            (b"CGC", b'R'),
            (b"CGA", b'R'),
            (b"CGG", b'R'),
            (b"ATT", b'I'),
            (b"ATC", b'I'),
            (b"ATA", b'I'),
            (b"ATG", b'M'),
            (b"ACT", b'T'),
            (b"ACC", b'T'),
            (b"ACA", b'T'),
            (b"ACG", b'T'),
            (b"AAT", b'N'),
            (b"AAC", b'N'),
            (b"AAA", b'K'),
            (b"AAG", b'K'),
            (b"AGT", b'S'),
            (b"AGC", b'S'),
            (b"AGA", b'R'),
            (b"AGG", b'R'),
            (b"GTT", b'V'),
            (b"GTC", b'V'),
            (b"GTA", b'V'),
            (b"GTG", b'V'),
            (b"GCT", b'A'),
            (b"GCC", b'A'),
            (b"GCA", b'A'),
            (b"GCG", b'A'),
            (b"GAT", b'D'),
            (b"GAC", b'D'),
            (b"GAA", b'E'),
            (b"GAG", b'E'),
            (b"GGT", b'G'),
            (b"GGC", b'G'),
            (b"GGA", b'G'),
            (b"GGG", b'G'),
        ];
        assert_eq!(expected.len(), 64);
        for &(codon, aa) in expected {
            assert_eq!(
                translate_codon(codon),
                aa,
                "codon {} should translate to {} but got {}",
                std::str::from_utf8(codon).unwrap(),
                aa as char,
                translate_codon(codon) as char,
            );
        }
    }

    /// The three standard stop codons should all translate to `*`.
    #[test]
    fn translate_stop_codons() {
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"TAG"), b'*');
        assert_eq!(translate_codon(b"TGA"), b'*');
    }

    /// Codons containing non-ACGT bases should translate to `X` (unknown).
    #[test]
    fn translate_ambiguous_codon() {
        assert_eq!(translate_codon(b"NNN"), b'X');
        assert_eq!(translate_codon(b"ANG"), b'X');
        assert_eq!(translate_codon(b"ATN"), b'X');
        // Lowercase is also non-ACGT in our convention (input must be uppercase)
        assert_eq!(translate_codon(b"atg"), b'X');
    }

    /// Verify the 4 codons that differ between the standard and mitochondrial
    /// genetic codes (NCBI table 2).
    #[test]
    fn translate_mitochondrial_differences() {
        // TGA: standard = stop, mito = Trp
        assert_eq!(translate_codon(b"TGA"), b'*');
        assert_eq!(translate_codon_mito(b"TGA"), b'W');

        // AGA: standard = Arg, mito = stop
        assert_eq!(translate_codon(b"AGA"), b'R');
        assert_eq!(translate_codon_mito(b"AGA"), b'*');

        // AGG: standard = Arg, mito = stop
        assert_eq!(translate_codon(b"AGG"), b'R');
        assert_eq!(translate_codon_mito(b"AGG"), b'*');

        // ATA: standard = Ile, mito = Met
        assert_eq!(translate_codon(b"ATA"), b'I');
        assert_eq!(translate_codon_mito(b"ATA"), b'M');

        // Verify dispatch helper
        assert_eq!(translate_codon_for_transcript(b"TGA", false), b'*');
        assert_eq!(translate_codon_for_transcript(b"TGA", true), b'W');
    }

    /// Watson-Crick complement: A↔T, C↔G, N passes through.
    #[test]
    fn complement_bases() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    /// Reverse complement of a 10-base sequence.
    #[test]
    fn reverse_complement_sequence() {
        let seq = b"ATCGATCGAT";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"ATCGATCGAT");

        let seq2 = b"AACCGGTT";
        let rc2 = reverse_complement(seq2);
        assert_eq!(rc2, b"AACCGGTT");

        // Asymmetric sequence
        let seq3 = b"AAACCCGGGT";
        let rc3 = reverse_complement(seq3);
        assert_eq!(rc3, b"ACCCGGGTTT");

        // complement_in_place
        let mut buf = b"ACGT".to_vec();
        complement_in_place(&mut buf);
        assert_eq!(buf, b"TGCA");
    }

    /// Verify all 20 standard amino acids + Ter + Xaa three-letter codes.
    #[test]
    fn aa_three_letter_all_20() {
        let cases: &[(u8, &str)] = &[
            (b'A', "Ala"),
            (b'C', "Cys"),
            (b'D', "Asp"),
            (b'E', "Glu"),
            (b'F', "Phe"),
            (b'G', "Gly"),
            (b'H', "His"),
            (b'I', "Ile"),
            (b'K', "Lys"),
            (b'L', "Leu"),
            (b'M', "Met"),
            (b'N', "Asn"),
            (b'P', "Pro"),
            (b'Q', "Gln"),
            (b'R', "Arg"),
            (b'S', "Ser"),
            (b'T', "Thr"),
            (b'V', "Val"),
            (b'W', "Trp"),
            (b'Y', "Tyr"),
            (b'*', "Ter"),
            (b'X', "Xaa"),
        ];
        for &(code, expected) in cases {
            assert_eq!(
                aa_three_letter(code),
                expected,
                "aa_three_letter({}) should be {}",
                code as char,
                expected,
            );
        }
    }

    /// VEP codon formatting: changed base at position 0 is uppercase.
    #[test]
    fn format_codons_position_0() {
        assert_eq!(format_codons(b"CGT", b"TGT", 0), "Cgt/Tgt");
    }

    /// VEP codon formatting: changed base at position 1 is uppercase.
    #[test]
    fn format_codons_position_1() {
        assert_eq!(format_codons(b"CGT", b"CAT", 1), "cGt/cAt");
    }

    /// VEP codon formatting: changed base at position 2 is uppercase.
    #[test]
    fn format_codons_position_2() {
        assert_eq!(format_codons(b"CGT", b"CGA", 2), "cgT/cgA");
    }
}
