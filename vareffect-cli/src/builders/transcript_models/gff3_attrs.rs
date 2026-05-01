//! GFF3 column-9 attribute extraction helpers.
//!
//! These helpers parse the semicolon-separated `key=value` attribute string
//! in GFF3 column 9. They handle URL-encoded values per the GFF3 spec and
//! extract cross-references from `Dbxref=` entries.

/// Extract a named attribute from a GFF3 column-9 attribute string.
///
/// GFF3 attributes are semicolon-separated `key=value` pairs with
/// URL-encoded values. The decoder covers the full set of reserved
/// characters per the GFF3 spec (tab, newline, carriage return, percent,
/// semicolon, equals, ampersand, comma) plus the space character which
/// upstream producers occasionally encode even though the spec allows
/// bare spaces inside values.
pub(super) fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    for pair in attrs.split(';') {
        let pair = pair.trim();
        if let Some((k, v)) = pair.split_once('=')
            && k.trim() == key
        {
            return Some(percent_decode_gff3(v.trim()));
        }
    }
    None
}

/// Decode the GFF3-reserved percent-encoded bytes in `raw`.
///
/// Covers the spec-required set (`%09 %0A %0D %25 %26 %2C %3B %3D`) plus
/// the space byte (`%20`) that some producers still emit. A trailing or
/// malformed `%` survives verbatim -- the decoder is intentionally lenient
/// so a bad byte sequence never makes a whole transcript row unparseable.
pub(super) fn percent_decode_gff3(raw: &str) -> String {
    // Fast path: no `%` anywhere -> nothing to decode, avoid the allocation.
    if !raw.contains('%') {
        return raw.to_string();
    }
    let bytes = raw.as_bytes();
    let mut out = String::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            // Decode a two-hex-digit escape. If parsing fails (e.g. `%ZZ`),
            // emit the literal `%` and keep walking.
            let hex = &raw[i + 1..i + 3];
            if let Ok(byte) = u8::from_str_radix(hex, 16) {
                out.push(byte as char);
                i += 3;
                continue;
            }
        }
        // SAFETY: indexing a byte is always valid; pushing a bare ASCII
        // byte is OK because we're reading through `as_bytes()` of a
        // valid `&str`.
        out.push(bytes[i] as char);
        i += 1;
    }
    out
}

/// Extract `HGNC:HGNC:xxxxx` from a `Dbxref=...` attribute value.
///
/// The Dbxref list is comma-separated after decoding. Gene rows typically
/// carry `Dbxref=GeneID:8831,HGNC:HGNC:11497,MIM:123456`.
pub(super) fn extract_hgnc_from_dbxref(attrs: &str) -> Option<String> {
    let dbxref = extract_attr(attrs, "Dbxref")?;
    for entry in dbxref.split(',') {
        let trimmed = entry.trim();
        if let Some(stripped) = trimmed.strip_prefix("HGNC:") {
            // `HGNC:HGNC:11497` -> `HGNC:11497`. The remaining `HGNC:`
            // prefix is the identifier scheme, not a nested duplicate.
            return Some(stripped.to_string());
        }
    }
    None
}

/// Extract `Ensembl:ENSTxxxxxxxxxxx.N` from a `Dbxref=...` attribute value.
pub(super) fn extract_ensembl_from_dbxref(attrs: &str) -> Option<String> {
    let dbxref = extract_attr(attrs, "Dbxref")?;
    for entry in dbxref.split(',') {
        let trimmed = entry.trim();
        if let Some(stripped) = trimmed.strip_prefix("Ensembl:") {
            return Some(stripped.to_string());
        }
    }
    None
}

/// Strip either `GenBank:` (the form NCBI ships in MANE v1.5) or
/// `Genbank:` (legacy older-release form) from a single Dbxref entry,
/// returning the identifier inside.
///
/// Does not lowercase the whole string -- the accessions (`NM_...`, `NP_...`)
/// are case-sensitive.
pub(super) fn strip_genbank_prefix(entry: &str) -> Option<&str> {
    entry
        .strip_prefix("GenBank:")
        .or_else(|| entry.strip_prefix("Genbank:"))
}

/// Extract a RefSeq transcript accession from the `Dbxref=...` list.
///
/// Looks for `GenBank:NM_...` or `GenBank:NR_...` entries. Used as a
/// fallback when the mRNA row's `Name` attribute is ambiguous.
pub(super) fn extract_refseq_acc_from_dbxref(attrs: &str) -> Option<String> {
    let dbxref = extract_attr(attrs, "Dbxref")?;
    for entry in dbxref.split(',') {
        let trimmed = entry.trim();
        if let Some(stripped) = strip_genbank_prefix(trimmed)
            && (stripped.starts_with("NM_") || stripped.starts_with("NR_"))
        {
            return Some(stripped.to_string());
        }
    }
    None
}

/// Fallback: derive a transcript accession from the GFF3 `ID` attribute
/// (e.g. `rna-NM_006772.2` -> `NM_006772.2`).
pub(super) fn extract_transcript_from_id(id: &str) -> Option<String> {
    let stripped = id.strip_prefix("rna-").unwrap_or(id);
    if stripped.starts_with("NM_") || stripped.starts_with("NR_") {
        Some(stripped.to_string())
    } else {
        None
    }
}

/// RefSeq protein-accession prefixes accepted by the CDS `protein_id`
/// helper: `NP_` (canonical), `YP_` (mitochondrial), `XP_` (predicted).
pub(super) const PROTEIN_ACCESSION_PREFIXES: &[&str] = &["NP_", "YP_", "XP_"];

/// Substring pattern NCBI uses in `Note=` attributes to flag transcripts
/// whose sequence differs from the reference assembly.
///
/// In the GRCh37.p13 GFF3 every divergence note is structured as one of:
///
/// > "The RefSeq transcript aligns at N% coverage compared to this genomic sequence"
/// > "The RefSeq transcript has N substitution(s) compared to this genomic sequence"
/// > "The RefSeq transcript has N frameshift(s) compared to this genomic sequence"
/// > "The RefSeq transcript has N non-frameshifting indel(s) compared to this genomic sequence"
/// > combinations of the above
///
/// All variants share the trailing phrase "compared to this genomic
/// sequence", and that phrase appears nowhere else in `Note=` attributes
/// on mRNA rows in the GRCh37 GFF3 (verified empirically: 6,402 / 6,402
/// of mRNA `Note=` attributes contain it; 0 don't). One pattern is
/// therefore enough.
///
/// Future drift in NCBI's wording surfaces as a ClinVar-concordance
/// regression rather than a silent miss: the ClinVar self-concordance
/// suite asserts the divergent-set fraction lands in [4 %, 6 %] of the
/// GRCh37 store, matching NCBI's published ~5 % figure for RefSeq Select.
pub(super) const DIVERGENCE_NOTE_PATTERN: &str = "compared to this genomic sequence";

/// Detect whether an mRNA or CDS row's `Note=` attribute flags the
/// transcript as having sequence divergence from the reference assembly.
///
/// Case-insensitive substring match against [`DIVERGENCE_NOTE_PATTERN`].
/// Returns `false` for absent or empty `Note=` attributes — divergence is
/// the rare case (~5 % on GRCh37, 0 % on GRCh38 MANE) and the default
/// must be permissive so non-divergent transcripts produce no warning.
pub(super) fn extract_divergence_flag(attrs: &str) -> bool {
    let Some(note) = extract_attr(attrs, "Note") else {
        return false;
    };
    if note.is_empty() {
        return false;
    }
    note.to_ascii_lowercase().contains(DIVERGENCE_NOTE_PATTERN)
}

/// Detect whether a CDS row's `transl_except=` attribute marks a
/// translational exception (selenocysteine / pyrrolysine readthrough).
///
/// Returns the residue identity. Distinct from divergence detection:
/// `transl_except` is a known biological mechanism that vareffect should
/// record but should NOT raise as a clinical warning, because the HGVS
/// position remains correct against the reference.
///
/// NCBI's qualifier format is
/// `transl_except=(pos:X..Y,aa:Sec)` (or `aa:Pyl`).
pub(super) fn extract_translational_exception(
    attrs: &str,
) -> Option<vareffect::TranslationalException> {
    let raw = extract_attr(attrs, "transl_except")?;
    // Look for `aa:Sec` / `aa:Pyl` / `aa:<other>` inside the qualifier.
    let lower = raw.to_ascii_lowercase();
    if lower.contains("aa:sec") || lower.contains("aa:selenocysteine") {
        Some(vareffect::TranslationalException::Selenocysteine)
    } else if lower.contains("aa:pyl") || lower.contains("aa:pyrrolysine") {
        Some(vareffect::TranslationalException::Pyrrolysine)
    } else if let Some(idx) = lower.find("aa:") {
        // Capture the rest of the value after `aa:` up to the next
        // delimiter, so an unknown amino acid round-trips as
        // `Other("Cys")` rather than getting silently dropped.
        let tail = &raw[idx + 3..];
        let stop = tail
            .find(|c: char| !c.is_ascii_alphabetic())
            .unwrap_or(tail.len());
        let aa = &tail[..stop];
        if aa.is_empty() {
            None
        } else {
            Some(vareffect::TranslationalException::Other(aa.to_string()))
        }
    } else {
        None
    }
}

/// Extract a protein accession from a CDS row's attribute column.
///
/// NCBI GFF3 places the protein accession in two possible locations:
/// (a) a direct `protein_id=NP_...` attribute, or (b) `Dbxref=...,
/// GenBank:NP_...,...`. We check `protein_id` first (cheaper) and fall
/// back to the Dbxref scan. Accepted prefixes are centralized in
/// [`PROTEIN_ACCESSION_PREFIXES`].
pub(super) fn extract_protein_id_from_cds_attrs(attrs: &str) -> Option<String> {
    // Primary: `protein_id=NP_000537.3`
    if let Some(pid) = extract_attr(attrs, "protein_id")
        && !pid.is_empty()
    {
        return Some(pid);
    }
    // Fallback: `Dbxref=...,GenBank:NP_000537.3,...`
    let dbxref = extract_attr(attrs, "Dbxref")?;
    for entry in dbxref.split(',') {
        let trimmed = entry.trim();
        if let Some(stripped) = strip_genbank_prefix(trimmed)
            && PROTEIN_ACCESSION_PREFIXES
                .iter()
                .any(|prefix| stripped.starts_with(prefix))
        {
            return Some(stripped.to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_divergence_flag_matches_real_ncbi_grch37_wording() {
        // Verbatim Note= attributes pulled from NCBI's
        // GCF_000001405.25_GRCh37.p13_genomic.gff.gz. All four shapes
        // appear on `tag=RefSeq Select` mRNAs and must be flagged.
        let cases = [
            "ID=rna-X;Note=The RefSeq transcript aligns at 99%25 coverage compared to this genomic sequence;tag=RefSeq Select",
            "ID=rna-X;Note=The RefSeq transcript has 1 substitution compared to this genomic sequence",
            "ID=rna-X;Note=The RefSeq transcript has 1 frameshift compared to this genomic sequence",
            "ID=rna-X;Note=The RefSeq transcript has 1 substitution%2C 1 non-frameshifting indel compared to this genomic sequence",
        ];
        for attrs in cases {
            assert!(
                extract_divergence_flag(attrs),
                "expected divergence flag on: {attrs}",
            );
        }
    }

    #[test]
    fn extract_divergence_flag_is_false_when_note_is_absent_or_unrelated() {
        assert!(!extract_divergence_flag("ID=rna-X;tag=RefSeq Select"));
        assert!(!extract_divergence_flag("ID=rna-X;Note="));
        // GRCh38 MANE rows carry no divergence Note. A benign Note must
        // not false-positive.
        assert!(!extract_divergence_flag(
            "ID=rna-X;Note=transcript variant 2;tag=MANE Select"
        ));
    }

    #[test]
    fn extract_attr_decodes_common_encodings() {
        assert_eq!(
            extract_attr("ID=abc;Name=Foo%2CBar", "Name").as_deref(),
            Some("Foo,Bar")
        );
        assert_eq!(
            extract_attr("ID=abc;Name=simple", "Name").as_deref(),
            Some("simple")
        );
        assert_eq!(extract_attr("ID=abc;Name=simple", "Missing"), None);
    }

    #[test]
    fn extract_hgnc_from_dbxref_works() {
        let attrs = "ID=g1;Dbxref=GeneID:8831,HGNC:HGNC:11497,MIM:603384";
        assert_eq!(
            extract_hgnc_from_dbxref(attrs).as_deref(),
            Some("HGNC:11497")
        );
    }

    #[test]
    fn extract_ensembl_from_dbxref_works() {
        let attrs = "ID=t1;Dbxref=Genbank:NM_006772.2,Ensembl:ENST00000418600.6";
        assert_eq!(
            extract_ensembl_from_dbxref(attrs).as_deref(),
            Some("ENST00000418600.6")
        );
    }

    #[test]
    fn extract_protein_id_prefers_direct_attr_over_dbxref() {
        let attrs = "ID=c1;protein_id=NP_006763.2;Dbxref=Genbank:NP_999999.9";
        assert_eq!(
            extract_protein_id_from_cds_attrs(attrs).as_deref(),
            Some("NP_006763.2")
        );
    }

    #[test]
    fn extract_protein_id_falls_back_to_dbxref() {
        let attrs = "ID=c1;Dbxref=Genbank:NP_000537.3";
        assert_eq!(
            extract_protein_id_from_cds_attrs(attrs).as_deref(),
            Some("NP_000537.3")
        );
    }

    #[test]
    fn extract_protein_id_from_dbxref_handles_capital_b() {
        // Real MANE v1.5 Dbxref entries use `GenBank:` with capital B.
        let attrs = "ID=c1;Dbxref=Ensembl:ENSP00000493376.2,GenBank:NP_006763.2";
        assert_eq!(
            extract_protein_id_from_cds_attrs(attrs).as_deref(),
            Some("NP_006763.2")
        );
    }

    #[test]
    fn extract_protein_id_accepts_predicted_xp_accession() {
        // `XP_*` isn't in today's MANE but the prefix table accepts it.
        let attrs = "ID=c1;Dbxref=GenBank:XP_123456.1";
        assert_eq!(
            extract_protein_id_from_cds_attrs(attrs).as_deref(),
            Some("XP_123456.1")
        );
    }

    #[test]
    fn extract_refseq_acc_from_dbxref_handles_capital_b() {
        let attrs = "ID=r1;Dbxref=GenBank:NM_006772.2";
        assert_eq!(
            extract_refseq_acc_from_dbxref(attrs).as_deref(),
            Some("NM_006772.2")
        );
    }

    #[test]
    fn extract_refseq_acc_from_dbxref_handles_legacy_lowercase_b() {
        // Keep compatibility with older NCBI releases that shipped
        // `Genbank:` (lowercase b).
        let attrs = "ID=r1;Dbxref=Genbank:NR_003051.4";
        assert_eq!(
            extract_refseq_acc_from_dbxref(attrs).as_deref(),
            Some("NR_003051.4")
        );
    }

    #[test]
    fn percent_decode_gff3_handles_full_reserved_set() {
        // Spec-required escapes.
        assert_eq!(percent_decode_gff3("a%09b"), "a\tb"); // tab
        assert_eq!(percent_decode_gff3("a%0Ab"), "a\nb"); // newline
        assert_eq!(percent_decode_gff3("a%0Db"), "a\rb"); // carriage return
        assert_eq!(percent_decode_gff3("a%25b"), "a%b"); // percent
        assert_eq!(percent_decode_gff3("a%26b"), "a&b"); // ampersand
        assert_eq!(percent_decode_gff3("a%2Cb"), "a,b"); // comma
        assert_eq!(percent_decode_gff3("a%3Bb"), "a;b"); // semicolon
        assert_eq!(percent_decode_gff3("a%3Db"), "a=b"); // equals
        assert_eq!(percent_decode_gff3("a%20b"), "a b"); // space
    }

    #[test]
    fn percent_decode_gff3_passes_through_malformed_escapes() {
        // A stray `%` or incomplete escape must not panic; the original
        // bytes survive verbatim.
        assert_eq!(percent_decode_gff3("a%ZZb"), "a%ZZb");
        assert_eq!(percent_decode_gff3("trailing%"), "trailing%");
    }
}
