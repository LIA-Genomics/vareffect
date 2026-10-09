"""Extract TranscriptModel JSON lines for named transcripts from GENCODE / RefSeq GFF3.

usage: python3 -I extract_fixture.py OUT.jsonl GFF_PATH ID [ID ...]
IDs are GFF transcript IDs (GENCODE `ENST...`, RefSeq `rna-NM_...`). Coordinates are
converted to 0-based half-open; exons and CDS segments are ordered 5'->3' on the
transcript; phase is GFF3 column 8. Mirrors vareffect-cli's builder conventions.

Store format 2 fields: `tier` from the MANE / RefSeq Select tag (`null` without
one); CDS completeness from GENCODE `cds_start_NF` / `cds_end_NF` tags or RefSeq
CDS `start_range` / `end_range` (strand-aware; `partial=true` without a range is
incomplete at both ends); `genome_mismatch` from a RefSeq `exception` or a
gapped `cDNA_match`; `cds_alignment_gap` when a plus-strand `cDNA_match` `Gap`
has an `I` op (transcript bases the genome lacks) strictly inside the CDS span.
Selenocysteine codons are not read (`[]`); a selenoprotein, or an `I` op on a
minus-strand alignment, stops the script.
"""
import gzip, json, re, sys

REFSEQ_CHROM = {f"NC_0000{n:02d}": f"chr{n}" for n in range(1, 23)}
REFSEQ_CHROM.update({"NC_000023": "chrX", "NC_000024": "chrY", "NC_012920": "chrM"})

def chrom_of(seqid):
    if seqid.startswith("chr"):
        return seqid
    return REFSEQ_CHROM[seqid.split(".")[0]]

def attr(attrs, key):
    m = re.search(r"(?:^|;)" + key + r"=([^;]*)", attrs)
    return m.group(1) if m else None

TIERS = [("MANE Select", "ManeSelect"), ("MANE Plus Clinical", "ManePlusClinical"), ("RefSeq Select", "RefSeqSelect")]

def insertion_points(start0, gap, target):
    """0-based genomic points of the `I` ops of a plus-strand cDNA_match `Gap`."""
    points, pos = [], start0
    for op in gap.split():
        kind, n = op[0], int(op[1:])
        if kind == "I":
            if target.split()[-1] != "+":
                sys.exit(f"{target}: I op on a minus-strand alignment is not supported")
            points.append(pos)
        else:
            pos += n
    return points

def main():
    out, path, wanted = sys.argv[1], sys.argv[2], set(sys.argv[3:])
    tx, exons, cds, ranges, gaps, gapped = {}, {}, {}, {}, {}, set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = f[8]
            if f[2] in ("exon", "CDS"):
                p = (attr(a, "Parent") or "").split(",")[0]
                if p in wanted:
                    (exons if f[2] == "exon" else cds).setdefault(p, []).append((int(f[3]) - 1, int(f[4]), f[7]))
                    if f[2] == "CDS":
                        ranges.setdefault(p, set()).update(k for k in ("start_range", "end_range") if attr(a, k))
            elif f[2] == "stop_codon_redefined_as_selenocysteine" and attr(a, "Parent") in wanted:
                sys.exit(f"{attr(a, 'Parent')}: selenocysteine transcripts are not supported by this extractor")
            elif f[2] == "cDNA_match" and attr(a, "Gap"):
                target = attr(a, "Target")
                p = "rna-" + target.split()[0]
                if p in wanted:
                    gapped.add(p)
                    gaps.setdefault(p, []).extend(insertion_points(int(f[3]) - 1, attr(a, "Gap"), target))
            else:
                i = attr(a, "ID")
                if i in wanted and f[2] not in ("gene",):
                    tx[i] = (f, a)
    with open(out, "w") as fo:
        for i in sorted(wanted):
            f, a = tx[i]
            plus = f[6] == "+"
            ex = sorted(exons[i], key=lambda e: e[0], reverse=not plus)
            exon_objs = [{"exon_number": n + 1, "genomic_start": s, "genomic_end": e} for n, (s, e, _) in enumerate(ex)]
            segs = []
            for s, e, ph in sorted(cds.get(i, []), key=lambda c: c[0], reverse=not plus):
                idx = next(k for k, (xs, xe, _) in enumerate(ex) if xs <= s and e <= xe)
                segs.append({"exon_index": idx, "genomic_start": s, "genomic_end": e, "phase": 0 if ph == "." else int(ph)})
            accession = attr(a, "transcript_id") or i
            if "Sec" in (attr(a, "transl_except") or ""):
                sys.exit(f"{accession}: selenocysteine transcripts are not supported by this extractor")
            raw_tags = (attr(a, "tag") or "").split(",")
            tags = [t.replace("_", " ") for t in raw_tags]
            tier = next((name for tag, name in TIERS if tag in tags), None)
            start_nf = "cds_start_NF" in raw_tags
            end_nf = "cds_end_NF" in raw_tags
            low, high = ("start_range" in ranges.get(i, ())), ("end_range" in ranges.get(i, ()))
            five, three = (low, high) if plus else (high, low)
            partial_without_range = attr(a, "partial") == "true" and not (low or high)
            start_complete = not (start_nf or five or partial_without_range)
            end_complete = not (end_nf or three or partial_without_range)
            cds_lo = min((s["genomic_start"] for s in segs), default=0)
            cds_hi = max((s["genomic_end"] for s in segs), default=0)
            cds_gap = any(cds_lo < g < cds_hi for g in gaps.get(i, []))
            biotype = attr(a, "transcript_type") or ("protein_coding" if segs else "non_coding_rna")
            model = {
                "accession": accession,
                "protein_accession": attr(a, "protein_id"),
                "gene_symbol": attr(a, "gene_name") or attr(a, "gene") or "",
                "hgnc_id": attr(a, "hgnc_id"),
                "ensembl_accession": None,
                "chrom": chrom_of(f[0]),
                "strand": "Plus" if plus else "Minus",
                "tx_start": int(f[3]) - 1,
                "tx_end": int(f[4]),
                "cds_genomic_start": min(s["genomic_start"] for s in segs) if segs else None,
                "cds_genomic_end": max(s["genomic_end"] for s in segs) if segs else None,
                "exons": exon_objs,
                "cds_segments": segs,
                "tier": tier,
                "biotype": biotype,
                "exon_count": len(exon_objs),
                "cds_start_complete": start_complete if segs else True,
                "cds_end_complete": end_complete if segs else True,
                "genome_mismatch": attr(a, "exception") is not None or i in gapped,
                "cds_alignment_gap": cds_gap,
                "selenocysteine_codons": [],
            }
            fo.write(json.dumps(model, separators=(",", ":")) + "\n")

main()
