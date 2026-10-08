"""Extract TranscriptModel JSON lines for named transcripts from GENCODE / RefSeq GFF3.

usage: python3 -I extract_fixture.py OUT.jsonl GFF_PATH ID [ID ...]
IDs are GFF transcript IDs (GENCODE `ENST...`, RefSeq `rna-NM_...`). Coordinates are
converted to 0-based half-open; exons and CDS segments are ordered 5'->3' on the
transcript; phase is GFF3 column 8. Mirrors vareffect-cli's builder conventions.
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

def main():
    out, path, wanted = sys.argv[1], sys.argv[2], set(sys.argv[3:])
    tx, exons, cds = {}, {}, {}
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
                "tier": "ManeSelect",
                "biotype": biotype,
                "exon_count": len(exon_objs),
            }
            fo.write(json.dumps(model, separators=(",", ":")) + "\n")

main()
