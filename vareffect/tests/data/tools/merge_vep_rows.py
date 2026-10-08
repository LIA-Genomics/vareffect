"""Merge the recorded VEP REST TSVs into one fixture with a fixed column set."""
import ast, sys
cols = ["label", "chrom", "pos", "ref", "alt", "transcript", "consequences", "impact", "hgvsc", "hgvsp",
        "cdna_start", "cds_start", "protein_start", "codons", "amino_acids", "exon", "intron"]
def norm(v):
    return "-" if v in (None, "None", "") else v
rows = []
def add(label, region, allele, hgvsc, consequences, impact, hgvsp, cdna, cds, prot, codons, aa, tx_id=None, exon="?", intron="?"):
    chrom, span = region.split(":")
    start = span.split("-")[0]
    ref, alt = allele.split("/")
    tx = hgvsc.split(":")[0] if hgvsc not in (None, "None") else (tx_id or "-")
    rows.append([label, "chr" + chrom if chrom != "MT" else "chrM", start, ref, alt, tx, consequences, impact,
                 norm(hgvsc), norm(hgvsp), norm(cdna), norm(cds), norm(prot), norm(codons), norm(aa),
                 exon if exon == "?" else norm(exon), intron if intron == "?" else norm(intron)])
padded, three, wp1 = sys.argv[1:4]
numbered = sys.argv[4:]
for line in open(padded):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    # strand phase tx region allele consequences hgvsc hgvsp 'cdna' cdna 'cds' cds 'prot' prot codons aa
    add(f"padded_{f[0]}_phase{f[1]}", f[3], f[4], f[6], ",".join(ast.literal_eval(f[5])), "-", f[7], f[9], f[11], f[13], f[14], f[15])
for line in open(three):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    add(f[0], f[1] + "-" + f[1].split(":")[1], f[2], f[6], f[4], f[5], f[7], f[8], f[9], f[10], f[11], f[12])
for line in open(wp1):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    add(f[0], f[1], f[2], f[6], f[4], f[5], f[7], f[8], f[9], f[10], f[11], f[12], f[3] if "." in f[3] else None)
for path in numbered:
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        add(f[0], f[1], f[2], f[6], f[4], f[5], f[7], f[8], f[9], f[10], f[11], f[12],
            f[3], f[13], f[14])
print("# Ensembl VEP REST release 116 (/info/software), GRCh38, fetched 2026-10-07; one row per (variant, transcript).")
print("# pos is 1-based; ref/alt as VEP's allele_string (`-` = empty). `-` in any other column = field absent in VEP's JSON;")
print("# `?` in exon/intron = not requested (rows queried without numbers=1).")
print("# Built by tools/merge_vep_rows.py (next to this file) from the per-query REST records; see tools/README.md.")
print("\t".join(cols))
# A row without HGVS carries the unversioned transcript id from the request;
# take the version from another row of the same transcript.
versioned = {r[5].split(".")[0]: r[5] for r in rows if "." in r[5]}
for r in rows:
    if r[5] == "-" or "." not in r[5]:
        r[5] = versioned.get(r[5].split(".")[0], r[5])
latest = {}
for r in rows:
    latest[(r[1], r[2], r[3], r[4], r[5])] = r
for r in latest.values():
    print("\t".join(r))
