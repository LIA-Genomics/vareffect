"""Query VEP REST for listed (label, chrom, start, end, alt, tx, extra) and print TSV rows."""
import json, subprocess, sys, time
def q(label, region, alt, tx, extra):
    filt = "" if tx.startswith(("NM_", "NR_")) else f"&transcript_id={tx}"
    url = f"https://rest.ensembl.org/vep/human/region/{region}/{alt}?content-type=application/json&hgvs=1{filt}{extra}"
    for attempt in range(3):
        r = subprocess.run(["curl", "-s", "-H", "Accept: application/json", url], capture_output=True, text=True).stdout
        try:
            d = json.loads(r); break
        except json.JSONDecodeError:
            time.sleep(2)
    else:
        print(label, region, alt, "NO_JSON", sep="\t"); return
    if isinstance(d, dict):
        print(label, region, alt, "ERR", d.get("error", "")[:80], sep="\t"); return
    for t in d[0].get("transcript_consequences", []):
        if not t.get("transcript_id", "").startswith(tx): continue
        print(label, region, d[0]["allele_string"], t.get("transcript_id"), ",".join(t.get("consequence_terms", [])), t.get("impact"),
              t.get("hgvsc"), t.get("hgvsp"), t.get("cdna_start"), t.get("cds_start"), t.get("protein_start"), t.get("codons"), t.get("amino_acids"),
              t.get("exon"), t.get("intron"), ",".join(t.get("flags", []) or []) or "-", sep="\t")
    time.sleep(0.15)
for line in sys.stdin:
    if not line.strip() or line.startswith("#"): continue
    label, region, alts, tx, extra = (line.rstrip("\n").split("\t") + [""])[:5]
    for alt in alts.split(","):
        q(label, region, alt, tx, extra)
