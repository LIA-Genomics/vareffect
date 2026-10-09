"""Fetch VEP REST rows for insertions and print fixture TSV lines.

Input (stdin, tab-separated): label, site, inserted length, transcript,
and optionally the plus-strand inserted bases (then used as given, e.g. to
query a ClinVar duplication).
`site` is either an HGVS c./n. string (only its position is used) or a
region `chrom:start-end` in VEP's insertion form (start = end + 1, the
insertion goes between `end` and `start`, 1-based). The inserted base is
the first of A, C, G, T that matches neither genomic flank, so the
insertion is never a duplication or a shiftable repeat.

Output: one TSV row per input line with label, chrom, left flank (1-based),
inserted plus-strand bases, transcript, VEP fields. Each raw JSON response
goes to `<raw_dir>/<label>.json` with the chosen site and bases; cases
already in `raw_dir` are not fetched again, so an interrupted run resumes.
Usage:

    python3 -I vep_insertion_rows.py RAW_DIR < cases.tsv > rows.tsv
"""
import json
import os
import subprocess
import sys
import time

SERVER = "https://rest.ensembl.org"
FIELDS = ["consequence_terms", "impact", "hgvsc", "hgvsp", "cdna_start", "cdna_end",
          "cds_start", "cds_end", "protein_start", "codons", "amino_acids", "exon", "intron"]


def get(path):
    for attempt in range(8):
        out = subprocess.run(["curl", "-s", "-H", "Accept: application/json", SERVER + path],
                             capture_output=True, text=True).stdout
        try:
            return json.loads(out)
        except json.JSONDecodeError:
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"no JSON from {path}")


def query_suffix(tx):
    refseq = "&refseq=1" if tx.startswith(("NM_", "NR_")) else ""
    return f"?content-type=application/json&hgvs=1&numbers=1{refseq}&transcript_id={tx}"


def fetch_case(raw_path, site, ins_len, tx, given):
    if site.split(":")[0].startswith(("NM_", "NR_", "ENST")):
        resolved = get(f"/vep/human/hgvs/{site}{query_suffix(tx)}")
        chrom, right = resolved[0]["seq_region_name"], resolved[0]["start"]
    else:
        chrom, span = site.split(":")
        right = int(span.split("-")[0])
    left = right - 1
    flanks = get(f"/sequence/region/human/{chrom}:{left}..{right}?content-type=application/json")["seq"]
    if given:
        inserted = given
    else:
        base = next(b for b in "ACGT" if b not in flanks)
        inserted = base * ins_len
    request = f"/vep/human/region/{chrom}:{right}-{left}/{inserted}{query_suffix(tx)}"
    response = get(request)
    with open(raw_path, "w") as raw:
        json.dump({"chrom": chrom, "left": left, "inserted": inserted, "request": request,
                   "response": response}, raw, indent=1)
    time.sleep(0.3)


def main():
    raw_dir = sys.argv[1]
    os.makedirs(raw_dir, exist_ok=True)
    release = get("/info/software?content-type=application/json")["release"]
    print(f"# Ensembl VEP REST release {release}, GRCh38, fetched {time.strftime('%Y-%m-%d')}, "
          f"region endpoint with hgvs=1&numbers=1 (refseq=1 for NM_/NR_) and transcript_id. "
          f"Tool: tests/data/tools/vep_insertion_rows.py")
    print("\t".join(["label", "chrom", "left_flank", "inserted", "transcript"] + FIELDS))
    for line in sys.stdin:
        if not line.strip() or line.startswith("#"):
            continue
        label, site, ins_len, tx, *given = line.rstrip("\n").split("\t")
        raw_path = os.path.join(raw_dir, f"{label}.json")
        if not os.path.exists(raw_path):
            fetch_case(raw_path, site, int(ins_len), tx, given[0] if given else "")
        with open(raw_path) as raw:
            case = json.load(raw)
        hit = next((t for t in case["response"][0].get("transcript_consequences", [])
                    if t["transcript_id"].startswith(tx)), {})
        values = []
        for field in FIELDS:
            value = hit.get(field)
            values.append(",".join(value) if isinstance(value, list) else str(value))
        print("\t".join([label, f"chr{case['chrom']}", str(case["left"]), case["inserted"],
                         hit.get("transcript_id", tx)] + values))


main()
