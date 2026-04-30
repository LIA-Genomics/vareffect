# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pysam>=0.22",
#     "httpx>=0.28",
# ]
# ///
"""Generate the VEP ground-truth TSV consumed by `vep_large_concordance.rs`
(GRCh38) or `vep_large_concordance_grch37.rs` (GRCh37).

Workflow:

    # GRCh38 (50k rows, ~1 hour wall time)
    uv run vareffect/scripts/generate_vep_ground_truth.py \\
      --assembly grch38 \\
      --input  data/clinvar/clinvar_grch38.vcf.gz \\
      --output vareffect/tests/data/vep_ground_truth.tsv

    # GRCh37
    uv run vareffect/scripts/generate_vep_ground_truth.py \\
      --assembly grch37 \\
      --input  data/clinvar/clinvar_grch37.vcf.gz \\
      --output vareffect/tests/data/vep_ground_truth_grch37.tsv

Algorithm:

1. Parse the ClinVar VCF, classifying every single-allele row into a
   stratum bucket (snv / del_1bp / del_multi / ins / delins).
2. Sample to per-stratum targets (defaults sum to 50,000; scale via
   `--sample-size`).
3. POST the sampled variants in batches of 200 to the assembly's VEP
   REST endpoint with `refseq=1&hgvs=1&shift_hgvs=1&numbers=1`.
4. For each variant, take the first `transcript_consequences` entry
   whose `transcript_id` starts with `NM_`/`NR_` (RefSeq), filtering out
   predicted `XM_`/`XR_` to keep the ground truth on the clinical-grade
   tier.
5. Write a TSV in the schema `chrom pos ref alt transcript_id
   consequence impact hgvsc hgvsp protein_start protein_end amino_acids
   codons exon intron strand biotype`. Both harnesses parse this layout
   identically.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import random
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import httpx
import pysam

# Ensembl REST endpoint per assembly. The main `rest.ensembl.org` host
# defaults to GRCh38; GRCh37 lives on a sibling subdomain.
REST_ENDPOINTS = {
    "grch38": "https://rest.ensembl.org",
    "grch37": "https://grch37.rest.ensembl.org",
}
ENDPOINT_POST = "/vep/human/region"
QUERY_PARAMS = {
    "refseq": 1,
    "hgvs": 1,
    "shift_hgvs": 1,
    "numbers": 1,
}
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
BATCH_SIZE = 200  # Ensembl REST max per POST.
ACGT = frozenset("ACGT")
DEFAULT_TARGETS = {
    "snv": 25_000,
    "del_1bp": 10_000,
    "del_multi": 7_500,
    "ins": 5_000,
    "delins": 2_500,
}
TSV_COLUMNS = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "transcript_id",
    "consequence",
    "impact",
    "hgvsc",
    "hgvsp",
    "protein_start",
    "protein_end",
    "amino_acids",
    "codons",
    "exon",
    "intron",
    "strand",
    "biotype",
)


@dataclass(frozen=True, slots=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str

    def stratum(self) -> str | None:
        if not self.ref or not self.alt:
            return None
        if not set(self.ref).issubset(ACGT) or not set(self.alt).issubset(ACGT):
            return None
        rl, al = len(self.ref), len(self.alt)
        if rl == al == 1:
            return "snv"
        if rl > al and al == 1 and rl == 2:
            return "del_1bp"
        if rl > al and al == 1:
            return "del_multi"
        if al > rl and rl == 1:
            return "ins"
        return "delins"

    def vcf_string(self) -> str:
        """VEP REST POST body wants whitespace-separated VCF rows.

        VEP accepts `<chrom> <pos> <id> <ref> <alt> . . .` — same as a
        VCF line, just space-separated. Bare-numeric chrom (`17`, not
        `chr17`).
        """
        bare = self.chrom.removeprefix("chr") or self.chrom
        return f"{bare} {self.pos} . {self.ref} {self.alt} . . ."


def scale_targets(base: dict[str, int], total: int) -> dict[str, int]:
    base_total = sum(base.values())
    if total >= base_total:
        return dict(base)
    factor = total / base_total
    return {k: max(1, int(v * factor)) for k, v in base.items()}


def index_vcf(path: Path) -> list[Variant]:
    pysam.set_verbosity(0)
    out: list[Variant] = []
    skipped: defaultdict[str, int] = defaultdict(int)
    seen = 0
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            seen += 1
            if seen % 500_000 == 0:
                print(f"  read {seen:>9,} / kept {len(out):>9,}", file=sys.stderr)
            if rec.alts is None or len(rec.alts) != 1:
                skipped["multi_allele"] += 1
                continue
            alt = rec.alts[0]
            if alt.startswith("<") or "*" in alt:
                skipped["symbolic"] += 1
                continue
            v = Variant(rec.chrom, rec.pos, (rec.ref or "").upper(), alt.upper())
            if v.stratum() is None:
                skipped["non_acgt_or_unclassified"] += 1
                continue
            out.append(v)
    print(
        f"  indexed {len(out):>9,} usable rows from {path.name}; skipped {dict(skipped)}",
        file=sys.stderr,
    )
    return out


def stratify(
    variants: list[Variant],
    targets: dict[str, int],
    seed: int,
) -> list[Variant]:
    bucket: dict[str, list[Variant]] = defaultdict(list)
    for v in variants:
        s = v.stratum()
        if s is not None:
            bucket[s].append(v)
    rng = random.Random(seed)
    out: list[Variant] = []
    for stratum, target in targets.items():
        pool = bucket.get(stratum, [])
        take = min(target, len(pool))
        if take < target:
            print(
                f"  warning: stratum {stratum!r} wanted {target}, only {take} available",
                file=sys.stderr,
            )
        out.extend(rng.sample(pool, take) if take < len(pool) else pool)
    return out


def post_batch(
    client: httpx.Client,
    rest_base: str,
    variants: list[Variant],
    *,
    max_retries: int = 4,
) -> list[dict] | None:
    body = {"variants": [v.vcf_string() for v in variants]}
    url = rest_base + ENDPOINT_POST
    for attempt in range(max_retries):
        resp = client.post(
            url,
            params=QUERY_PARAMS,
            json=body,
            headers=HEADERS,
            timeout=120.0,
        )
        if resp.status_code == 429:
            wait = float(resp.headers.get("Retry-After", "2"))
            time.sleep(wait + 0.5 * attempt)
            continue
        if resp.status_code >= 500:
            time.sleep(2.0 * (attempt + 1))
            continue
        if resp.status_code >= 400:
            print(
                f"  HTTP {resp.status_code}: {resp.text[:200]}",
                file=sys.stderr,
            )
            return None
        return resp.json()
    print("  exhausted retries on batch", file=sys.stderr)
    return None


def pick_refseq_consequence(entry: dict) -> dict | None:
    """Pick the first NM_/NR_ transcript_consequences row from a VEP entry.

    Filters out XM_/XR_ predicted transcripts to keep the ground truth on
    the clinical-grade tier.
    """
    for tc in entry.get("transcript_consequences", []) or []:
        tid = tc.get("transcript_id") or ""
        if tid.startswith(("NM_", "NR_")):
            return tc
    return None


def opt(value: object) -> str:
    return "" if value is None else str(value)


def write_tsv_row(
    fh,
    variant: Variant,
    tc: dict,
) -> None:
    csq = "|".join(tc.get("consequence_terms") or [])
    strand_int = tc.get("strand")
    strand = "+" if strand_int == 1 else "-" if strand_int == -1 else ""
    bare_chrom = variant.chrom.removeprefix("chr") or variant.chrom
    cols = [
        bare_chrom,
        str(variant.pos),
        variant.ref,
        variant.alt,
        opt(tc.get("transcript_id")),
        csq,
        opt(tc.get("impact")),
        opt(tc.get("hgvsc")),
        opt(tc.get("hgvsp")),
        opt(tc.get("protein_start")),
        opt(tc.get("protein_end")),
        opt(tc.get("amino_acids")),
        opt(tc.get("codons")),
        opt(tc.get("exon")),
        opt(tc.get("intron")),
        strand,
        opt(tc.get("biotype")),
    ]
    fh.write("\t".join(cols) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--assembly",
        choices=tuple(REST_ENDPOINTS),
        required=True,
        help="Target assembly. Selects the VEP REST endpoint and the harness "
        "the output TSV is paired with.",
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--sample-size",
        type=int,
        default=sum(DEFAULT_TARGETS.values()),
        help="Total target row count. Default: 50,000.",
    )
    parser.add_argument("--seed", type=int, default=20260430)
    args = parser.parse_args()

    if not args.input.exists():
        print(f"error: {args.input} not found", file=sys.stderr)
        return 2

    rest_base = REST_ENDPOINTS[args.assembly]
    targets = scale_targets(DEFAULT_TARGETS, args.sample_size)
    print(
        f"Targeting {sum(targets.values()):,} rows for {args.assembly.upper()}: {targets}",
        file=sys.stderr,
    )
    print(f"VEP REST endpoint: {rest_base}{ENDPOINT_POST}", file=sys.stderr)

    print(f"Reading {args.input}", file=sys.stderr)
    variants = index_vcf(args.input)
    sampled = stratify(variants, targets, args.seed)
    print(f"Sampled {len(sampled):,} variants for VEP query", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    today = _dt.date.today().isoformat()
    rows_written = 0
    rows_skipped_no_refseq = 0
    harness = (
        "vep_large_concordance.rs"
        if args.assembly == "grch38"
        else "vep_large_concordance_grch37.rs"
    )

    with httpx.Client() as client, args.output.open("w", encoding="utf-8") as fh:
        fh.write(f"# VEP ground truth for {harness}\n")
        fh.write(f"# Assembly: {args.assembly.upper()}\n")
        fh.write(f"# Generated: {today}\n")
        fh.write(f"# Input VCF: {args.input}\n")
        fh.write(f"# VEP endpoint: {rest_base}{ENDPOINT_POST}\n")
        fh.write(f"# VEP parameters: {QUERY_PARAMS}\n")
        fh.write(f"# Stratum targets: {targets}\n")
        fh.write("# Format: tab-separated, see TSV_COLUMNS in the generator.\n")
        fh.write("#\n")
        fh.write("\t".join(TSV_COLUMNS) + "\n")

        total_batches = (len(sampled) + BATCH_SIZE - 1) // BATCH_SIZE
        for batch_idx in range(0, len(sampled), BATCH_SIZE):
            batch = sampled[batch_idx : batch_idx + BATCH_SIZE]
            batch_num = batch_idx // BATCH_SIZE + 1
            print(
                f"  batch {batch_num:>4}/{total_batches} ({len(batch)} variants, "
                f"{rows_written:>6} rows written so far)",
                file=sys.stderr,
            )
            response = post_batch(client, rest_base, batch)
            if response is None:
                continue

            by_input: dict[str, dict] = {}
            for entry in response:
                key = entry.get("input")
                if key:
                    by_input[key] = entry

            for v in batch:
                entry = by_input.get(v.vcf_string())
                if entry is None:
                    rows_skipped_no_refseq += 1
                    continue
                tc = pick_refseq_consequence(entry)
                if tc is None:
                    rows_skipped_no_refseq += 1
                    continue
                write_tsv_row(fh, v, tc)
                rows_written += 1

    print(
        f"\nWrote {rows_written:,} rows to {args.output} "
        f"({rows_skipped_no_refseq:,} variants skipped — no NM_/NR_ consequence)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
