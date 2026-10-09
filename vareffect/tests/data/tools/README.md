# Incomplete-CDS concordance fixtures

`../incomplete_cds_transcripts.jsonl` and `../incomplete_cds_vep116.tsv` back
`tests/incomplete_cds_concordance.rs`. They are generated, not hand-edited.

## Transcript models

`extract_fixture.py` turns named GFF3 transcripts into `TranscriptModel` JSON lines
(0-based half-open coordinates, exons and CDS segments in transcript order,
GFF3 phase kept). Sources:

| File | SHA-256 | Transcripts |
|---|---|---|
| GENCODE v50 `gencode.v50.annotation.gff3.gz` (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/) | `2aaf245c91ed00e80920953add6cfaffcccc876dc0aceeb6ca0c86d15875899a` | ENST00000518351.2, ENST00000616125.5, ENST00000361381.2, ENST00000379099.3, ENST00000400999.7, ENST00000403997.2, ENST00000420515.1, ENST00000436261.7, ENST00000455979.1, ENST00000511286.1, ENST00000609389.5, ENST00000654733.1, ENST00000428369.7 |
| NCBI RefSeq RS_2025_08 `GCF_000001405.40_GRCh38.p14_genomic.gff.gz`, NC_ rows only | `4920f0eae7e2197c50b67a201e06d657387137b49dd60f474b4f1d5b29334051` | NM_001304717.5 (ID `rna-NM_001304717.5`) |
| MANE 1.5 `MANE.GRCh38.v1.5.refseq_genomic.gff.gz` | `040f0d4056de2e9a416cd52bc20ff07ef403baadf5f5968faf188907893f6002` | NM_001424184.1 (ID `rna-NM_001424184.1`) |

```bash
python3 -I extract_fixture.py out_gencode.jsonl gencode.v50.annotation.gff3.gz ENST00000361381.2 ...
python3 -I extract_fixture.py out_gencode_r4.jsonl gencode.v50.annotation.gff3.gz ENST00000654733.1 ENST00000428369.7
python3 -I extract_fixture.py out_gencode_wp2.jsonl gencode.v50.annotation.gff3.gz ENST00000616125.5
python3 -I extract_fixture.py out_refseq.jsonl refseq_nc_only.gff rna-NM_001304717.5
python3 -I extract_fixture.py out_mane.jsonl MANE.GRCh38.v1.5.refseq_genomic.gff.gz rna-NM_001424184.1
cat out_gencode.jsonl out_refseq.jsonl out_mane.jsonl out_gencode_r4.jsonl out_gencode_wp2.jsonl \
  > ../incomplete_cds_transcripts.jsonl
```

The models carry the store format 2 fields (see the extractor's docstring):
`tier` from the MANE tag (`null` otherwise), CDS completeness from the GENCODE
`cds_start_NF` / `cds_end_NF` tags or the RefSeq ranges, and the genome-mismatch
flags. NM_001304717.5 has a `D` gap (a 1-base frameshift intron), so
`genome_mismatch` is set and `cds_alignment_gap` is not. Lines keep the order of
the earlier fixture; ENST00000616125.5 (`cds_start_NF` with a whole ATG codon 1)
is appended.

## VEP rows

`vep_rows.py` queries Ensembl VEP REST (release 116 at the time of capture, GRCh38)
for tab-separated `label, region, alt alleles, transcript, extra query` lines on
stdin. The raw records are in `vep_rest/` (each file names its release, date and
request). `merge_vep_rows.py` folds them into the fixture:

```bash
python3 -I merge_vep_rows.py vep_rest/vep_rest_padded_cds.tsv vep_rest/vep_rest_3prime_partial.tsv \
  vep_rest/vep_rest_wp1.tsv vep_rest/vep_rest_wp1_review.tsv vep_rest/vep_rest_wp1_r1.tsv vep_rest/vep_rest_wp1_r2.tsv vep_rest/vep_rest_wp1_r3.tsv vep_rest/vep_rest_wp1_r4.tsv vep_rest/vep_rest_wp1_r5.tsv \
  vep_rest/vep_rest_wp2.tsv > ../incomplete_cds_vep116.tsv
```

## Insertion-point fixture

`../insertion_point_transcripts.jsonl` (MANE 1.5 NM_000546.6, NM_000059.4,
NM_006231.4, NM_000548.5, NM_003036.4, NM_206933.4, NM_001232.4, NM_213653.4,
NM_001854.4; GENCODE v50 ENST00000424518.6, ENST00000617316.2,
ENST00000372299.4, ENST00000336338.10, extracted with `extract_fixture.py`) and
`../insertion_point_vep116.tsv` back `tests/insertion_point_concordance.rs`:

```bash
python3 -I vep_insertion_rows.py RAW_DIR < vep_rest/insertion_point_cases.tsv > ../insertion_point_vep116.tsv
```

## Built-store rows

`../source_stores_vep116.tsv` adds VEP 116 rows for the stores built from
RefSeq RS_2025_08 and GENCODE v50 (split stop codons on both strands, the
chrY PAR copy of NM_178129.5, the `I`-gapped NM_002111.8, the non-MANE
NM_001126113.3, and SHANK3 NM_001372044.2 / NBPF1 NM_001405693.2 at their
coverage breaks). `tests/incomplete_cds_concordance.rs`
(`built_stores_match_vep_116`) runs it and `../incomplete_cds_vep116.tsv`
against `data/vareffect/transcript_models_{refseq,gencode}.bin`:

```bash
python3 -I vep_rows.py < vep_rest/wp3_cases.tsv   # recorded as vep_rest/vep_rest_wp3.tsv (2026-10-09)
python3 -I merge_vep_rows.py /dev/null /dev/null /dev/null vep_rest/vep_rest_wp3.tsv \
  > ../source_stores_vep116.tsv
```
