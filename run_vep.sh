#!/usr/bin/env bash
# Annotate the proband VCF with Ensembl VEP (offline cache).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VCF_NAME=WGS_EX2312012_HGWCNDSX7.vcf.gz
OUT=${1:-vep_annotated.tsv.gz}

docker run --rm \
  -v "$HERE/data/vep:/opt/vep/.vep" \
  -v "$HERE/data:/work" \
  ensemblorg/ensembl-vep:release_115.2 \
  vep --offline --cache --dir_cache /opt/vep/.vep \
      --assembly GRCh38 --species homo_sapiens --cache_version 115 \
      --fasta /opt/vep/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
      --input_file "/work/$VCF_NAME" \
      --output_file "/work/$OUT" --compress_output gzip --tab \
      --fork 10 --buffer_size 20000 \
      --hgvs --symbol --canonical --mane --pick_allele_gene \
      --check_existing --af_gnomadg --af_gnomade \
      --sift b --polyphen b --numbers --biotype \
      --stats_file /work/vep_summary.html \
      --force_overwrite --no_escape
