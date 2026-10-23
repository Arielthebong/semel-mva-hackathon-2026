#!/usr/bin/env bash
# Run LIRICAL 2.4.1 (genotype-aware, hg38) on the proband, reusing the Exomiser 2512 databases.
#   ./run_lirical.sh                 # all HPO terms in data/hpo_terms.txt
#   ./run_lirical.sh HP:0200067      # same, but leave out the listed term(s) (comma-separated)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAVA=/opt/homebrew/opt/openjdk@21/bin/java
LIR="$HERE/tools/lirical-cli-2.4.1"
E="$HERE/data/exomiser/data/2512_hg38"
DROP=${1:-}
PREFIX=proband${DROP:+_minus_${DROP//[:,]/_}}

# Data fixes for the current HGNC release (LIRICAL 2.4.1 reads it by fixed column position):
#  1. HGNC dropped the `location_sortable` column, shifting entrez_id onto ensembl_gene_id.
#     Without this, no disease links to a gene and genotype evidence is silently ignored.
#  2. Five lncRNA/pseudogene rows share an Ensembl or Entrez ID, which crashes the loader.
python3 - "$LIR/data/hgnc_complete_set.txt" <<'PY'
import csv, sys
p = sys.argv[1]
rows = list(csv.reader(open(p), delimiter="\t"))
h = rows[0]
if "location_sortable" not in h:
    i = h.index("location") + 1
    for r in rows:
        r.insert(i, "location_sortable" if r is rows[0] else "")
    h = rows[0]
cols = [h.index("entrez_id"), h.index("ensembl_gene_id")]
seen = {c: set() for c in cols}
for r in rows[1:]:
    for c in cols:
        if len(r) > c and r[c]:
            if r[c] in seen[c]: r[c] = ""
            else: seen[c].add(r[c])
csv.writer(open(p, "w", newline=""), delimiter="\t", lineterminator="\n").writerows(rows)
PY

TERMS=$(cut -f1 "$HERE/data/hpo_terms.txt" | grep '^HP:' | grep -vxF -f <(tr ',' '\n' <<<"$DROP") | paste -sd, -)
echo "HPO terms: $TERMS"
mkdir -p "$HERE/results/lirical"
"$JAVA" -Xmx16g -jar "$LIR/lirical-cli-2.4.1.jar" prioritize \
  -d "$LIR/data" -e38 "$E/2512_hg38_variants.mv.db" -c38 "$E/2512_hg38_clinvar.mv.db" \
  --assembly hg38 --vcf "$HERE/data/WGS_EX2312012_HGWCNDSX7.vcf.gz" \
  --sample-id WGS_EX2312012 --sex MALE -p "$TERMS" --parallelism 8 \
  -o "$HERE/results/lirical" -x "$PREFIX" -f tsv -f html -f json > "$HERE/results/lirical/$PREFIX.log" 2>&1
grep -v '^!' "$HERE/results/lirical/$PREFIX.tsv" | awk -F'\t' 'NR>1 && NR<=6 {printf "%s\t%-60s\tpost-test %s\n",$1,$2,$5}'
