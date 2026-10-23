#!/usr/bin/env bash
# Unpack Exomiser data (once) and run the proband analysis.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAVA=/opt/homebrew/opt/openjdk@21/bin/java
EXO="$HERE/tools/exomiser-cli-15.1.0"
DATA="$HERE/data/exomiser/data"

mkdir -p "$DATA"
for z in 2512_hg38 2512_phenotype; do
  if [[ ! -d "$DATA/$z" ]]; then
    echo "unpacking $z.zip ..."
    unzip -q -o "$HERE/data/exomiser/$z.zip" -d "$DATA"
  fi
done
ls -d "$DATA"/*/ 2>/dev/null

echo "running exomiser ..."
cd "$EXO"
"$JAVA" -Xmx16g \
  -Dexomiser.data-directory="$DATA" \
  -Dexomiser.hg38.data-version=2512 \
  -Dexomiser.phenotype.data-version=2512 \
  -jar exomiser-cli-15.1.0.jar \
  analyse --analysis "$HERE/analysis/${1:-proband_job.yml}"
