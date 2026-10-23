#!/usr/bin/env bash
# Local (macOS, Apple Silicon) version of cloud/run_cloud.sh:
#   reads + REMM download -> Exomiser genome rerun WITH REMM -> minimap2 alignment ->
#   Manta (Docker) + CNVpytor -> Exomiser on the SV calls -> results/sv/SUMMARY.txt
#
#   nohup caffeinate -dimsu ./run_local_sv.sh > results/sv/run.log 2>&1 &
#
# Steps leave markers in data/sv/done/, so rerunning resumes where it stopped.
# Steps that use a lot of memory run one at a time (32 GB RAM).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
D="$HERE/data/sv"; R="$HERE/results/sv"
mkdir -p "$D"/{reads,bam,manta,cnv,done} "$R"
SAMPLE=WGS_EX2312012
REF="$HERE/data/ref_primary.fa"           # GRCh38 primary, contigs 1..22,X,Y,MT (same names as the VCF)
JAVA=/opt/homebrew/opt/openjdk@21/bin/java
EXO="$HERE/tools/exomiser-cli-15.1.0"
REMM="$HERE/data/exomiser/data/remm/ReMM.v0.4.hg38.tsv.gz"
MANTA_IMG=quay.io/biocontainers/manta:1.6.0--h9ee0642_3
T=11

log()  { echo "[$(date '+%F %T')] $*"; }
step() { [[ -e "$D/done/$1" ]] && { log "skip $1 (done)"; return 1; }; log "start $1"; }
ok()   { touch "$D/done/$1"; log "done  $1"; }

# ------------------------------------------------------------------ downloads (in parallel)
if step downloads; then
  hf download SageBio/mva-hackathon-2026-data --repo-type dataset --local-dir "$D/reads" \
     --include "*_R[12]_001.fastq.gz" > "$D/reads/dl.log" 2>&1 &
  P1=$!
  ( cd "$(dirname "$REMM")" && for i in $(seq 1 50); do
      curl -sS -L -C - -o ReMM.v0.4.hg38.tsv.gz https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz && break
      sleep 10; done ) &
  P2=$!
  wait $P1 $P2
  [[ $(ls "$D"/reads/*_R[12]_001.fastq.gz | wc -l) -eq 8 ]] || { echo "expected 8 FASTQs"; exit 1; }
  [[ $(stat -f %z "$REMM") -eq 15967302704 ]] || { echo "REMM incomplete: $(stat -f %z "$REMM") bytes"; exit 1; }
  ok downloads
fi

# ------------------------------------------------------------------ Exomiser genome rerun WITH REMM (minutes)
exo_job() {   # $1 = output name, $2 = vcf, $3 = "drop-abortion" to leave out HP:0200067
  sed -e 's/pathogenicitySources: \[ REVEL, MVP, ALPHA_MISSENSE, SPLICE_AI \]/pathogenicitySources: [ REVEL, MVP, ALPHA_MISSENSE, REMM, SPLICE_AI ]/' \
      -e "s#outputFileName: .*#outputFileName: $1#" \
      -e "s#uri: .*#uri: $2#" \
      "$HERE/analysis/proband_genome.yml" | { if [[ "${3:-}" == drop-abortion ]]; then grep -v 'HP:0200067'; else cat; fi; }
}
exo_run() {   # $1 = yml path
  ( cd "$EXO" && "$JAVA" -Xmx16g -Dexomiser.hg38.remm-path="$REMM" -jar exomiser-cli-15.1.0.jar analyse --analysis "$1" )
}
if step exomiser_remm; then
  VCF="$HERE/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
  exo_job proband-genome-remm             "$VCF"               > "$HERE/analysis/proband_genome_remm.yml"
  exo_job proband-genome-remm-no-abortion "$VCF" drop-abortion > "$HERE/analysis/proband_genome_remm_no_abortion.yml"
  grep -q REMM "$HERE/analysis/proband_genome_remm.yml" || { echo "REMM not added to yml"; exit 1; }
  exo_run "$HERE/analysis/proband_genome_remm.yml"             > "$R/exomiser_remm.log" 2>&1
  exo_run "$HERE/analysis/proband_genome_remm_no_abortion.yml" > "$R/exomiser_remm_no_abortion.log" 2>&1
  ok exomiser_remm
fi

# ------------------------------------------------------------------ alignment (lanes one at a time, ~11 GB each)
BAM="$D/bam/$SAMPLE.bam"
if step align; then
  [[ -s "$REF.fai" ]] || samtools faidx "$REF"
  [[ -s "$D/ref.sr.mmi" ]] || minimap2 -t $T -x sr -d "$D/ref.sr.mmi" "$REF"
  for L in L001 L002 L003 L004; do
    [[ -e "$D/done/align_$L" ]] && continue
    log "aligning $L"
    minimap2 -t $T -ax sr -R "@RG\tID:$SAMPLE.$L\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA\tPU:$L" "$D/ref.sr.mmi" \
        "$D"/reads/*_${L}_R1_001.fastq.gz "$D"/reads/*_${L}_R2_001.fastq.gz 2> "$D/bam/$L.mm2.log" \
      | samtools fixmate -u -m - - \
      | samtools sort -u -@3 -m 2G -T "$D/bam/tmp_$L" -o "$D/bam/$L.bam" -
    samtools quickcheck "$D/bam/$L.bam"
    touch "$D/done/align_$L"
  done
  samtools merge -@4 -u -o - "$D"/bam/L00?.bam | samtools markdup -@6 --write-index - "$BAM##idx##$BAM.bai"
  samtools flagstat -@8 "$BAM" > "$R/flagstat.txt"
  rm -f "$D"/bam/L00?.bam
  ok align
fi

# ------------------------------------------------------------------ Manta (Linux-only -> Docker) and CNVpytor, in parallel
if step manta; then
  ( rm -rf "$D/manta/run"
    docker run --rm --platform linux/amd64 -v "$D:$D" -v "$(dirname "$REF"):$(dirname "$REF")" "$MANTA_IMG" \
      bash -c "configManta.py --bam $BAM --referenceFasta $REF --runDir $D/manta/run && $D/manta/run/runWorkflow.py -j 8 -g 7" \
      > "$R/manta.log" 2>&1
    cp "$D/manta/run/results/variants/diploidSV.vcf.gz"* "$R/" && touch "$D/done/manta" ) &
  PM=$!
fi
if step cnv; then
  ( P="$D/cnv/$SAMPLE.pytor"
    # cnvpytor 1.3.2's -download crashes (PosixPath.split); 1.3.1 works but needs pkg_resources
    # macOS defaults to 'spawn' multiprocessing, which cnvpytor's local worker functions can't use: force 'fork'
    CNV() { OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES uv run -q --python 3.12 --with 'cnvpytor==1.3.1' --with 'setuptools<70' --with 'numpy<2' \
              python -c "import sys, multiprocessing as m; m.set_start_method('fork'); sys.argv[0]='cnvpytor'; from cnvpytor.__main__ import main; main()" "$@"; }
    CNV -download
    CNV -root "$P" -rd "$BAM" && CNV -root "$P" -his 10000 100000 && CNV -root "$P" -partition 10000 100000 \
      && CNV -root "$P" -call 10000 > "$R/cnvpytor_calls_10000.tsv" && CNV -root "$P" -call 100000 > "$R/cnvpytor_calls_100000.tsv" \
      && touch "$D/done/cnv" ) > "$R/cnv.log" 2>&1 &
  PC=$!
fi
wait ${PM:-} ${PC:-} || true
[[ -e "$D/done/manta" ]] || { echo "Manta failed, see $R/manta.log"; exit 1; }
[[ -e "$D/done/cnv" ]]   || log "WARNING: CNVpytor failed (see $R/cnv.log); continuing with Manta only"

# ------------------------------------------------------------------ Exomiser on the SV calls
if step exomiser_sv; then
  exo_job proband-sv "$R/diploidSV.vcf.gz" > "$HERE/analysis/proband_sv.yml"
  exo_run "$HERE/analysis/proband_sv.yml" > "$R/exomiser_sv.log" 2>&1 || log "WARNING: Exomiser SV run failed (see $R/exomiser_sv.log)"
  ok exomiser_sv
fi

# ------------------------------------------------------------------ focused summary
{
  echo "# SV / CNV calls near the MVA genes (+/- ~100 kb)"
  for g in "BUB1B 15:40100000-40330000" "CEP57 11:95650000-95910000" "TRIP13 5:790000-1000000"; do
    set -- $g; echo "## $1 ($2)"
    echo "Manta:"; bcftools view -H -r "$2" "$R/diploidSV.vcf.gz" | cut -f1-7 || true
    if [[ -s "$R/cnvpytor_calls_10000.tsv" ]]; then
      echo "CNVpytor (10 kb bins):"
      awk -v r="$2" 'BEGIN{split(r,a,"[:-]")} {split($2,b,"[:-]"); if (b[1]==a[1] && b[3]>=a[2] && b[2]<=a[3]) print}' "$R/cnvpytor_calls_10000.tsv"
    fi
  done
  echo; echo "# Manta PASS calls by type"
  bcftools view -H -f PASS "$R/diploidSV.vcf.gz" | grep -o 'SVTYPE=[A-Z]*' | sort | uniq -c
  for f in proband-genome-remm proband-genome-remm-no-abortion proband-sv; do
    echo; echo "# Exomiser top 10: $f"; head -11 "$HERE/results/$f.genes.tsv" 2>/dev/null | cut -f1-6 || echo "(missing)"
  done
} > "$R/SUMMARY.txt" 2>&1
log "ALL DONE"; cat "$R/SUMMARY.txt"
