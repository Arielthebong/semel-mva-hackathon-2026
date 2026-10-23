#!/usr/bin/env bash
# Structural-variant calling + REMM Exomiser rerun on an AWS EC2 VM (Ubuntu 24.04, x86_64).
#
#   read -rs HF_TOKEN && export HF_TOKEN      # paste your Hugging Face token (not echoed, not in history)
#   nohup bash run_cloud.sh > run.log 2>&1 &  # survives SSH disconnects
#   tail -f run.log
#
# Every step writes a marker in $W/done/, so rerunning after a failure resumes where it stopped.
# Output: ~/mva_results.tar.gz (a few MB). Raw data stays on the VM and is destroyed with it.
set -euo pipefail

W=${W:-$HOME/data}                 # needs ~350 GB free: launch with a 600 GB root volume
T=$(nproc)
SAMPLE=WGS_EX2312012
DS=SageBio/mva-hackathon-2026-data
EXO_VER=15.1.0
DATA_VER=2512
mkdir -p "$W"/{ref,reads,bam,manta,cnv,exomiser,results,done}
cd "$W"

log()  { echo "[$(date '+%F %T')] $*"; }
step() { [[ -e "$W/done/$1" ]] && { log "skip $1 (done)"; return 1; }; log "start $1"; }
ok()   { touch "$W/done/$1"; log "done  $1"; }

: "${HF_TOKEN:?export HF_TOKEN first (read -rs HF_TOKEN && export HF_TOKEN)}"

# ---------------------------------------------------------------- tools
if step tools; then
  sudo apt-get update -qq
  sudo apt-get install -y -qq openjdk-21-jre-headless pigz unzip aria2 python3-venv >/dev/null
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME" bin/micromamba
  MM="$HOME/bin/micromamba"
  "$MM" create -y -q -p "$HOME/env/bio" -c conda-forge -c bioconda minimap2 samtools bcftools mosdepth
  "$MM" create -y -q -p "$HOME/env/manta" -c conda-forge -c bioconda manta
  python3 -m venv "$HOME/env/py"
  "$HOME/env/py/bin/pip" install -q "huggingface_hub[cli,hf_transfer]" cnvpytor
  "$HOME/env/py/bin/cnvpytor" -download
  ok tools
fi
export PATH="$HOME/env/bio/bin:$HOME/env/py/bin:$PATH"
MANTA="$HOME/env/manta/bin"

# ---------------------------------------------------------------- REMM: background, only the Exomiser step waits for it
# Its server has been slow, so it must never hold up alignment. aria2c -c resumes after a restart.
if [[ ! -e "$W/done/remm" ]]; then
  ( cd "$W/exomiser" && aria2c -q -c -x8 -s8 --max-tries=0 --retry-wait=30 \
       https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz.tbi \
       https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz \
    && touch "$W/done/remm" ) > "$W/exomiser/remm_dl.log" 2>&1 &
  log "REMM download running in background"
fi

# ---------------------------------------------------------------- downloads (in parallel)
if step downloads; then
  export HF_HUB_ENABLE_HF_TRANSFER=1
  hf download "$DS" --repo-type dataset --local-dir "$W/reads" \
     --include "*.fastq.gz" --include "*.vcf.gz*" --token "$HF_TOKEN" > "$W/reads/dl.log" 2>&1 &
  P1=$!
  ( cd "$W/ref" && aria2c -q -x8 -s8 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    && pigz -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38.fa && samtools faidx GRCh38.fa ) &
  P2=$!
  ( cd "$W/exomiser" && aria2c -q -x8 -s8 \
       "https://github.com/exomiser/Exomiser/releases/download/$EXO_VER/exomiser-cli-$EXO_VER-distribution.zip" \
       "https://data.monarchinitiative.org/exomiser/data/${DATA_VER}_hg38.zip" \
       "https://data.monarchinitiative.org/exomiser/data/${DATA_VER}_phenotype.zip" ) &
  P3=$!
  wait $P1 $P2 $P3
  ls -la "$W/reads" | grep -c fastq.gz | grep -qx 8 || { echo "expected 8 FASTQs"; exit 1; }
  ok downloads
fi
REF="$W/ref/GRCh38.fa"

# ---------------------------------------------------------------- Exomiser data unpack (background while aligning)
if step exo_unpack; then
  ( cd "$W/exomiser" && unzip -q -o "exomiser-cli-$EXO_VER-distribution.zip" \
    && mkdir -p data && unzip -q -o "${DATA_VER}_hg38.zip" -d data && unzip -q -o "${DATA_VER}_phenotype.zip" -d data \
    && touch "$W/done/exo_unpack" ) &
  EXO_UNPACK=$!
fi

# ---------------------------------------------------------------- alignment: 4 lanes in parallel, then merge
BAM="$W/bam/$SAMPLE.bam"
if step align; then
  [[ -s "$W/ref/GRCh38.sr.mmi" ]] || minimap2 -t "$T" -d "$W/ref/GRCh38.sr.mmi" -x sr "$REF"
  LT=$(( (T - 8) / 4 ))                       # threads per lane, leaving room for sort/merge
  for L in L001 L002 L003 L004; do
    ( minimap2 -t "$LT" -ax sr -R "@RG\tID:$SAMPLE.$L\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA\tPU:$L" \
          "$W/ref/GRCh38.sr.mmi" "$W"/reads/*_${L}_R1_001.fastq.gz "$W"/reads/*_${L}_R2_001.fastq.gz \
        | samtools fixmate -u -m - - \
        | samtools sort -u -@2 -m 3G -T "$W/bam/tmp_$L" -o "$W/bam/$L.bam" - ) 2> "$W/bam/$L.log" &
  done
  wait
  for L in L001 L002 L003 L004; do samtools quickcheck "$W/bam/$L.bam" || { echo "lane $L failed, see $W/bam/$L.log"; exit 1; }; done
  samtools merge -@8 -u -o - "$W"/bam/L00?.bam \
    | samtools markdup -@8 --write-index - "$BAM##idx##$BAM.bai"
  rm -f "$W"/bam/L00?.bam
  samtools flagstat -@8 "$BAM" > "$W/results/flagstat.txt"
  mosdepth -t 8 -n --by 1000000 "$W/results/depth_1mb" "$BAM"
  ok align
fi

# ---------------------------------------------------------------- Manta (SVs: deletions, duplications, inversions, translocations)
if step manta; then
  rm -rf "$W/manta/run"
  "$MANTA/configManta.py" --bam "$BAM" --referenceFasta "$REF" --runDir "$W/manta/run"
  "$W/manta/run/runWorkflow.py" -j "$T" -g 100
  cp "$W/manta/run/results/variants/diploidSV.vcf.gz"* "$W/results/"
  ok manta
fi

# ---------------------------------------------------------------- CNVpytor (read-depth copy-number)
if step cnv; then
  P="$W/cnv/$SAMPLE.pytor"
  cnvpytor -root "$P" -rd "$BAM" -T "$REF"
  cnvpytor -root "$P" -his 1000 10000 100000
  cnvpytor -root "$P" -partition 1000 10000 100000
  for b in 1000 10000 100000; do cnvpytor -root "$P" -call "$b" > "$W/results/cnvpytor_calls_${b}.tsv"; done
  ok cnv
fi

# ---------------------------------------------------------------- Exomiser: SNV/indel rerun with REMM, and SV run
[[ -n "${EXO_UNPACK:-}" ]] && wait "$EXO_UNPACK"
[[ -e "$W/done/exo_unpack" ]] || { echo "Exomiser unpack failed"; exit 1; }
EXO="$W/exomiser/exomiser-cli-$EXO_VER"
cat > "$EXO/application.properties" <<EOF
exomiser.data-directory=$W/exomiser/data
exomiser.hg38.data-version=$DATA_VER
exomiser.phenotype.data-version=$DATA_VER
exomiser.hg38.remm-path=$W/exomiser/ReMM.v0.4.hg38.tsv.gz
remm.version=0.4
logging.level.com.zaxxer.hikari=ERROR
EOF

job() {   # $1 = vcf, $2 = output name
cat <<EOF
---
phenopacket:
  id: $SAMPLE
  subject: { id: $SAMPLE, sex: MALE }
  phenotypicFeatures:
    - type: { id: HP:0002859, label: Rhabdomyosarcoma }
    - type: { id: HP:0000121, label: Nephrocalcinosis }
    - type: { id: HP:0004322, label: Short stature }
    - type: { id: HP:0001508, label: Failure to thrive }
    - type: { id: HP:0003202, label: Skeletal muscle atrophy }
    - type: { id: HP:0001622, label: Premature birth }
    - type: { id: HP:0001518, label: Small for gestational age }
    - type: { id: HP:0200067, label: Recurrent spontaneous abortion }
  htsFiles:
    - { uri: $1, htsFormat: VCF, genomeAssembly: hg38 }
  metaData: { created: '2026-10-15T00:00:00Z', createdBy: mva-hackathon-2026, phenopacketSchemaVersion: 1.0 }
analysis:
  # examples/preset-genome-analysis.yml (exomiser_model_4), unchanged
  analysisMode: PASS_ONLY
  inheritanceModes: { AUTOSOMAL_DOMINANT: 0.1, AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1, AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
                      X_DOMINANT: 0.1, X_RECESSIVE_HOM_ALT: 0.1, X_RECESSIVE_COMP_HET: 2.0, MITOCHONDRIAL: 0.2 }
  frequencySources: [ UK10K, GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_EAS, GNOMAD_E_NFE, GNOMAD_E_SAS,
                      GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_EAS, GNOMAD_G_NFE, GNOMAD_G_SAS ]
  pathogenicitySources: [ REVEL, MVP, ALPHA_MISSENSE, REMM, SPLICE_AI ]
  steps: [
    hiPhivePrioritiser: { },
    priorityScoreFilter: { priorityType: HIPHIVE_PRIORITY, minPriorityScore: 0.501 },
    geneBlacklistFilter: { },
    failedVariantFilter: { },
    alleleBalanceFilter: { },
    regulatoryFeatureFilter: { },
    frequencyFilter: { maxFrequency: 2.0 },
    pathogenicityFilter: { keepNonPathogenic: false, target: NON_CODING },
    inheritanceFilter: { },
    omimPrioritiser: { }
  ]
outputOptions:
  outputContributingVariantsOnly: false
  numGenes: 0
  outputDirectory: $W/results
  outputFileName: $2
  outputFormats: [ HTML, JSON, TSV_GENE, TSV_VARIANT ]
EOF
}

if step exomiser; then
  until [[ -e "$W/done/remm" ]]; do log "waiting for REMM download ($(du -h "$W/exomiser/ReMM.v0.4.hg38.tsv.gz" 2>/dev/null | cut -f1) of ~16G)"; sleep 300; done
  job "$W/reads/$SAMPLE"_HGWCNDSX7.vcf.gz proband-genome-remm > "$W/exomiser/snv.yml"
  job "$W/results/diploidSV.vcf.gz"       proband-sv          > "$W/exomiser/sv.yml"
  cd "$EXO"
  for j in snv sv; do
    java -Xmx48g -jar "exomiser-cli-$EXO_VER.jar" analyse --analysis "$W/exomiser/$j.yml" > "$W/results/exomiser_$j.log" 2>&1
  done
  cd "$W"
  ok exomiser
fi

# ---------------------------------------------------------------- focused summary for the candidate genes
if step summary; then
  {
    echo "# SV / CNV calls overlapping MVA genes (+/- 100 kb)"
    for g in "BUB1B chr15:40100000-40330000" "CEP57 chr11:95650000-95910000" "TRIP13 chr5:790000-1000000"; do
      set -- $g; echo "## $1 ($2)"
      bcftools view -H -r "$2" "$W/results/diploidSV.vcf.gz" | cut -f1-7 || true
      echo "cnvpytor (10 kb bins):"
      awk -v r="$2" 'BEGIN{split(r,a,"[:-]")} {split($2,b,"[:-]"); if (b[1]==a[1] && b[3]>=a[2] && b[2]<=a[3]) print}' \
          "$W/results/cnvpytor_calls_10000.tsv" || true
    done
    echo; echo "# Manta PASS SV counts by type"
    bcftools view -H -f PASS "$W/results/diploidSV.vcf.gz" | grep -o 'SVTYPE=[A-Z]*' | sort | uniq -c
    echo; echo "# Exomiser top 10 (SNV/indel with REMM)"; head -11 "$W/results/proband-genome-remm.genes.tsv" | cut -f1-6
    echo; echo "# Exomiser top 10 (SVs)";               head -11 "$W/results/proband-sv.genes.tsv" | cut -f1-6
  } > "$W/results/SUMMARY.txt" 2>&1
  ok summary
fi

tar -czf "$HOME/mva_results.tar.gz" -C "$W" results
log "ALL DONE -> $HOME/mva_results.tar.gz ($(du -h "$HOME/mva_results.tar.gz" | cut -f1))"
cat "$W/results/SUMMARY.txt"
