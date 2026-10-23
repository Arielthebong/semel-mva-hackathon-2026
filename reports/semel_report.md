# Track 1 report — team semel

**Rare Disease, Real Kid: The MVA Hackathon 2026 · Track 1 (variant prediction)**
Proband WGS_EX2312012 (GRCh38, singleton, male) · Predictions file: `semel_predictions.csv` · Code: https://github.com/Arielthebong/semel-mva-hackathon-2026

**AI disclosure.** Anthropic Claude Code (model Claude Opus 5.5), Claude Max plan, training opt-out enabled, no feedback provided. Claude was used as a coding and analysis assistant: writing the pipeline scripts, querying annotation output, drafting curation notes and this report. Every submitted variant and every number in this report was checked against the VCF, the annotation output or the tool logs.

## Answer

Our primary candidate is a **compound-heterozygous pair in BUB1B** (NM_001211.6), causing **mosaic variegated aneuploidy type 1 (MVA1, OMIM 257300)**. We estimate the probability that this causal relationship is correct (epcr) at **0.72**.

| | Variant 1 | Variant 2 |
|---|---|---|
| Position (GRCh38) | chr15:40209701 T>G | chr15:40220612 T>G |
| HGVS | c.2210T>G p.(Leu737*) | c.3006T>G p.(Asn1002Lys) |
| Effect | stop-gain, exon 17/23 | missense, last exon, kinase domain |
| Population frequency | gnomAD exomes 7.9e-05 | absent from gnomAD |
| ClinVar / predictors | pathogenic / likely pathogenic | AlphaMissense 0.92, SIFT deleterious 0.01, PolyPhen-2 0.997 |
| Genotype (AD ref/alt) | 0/1 (21/25), PASS | 0/1 (15/13), PASS |

A truncating allele combined with a hypomorphic kinase-domain missense is the classic MVA1 genotype. MVA1 explains the embryonal rhabdomyosarcoma, the pre- and postnatal growth failure (SGA, short stature, failure to thrive) and is compatible with the family's reproductive history.

BUB1B ranked first in every analysis we ran:

| Analysis | BUB1B result | Next best |
|---|---|---|
| Exomiser 15.1.0, coding (hiPHIVE) | rank 1, p = 0.0005 | FANCD2, p = 0.015 |
| Exomiser, whole genome + SpliceAI | rank 1, p = 0.0004 | FANCD2 |
| Exomiser, whole genome + SpliceAI + REMM | rank 1, p = 0.0004 | FANCD2 |
| Same, without HP:0200067 (recurrent miscarriage) | rank 1, p = 0.0004 | FANCD2 |
| LIRICAL 2.4.1, all 8 HPO terms | MVA1 is disease #1, post-test probability 18.4%, genotype LR 1240 | — |
| LIRICAL, without HP:0200067 | MVA1 #1, post-test probability 95.8% | MDDS13/FBXL4, 12.6% |

HP:0200067 (recurrent spontaneous abortion) describes the mother's history rather than the child's phenotype. Removing it raises LIRICAL's confidence sharply, because neither it nor muscle atrophy is annotated to MVA1. We kept all eight terms in the scored runs and report both.

## Data

- Input: proband short-read WGS VCF (`WGS_EX2312012_HGWCNDSX7.vcf.gz`), the raw FASTQ reads (4 lanes), and eight HPO terms: HP:0002859 rhabdomyosarcoma, HP:0000121 nephrocalcinosis, HP:0004322 short stature, HP:0001508 failure to thrive, HP:0003202 muscle atrophy, HP:0001622 premature birth, HP:0001518 small for gestational age, HP:0200067 recurrent spontaneous abortion.
- Only public resources were used: GRCh38; Ensembl VEP 115 cache (MANE, SIFT, PolyPhen-2, gnomAD v4, ClinVar, dbSNP); Exomiser data release 2512 (gnomAD/UK10K frequencies, REVEL, MVP, AlphaMissense, SpliceAI, ClinVar, OMIM/Orphanet, HPO, model-organism phenotypes); REMM v0.4; LIRICAL 2.4.1 data; gnomAD-SV v4; ACMG SF v3.2; published literature.

## Pipeline

1. **QC.** Per-chromosome depth, B-allele-frequency (BAF) histograms, runs of homozygosity (bcftools roh), sex check (male from X/Y).
2. **Aneuploidy screen.** GC-corrected depth ratio was 1.00 on every autosome, and BAF deviation was 6–8% everywhere once centromeric and acrocentric windows were masked. No clonal whole-chromosome gain or loss. Variegated aneuploidy averages out in bulk blood WGS, so this neither supports nor refutes MVA.
3. **Annotation.** Ensembl VEP 115.2 (offline cache, GRCh38) on all ~5.0 M records.
4. **Phenotype-driven prioritisation.** Exomiser 15.1.0 (data 2512, hg38) was run in three configurations (`analysis/*.yml`):
    - coding only;
    - whole genome with the allele-balance filter and SpliceAI;
    - whole genome with REMM regulatory scores added.
5. **Independent likelihood-ratio ranking.** LIRICAL 2.4.1, reusing the Exomiser 2512 variant and ClinVar databases. The current HGNC file breaks LIRICAL 2.4.1 silently: every disease gets a null genotype likelihood ratio. `run_lirical.sh` repairs the file before each run.
6. **Structural and copy-number variants.** The raw reads were aligned with minimap2 2.31 (lanes L001–L002, about 22× after duplicate marking; 98.4% mapped). Structural variants were called with Manta 1.6.0 and copy-number variants with CNVpytor 1.3.1 (10 kb and 100 kb bins). Exomiser was then run on the Manta calls.
7. **Targeted review.** Every rare (gnomAD < 1%) variant was screened for:
    - ClinVar pathogenic or likely pathogenic status;
    - HIGH or MODERATE impact in phenotype panels (MVA and spindle-checkpoint genes, rhabdomyosarcoma predisposition and RASopathies, nephrocalcinosis, growth and DNA repair);
    - presence in an ACMG SF v3.2 gene.
8. **Curation and epcr** (next section).

## How epcr is computed

epcr = P_variant × P_gene × P_inheritance, computed by `submission/epcr.py`. Every row's inputs are listed in `submission/epcr_audit.tsv`.

- **P_variant** is the Bayesian ACMG/AMP posterior (Tavtigian et al. 2018, *Genet Med*; 2020, *Hum Mutat*). The prior is 0.10 and each point multiplies the odds by 2.08. The scale gives Very strong = 8, Strong = 4, Moderate = 2, Supporting = 1, and benign criteria are negative. For a pair, P_variant is the product of both variants.
- **P_gene** is the ClinGen-style gene–disease validity: definitive 1.0, strong 0.9, moderate 0.7, limited 0.3, none 0.02.
- **P_inheritance** is the probability that the inheritance requirement is met:
    - trans phase for a pair: 0.9 for BUB1B (a matching phenotype raises the 0.5 prior);
    - an undetected second hit for a lone heterozygote in a recessive gene: 0.05;
    - penetrance for dominant secondary findings: LZTR1 0.4, FLG 0.15.
    PM3 is not also scored, to avoid counting phase twice.

**Worked example (BUB1B):**
- p.Leu737* scores PVS1 + PM2_Supporting = 9 points, giving 0.988.
- p.Asn1002Lys scores PM2_Supporting + PM1 + PP3 + PP4 = 5 points, giving 0.812.
- epcr = 0.988 × 0.812 × 1.0 × 0.9 = **0.72**.

The novel missense is the limiting term. Without PM1, the pair would score 0.44. Confirming trans phase would raise it to about 0.80, and a functional assay on Asn1002Lys could take it above 0.9.

## Submitted rows

| Rank | Gene | Variant | Type | epcr | Rationale |
|---|---|---|---|---|---|
| 1 | BUB1B | p.Leu737* + p.Asn1002Lys | primary | 0.722 | Biallelic loss → MVA1 |
| 2 | CACNB1 | chr17:39175193 T>C, stop-lost | primary | 0.010 | Weak phenotype match; no established disease |
| 3 | DDX28 | chr16:68021897 T>A, p.Arg436Trp | primary | 0.006 | Rare missense; no established disease |
| 4 | TRDMT1 | chr10:17149087 G>A, p.Leu377Phe | primary | 0.006 | Rare missense; no established disease |
| 5 | SLC34A1 | chr5:177386432 C>T, p.Ala133Val | primary | 0.004 | Could contribute to nephrocalcinosis only; gnomAD 0.4% |
| 6 | MAD1L1 | chr7:2225526 G>A, p.Arg59Cys | primary | 0.002 | Same checkpoint pathway; single heterozygote |
| 7 | FANCD2 | chr3:10046723 AG>A, c.1278+1del | primary | 0.001 | Single hit (see below); pseudogene region; ClinVar benign |
| — | LZTR1 | chr22:20996720 C>G, p.Tyr748* | secondary | 0.395 | Schwannomatosis predisposition; Noonan 2 carrier |
| — | FLG | chr1:152307547 G>A, p.Arg2447* | secondary | 0.146 | Ichthyosis vulgaris / atopy risk allele |
| — | GNRHR | chr4:67753920 C>T, p.Arg139His | secondary | 0.018 | Carrier of recessive hypogonadotropic hypogonadism |

No reportable variant was found in an ACMG SF v3.2 gene.

## False positives caught by curation

Each of these looked like a strong candidate in automated output and was rejected on evidence:

- **FANCD2 "compound het"** (Exomiser rank 3). The two deletions share GATK phase set 10046720, so they are on the same copy: one hit, not two. The same check rejected apparent pairs in SERPINA1, MTCH2 and ADAMTS1.
- **KCNJ1** (rank 4 once REMM is added). The variant is a homozygous 10-bp insertion (rs754659348) downstream of the gene, inside a CCCACCCCCC repeat, with REMM 0.98. Exomiser's frequency database lacks it, but gnomAD genomes put it at **69.8%** (49–90% across ancestries). Lesson: cross-check the frequency of every REMM-driven non-coding hit.
- **PRKRA** (the only gene from Exomiser on structural variants). Three heterozygous "deletions" whose breakpoints match PRKRA exon boundaries exactly. This is the signature of a polymorphic processed retrocopy of PRKRA (reads span exon–exon junctions), not a deletion; one carries dbVar AF 6.2%.
- **6q27 duplication.** A 264 kb heterozygous tandem duplication (AFDN/KIF25/FRMD1) was seen by both Manta and CNVpytor. gnomAD-SV v4 lists near-identical duplications at about 0.9% (over 1,000 alleles), so it is benign.
- **Large CNVpytor calls.** All fall in centromeres, heterochromatin, acrocentric p-arms or segmental duplications.

Structural-variant and copy-number calling found **no call within 100 kb of BUB1B, CEP57 or TRIP13**, and no intronic or splice variant in those genes reached SpliceAI ≥ 0.5.

## Limitations

- Singleton data: the BUB1B variants are about 11 kb apart, so trans phase cannot be proven from short reads. Parental testing or long-read sequencing is needed.
- Asn1002Lys has no functional data. PM1 (critical domain) is the one judgement that most moves epcr.
- Variegated aneuploidy cannot be seen in bulk WGS. A karyotype showing a high fraction of aneuploid cells (typically above 25%) would confirm MVA.
- Structural-variant calling used 2 of 4 lanes (about 22×). That is adequate for large deletions and duplications, but less sensitive to small or low-allele-fraction events. Repeat expansions and mitochondrial variants were not assessed.
- Ranks 2–7 complete the top 10 and carry very low epcr. Assigning ACMG criteria, validity tiers and penetrance is a curated step, so the pipeline is semi-automated.

## Reproducibility, run time and cost

All scripts and configs are in the repository (`run_vep.sh`, `run_exomiser.sh <yml>`, `run_lirical.sh`, `run_local_sv.sh`, `analysis/*.yml`, `submission/epcr.py`). Patient data are not included.

Everything ran on a laptop (Apple M2 Max, 12 cores, 32 GB RAM), with no cloud cost:

| Step | Run time |
|---|---|
| VEP, whole genome (Docker, unoptimised) | 36.6 h |
| Each Exomiser run | ~1 min |
| LIRICAL | ~3.5 min |
| minimap2 alignment, 2 lanes (memory-bound because the machine was swapping) | ~10 h |
| Manta | 26 min |
| CNVpytor | 5 min |

A cloud kit for the alignment and SV steps (`cloud/`, c7i.24xlarge) is estimated at about 3.5 h and about $15 per run.

**Data handling.** The data were processed only on the team's own laptop. All data and derived files will be deleted within 30 days of the hackathon's close, and deletion will be confirmed to cnb@sagebase.org.
