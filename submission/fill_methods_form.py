# Fill the Track 1 methods form (sheet "Track 1 methods", answers in column B).
# Run: uv run --python 3.12 --with openpyxl python submission/fill_methods_form.py
from pathlib import Path

import openpyxl

HERE = Path(__file__).parent
SRC = HERE / "methods_description_form.xlsx"
OUT = HERE / "methods_description_form_filled.xlsx"

APPROACH = """\
1. Input: proband short-read WGS VCF (WGS_EX2312012, GRCh38, singleton, no parental data) and the 8 HPO terms from the clinical phenotype (HP:0002859, HP:0000121, HP:0004322, HP:0001508, HP:0003202, HP:0001622, HP:0001518, HP:0200067).
2. QC: per-chromosome depth, B-allele-frequency histogram, runs of homozygosity (bcftools roh) and 1 Mb read-depth bins from a re-aligned read subsample (minimap2 2.31). Sex verified as male from X/Y genotypes. Bulk aneuploidy screen: per-chromosome GC-corrected depth ratio (1.00 on all autosomes) and B-allele-frequency deviation of heterozygous PASS SNPs with centromeric/acrocentric windows masked (6-8% off-centre on every autosome). No clonal whole-chromosome gain or loss is detectable in blood WGS; note that variegated (cell-to-cell random) aneuploidy averages out in bulk sequencing and is not expected to be visible here. Used as sanity checks; not used to score variants.
3. Annotation: Ensembl VEP 115.2 (Docker, offline cache 115, GRCh38) with HGVS, MANE/canonical, SIFT, PolyPhen, gnomAD exome and genome allele frequencies, and ClinVar via --check_existing.
4. Prioritisation: Exomiser 15.1.0 (data release 2512, hg38), PASS-only, removing UTR/intronic/intergenic/regulatory effects; maximum population frequency 2%, with mode-specific cutoffs of 0.1% for dominant/hom-alt and 2% for compound het; pathogenicity sources REVEL, MVP, AlphaMissense and SpliceAI; OMIM and hiPHIVE phenotype prioritisers. Exomiser scores all inheritance modes, including compound-heterozygous pairs.
5. Targeted review: every rare (max gnomAD AF <1%) MANE/canonical variant from the VEP output was screened for (a) ClinVar pathogenic/likely pathogenic, (b) HIGH/MODERATE impact in a phenotype-driven gene panel (spindle assembly checkpoint / MVA genes, rhabdomyosarcoma predisposition and RASopathy genes, nephrocalcinosis genes, primordial dwarfism and DNA repair genes) and (c) ACMG SF v3.2 genes for secondary findings.
6. Curation and ranking: candidates were kept or dropped by genotype, read support, allele frequency, predicted effect, gene-disease fit and inheritance. Apparent compound-het pairs whose two variants lie a few bp apart (FANCD2, SERPINA1, MTCH2, ADAMTS1) were treated as cis/single events, not two hits. Recurrent artefact-prone loci (HLA, MUC4, KRTAP, olfactory receptors) were excluded. Every call was checked against the VCF (PASS, genotype, allele depths). epcr is a manually assigned probability that the row is causal (primary) or a true reportable finding (secondary). Rows are listed primary first, by decreasing epcr, then secondary."""

AI = ("Anthropic Claude Code CLI, model Claude Opus 5.5, [TODO: plan/tier, e.g. Claude Max / Team / API], "
      "[TODO: data setting, e.g. commercial terms, no training on customer content]. Used as a coding and analysis assistant "
      "(writing pipeline scripts, querying annotation output, drafting curation notes and this form). "
      "All variant calls and final rankings were checked against the underlying VCF and annotations.")

AUTOMATED = ("Semi-automated. Candidates were generated automatically (VEP + Exomiser + scripted filters), "
             "then the submission file was produced after downstream manual review and curation.")

MANUAL = ("Manual review of the Exomiser ranking and filtered VEP candidates: (1) removed artefact-prone loci (HLA class I/II, MUC4, KRTAP, "
          "olfactory receptors); (2) re-evaluated Exomiser compound-het calls and rejected pairs whose variants lie within a few bp "
          "(likely the same haplotype/complex event); (3) checked gene-disease fit against the phenotype using OMIM/ClinVar knowledge "
          "(BUB1B biallelic loss -> mosaic variegated aneuploidy type 1 with rhabdomyosarcoma and growth failure); (4) confirmed every "
          "submitted variant in the VCF (PASS, 0/1, allele depths about 50%); (5) assigned epcr values by hand and separated primary from "
          "secondary findings. Rows 5-7 are the next-best non-artefact Exomiser candidates, kept at low epcr to complete the top 10.")

DATA_PUBLIC = "Only publicly available data were used."

PUBLIC = ("GRCh38 reference; Ensembl VEP 115 cache (Ensembl/GENCODE gene models incl. MANE Select, SIFT, PolyPhen-2, "
          "gnomAD v4 exome and genome allele frequencies, ClinVar and dbSNP identifiers); Exomiser 2512 data release "
          "(gnomAD, UK10K frequencies; REVEL, MVP, AlphaMissense and SpliceAI scores; ClinVar; OMIM/Orphanet disease-gene associations; "
          "HPO; mouse and zebrafish phenotype data for hiPHIVE); gnomAD gene constraint (LOEUF); ACMG SF v3.2 gene list; "
          "published literature on BUB1B/MVA1 and spindle assembly checkpoint disorders.")

PROPRIETARY = "None."

COMPHET = ("Yes. Our approach outputs compound-heterozygous candidate pairs (Exomiser AR comp-het mode plus manual review). The top-ranked "
           "row is a BUB1B pair. Limitation: with singleton short-read data and variants ~11 kb apart, phase (trans) cannot be confirmed; "
           "pairs within a few bp were treated as cis and not submitted as pairs.")

SECONDARY = ("Secondary/incidental findings were identified by screening rare variants for ClinVar pathogenic/likely pathogenic calls "
             "and by checking ACMG SF v3.2 genes (none reportable there). Heterozygous P/LP variants unrelated to the primary phenotype were "
             "submitted as finding_type=secondary: LZTR1 p.Tyr748Ter (schwannomatosis predisposition / Noonan 2 carrier), GNRHR p.Arg139His "
             "(AR hypogonadotropic hypogonadism carrier) and FLG p.Arg2447Ter (ichthyosis vulgaris risk). These are listed after all primary candidates.")

RUNTIME = ("Local laptop (Apple M2 Max, 12 cores, 32 GB RAM); no cloud cost. Exomiser: ~42 s after data unpacking. "
           "VEP whole-genome annotation of ~5.0 M records: ~36 h wall time (Docker, --fork 10). Data/cache downloads: several hours. "
           "Manual curation and review: ~[TODO: hours] person-hours.")

ABSTRACT = """\
We ranked candidate variants for a single proband with short-read whole-genome sequencing and a phenotype of rhabdomyosarcoma, nephrocalcinosis, short stature, failure to thrive, muscle atrophy, prematurity/small for gestational age and a family history of recurrent miscarriage. The pipeline pairs standard open-source tools with targeted expert-style curation.

Methods. The VCF (GRCh38) was annotated with Ensembl VEP 115 (MANE/canonical transcripts, HGVS, SIFT, PolyPhen-2, gnomAD v4 exome/genome frequencies, ClinVar). Genes were prioritised with Exomiser 15.1.0 (data release 2512) using the eight HPO terms, PASS-only calls, coding/splice effects, inheritance-specific frequency cutoffs, REVEL/MVP/AlphaMissense/SpliceAI pathogenicity and the OMIM and hiPHIVE phenotype prioritisers. Rare variants (<1% in gnomAD) were then screened for ClinVar pathogenic calls, for HIGH/MODERATE impact in phenotype-driven gene panels (spindle assembly checkpoint/MVA, rhabdomyosarcoma predisposition and RASopathy, nephrocalcinosis, growth/DNA repair) and in ACMG SF v3.2 genes. Candidates were curated for genotype quality, allele balance, frequency, predicted effect, gene-disease fit and inheritance plausibility. A manually assigned epcr reflects our estimated probability of causality.

Result. The top candidate is a presumed compound-heterozygous pair in BUB1B: c.2210T>G p.(Leu737*) (ClinVar pathogenic/likely pathogenic, gnomAD 7.9e-05) and c.3006T>G p.(Asn1002Lys), a novel missense in the kinase domain (absent from gnomAD, AlphaMissense 0.92, SIFT deleterious, PolyPhen 0.997). Both are high-quality heterozygous PASS calls with balanced allele depths. Biallelic BUB1B loss causes mosaic variegated aneuploidy type 1 (MVA1), which fits the embryonal rhabdomyosarcoma, pre- and postnatal growth failure and the reproductive history. A truncating allele plus a hypomorphic kinase-domain missense is the classic MVA1 genotype. Exomiser ranked BUB1B first (p = 0.0005, ~30-fold ahead of the next gene). Lower-ranked primary rows include a possible nephrocalcinosis contributor (SLC34A1 p.Ala133Val), a same-pathway heterozygous variant (MAD1L1 p.Arg59Cys) and the best remaining Exomiser candidates. Three heterozygous P/LP variants unrelated to the phenotype are reported as secondary findings (LZTR1, GNRHR, FLG).

Strengths. Uses only public data and open-source tools and is reproducible on a laptop. It combines phenotype-driven whole-genome ranking with curation that removes common false positives: artefact-prone loci (HLA, MUC4, KRTAP, olfactory receptors) and pseudo-compound-het pairs whose variants lie a few bp apart and are almost certainly in cis. Compound-heterozygous pairs are output natively. Every submitted variant was checked in the VCF.

Limitations. Singleton data: the BUB1B variants are ~11 kb apart, so trans phase cannot be established from short reads. Parental testing or long-read sequencing is needed. The missense allele has no functional data. Structural variants, copy-number changes, repeat expansions, mitochondrial and deep-intronic/regulatory variants were not systematically assessed, since the Exomiser run excluded intronic/UTR effects. The ranking depends on curated knowledge and manually assigned epcr values, so it is not fully automated. Lower-ranked rows have low prior probability and mainly complete the top-10 list. A bulk-WGS screen (GC-corrected depth and allele-balance per chromosome, centromeres masked) found no clonal whole-chromosome gain or loss, but variegated aneuploidy averages out in bulk data, so this neither supports nor refutes MVA; karyotyping is needed."""

ANSWERS = {
    7: "[TODO: team name]",
    8: 1,
    9: APPROACH,
    10: AI,
    11: AUTOMATED,
    12: MANUAL,
    13: DATA_PUBLIC,
    14: PUBLIC,
    15: PROPRIETARY,
    16: COMPHET,
    17: SECONDARY,
    18: RUNTIME,
    19: ABSTRACT,
}

wb = openpyxl.load_workbook(SRC)
ws = wb["Track 1 methods"]
for row, value in ANSWERS.items():
    cell = ws.cell(row=row, column=2, value=value)
    cell.alignment = openpyxl.styles.Alignment(wrap_text=True, vertical="top")
wb.save(OUT)
print(f"wrote {OUT}; abstract words = {len(ABSTRACT.split())}")
