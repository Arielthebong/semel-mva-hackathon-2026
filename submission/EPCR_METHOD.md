# How epcr is computed

`epcr` in `track1_submission.csv` is the estimated probability of a causal relationship between the row's variant(s) and the child's phenotype. It is **computed by `epcr.py`**, not typed in, so that every input can be inspected and changed.

## The formula

```
epcr = P_variant × P_gene × P_inheritance
```

Three questions have to be answered "yes" for a variant to be causal, and each factor is the probability of one of them.

### 1. P_variant — is the variant actually pathogenic?

We use the Bayesian form of the ACMG/AMP classification rules (Tavtigian et al. 2018, *Genet Med* 20:1054; point system in Tavtigian et al. 2020, *Hum Mutat* 41:1734). Each ACMG criterion contributes points by strength:

| Strength | Points | Odds multiplier |
|---|---|---|
| Very strong (PVS1) | +8 | ×350 |
| Strong (PS) | +4 | ×18.7 |
| Moderate (PM) | +2 | ×4.33 |
| Supporting (PP, PM2_Supporting) | +1 | ×2.08 |
| Benign supporting (BP) | −1 | ÷2.08 |
| Benign strong (BS) | −4 | ÷18.7 |

Starting from the framework's prior of 0.10 that a rare variant is pathogenic:

```
odds      = (0.10 / 0.90) × 2.08^points
P_variant = odds / (1 + odds)
```

This reproduces the ACMG categories: ≥10 points → ≥0.99 (Pathogenic), 6–9 → 0.90–0.99 (Likely pathogenic), 0–5 → 0.10–0.81 (VUS), ≤−7 → <0.001 (Benign). For a compound-heterozygous row, P_variant is the product for both variants, since both must be pathogenic.

### 2. P_gene — does this gene cause this phenotype?

ClinGen-style gene–disease validity, tiered: definitive 1.0, strong 0.9, moderate 0.7, limited 0.3, no established association 0.02. Exomiser's combined score is *not* used here: it is a logistic-regression score trained to separate pathogenic from benign variants, not a probability that the gene is this patient's diagnosis.

### 3. P_inheritance — is the inheritance requirement met in this child?

- Compound-het pair, phase unknown: P(the two variants are on different copies). The prior is 0.5 for two independent rare alleles; observing a phenotype that matches the gene updates this strongly (if they were in cis the child would be an unaffected carrier and the phenotype would have to be a coincidence). We use **0.9**. PM3 ("in trans with a pathogenic variant") is deliberately *not* also scored as an ACMG criterion, to avoid counting phase twice.
- Single het in a recessive gene: probability of an undetected second hit (deep-intronic, structural) — **0.05**; for a known carrier state reported as secondary, **0.02**.
- Het in a dominant gene: 1.0 for primary candidates; for secondary findings this term is the **penetrance** — the probability the child actually develops the associated condition (LZTR1 schwannomatosis 0.4, one FLG null allele 0.15; both flagged for literature refinement).

## Worked example: the BUB1B pair

| Variant | Criteria | Points | P |
|---|---|---|---|
| p.Leu737* | PVS1 (+8), PM2_Supporting (+1) | 9 | 0.988 |
| p.Asn1002Lys | PM2_Supporting (+1), PM1 kinase domain (+2), PP3 (+1), PP4 (+1) | 5 | 0.812 |

P_variant = 0.988 × 0.812 = 0.802; P_gene = 1.0 (BUB1B–MVA1 definitive); P_inheritance = 0.9.
**epcr = 0.72.** The number is limited almost entirely by the novel missense, which is a strong VUS rather than likely pathogenic under strict ACMG scoring. Confirming trans phase (parental testing or long reads) would raise the pair to ~0.80; a functional assay on Asn1002Lys would push it above 0.9.

## Where judgement still enters

The framework does not remove judgement; it moves it into named, standardized decisions a reader can check one by one: which ACMG criteria apply (e.g. whether PM1 is justified for the BUB1B kinase domain), which validity tier a gene–disease pair earns, and the penetrance and phase values. All of these are literal constants in `epcr.py`, each with a comment. Items marked `TWEAK` are the ones we consider least settled.

## Files

- `epcr.py` — inputs and computation; rewrites the `epcr` column of `track1_submission.csv`.
- `epcr_audit.tsv` — one line per row with every factor, for the judges.
