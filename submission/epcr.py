#!/usr/bin/env python3
"""Compute epcr (estimated probability of causal relationship) for each submission row.

    epcr = P_variant * P_gene * P_inheritance

P_variant     Tavtigian et al. 2018/2020 Bayesian ACMG framework: prior 0.10, each
              point multiplies the odds by 350**(1/8); product over both variants of a pair.
P_gene        gene-disease validity for the phenotype in question (ClinGen-style tier).
P_inheritance probability the inheritance requirement is met in this child
              (trans phase for a pair; second-hit chance for a lone het in a recessive
              gene; penetrance for a dominant secondary finding).

Every number below is an explicit, tweakable assumption. Rerun after editing:
    python3 submission/epcr.py            # prints audit table, rewrites CSV epcr column
See EPCR_METHOD.md for the reasoning behind each choice.
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
CSV = HERE / "track1_submission.csv"
AUDIT = HERE / "epcr_audit.tsv"

# --- Tavtigian point system (Hum Mutat 2020) -------------------------------
PRIOR = 0.10
ODDS_PER_POINT = 350 ** (1 / 8)        # very strong = 8 points = odds x350
POINTS = {
    "PVS1": 8, "PS": 4, "PM": 2, "PP": 1,   # pathogenic: very strong / strong / moderate / supporting
    "BS": -4, "BP": -1,                     # benign: strong / supporting
}


def p_variant(criteria):
    """criteria: list of (code, strength) e.g. [("PVS1","PVS1"),("PM2","PP")]."""
    pts = sum(POINTS[s] for _, s in criteria)
    odds = PRIOR / (1 - PRIOR) * ODDS_PER_POINT ** pts
    return pts, odds / (1 + odds)


# --- gene-disease validity tiers (ClinGen-style) ---------------------------
GENE = {"definitive": 1.0, "strong": 0.9, "moderate": 0.7, "limited": 0.3, "none": 0.02}

# --- rows: keyed by chrom_1:pos_1 as in the CSV ------------------------------
# criteria strengths: PVS1 | PS | PM | PP | BS | BP
ROWS = {
    "chr15:40209701": dict(
        label="BUB1B p.Leu737* + p.Asn1002Lys (pair)",
        variants=[
            [("PVS1", "PVS1"),          # nonsense, exon 17/23, NMD expected, LoF = known mechanism
             ("PM2", "PP")],            # gnomAD 7.9e-5 (very rare; supporting per ClinGen SVI)
            [("PM2", "PP"),             # absent from gnomAD
             ("PM1", "PM"),             # kinase domain (TWEAK: needs hotspot evidence to keep moderate)
             ("PP3", "PP"),             # AlphaMissense 0.92, SIFT 0.01, PolyPhen 0.997
             ("PP4", "PP")],            # phenotype specific for MVA (rhabdomyosarcoma + growth failure)
        ],
        gene="definitive",              # BUB1B - MVA type 1, OMIM 257300
        inheritance=0.9,                # P(trans). Prior 0.5; phenotype match updates to ~0.9-0.95.
        note="PM3 deliberately NOT used: phase is folded into the inheritance term instead",
    ),
    "chr5:177386432": dict(
        label="SLC34A1 p.Ala133Val (het)",
        variants=[[("BS1", "BS"),       # gnomAD 0.4%: too common for a dominant rare disease
                   ("PP3", "PP")]],     # SIFT 0, PolyPhen 0.992
        gene="limited",                 # heterozygous SLC34A1 -> AD nephrolithiasis/nephrocalcinosis is disputed
        inheritance=1.0,                # het satisfies AD
    ),
    "chr7:2225526": dict(
        label="MAD1L1 p.Arg59Cys (het)",
        variants=[[]],                  # 0.46% (acceptable for AR), predictors split, ClinVar P/B conflict -> 0 pts
        gene="limited",                 # biallelic MAD1L1 MVA-like: case reports only
        inheritance=0.05,               # lone het in a recessive gene: needs an undetected second hit
    ),
    "chr3:10046723": dict(
        label="FANCD2 c.1278+1del (het)",
        variants=[[("BS", "BS")]],      # ClinVar benign; 30% allele fraction; pseudogene region
        gene="definitive",              # FANCD2 - Fanconi anaemia
        inheritance=0.05,               # het only: adjacent indel is in cis (GATK phase set 10046720)
    ),
    "chr17:39175193": dict(
        label="CACNB1 stop-lost (het)",
        variants=[[("PM2", "PP"), ("PM4", "PM")]],   # rare; protein-length change (stop-loss is not PVS1)
        gene="none", inheritance=1.0,
    ),
    "chr16:68021897": dict(
        label="DDX28 p.Arg436Trp (het)",
        variants=[[("PM2", "PP"), ("PP3", "PP")]],
        gene="none", inheritance=1.0,
    ),
    "chr10:17149087": dict(
        label="TRDMT1 p.Leu377Phe (het)",
        variants=[[("PM2", "PP"), ("PP3", "PP")]],
        gene="none", inheritance=1.0,
    ),
    # --- secondary findings: inheritance term = P(phenotype expressed in this child) ---
    "chr22:20996720": dict(
        label="LZTR1 p.Tyr748* (het) - secondary",
        variants=[[("PVS1", "PVS1"), ("PM2", "PP")]],   # gnomAD 1.4e-6; ClinVar P/LP
        gene="definitive",              # LZTR1 - schwannomatosis (AD)
        inheritance=0.4,                # TWEAK: penetrance of LZTR1 LoF for schwannomatosis (literature: incomplete)
    ),
    "chr4:67753920": dict(
        label="GNRHR p.Arg139His (het) - secondary",
        variants=[[("PS3", "PS"), ("PM2", "PP"), ("PP3", "PP")]],   # functionally characterised; 2e-4; MVP 0.985
        gene="definitive",              # GNRHR - AR hypogonadotropic hypogonadism
        inheritance=0.02,               # carrier: phenotype only if an undetected second allele exists
    ),
    "chr1:152307547": dict(
        label="FLG p.Arg2447* (het) - secondary",
        variants=[[("PVS1", "PVS1")]],  # 0.4% but FLG null alleles are common; no frequency penalty either way
        gene="definitive",              # FLG - ichthyosis vulgaris (semi-dominant)
        inheritance=0.15,               # TWEAK: penetrance of one FLG null allele
    ),
}


def compute():
    out = {}
    for key, r in ROWS.items():
        pv = 1.0
        detail = []
        for crit in r["variants"]:
            pts, p = p_variant(crit)
            pv *= p
            detail.append(f"{'+'.join(c for c, _ in crit) or 'none'}={pts}pts->{p:.3f}")
        pg = GENE[r["gene"]]
        epcr = pv * pg * r["inheritance"]
        out[key] = dict(label=r["label"], variant=" ; ".join(detail), p_variant=pv,
                        gene=f"{r['gene']}={pg}", p_inh=r["inheritance"],
                        epcr=max(round(epcr, 3), 0.001))
    return out


def main():
    res = compute()
    rows = list(csv.reader(open(CSV)))
    hdr = rows[0]
    ic, ip, ie = hdr.index("chrom_1"), hdr.index("pos_1"), hdr.index("epcr")
    for r in rows[1:]:
        key = f"{r[ic]}:{r[ip]}"
        if key in res:
            r[ie] = f"{res[key]['epcr']:.3f}"
    csv.writer(open(CSV, "w", newline=""), lineterminator="\n").writerows(rows)

    with open(AUDIT, "w") as f:
        f.write("row\tvariant_evidence\tP_variant\tP_gene\tP_inheritance\tepcr\n")
        for key, v in res.items():
            f.write(f"{v['label']}\t{v['variant']}\t{v['p_variant']:.3f}\t{v['gene']}\t{v['p_inh']}\t{v['epcr']:.3f}\n")
    for key, v in res.items():
        print(f"{v['epcr']:.3f}  {v['label']:42s} Pvar={v['p_variant']:.3f} ({v['variant']})  gene {v['gene']}  inh {v['p_inh']}")
    print(f"wrote {CSV.name} and {AUDIT.name}")


if __name__ == "__main__":
    main()
