# semel — Rare Disease, Real Kid: MVA Hackathon 2026, Track 1

**Result:** compound-heterozygous BUB1B c.2210T>G p.(Leu737*) + c.3006T>G p.(Asn1002Lys) → mosaic variegated aneuploidy type 1 (epcr 0.72).

| File | What it is |
|---|---|
| `submission/semel_exomiser-lirical-bayesian-acmg.csv` | Submitted predictions (template format, 10 rows) |
| `submission/semel_track1_report.md` | Submitted report: methods, evidence, limitations, AI disclosure |
| `submission/epcr.py`, `EPCR_METHOD.md`, `epcr_audit.tsv`, `epcr_method.pdf` | How epcr is computed (Bayesian ACMG × gene validity × inheritance), with per-row audit |
| `submission/track1_submission.csv` | Working copy with extra `gene` columns (not the submitted file) |
| `submission/fill_methods_form.py` → `methods_description_form_filled.xlsx` | Organisers' methods template, filled (earlier draft; the report supersedes it) |
| `run_vep.sh` | Ensembl VEP 115 annotation (Docker) |
| `run_exomiser.sh <yml>` + `analysis/*.yml` | Exomiser 15.1.0 runs: coding, whole-genome + SpliceAI, + REMM, SV |
| `run_lirical.sh [terms-to-drop]` | LIRICAL 2.4.1 (includes the HGNC-file repair it needs) |
| `run_local_sv.sh` | Reads → minimap2 → Manta + CNVpytor → Exomiser on SVs (resumable; macOS) |
| `cloud/` | Same SV pipeline for one AWS VM |
| `vq` | Query the VCF by gene or region |

## Reproducing
Request access to `SageBio/mva-hackathon-2026-data` on Hugging Face and place the VCF in `data/`. Tools and databases (VEP cache, Exomiser 15.1.0 + data 2512, LIRICAL 2.4.1, REMM v0.4) go in `tools/` and `data/` as described at the top of each script. Then run, in order: `run_vep.sh`, `run_exomiser.sh proband_job.yml`, `run_exomiser.sh proband_genome.yml`, `run_lirical.sh`, `run_local_sv.sh`, and `submission/epcr.py`.

No patient data is in this repository. Code and reports are released under CC BY 4.0, per the hackathon rules.
