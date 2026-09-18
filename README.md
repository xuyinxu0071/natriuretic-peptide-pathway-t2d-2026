# Natriuretic Peptide Pathway as a Genetically Prioritized Entry Point — Analysis Code

**Manuscript:** Barometer or Driver? Residual Vascular Risk in Glycemic-Controlled Type 2 Diabetes and the Natriuretic Peptide Pathway as a Genetically Prioritized Entry Point  
**Journal:** *Cardiovascular Diabetology* — Special Issue "Vascular Drivers of Cardio-Metabolic Disease: Mechanisms and Clinical Implications"  
**Authors:** Yin Xu, Xinmei Wang, Guofeng Wang, Wei Wei\*, Ning Li\*  
**Department:** Department of Geriatrics, Shifan Road, Jinan, China

## What this repository contains

This repository deposits all analysis code underlying the manuscript, including the **ported SMR/HEIDI implementation** and its cross-validation against the official SMR 1.3.1 binary. It supports full reproducibility of the observational (HRS + CHARLS) and Mendelian-randomization (natriuretic-peptide-pathway drug-target) components.

| Folder | Content |
|---|---|
| `observational/` | HRS/CHARLS cohort analyses: penalized/reduced-model Cox, restricted cubic splines (RCS), KM/CIF, fragility check, covariate-audit sensitivity, figure generation |
| `mr/` | Two-sample MR for NPPA/NPPB/NPR1–3/MME/CORIN with IL6R+/CRP− dual calibration; replication; cis-MR; two-step mediation MR; virtual knockout |
| `rs5068_reaudit/` | Seven-dataset functional-variant re-audit of rs5068 and the three-family (CAD/MI, HF, AF) forest plot |
| `colocalization/` | Colocalization (coloc/SuSiE) and the in-house **SMR/HEIDI port** (`step13q_smr_heidi.py`) |
| `verification/` | Cross-validation of the SMR port against the official `smr-1.3.1` binary (`step18d_official_recalc.py`) |
| `figures/` | Manuscript figures (Fig1–5, graphical abstract, key supplementary figures) |
| `results/` | Summary CSVs: Table 1–3, full MR results, sample-overlap matrix |

## Reproducibility notes

- Observational analyses use R/Python with fixed random seeds; MR uses OpenGWAS/IEU summary statistics (dataset inventory in Supplement S9). Raw HRS/CHARLS individual-level data are **not** included (access-restricted) — only analysis scripts and summary outputs are deposited.
- **OpenGWAS API token:** scripts that query the OpenGWAS/IEU API read the token from the environment variable `OPENGWAS_TOKEN`, or from a local file `.opengwas_token.txt` placed next to the script. That file is gitignored and never committed. Obtain your own token at https://api.opengwas.io (free registration). Without a token, scripts that call the API will raise a clear `RuntimeError`; all offline/post-download steps run without it.
- The SMR/HEIDI implementation is a port validated against the official `smr-1.3.1` Linux binary (see `verification/`); results are reported as *suggestive/prioritized*, not as causal validation, consistent with the manuscript's graded-evidence language.

## License

Code: MIT. Summary tables and figures: CC-BY-4.0 unless otherwise noted.

## Citation

**Archived DOI:**&#8203; [10.5281/zenodo.22815831](https://doi.org/10.5281/zenodo.22815831) (Zenodo; full metadata in `CITATION.cff`).
