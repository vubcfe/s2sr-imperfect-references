# Reproducing the study

All commands run from the repository root with the package installed (`pip install -e .`). Steps that train
networks need a GPU. Each step writes its outputs next to the paths given; the analysis summaries shipped in
`results/` and `runs/` are the outputs used in the paper.

## 0. Data

| Step | Command | Output |
|---|---|---|
| Download SEN2NEON | `python scripts/download_sen2neon.py --out data/sen2neon --workers 8` | `data/sen2neon/` (~30 GB) |

## 1. Reference error analysis (E0) — Sections 4 and 6, Figs. 2–3, Tables 1 and 3

| Step | Command | Output |
|---|---|---|
| Error analysis of all tiles | `python scripts/p6_e0_reference_errors.py --root data/sen2neon --out results/p6_e0` | `results/p6_e0/{tiles.csv,bands.csv,summary.md,harmonisation_regression_check.csv}` |
| MTF sensitivity | `python scripts/p6_e0_reference_errors.py --mtf-delta 0.05 --every 5 --out results/p6_e0_mtf0.05` (and `-0.05`) | `results/p6_e0_mtf±0.05/` |
| Site split | `python scripts/p6_make_split.py` | `configs/p6_split.yaml`, `configs/p6_split.sites.csv` |

## 2. Domains

| Step | Command | Output |
|---|---|---|
| Domain S (simulated) | `python scripts/p6_build_domain_s.py` (config `configs/p6_domain_s.yaml`) | `data/p6_domain_s/` |
| Domain R (real, corrected references) | `python scripts/p6_build_domain_r.py --sen2neon data/sen2neon --out data/p6_domain_r` | `data/p6_domain_r/` (also on Zenodo) |
| E1 corrupted targets | `python scripts/p6_make_e1_targets.py` | `data/p6_e1/` |
| E1-R harmonised targets | `python scripts/p6_make_e1r_targets.py` | `data/p6_e1r/` |

## 3. Experiments

Each job list holds one command per line; `bash scripts/p6_queue.sh <jobs file> <state dir> "0 1"` runs it on GPUs
0 and 1 and can be restarted after an interruption (finished runs are skipped, interrupted ones resume).

| Experiment | Configuration | Job lists | Analysis | Paper |
|---|---|---|---|---|
| E1, misregistration (domain S) | `configs/p6_e1.yaml` | `jobs/p6_e1_*.txt` | `python scripts/p6_e1_analyze.py --runs runs/p6_e1 --split eval` → `runs/p6_e1/summary_eval.md` | Table 4 |
| E1-R, harmonisation (domain S) | `configs/p6_e1r.yaml` | `jobs/p6_e1r.txt`, `jobs/p6_e1r_fix.txt` | `python scripts/p6_e1r_analyze.py`; `python scripts/p6_e1r_tau_ci.py` → `runs/p6_e1r/analysis/` | Table 5, Fig. 4 |
| E2, real data | `configs/p6_e2.yaml` | `jobs/p6_e2.txt`, `jobs/p6_e2_seeds12.txt` | `python scripts/p6_e2_analyze.py --runs runs/p6_e2 --extra runs/p6_e2_posthoc` → `runs/p6_e2/analysis/summary.md` | Table 6, Fig. 5 |
| E2 post hoc (Reg-HP Wald weight) | `configs/p6_e2_posthoc.yaml` | `jobs/p6_e2_posthoc*.txt` | included above | Section 8 |

## 4. Additional analyses

| Analysis | Command | Output | Paper |
|---|---|---|---|
| Scoring against the reference as released | `python scripts/p6_e2_eval_raw.py` | `runs/p6_e2/*/seed*/test_r_raw.csv` | Table 6 |
| Gain at several scales | `python scripts/p6_gain_scales.py` | `runs/p6_e2/analysis/gain_scales.csv` | Section 6 |
| Full-resolution consistency and bicubic baseline | `python scripts/p6_e2_fullres.py` | `runs/p6_e2/analysis/{fullres_consistency.csv,bicubic_*.csv}` | Table 6 |
| 10 m positive control | `python scripts/p6_guide_control.py` | `runs/p6_e2/analysis/guide_control.csv` | Section 6 |
| Site-level paired tests | `python scripts/p6_e2_sitepairs.py`; `python scripts/p6_revision_gpt2.py` | `runs/p6_e2/analysis/{site_pairs.md,gpt2_revision.md}` | Section 8 |
| Summaries of the above | `python scripts/p6_revision_gpt_summary.py` | `runs/p6_e2/analysis/gpt_revision_summary.md` | — |
| Qualitative crop | `python scripts/p6_qualitative.py` | `runs/p6_e2/analysis/qualitative.{npz,json}` | Fig. 6 |
| All remaining numbers | `python scripts/p6_paper_facts.py` | `results/p6_paper_facts/facts.md` | throughout |

## 5. Figures and paper

```bash
python paper/make_figures.py       # reads results/ and runs/; writes paper/figures/*.pdf
cd paper && latexmk -pdf main.tex  # needs the Springer Nature LaTeX template (sn-jnl.cls, sn-basic.bst)
```
