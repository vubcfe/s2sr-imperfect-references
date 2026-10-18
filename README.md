# Imperfect cross-sensor references for Sentinel-2 super-resolution

Code, configurations, results and derived data for the paper

> *Imperfect cross-sensor references for Sentinel-2 super-resolution: hidden benchmark errors and their
> consequences for training and evaluating deep networks* (submitted to PFG – Journal of Photogrammetry, Remote
> Sensing and Geoinformation Science).

The study analyses the SEN2NEON benchmark (Donike et al., IEEE GRSL 2026), which provides airborne (AVIRIS-NG)
references for all Sentinel-2 bands. It finds a per-tile misregistration, a contrast compression caused by the
direction of the harmonisation regression, and 131 tiles with resampled 20 m inputs; it compensates the first two
with a per-tile registration and one compression factor per tile and band, and studies their consequences for deep
super-resolution training and evaluation in simulation (domain S) and on nine real test sites (domain R).

## Contents

| Path | What it holds |
|---|---|
| `src/s2sr/` | Python package: degradation model (`degrade.py`), registration (`registration.py`), SEN2NEON reader, networks, and the study code in `p6/` (losses, metrics, training, corruption models) |
| `scripts/` | One script per step of the study (see [REPRODUCE.md](REPRODUCE.md)) |
| `configs/` | Frozen experiment configurations with their dates and dated amendments, MTF values, site split |
| `jobs/` | The job lists that were run (one command per line; `p6_queue.sh` runs them on two GPUs) |
| `tests/` | Unit tests (registration, metrics, losses, training/resume, data preparation) |
| `results/` | Reference-error analysis (E0) and every number of the paper that is not in an experiment summary (`results/p6_paper_facts/facts.md`) |
| `runs/` | Per-run metrics (CSV/JSON) and analysis summaries of E1, E1-R, E2 and the post-hoc runs; `*_superseded/` holds runs invalidated by the Reg-HP implementation fix (kept for transparency) |
| `paper/` | LaTeX source, bibliography, figures and `make_figures.py` |

Model weights and the prepared real-data domain (corrected references) are too large for git and are archived on
Zenodo (see below).

## Installation

Python ≥ 3.10 and PyTorch ≥ 2.1 (CUDA recommended for training).

```bash
pip install -e .            # or: pip install -r requirements.txt
pytest -q                  # 67 unit tests, CPU only, about three minutes
```

## Data

SEN2NEON (CC-BY-4.0) is downloaded from Hugging Face (`isp-uv-es/SEN2NEON`):

```bash
python scripts/download_sen2neon.py --out data/sen2neon --workers 8
```

The derived products of this study (tile shifts, compression factors, site split, corrected references) are released
under CC-BY-4.0 with attribution to SEN2NEON (see [LICENSE-DATA.md](LICENSE-DATA.md)):

- in this repository: per-tile shifts and peaks (`results/p6_e0/tiles.csv`, columns `g_dy, g_dx` for the 10 m bands
  and `t_dy, t_dx` for the 20 m target bands, in 10 m pixels; the reference is shifted by the negative), per-band
  statistics (`results/p6_e0/bands.csv`), compression factors (`results/p6_paper_facts/domain_r_gains.csv`), site
  split (`configs/p6_split.yaml`, `configs/p6_split.sites.csv`) and the list of excluded tiles (status column of
  `tiles.csv`);
- on Zenodo: `domain_r/{train,val,test,train_dt40}.h5` with, per tile, the native 20 m bands `x`, the 10 m bands `g`,
  the reference as released `r_raw`, registered `r_reg` and corrected `y` (= y_c), and validity masks (uint16 ×
  10⁻⁴ reflectance), plus `gains.csv` and `build_report.json`; and the final weights of every reported model
  (`weights/<experiment>/<method>/seed<k>/student.pt`).

## Reproducing the study

[REPRODUCE.md](REPRODUCE.md) lists every step in order, with the command, its inputs and outputs, and the tables and
figures it produces. The experiments were run on two RTX 4090 GPUs (about 20,000 steps, 30–60 minutes per run).

## Notes on the analysis

- Hypotheses H1–H3 and selection rules were fixed in dated configuration files before the corresponding runs
  (`configs/p6_e1.yaml`, `configs/p6_e2.yaml`, `configs/p6_e2_posthoc.yaml`); amendments are dated in the same files.
- An implementation error in the first version of the Reg-HP loss was fixed on 15 October 2026 before all reported
  Reg-HP runs; the invalid runs are kept in `runs/*_superseded/`.
- The post-hoc Reg-HP weight (`λ_W = 30`) was chosen on validation sites by a rule recorded before evaluation.

## Licences

- Code: MIT ([LICENSE](LICENSE)).
- Derived data and results: CC-BY-4.0, with attribution to SEN2NEON ([LICENSE-DATA.md](LICENSE-DATA.md)).

## Citation

See [CITATION.cff](CITATION.cff). Please also cite SEN2NEON:
Donike, S., Aybar, C., Contreras, J., Gómez-Chova, L. (2026). SEN2NEON: Enabling quantitative benchmarking of
Sentinel-2 superresolution for all multispectral bands. IEEE Geoscience and Remote Sensing Letters 23, 6013905.
https://doi.org/10.1109/LGRS.2026.3703947
