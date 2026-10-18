"""Uncertainty of the E1-R ranking reversal (revision of the P6 paper).

Kendall tau between the rankings of models by an endpoint scored against the true target y and against the
harmonised y_h (and the inverted y_inv), from runs/p6_e1r/analysis/scores_by_reference.csv. Seeds are averaged per
model before ranking. CI: 2,000 bootstrap resamples of tiles (the two evaluation sites are too few to resample);
p: one-sided permutation test of tau > 0. Two model sets: all models, and one model per method (E1-R methods plus
the E1 S0 methods, i.e. without the S2 duplicates of the same methods).

  python scripts/p6_e1r_tau_ci.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = Path(__file__).resolve().parents[1]
N_BOOT, N_PERM, SEED = 2000, 5000, 2026


def main():
    s = pd.read_csv(ROOT / "runs/p6_e1r/analysis/scores_by_reference.csv")
    rng = np.random.default_rng(SEED)
    md = ["# E1-R ranking reversal: Kendall tau with uncertainty", "",
          "| model set | n models | endpoint | refs | tau | CI (tile bootstrap) | p (permutation) |",
          "|---|---|---|---|---|---|---|"]
    sets = {"all": sorted(s.model.unique()),
            "one per method": sorted(m for m in s.model.unique() if not m.startswith("S2_realistic/"))}
    for set_name, models in sets.items():
        sub = s[s.model.isin(models)]
        for k in ("hp_psnr", "hp_corr"):
            # tile x model tables (band mean per tile, then seed mean per model)
            t = {r: sub[sub.ref == r].groupby(["tile", "model", "seed"])[k].mean().groupby(["tile", "model"]).mean()
                 .unstack("model")[models] for r in ("y", "y_h", "y_inv")}
            tiles = t["y"].index
            for r2 in ("y_h", "y_inv"):
                tau = kendalltau(t["y"].mean(), t[r2].mean())[0]
                bs = []
                for _ in range(N_BOOT):
                    idx = rng.integers(0, len(tiles), len(tiles))
                    bs.append(kendalltau(t["y"].iloc[idx].mean(), t[r2].iloc[idx].mean())[0])
                perm = np.array([kendalltau(t["y"].mean().to_numpy(), rng.permutation(t[r2].mean().to_numpy()))[0]
                                 for _ in range(N_PERM)])
                lo, hi = np.nanquantile(bs, [0.025, 0.975])
                md.append(f"| {set_name} | {len(models)} | {k} | y vs {r2} | {tau:.3f} | [{lo:.3f}, {hi:.3f}] | "
                          f"{(perm >= tau).mean():.4f} |")
    # ranks by HP-PSNR (tile means, seeds averaged) under y and y_h, all models (as in Fig. 4 of the paper)
    t = {r: s[s.ref == r].groupby(["tile", "model", "seed"]).hp_psnr.mean().groupby(["tile", "model"]).mean()
         .unstack("model").mean() for r in ("y", "y_h")}
    ry, rh = t["y"].rank(ascending=False).astype(int), t["y_h"].rank(ascending=False).astype(int)
    md += ["", "Ranks by HP-PSNR (tile means, seeds averaged; all models):", "", "| model | rank vs y | rank vs y_h |",
           "|---|---|---|"]
    md += [f"| {m} | {ry[m]} | {rh[m]} |" for m in ry.sort_values().index]
    out = ROOT / "runs/p6_e1r/analysis/tau_ci.md"
    out.write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
