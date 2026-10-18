"""Site-level paired evidence for the key E2 comparisons (revision after the GPT review): the nine seed-averaged
site differences, sign counts, a paired t-interval and an exact sign-flip permutation p-value (2^9 = 512 flips).

  python scripts/p6_e2_sitepairs.py   -> runs/p6_e2/analysis/site_pairs.md
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PAIRS = [("reg_l1", "naive_l1"), ("reg_l1_nogain", "naive_l1"), ("reg_l1", "reg_l1_nogain"), ("reg_hp", "reg_l1"),
         ("shift_l1", "reg_l1"), ("p6", "reg_l1"), ("p6", "b_sim_ft"), ("p6", "reg_hp"), ("b_sim_ft", "shift_l1"),
         ("reg_hp_w30", "reg_l1")]


def main():
    df = pd.read_csv(ROOT / "runs/p6_e2/analysis/site_means.csv")
    md = ["# E2 site-level paired differences (seed-averaged site means; 9 test sites)", "",
          "t: paired t-interval (95 %); p: exact two-sided sign-flip permutation test over the 9 sites.", ""]
    for ref, k in (("y", "hp_corr"), ("r_raw", "hp_corr"), ("y", "hp_psnr"), ("y", "edge_error")):
        t = df[df.ref == ref].groupby(["site", "method"])[k].mean().unstack()
        md += [f"## {k} against {ref}", "", "| comparison | " + " | ".join(t.index) + " | + / − | mean | t-interval | p |",
               "|---|" + "---|" * len(t.index) + "---|---|---|---|"]
        for a, b in PAIRS:
            if a not in t or b not in t:
                continue
            d = (t[a] - t[b]).to_numpy()
            n = len(d)
            lo, hi = stats.t.interval(0.95, n - 1, loc=d.mean(), scale=d.std(ddof=1) / np.sqrt(n))
            flips = np.array([np.mean(d * np.array(s)) for s in itertools.product([1, -1], repeat=n)])
            p = float(np.mean(np.abs(flips) >= abs(d.mean()) - 1e-15))
            md.append(f"| {a} − {b} | " + " | ".join(f"{v:+.4f}" for v in d) +
                      f" | {(d > 0).sum()} / {(d < 0).sum()} | {d.mean():+.4f} | [{lo:+.4f}, {hi:+.4f}] | {p:.3f} |")
        md.append("")
    out = ROOT / "runs/p6_e2/analysis/site_pairs.md"
    out.write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
