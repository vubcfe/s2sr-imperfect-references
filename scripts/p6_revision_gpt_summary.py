"""Summary of the analyses added after the GPT-5.5 review: full-resolution consistency of the E2 models, the
bicubic baseline and the positive control of the correction on the 10 m bands.

  python scripts/p6_revision_gpt_summary.py  -> runs/p6_e2/analysis/gpt_revision_summary.md
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "runs/p6_e2/analysis"


def main():
    md = ["# Analyses added after the GPT-5.5 review (18/10/2026)", ""]
    d = pd.read_csv(A / "fullres_consistency.csv")
    s = d.groupby(["method", "seed", "site"])[["psnr", "rmse", "corr", "contrast"]].mean().groupby(["method", "seed"]).mean()
    g = s.groupby("method").agg(["mean", "std"])
    md += ["## Full-resolution consistency: D(y_hat) vs observed 20 m bands x (test sites; mean of site means; "
           "mean and sd over seeds)", "", "| method | PSNR (dB) | sd | contrast sd(D(y_hat))/sd(x) | corr |",
           "|---|---|---|---|---|"]
    for m, r in g.sort_values(("psnr", "mean"), ascending=False).iterrows():
        md.append(f"| {m} | {r[('psnr', 'mean')]:.2f} | {r[('psnr', 'std')]:.2f} | {r[('contrast', 'mean')]:.3f} | "
                  f"{r[('corr', 'mean')]:.4f} |")
    t = d.groupby(["method", "site"]).psnr.mean().unstack(0)
    md += ["", "Site-level paired differences of full-resolution PSNR (dB):", ""]
    for a, b in [("reg_hp", "reg_l1"), ("reg_l1", "naive_l1"), ("reg_l1", "reg_l1_nogain"), ("reg_hp_w30", "reg_l1")]:
        x = t[a] - t[b]
        md.append(f"- {a} − {b}: mean {x.mean():+.2f}, sites + / −: {(x > 0).sum()} / {(x < 0).sum()}")
    md += ["", "## Bicubic baseline (y_hat = bicubic upsampling of x), high-pass metrics", ""]
    for ref in ("y", "r_reg", "r_raw"):
        b = pd.read_csv(A / f"bicubic_{ref}.csv").groupby(["site", "tile"])[["hp_psnr", "hp_corr", "edge_error",
                                                                              "hp_gain"]].mean().groupby("site").mean().mean()
        md.append(f"- against {ref}: " + ", ".join(f"{k} {v:.4f}" for k, v in b.items()))
    c = pd.read_csv(A / "guide_control.csv")
    md += ["", f"## Positive control on the 10 m bands ({c.tile.nunique()} tiles, {c.site.nunique()} sites): "
           "reference versions vs observed Sentinel-2 10 m (medians over tile-bands)", "",
           "| reference | HP-gain | HP-corr | HP-PSNR | edge error |", "|---|---|---|---|---|"]
    for ref, r in c.groupby("ref")[["hp_gain", "hp_corr", "hp_psnr", "edge_error"]].median().iterrows():
        md.append(f"| {ref} | {r.hp_gain:.3f} | {r.hp_corr:.3f} | {r.hp_psnr:.2f} | {r.edge_error:.3f} |")
    md += ["", "Per band (median HP-gain r -> y_c; a_hat): " + ", ".join(
        f"{b} {c[(c.band == b) & (c.ref == 'r')].hp_gain.median():.2f}->{c[(c.band == b) & (c.ref == 'y_c')].hp_gain.median():.2f}"
        f" (a_hat {c[(c.band == b) & (c.ref == 'y_c')].a_hat.median():.2f})" for b in sorted(c.band.unique()))]
    site = c.groupby(["ref", "site"])[["hp_gain", "hp_psnr"]].median().unstack(0)
    md += ["", f"Sites where y_c has higher median HP-PSNR than r_reg: "
           f"{(site[('hp_psnr', 'y_c')] > site[('hp_psnr', 'r_reg')]).sum()} of {len(site)}"]
    (A / "gpt_revision_summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
