"""Analyses requested in the second GPT-5.5 review round (18/10/2026): exact paired sign-flip tests over the nine
test sites for H1/H3, paired site intervals for the full-resolution consistency, and site-level intervals for the
10 m positive control.

  python scripts/p6_revision_gpt2.py  -> runs/p6_e2/analysis/gpt2_revision.md
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "runs/p6_e2/analysis"
FLIPS = np.array(list(itertools.product([1, -1], repeat=9)))


def paired(d):
    d = np.asarray(d, float)
    lo, hi = stats.t.interval(0.95, len(d) - 1, loc=d.mean(), scale=d.std(ddof=1) / np.sqrt(len(d)))
    p = float(np.mean(np.abs((FLIPS[:, :len(d)] * d).mean(1)) >= abs(d.mean()) - 1e-15))
    return d.mean(), lo, hi, (d > 0).sum(), (d < 0).sum(), p


def holm(ps):
    o = np.argsort(ps); out = np.empty(len(ps)); run = 0
    for k, i in enumerate(o):
        run = max(run, (len(ps) - k) * ps[i]); out[i] = min(1, run)
    return out


def main():
    md = ["# Second GPT-review round: exact site-level tests (nine test sites, seed-averaged site means)", "",
          "Columns: mean difference, paired t-interval (95 %), sites + / −, exact two-sided sign-flip p (512 flips).", ""]
    sm = pd.read_csv(A / "site_means.csv")
    tab = lambda k: sm[sm.ref == "y"].groupby(["site", "method"])[k].mean().unstack()  # noqa: E731
    md += ["## H1 and H3 (reference y_c)", "", "| test | endpoint | diff | t-interval | + / − | p | Holm p |",
           "|---|---|---|---|---|---|---|"]
    for name, a, b, ks in (("H1", "reg_hp", "wald", ("hp_psnr", "hp_corr", "wald_psnr")),
                           ("H3 step 1", "p6", "b_sim_ft", ("hp_psnr", "hp_corr")),
                           ("H3 step 2", "p6", "reg_hp", ("hp_psnr", "hp_corr"))):
        res = [paired(tab(k)[a] - tab(k)[b]) for k in ks]
        hp = holm([r[5] for r in res[:2]])
        for i, (k, r) in enumerate(zip(ks, res)):
            md.append(f"| {name} | {k} | {r[0]:+.4f} | [{r[1]:+.4f}, {r[2]:+.4f}] | {r[3]} / {r[4]} | {r[5]:.3f} | "
                      f"{hp[i]:.3f} |" if i < 2 else
                      f"| {name} | {k} | {r[0]:+.4f} | [{r[1]:+.4f}, {r[2]:+.4f}] | {r[3]} / {r[4]} | {r[5]:.3f} | -- |")
    fr = pd.read_csv(A / "fullres_consistency.csv").groupby(["method", "site"]).psnr.mean().unstack(0)
    md += ["", "## Full-resolution consistency PSNR(D(y_hat), x), paired site differences (dB)", "",
           "| comparison | diff | t-interval | + / − | p |", "|---|---|---|---|---|"]
    for a, b in (("reg_hp", "reg_l1"), ("reg_hp", "p6"), ("reg_hp", "shift_l1"), ("reg_hp_w30", "reg_l1"),
                 ("reg_l1", "naive_l1"), ("reg_l1", "reg_l1_nogain")):
        r = paired(fr[a] - fr[b])
        md.append(f"| {a} − {b} | {r[0]:+.2f} | [{r[1]:+.2f}, {r[2]:+.2f}] | {r[3]} / {r[4]} | {r[5]:.3f} |")
    c = pd.read_csv(A / "guide_control.csv")
    s = c.groupby(["ref", "site"])[["hp_gain", "hp_corr", "hp_psnr", "edge_error"]].median()
    rng = np.random.default_rng(2026)
    sites = s.loc["r"].index.to_numpy()
    md += ["", f"## 10 m positive control: medians of site medians ({len(sites)} sites) with site-bootstrap 95 % CI", "",
           "| reference | HP-gain | HP-corr | HP-PSNR | edge error |", "|---|---|---|---|---|"]
    for ref in ("r", "r_reg", "y_c"):
        cells = []
        for k in ("hp_gain", "hp_corr", "hp_psnr", "edge_error"):
            v = s.loc[ref][k].to_numpy()
            bs = [np.median(rng.choice(v, len(v))) for _ in range(2000)]
            cells.append(f"{np.median(v):.3f} [{np.quantile(bs, .025):.3f}, {np.quantile(bs, .975):.3f}]")
        md.append(f"| {ref} | " + " | ".join(cells) + " |")
    d = (s.loc["y_c"].hp_gain - s.loc["r_reg"].hp_gain)
    md += ["", f"Compensation effect on HP-gain per site: higher at {(d > 0).sum()} of {len(d)} sites."]
    # frequency separation of the gain estimation (low-pass sigma 60 m) and the evaluation (high-pass sigma 20 m)
    f50 = np.sqrt(np.log(2) / (2 * np.pi ** 2 * 20.0 ** 2))
    lp = np.exp(-2 * np.pi ** 2 * 60.0 ** 2 * f50 ** 2)
    md += ["", f"Frequency separation: the 20 m high-pass passes half of the signal at {f50:.4f} cycles/m, where the "
           f"60 m low-pass used for a_hat transmits {lp:.4f}."]
    (A / "gpt2_revision.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
