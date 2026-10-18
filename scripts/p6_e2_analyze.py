"""P6 E2 analysis (configs/p6_e2.yaml, rules fixed before the runs).

Unit of analysis: the test site (9 sites). Per run: tile/band values -> band mean per tile -> site mean; core
methods are averaged over their 3 seeds per site before the paired site bootstrap (2,000 resamples, percentile CI).
  H1: Reg-HP - Wald > 0 on HP-PSNR and HP-corr (both CIs exclude 0, reference y) and
      Wald self-consistency PSNR(Reg-HP) - PSNR(Wald) > -0.10 dB (lower CI bound).
  H3: step 1 P6 - B-sim+FT > 0 on both endpoints; step 2 (only if step 1 holds) P6 - Reg-HP, Holm over the two
      endpoints.
Secondary (descriptive, no adjustment): pairwise differences of interest, reference r_reg (contrast as delivered),
per-band HP-corr, ranking agreement between references.

  python scripts/p6_e2_analyze.py --runs runs/p6_e2
Outputs <runs>/analysis/: site_means.csv, summary.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

METRICS = ["hp_psnr", "hp_corr", "hp_gain", "edge_error"]
ORDER = ["wald", "naive_l1", "reg_l1_nogain", "reg_l1", "shift_l1", "reg_hp", "b_sim", "b_sim_ft", "p6",
         "reg_hp_w30"]  # the last one only with --extra runs/p6_e2_posthoc (exploratory)
N_BOOT, SEED, MARGIN = 2000, 2026, 0.10


def load(runs: Path) -> pd.DataFrame:
    rows = []
    for d in sorted(runs.glob("*/seed*")):
        m, s = d.parts[-2], d.parts[-1]
        if m.startswith("_") or not (d / "run.json").exists():
            continue
        for ref in ("y", "r_reg", "r_raw"):
            if not (d / f"test_{ref}.csv").exists():
                continue
            t = pd.read_csv(d / f"test_{ref}.csv")
            t = t.groupby(["site", "tile"])[METRICS].mean().groupby("site").mean().reset_index()
            t["ref"] = ref
            w = pd.read_csv(d / "wald_test.csv").groupby("site").wald_psnr.mean()
            t["wald_psnr"] = t.site.map(w)
            t["method"], t["seed"] = m, s
            rows.append(t)
    return pd.concat(rows, ignore_index=True)


def site_table(df, ref, metric):
    """sites x methods, seed-averaged."""
    return df[df.ref == ref].groupby(["site", "method"])[metric].mean().unstack()


def boot(diff: np.ndarray, rng) -> tuple[float, float, float]:
    idx = rng.integers(0, diff.size, (N_BOOT, diff.size))
    bs = diff[idx].mean(1)
    return float(diff.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def boot_p(diff: np.ndarray, rng) -> float:
    """Two-sided bootstrap p-value (centred percentile)."""
    idx = rng.integers(0, diff.size, (N_BOOT, diff.size))
    bs = diff[idx].mean(1)
    return float(min(1.0, 2 * min((bs <= 0).mean(), (bs >= 0).mean())))


def holm(ps):
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run)
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs/p6_e2")
    ap.add_argument("--extra", nargs="*", default=[], help="further run directories (e.g. runs/p6_e2_posthoc)")
    a = ap.parse_args()
    runs = Path(a.runs)
    out = runs / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.concat([load(runs)] + [load(Path(e)) for e in a.extra], ignore_index=True)
    df = df[df.method.isin(ORDER)]
    df.to_csv(out / "site_means.csv", index=False)
    methods = [m for m in ORDER if m in set(df.method)]
    def rng_for(section: str):
        """Independent stream per section, so adding or reordering sections never changes other intervals."""
        return np.random.default_rng([SEED, sum(map(ord, section))])

    rng = rng_for("pre-specified")
    md = ["# E2 analysis (test sites: %d, tiles per site averaged)" % df.site.nunique(), ""]

    # main table
    for ref in [r for r in ("y", "r_reg", "r_raw") if r in set(df.ref)]:
        md += [f"## Results, reference {ref} (mean of site means; ± sd over seeds; n seeds)", "",
               "| method | n | HP-PSNR | HP-corr | HP-gain | edge error | Wald PSNR |", "|---|---|---|---|---|---|---|"]
        sub = df[df.ref == ref]
        per_seed = sub.groupby(["method", "seed"])[METRICS + ["wald_psnr"]].mean()
        for m in methods:
            p = per_seed.loc[m]
            cell = lambda k: f"{p[k].mean():.4f}" + (f" ± {p[k].std():.4f}" if len(p) > 1 else "")  # noqa: E731
            md.append(f"| {m} | {len(p)} | " + " | ".join(cell(k) for k in METRICS + ["wald_psnr"]) + " |")
        md.append("")

    # H1
    md += ["## H1: Reg-HP vs Wald (reference y; site bootstrap, 95 % CI)", "", "| endpoint | diff | CI | holds |",
           "|---|---|---|---|"]
    h1 = []
    for k in ("hp_psnr", "hp_corr"):
        t = site_table(df, "y", k)
        mu, lo, hi = boot((t["reg_hp"] - t["wald"]).to_numpy(), rng)
        h1.append(lo > 0)
        md.append(f"| {k} | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {lo > 0} |")
    t = site_table(df, "y", "wald_psnr")
    mu, lo, hi = boot((t["reg_hp"] - t["wald"]).to_numpy(), rng)
    h1.append(lo > -MARGIN)
    md += [f"| Wald PSNR (non-inferiority, margin {MARGIN} dB) | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {lo > -MARGIN} |",
           "", f"**H1 {'holds' if all(h1) else 'does not hold'}.**", ""]

    # H3
    md += ["## H3 (hierarchical)", "", "Step 1: P6 vs B-sim+FT", "", "| endpoint | diff | CI | holds |",
           "|---|---|---|---|"]
    s1 = []
    for k in ("hp_psnr", "hp_corr"):
        t = site_table(df, "y", k)
        mu, lo, hi = boot((t["p6"] - t["b_sim_ft"]).to_numpy(), rng)
        s1.append(lo > 0)
        md.append(f"| {k} | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {lo > 0} |")
    md += ["", f"Step 1 {'holds' if all(s1) else 'does not hold'}."]
    ps, rows = [], []
    for k in ("hp_psnr", "hp_corr"):
        t = site_table(df, "y", k)
        diff = (t["p6"] - t["reg_hp"]).to_numpy()
        rows.append((k, *boot(diff, rng)))
        ps.append(boot_p(diff, rng))
    adj = holm(np.array(ps))
    md += ["", "Step 2: P6 vs Reg-HP (" + ("confirmatory" if all(s1) else "descriptive only, step 1 failed") + ")", "",
           "| endpoint | diff | CI | Holm p |", "|---|---|---|---|"]
    md += [f"| {k} | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {p:.3g} |" for (k, mu, lo, hi), p in zip(rows, adj)]
    md += ["", f"**H3 {'holds' if all(s1) and all(r[2] > 0 for r in rows) and (adj < 0.05).all() else 'does not hold'}.**",
           ""]

    # secondary comparisons
    rng = rng_for("secondary")
    pairs = [("shift_l1", "reg_l1"), ("reg_l1", "naive_l1"), ("reg_l1", "reg_l1_nogain"), ("reg_hp", "reg_l1"),
             ("reg_hp", "shift_l1"), ("p6", "reg_l1"), ("b_sim_ft", "b_sim"), ("reg_l1", "b_sim_ft"),
             ("reg_l1", "wald")]
    md += ["## Secondary paired differences (site bootstrap, 95 % CI, unadjusted)", "",
           "| comparison | ref | HP-PSNR | HP-corr | HP-gain | edge error | Wald PSNR |", "|---|---|---|---|---|---|---|"]
    for ref in [r for r in ("y", "r_reg", "r_raw") if r in set(df.ref)]:
        for x, z in pairs:
            if x not in methods or z not in methods:
                continue
            cells = []
            for k in METRICS + ["wald_psnr"]:
                t = site_table(df, ref, k)
                mu, lo, hi = boot((t[x] - t[z]).to_numpy(), rng)
                star = "*" if lo > 0 or hi < 0 else ""
                cells.append(f"{mu:+.4f} [{lo:+.4f}, {hi:+.4f}]{star}")
            md.append(f"| {x} − {z} | {ref} | " + " | ".join(cells) + " |")
    md += ["", "\\* CI excludes 0. Seed-averaged site means (three seeds per method).", ""]

    # per band HP-corr (reference y)
    rows = []
    for d in sorted(runs.glob("*/seed*")):
        if d.parts[-2].startswith("_") or not (d / "test_y.csv").exists():
            continue
        t = pd.read_csv(d / "test_y.csv")
        t = t.groupby(["band", "site"]).hp_corr.mean().groupby("band").mean().rename(d.parts[-2])
        rows.append(t)
    band = pd.concat(rows, axis=1).T.groupby(level=0).mean()
    band = band.loc[[m for m in methods if m in band.index]]
    md += ["## HP-corr per band (reference y, seed means)", "", "| method | " + " | ".join(band.columns) + " |",
           "|---|" + "---|" * band.shape[1]]
    md += [f"| {m} | " + " | ".join(f"{v:.4f}" for v in r) + " |" for m, r in band.iterrows()]

    # ranking agreement between references, with site-bootstrap CI and a permutation p-value
    rng = rng_for("ranking")
    md += ["", "## Ranking agreement between references (Kendall tau over methods; seed-averaged site means)", "",
           "CI: 2,000 site-bootstrap resamples; p: one-sided permutation test of tau > 0 (5,000 permutations of "
           "the method labels of one ranking).", "", "| refs | endpoint | tau | CI | p |", "|---|---|---|---|---|"]
    for r1, r2 in (("y", "r_reg"), ("y", "r_raw")):
        if r2 not in set(df.ref):
            continue
        conf = [m for m in methods if m != "reg_hp_w30"]  # rankings over the pre-specified strategies only
        for k in METRICS:
            t1, t2 = site_table(df, r1, k)[conf], site_table(df, r2, k)[conf]
            tau = kendalltau(t1.mean(), t2.mean())[0]
            bs = []
            for _ in range(N_BOOT):
                idx = rng.integers(0, len(t1), len(t1))
                bs.append(kendalltau(t1.iloc[idx].mean(), t2.iloc[idx].mean())[0])
            perm = np.array([kendalltau(t1.mean().to_numpy(), rng.permutation(t2.mean().to_numpy()))[0]
                             for _ in range(5000)])
            lo, hi = np.nanquantile(bs, [0.025, 0.975])
            md.append(f"| {r1} vs {r2} | {k} | {tau:.3f} | [{lo:.3f}, {hi:.3f}] | {(perm >= tau).mean():.4f} |")

    # hierarchical bootstrap (sites, then seeds within method) for the key comparisons
    rng = rng_for("hierarchical")
    pairs_h = [("b_sim_ft", "shift_l1"), ("reg_l1", "naive_l1"), ("reg_l1", "reg_l1_nogain"), ("reg_hp", "reg_l1"), ("shift_l1", "reg_l1"),
               ("p6", "reg_l1"), ("p6", "b_sim_ft"), ("reg_l1", "b_sim_ft"), ("reg_l1", "b_sim"), ("reg_hp", "wald"),
               ("p6", "reg_hp"), ("reg_hp_w30", "reg_l1"), ("reg_hp_w30", "reg_hp")]
    md += ["", "## Hierarchical bootstrap (sites, then seeds within each method), reference y", "",
           "| comparison | endpoint | diff | CI |", "|---|---|---|---|"]
    sub = df[df.ref == "y"]
    cube = {k: sub.pivot_table(index="site", columns=["method", "seed"], values=k)
            for k in ("hp_corr", "hp_psnr", "edge_error", "hp_gain", "wald_psnr")}
    sites = cube["hp_corr"].index.to_numpy()
    for x_, z_ in pairs_h:
        if x_ not in methods or z_ not in methods:
            continue
        for k, cb in cube.items():
            sx, sz = cb[x_], cb[z_]
            bs = []
            for _ in range(N_BOOT):
                si = rng.integers(0, len(sites), len(sites))
                cx = rng.integers(0, sx.shape[1], sx.shape[1])
                cz = rng.integers(0, sz.shape[1], sz.shape[1])
                bs.append(sx.iloc[si, cx].to_numpy().mean() - sz.iloc[si, cz].to_numpy().mean())
            mu = sx.to_numpy().mean() - sz.to_numpy().mean()
            lo, hi = np.quantile(bs, [0.025, 0.975])
            md.append(f"| {x_} − {z_} | {k} | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}]{'*' if lo > 0 or hi < 0 else ''} |")
    if "r_raw" in set(df.ref):
        rng = rng_for("hierarchical r_raw")
        md += ["", "## Hierarchical bootstrap, reference as released (r_raw), HP-corr", "",
               "| comparison | diff | CI |", "|---|---|---|"]
        cr = df[df.ref == "r_raw"].pivot_table(index="site", columns=["method", "seed"], values="hp_corr")
        for x_, z_ in pairs_h:
            if x_ not in cr.columns.get_level_values(0) or z_ not in cr.columns.get_level_values(0):
                continue
            sx, sz = cr[x_], cr[z_]
            bs = []
            for _ in range(N_BOOT):
                si = rng.integers(0, len(sx), len(sx))
                bs.append(sx.iloc[si, rng.integers(0, sx.shape[1], sx.shape[1])].to_numpy().mean()
                          - sz.iloc[si, rng.integers(0, sz.shape[1], sz.shape[1])].to_numpy().mean())
            lo, hi = np.quantile(bs, [0.025, 0.975])
            md.append(f"| {x_} − {z_} | {sx.to_numpy().mean() - sz.to_numpy().mean():+.4f} | [{lo:+.4f}, {hi:+.4f}] |")
    (out / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
