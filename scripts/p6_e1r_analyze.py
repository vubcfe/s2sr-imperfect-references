"""P6 E1-R analysis (configs/p6_e1r.yaml).

Part A (training consequence): the E1-R students scored against the true y on the E1 `evaluate` sites
  (eval.csv written by p6_e1.py): per-seed means of site means, seed mean/sd, per-band HP-gain, tile-bootstrap CIs of
  the paired differences to reg_l1_clean.
Part B (evaluation consequence): the same predictions scored against three references of the same tiles, the true y,
  the harmonised y_h (as SEN2NEON is used today) and the inverted y_inv. Models: the E1-R students and the E1 S0/S2
  students. The ranking of the methods under y_h vs y (Kendall tau) shows whether the harmonised reference changes
  conclusions.

  python scripts/p6_e1r_analyze.py --gpu 1
Outputs runs/p6_e1r/analysis/: scores_by_reference.csv, summary.md
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import yaml
from scipy.stats import kendalltau

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.p6 import metrics as M  # noqa: E402
from s2sr.p6.models import build_student  # noqa: E402

Q = 1e-4
METRICS = ["hp_psnr", "hp_corr", "hp_gain", "edge_error"]


def site_seed_table(runs: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(runs.glob("*/seed*/eval.csv")):
        d = pd.read_csv(f)
        d["method"], d["seed"] = f.parts[-3], f.parts[-2]
        rows.append(d)
    return pd.concat(rows)


def boot_diff(d: pd.DataFrame, a: str, b: str, metric: str, n=2000, seed=2026):
    """Tile-bootstrap CI of mean(a - b) (seed-averaged per tile and band, then band-averaged per tile)."""
    t = d.groupby(["method", "tile", "band"])[metric].mean().groupby(["method", "tile"]).mean().unstack(0)
    diff = (t[a] - t[b]).dropna().to_numpy()
    rng = np.random.default_rng(seed)
    bs = rng.choice(diff, (n, diff.size)).mean(1)
    return diff.mean(), np.quantile(bs, 0.025), np.quantile(bs, 0.975)


@torch.no_grad()
def score_refs(ckpt, mtf, s_val, t_val, tiles, device):
    model = build_student(mtf).to(device).eval()
    model.load_state_dict(torch.load(ckpt, map_location=device))
    rows = []
    with h5py.File(s_val, "r") as fs, h5py.File(t_val, "r") as ft:
        for tile in tiles:
            gs, gt = fs[tile], ft[tile]
            x = torch.from_numpy(gs["x"][:] * Q).float()[None].to(device)
            g = torch.from_numpy(gs["g"][:] * Q).float()[None].to(device)
            pred = model(x, g)[0].cpu().double().numpy()
            v = gs["valid10"][:] > 0
            refs = {"y": gs["y"][:] * Q, "y_h": gt["y_h"][:] * Q, "y_inv": gt["y_inv"][:] * Q}
            for rk, r in refs.items():
                vals = dict(hp_psnr=M.hp_psnr(pred, r, v), hp_corr=M.hp_corr(pred, r, v),
                            hp_gain=M.hp_gain(pred, r, v), edge_error=M.edge_error(pred, r, v))
                for k, band in enumerate(TARGET_BANDS):
                    rows.append(dict(tile=tile, site=tile.split("_")[1], band=band, ref=rk,
                                     **{m: vals[m][k] for m in METRICS}))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/p6_e1r.yaml")
    ap.add_argument("--e1_runs", default="runs/p6_e1")
    ap.add_argument("--e1_scenarios", nargs="+", default=["S0_clean", "S2_realistic"])
    ap.add_argument("--gpu", type=int, default=0)
    a = ap.parse_args()
    cfg = yaml.safe_load((ROOT / a.config).read_text())
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    device = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    runs = ROOT / cfg["runs_dir"] / cfg["scenarios"][0]["name"]
    out = ROOT / cfg["runs_dir"] / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    md = ["# E1-R analysis", ""]

    # Part A
    d = site_seed_table(runs)
    seed = d.groupby(["method", "seed", "site"])[METRICS].mean().groupby(["method", "seed"]).mean()
    agg = seed.groupby("method").agg(["mean", "std"])
    md += ["## A. Training on the harmonised reference (scored against the true y; mean of site means; 3 seeds)", "",
           "| method | HP-PSNR | HP-corr | HP-gain | edge error |", "|---|---|---|---|---|"]
    for m, r in agg.iterrows():
        md.append(f"| {m} | " + " | ".join(f"{r[(k, 'mean')]:.4f} ± {r[(k, 'std')]:.4f}" for k in METRICS) + " |")
    band = d.groupby(["method", "band"]).hp_gain.mean().unstack()[list(TARGET_BANDS)]
    md += ["", "HP-gain per band:", "", "| method | " + " | ".join(band.columns) + " |",
           "|---|" + "---|" * band.shape[1]]
    md += [f"| {m} | " + " | ".join(f"{v:.3f}" for v in r) + " |" for m, r in band.iterrows()]
    md += ["", "Paired differences to reg_l1_clean (tile bootstrap, 95 % CI):", "",
           "| method | metric | diff | CI |", "|---|---|---|---|"]
    for m in sorted(set(d.method) - {"reg_l1_clean"}):
        for k in METRICS:
            mu, lo, hi = boot_diff(d, m, "reg_l1_clean", k)
            md.append(f"| {m} | {k} | {mu:+.4f} | [{lo:+.4f}, {hi:+.4f}] |")

    # Part B
    s_val = ROOT / cfg["domain_s"] / "val.h5"
    t_val = ROOT / cfg["targets_dir"] / "targets_val.h5"
    with h5py.File(s_val, "r") as f:
        tiles = [t for t in sorted(f) if t.split("_")[1] in cfg["val_sites"]["evaluate"]]
    models = [(f"E1R/{p.parts[-3]}", p.parts[-2], p) for p in sorted(runs.glob("*/seed*/student.pt"))]
    for sc in a.e1_scenarios:
        models += [(f"{sc}/{p.parts[-3]}", p.parts[-2], p)
                   for p in sorted((ROOT / a.e1_runs / sc).glob("*/seed*/student.pt")) if p.parts[-3] != "_teacher"]
    cache = out / "scores_by_reference.csv"
    done = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=["model", "seed"])
    have = set(zip(done.model, done.seed))
    parts = [done]
    for name, sd, ck in models:
        if (name, sd) in have:
            continue
        s = score_refs(ck, mtf, s_val, t_val, tiles, device)
        s["model"], s["seed"] = name, sd
        parts.append(s)
        pd.concat(parts).to_csv(cache, index=False)
        print("scored", name, sd, flush=True)
    sc = pd.concat(parts)
    tab = (sc.groupby(["model", "seed", "ref", "site"])[METRICS].mean().groupby(["model", "seed", "ref"]).mean()
           .groupby(["model", "ref"]).mean())
    md += ["", "## B. Scoring the same predictions against y, y_h and y_inv (seed means)", ""]
    for k in ["hp_psnr", "hp_corr", "hp_gain"]:
        t = tab[k].unstack("ref")[["y", "y_inv", "y_h"]]
        md += [f"### {k}", "", "| model | y | y_inv | y_h | rank y | rank y_h |", "|---|---|---|---|---|---|"]
        asc = k == "edge_error"
        ry, rh = t["y"].rank(ascending=asc), t["y_h"].rank(ascending=asc)
        for m, r in t.sort_values("y", ascending=asc).iterrows():
            md.append(f"| {m} | {r.y:.4f} | {r.y_inv:.4f} | {r.y_h:.4f} | {ry[m]:.0f} | {rh[m]:.0f} |")
        tau, p = kendalltau(t["y"], t["y_h"])
        tau_i, _ = kendalltau(t["y"], t["y_inv"])
        md += ["", f"Kendall tau (y vs y_h) = {tau:.3f} (p = {p:.3g}); (y vs y_inv) = {tau_i:.3f}", ""]
    (out / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
