"""Full-resolution consistency of every E2 model and a bicubic baseline (revision after the GPT review).

For each run (seed) and test tile: the prediction y_hat at 10 m, degraded with the MTF-matched D, is compared with
the observed 20 m bands x on valid 20 m pixels: PSNR (peak 1), RMSE, correlation and the contrast ratio
sd(D(y_hat))/sd(x). This is the low-frequency consistency the LFC is meant to enforce; it differs from the
reduced-scale Wald fidelity PSNR(f(D(x), D(g)), x). The bicubic baseline (y_hat = U(x)) is also scored with the
paper's high-pass metrics against y_c, r_reg and r (the 10 m aggregate of the released reference).

  python scripts/p6_e2_fullres.py --gpu 0
Writes runs/p6_e2/analysis/fullres_consistency.csv and bicubic_<ref>.csv (tile x band rows).
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.degrade import band_sigmas, degrade, upsample  # noqa: E402
from s2sr.p6 import metrics as M  # noqa: E402
from s2sr.p6.models import build_student  # noqa: E402

Q = 1e-4


def consistency(pred10: torch.Tensor, x20: np.ndarray, v20: np.ndarray, sig) -> list[dict]:
    d = degrade(pred10[None], sig, 2)[0].double().cpu().numpy()
    out = []
    for b, band in enumerate(TARGET_BANDS):
        a, t = d[b][v20], x20[b][v20]
        mse = float(((a - t) ** 2).mean())
        out.append(dict(band=band, psnr=10 * np.log10(1 / max(mse, 1e-20)), rmse=mse ** 0.5,
                        corr=float(np.corrcoef(a, t)[0, 1]), contrast=float(a.std() / max(t.std(), 1e-12))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", type=int, default=0)
    a = ap.parse_args()
    dev = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    sig = band_sigmas(mtf, TARGET_BANDS, 2)
    test_h5 = ROOT / "data/p6_domain_r/test.h5"
    with h5py.File(test_h5, "r") as f:
        tiles = sorted(f)
        data = {t: {k: f[t][k][:] for k in ("x", "g", "y", "r_reg", "r_raw", "valid10", "valid20")} for t in tiles}
    runs = [d for rd in ("runs/p6_e2", "runs/p6_e2_posthoc") for d in sorted((ROOT / rd).glob("*/seed*"))
            if not d.parts[-2].startswith("_") and (d / "student.pt").exists() and (d / "test_y.csv").exists()]
    rows, bic = [], {k: [] for k in ("y", "r_reg", "r_raw")}
    for t in tiles:  # bicubic baseline
        x = torch.from_numpy(data[t]["x"] * Q).float()[None].to(dev)
        up = upsample(x, 2)[0]
        v10, v20 = data[t]["valid10"] > 0, data[t]["valid20"] > 0
        for r in consistency(up, data[t]["x"] * Q, v20, sig):
            rows.append(dict(method="bicubic", seed="-", tile=t, site=t.split("_")[1], **r))
        p = up.double().cpu().numpy()
        for ref in bic:
            yr = data[t][ref] * Q
            hc, hps, ee, hg = M.hp_corr(p, yr, v10), M.hp_psnr(p, yr, v10), M.edge_error(p, yr, v10), M.hp_gain(p, yr, v10)
            bic[ref] += [dict(tile=t, site=t.split("_")[1], band=b, hp_corr=hc[k], hp_psnr=hps[k], edge_error=ee[k],
                              hp_gain=hg[k]) for k, b in enumerate(TARGET_BANDS)]
    for d in runs:
        model = build_student(mtf).to(dev).eval()
        model.load_state_dict(torch.load(d / "student.pt", map_location=dev))
        with torch.no_grad():
            for t in tiles:
                x = torch.from_numpy(data[t]["x"] * Q).float()[None].to(dev)
                g = torch.from_numpy(data[t]["g"] * Q).float()[None].to(dev)
                pred = model(x, g)[0]
                for r in consistency(pred, data[t]["x"] * Q, data[t]["valid20"] > 0, sig):
                    rows.append(dict(method=d.parts[-2], seed=d.parts[-1], tile=t, site=t.split("_")[1], **r))
        print("done", d.relative_to(ROOT), flush=True)
    out = ROOT / "runs/p6_e2/analysis"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "fullres_consistency.csv", index=False)
    for ref, rr in bic.items():
        pd.DataFrame(rr).to_csv(out / f"bicubic_{ref}.csv", index=False)
    s = (df.groupby(["method", "seed", "site"])[["psnr", "rmse", "corr", "contrast"]].mean()
         .groupby(["method", "seed"]).mean().groupby("method").mean())
    print(s.round(4).to_string())
    for ref, rr in bic.items():
        b = pd.DataFrame(rr).groupby(["site", "tile"])[["hp_corr", "hp_psnr", "edge_error", "hp_gain"]].mean() \
            .groupby("site").mean().mean()
        print("bicubic vs", ref, b.round(4).to_dict())


if __name__ == "__main__":
    main()
