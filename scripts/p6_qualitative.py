"""Qualitative crops for the P6 paper (E2 test sites).

Tile choice, fixed before looking at any image: among the test tiles whose corrected reference has a high-pass
standard deviation (mean over the six bands) in the top quartile, take the tile whose Reg-L1 (seed 0) HP-corr is
closest to the median of that subset. Crop: the 64 x 64 pixel (640 m) window with the largest reference high-pass
energy, restricted to valid pixels.

  python scripts/p6_qualitative.py --gpu 0
Writes runs/p6_e2/analysis/qualitative.npz (crops of x upsampled, r_raw, y, and the predictions) and the choice.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import yaml
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.p6.metrics import highpass, stat_mask  # noqa: E402
from s2sr.p6.models import build_student  # noqa: E402

Q = 1e-4
CROP = 64
MODELS = {"wald": "runs/p6_e2/wald/seed0", "naive_l1": "runs/p6_e2/naive_l1/seed0",
          "reg_l1": "runs/p6_e2/reg_l1/seed0", "reg_hp_w30": "runs/p6_e2_posthoc/reg_hp_w30/seed0"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", type=int, default=0)
    a = ap.parse_args()
    dev = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    test_h5 = ROOT / "data" / "p6_domain_r" / "test.h5"
    corr = pd.read_csv(ROOT / MODELS["reg_l1"] / "test_y.csv").groupby("tile").hp_corr.mean()
    rows = []
    with h5py.File(test_h5, "r") as f:
        for t in f:
            y, v = f[t]["y"][:] * Q, f[t]["valid10"][:] > 0
            m = stat_mask(v)
            rows.append(dict(tile=t, hp_sd=float(highpass(y, v)[:, m].std(1).mean()), hp_corr=corr.get(t, np.nan)))
    d = pd.DataFrame(rows).dropna()
    top = d[d.hp_sd >= d.hp_sd.quantile(0.75)]
    tile = top.iloc[(top.hp_corr - top.hp_corr.median()).abs().argmin()].tile
    with h5py.File(test_h5, "r") as f:
        g = {k: f[tile][k][:] for k in ("x", "g", "r_raw", "y", "valid10")}
    x, gg = g["x"] * Q, g["g"] * Q
    y, rraw, v = g["y"] * Q, g["r_raw"] * Q, g["valid10"] > 0
    hp = highpass(y, v)
    energy = ndimage.uniform_filter((hp ** 2).mean(0) * stat_mask(v), CROP)
    half = CROP // 2
    energy[:half], energy[-half:], energy[:, :half], energy[:, -half:] = -1, -1, -1, -1
    cy, cx = np.unravel_index(np.argmax(energy), energy.shape)
    sl = (slice(None), slice(cy - half, cy + half), slice(cx - half, cx + half))
    xt = torch.from_numpy(x).float()[None].to(dev)
    gt = torch.from_numpy(gg).float()[None].to(dev)
    out = {"x_bicubic": F.interpolate(xt, scale_factor=2, mode="bicubic", align_corners=False)[0].cpu().numpy()[sl],
           "r_raw": rraw[sl], "y": y[sl]}
    for name, rd in MODELS.items():
        model = build_student(mtf).to(dev).eval()
        model.load_state_dict(torch.load(ROOT / rd / "student.pt", map_location=dev))
        with torch.no_grad():
            out[name] = model(xt, gt)[0].cpu().numpy()[sl]
    dst = ROOT / "runs" / "p6_e2" / "analysis"
    dst.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dst / "qualitative.npz", **out)
    info = dict(tile=tile, crop_centre=[int(cy), int(cx)], crop_px=CROP, n_top=len(top),
                reg_l1_hp_corr=float(top.set_index("tile").hp_corr[tile]), subset_median=float(top.hp_corr.median()))
    (dst / "qualitative.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
