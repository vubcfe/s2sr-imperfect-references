"""Contrast gain of a reference against Sentinel-2 at several spatial scales (revision of the P6 paper).

Per tile and band, on the 20 m grid: ratio sd(B(D(ref))) / sd(B(x)), where B is a band-pass (difference of
Gaussians, sigma in metres) or the 60 m low-pass used by the gain inversion (Eq. 3 of the paper), on pixels whose
whole filter footprint is valid.
  domain R (real):      ref = r_reg (registered, contrast as delivered), x = native Sentinel-2 20 m bands
  domain S (simulated): ref = y (true, gain 1) and ref = y_h (harmonised, known gain m), x = simulated 20 m input
In domain S the true gain is scale-independent by construction (1 for y, m for y_h), so any drop of the measured
ratio at fine scales there is a measurement effect (noise, aliasing, MTF mismatch of x), not a property of the
reference. Comparing that drop with the one on real data tells whether the real fine-scale gain differs from the
low-pass gain beyond the measurement effect.

  python scripts/p6_gain_scales.py --out runs/p6_e2/analysis/gain_scales.csv
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402

Q = 1e-4
GRID = 20.0
SCALES = {"bp_20_40": (20, 40), "bp_40_80": (40, 80), "bp_80_160": (80, 160),  # larger scales leave no valid pixel on a 2.56 km tile
          "lp_60": (None, 60)}


def _smooth(z, s, v):
    return ndimage.gaussian_filter(np.where(v, z, z[v].mean()), s)


def ratios(x20, r20, v20) -> dict:
    out = {}
    for name, (s1, s2) in SCALES.items():
        reach = int(math.ceil(4 * s2 / GRID))
        m = ndimage.binary_erosion(v20, iterations=reach, border_value=0)
        if m.sum() < 200:
            continue
        vals = []
        for b in range(x20.shape[0]):
            if s1 is None:
                fx, fr = _smooth(x20[b], s2 / GRID, v20), _smooth(r20[b], s2 / GRID, v20)
            else:
                fx = _smooth(x20[b], s1 / GRID, v20) - _smooth(x20[b], s2 / GRID, v20)
                fr = _smooth(r20[b], s1 / GRID, v20) - _smooth(r20[b], s2 / GRID, v20)
            vals.append(fr[m].std() / max(fx[m].std(), 1e-12))
        out[name] = vals
    return out


def _job(args):
    kind, path, tpath, tile, sig = args
    torch.set_num_threads(1)
    rows = []
    with h5py.File(path, "r") as f:
        g = f[tile]
        x = g["x"][:] * Q
        if kind == "R":
            refs = {"r_reg": g["r_reg"][:] * Q}
            v20 = g["valid20"][:] > 0
        else:
            refs = {"y": g["y"][:] * Q}
            v10 = g["valid10"][:] > 0
            v20 = v10.reshape(128, 2, 128, 2).all(axis=(1, 3))
            with h5py.File(tpath, "r") as ft:
                refs["y_h"] = ft[tile]["y_h"][:] * Q
                m_true = np.array(ft[tile].attrs["m"])
    for rk, r in refs.items():
        with torch.no_grad():
            d = degrade(torch.from_numpy(r).float()[None], sig, 2)[0].double().numpy()
        res = ratios(x.astype(np.float64), d, v20)
        for scale, vals in res.items():
            for b, val in zip(TARGET_BANDS, vals):
                row = dict(domain=kind, tile=tile, ref=rk, scale=scale, band=b, ratio=val)
                if kind == "S":
                    row["m_true"] = 1.0 if rk == "y" else float(m_true[TARGET_BANDS.index(b)])
                rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/p6_e2/analysis/gain_scales.csv")
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    sig = band_sigmas(mtf, TARGET_BANDS, 2)
    jobs = []
    for role in ("train", "val", "test"):
        p = ROOT / "data/p6_domain_r" / f"{role}.h5"
        with h5py.File(p, "r") as f:
            jobs += [("R", p, None, t, sig) for t in f]
    ps, pt = ROOT / "data/p6_domain_s/val.h5", ROOT / "data/p6_e1r/targets_val.h5"
    with h5py.File(ps, "r") as f:
        jobs += [("S", ps, pt, t, sig) for t in f]
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for r in ex.map(_job, jobs, chunksize=4):
            rows += r
    d = pd.DataFrame(rows)
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(out, index=False)
    med = d.groupby(["domain", "ref", "scale"]).ratio.median().unstack("scale")[list(SCALES)]
    print(med.round(3).to_string())
    if "m_true" in d:
        s = d[d.domain == "S"].copy()
        s["rel"] = s.ratio / s.m_true
        print("domain S, ratio / true gain:")
        print(s.groupby(["ref", "scale"]).rel.median().unstack("scale")[list(SCALES)].round(3).to_string())
    r = d[d.domain == "R"].pivot_table(index=["tile", "band"], columns="scale", values="ratio")
    rel = r.div(r["lp_60"], axis=0)
    print("domain R, ratio / low-pass gain (median, IQR):")
    print(rel.quantile([0.25, 0.5, 0.75]).round(3).to_string())


if __name__ == "__main__":
    main()
