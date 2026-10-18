"""Positive control of the registration + gain correction on the 10 m bands (revision after the GPT review).

For B2, B3, B4 and B8 Sentinel-2 observes the 10 m truth, so the correction of the SEN2NEON reference can be
checked directly. Per tile (the 521 domain-R tiles of the site split) and band:
  r     = 10 m aggregate (MTF-matched, ratio 4) of the released 2.5 m reference
  r_reg = r shifted by minus the guide-band tile shift of E0
  y_c   = mean(g) + (r_reg - mean(r_reg)) / a_hat, a_hat = sd(L r_reg) / sd(L g), L = Gaussian low-pass sigma 60 m
and each version is compared with the observed Sentinel-2 band g at 10 m with the paper's high-pass metrics
(high-pass sigma 20 m), i.e. at frequencies above the band used to estimate the gain.

  python scripts/p6_guide_control.py --workers 16
Writes runs/p6_e2/analysis/guide_control.csv and prints medians.
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
import rasterio  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import GUIDE_BANDS  # noqa: E402
from s2sr.data.sen2neon import GUIDE_IDX, HR_NODATA, LR_NODATA, SCALE, reference_valid  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.p6 import metrics as M  # noqa: E402
from s2sr.registration import subpixel_shift  # noqa: E402

LOWPASS_PX = 60.0 / 10.0


def _job(args):
    root, tile, lr_rel, hr_rel, shift, sig4 = args
    torch.set_num_threads(1)
    with rasterio.open(root / lr_rel) as d:
        lr = d.read()
    with rasterio.open(root / hr_rel) as d:
        hr = d.read()
    g = lr[GUIDE_IDX].astype(np.float64) / SCALE
    lr_ok = np.all(lr[GUIDE_IDX] != LR_NODATA, axis=0)
    z = hr[GUIDE_IDX].astype(np.float32) / SCALE
    hr_ok = np.all(hr[GUIDE_IDX] != HR_NODATA, axis=0)
    for b in range(z.shape[0]):
        if (~hr_ok).any() and hr_ok.any():
            z[b][~hr_ok] = z[b][hr_ok].mean()
    with torch.no_grad():
        r = degrade(torch.from_numpy(z)[None], sig4, 4)[0].double().numpy()
    v = lr_ok & reference_valid(hr_ok, sig4, 4)
    dy, dx = shift
    v = ndimage.binary_erosion(v, iterations=math.ceil(max(abs(dy), abs(dx))) + 3, border_value=1)
    if v.mean() < 0.3:
        return []
    r_reg = subpixel_shift(r, -dy, -dx)
    m = ndimage.binary_erosion(v, iterations=int(math.ceil(4 * LOWPASS_PX)), border_value=0)
    if m.sum() < 500:
        return []
    a_hat = np.empty(4)
    for b in range(4):
        lg = ndimage.gaussian_filter(np.where(v, g[b], g[b][v].mean()), LOWPASS_PX)[m]
        lr_ = ndimage.gaussian_filter(np.where(v, r_reg[b], r_reg[b][v].mean()), LOWPASS_PX)[m]
        a_hat[b] = lr_.std() / max(lg.std(), 1e-12)
    gbar = np.array([g[b][v].mean() for b in range(4)])
    rbar = np.array([r_reg[b][v].mean() for b in range(4)])
    y_c = gbar[:, None, None] + (r_reg - rbar[:, None, None]) / a_hat[:, None, None]
    rows = []
    for name, ref in (("r", r), ("r_reg", r_reg), ("y_c", y_c)):
        hc, hps, ee = M.hp_corr(ref, g, v), M.hp_psnr(ref, g, v), M.edge_error(ref, g, v)
        hg = M.hp_gain(ref, g, v)  # sd HP(ref) / sd HP(g): 1 = the reference has the Sentinel-2 detail amplitude
        rows += [dict(tile=tile, site=tile.split("_")[1], band=b, ref=name, a_hat=a_hat[k], hp_corr=hc[k],
                      hp_psnr=hps[k], edge_error=ee[k], hp_gain=hg[k]) for k, b in enumerate(GUIDE_BANDS)]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sen2neon", default="data/sen2neon")
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    root = ROOT / a.sen2neon
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    sig4 = band_sigmas(mtf, GUIDE_BANDS, ratio=4)
    e0 = pd.read_csv(ROOT / "results/p6_e0/tiles.csv").set_index("tile")
    meta = pd.read_csv(root / "metadata.csv").set_index("id")
    tiles = []
    for role in ("train", "val", "test"):
        with h5py.File(ROOT / "data/p6_domain_r" / f"{role}.h5", "r") as f:
            tiles += list(f)
    jobs = [(root, t, meta.loc[t, "lr"], meta.loc[t, "hr_2_5m_path"],
             (float(e0.loc[t, "g_dy"]), float(e0.loc[t, "g_dx"])), sig4) for t in tiles]
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for r in ex.map(_job, jobs, chunksize=4):
            rows += r
    d = pd.DataFrame(rows)
    out = ROOT / "runs/p6_e2/analysis/guide_control.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(out, index=False)
    print("tiles", d.tile.nunique(), "sites", d.site.nunique())
    print(d.groupby(["ref", "band"])[["hp_gain", "hp_corr", "hp_psnr", "edge_error"]].median().round(3).to_string())
    print(d.groupby("ref")[["hp_gain", "hp_corr", "hp_psnr", "edge_error"]].median().round(3).to_string())
    print("a_hat median per band:", d[d.ref == "y_c"].groupby("band").a_hat.median().round(3).to_dict())


if __name__ == "__main__":
    main()
