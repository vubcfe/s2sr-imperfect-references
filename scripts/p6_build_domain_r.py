"""P6 domain R: real Sentinel-2 inputs paired with the SEN2NEON reference, prepared for E2 (plan §4.1, §7).

For every clean tile of configs/p6_split.yaml (train: |dt| <= 15 days of train sites; val / test: |dt| <= 7 days):
  x      (6, 128, 128)  native 20 m target bands (2x2 replicas undone; non-replicated blocks masked)
  g      (4, 256, 256)  10 m guide bands B2, B3, B4, B8
  r_raw  (6, 256, 256)  SEN2NEON reference aggregated to 10 m with the MTF (ratio 4), as delivered
  r_reg  r_raw shifted by minus the joint target-band tile shift of E0 (target_shift = t_dy, t_dx)
  y      = r_corr: r_reg with its contrast restored by one frequency-independent gain per band,
           y = mean(x) + (r_reg - mean(r_reg)) / a_hat, a_hat = low-pass (sigma 60 m) sd ratio of D(r_reg) vs x
           (inverts the OLS harmonisation of SEN2NEON, which compresses contrast by about R; notes/baselines_fulltext.md)
  valid10 / valid20 masks
Layout matches domain S (`y` is the default training/evaluation target), so the same training code reads it; the
other references are separate datasets of the same group. Stored uint16 x 1e4.

  python scripts/p6_build_domain_r.py --sen2neon data/sen2neon --out data/p6_domain_r
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
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

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.data.sen2neon import (GUIDE_IDX, HR_NODATA, LR_NODATA, SCALE, TARGET_IDX,  # noqa: E402
                                reference_valid)
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.registration import subpixel_shift  # noqa: E402

Q = 1e-4
LOWPASS_M = 60.0
MAX_RESAMPLED = 0.01


def lowpass_gain(x20, r20, valid20):
    s = LOWPASS_M / 20.0
    m = ndimage.binary_erosion(valid20, iterations=int(math.ceil(4 * s)), border_value=0)
    out = []
    for b in range(x20.shape[0]):
        lx = ndimage.gaussian_filter(np.where(valid20, x20[b], x20[b][valid20].mean()), s)[m]
        lr = ndimage.gaussian_filter(np.where(valid20, r20[b], r20[b][valid20].mean()), s)[m]
        out.append(lr.std() / max(lx.std(), 1e-12))
    return np.array(out)


def prepare(lr: np.ndarray, hr: np.ndarray, shift10, mtf: dict) -> dict:
    lr_ok = np.all(lr != LR_NODATA, axis=0)
    t = lr[TARGET_IDX]
    blk = t.reshape(6, 128, 2, 128, 2)
    rep_ok = ~(blk.max(axis=(2, 4)) != blk.min(axis=(2, 4))).any(axis=0)
    if (~rep_ok).mean() > MAX_RESAMPLED:
        raise ValueError("20 m bands resampled")
    x = t[:, ::2, ::2].astype(np.float64) / SCALE
    g = lr[GUIDE_IDX].astype(np.float64) / SCALE
    z = hr[TARGET_IDX].astype(np.float32) / SCALE
    hr_ok = np.all(hr[TARGET_IDX] != HR_NODATA, axis=0)
    for b in range(6):
        if (~hr_ok).any() and hr_ok.any():
            z[b][~hr_ok] = z[b][hr_ok].mean()
    sig4 = band_sigmas(mtf, TARGET_BANDS, ratio=4)
    with torch.no_grad():
        r_raw = degrade(torch.from_numpy(z)[None], sig4, 4)[0].double().numpy()
    v10 = lr_ok & reference_valid(hr_ok, sig4, 4)
    dy, dx = shift10
    border = math.ceil(max(abs(dy), abs(dx))) + 3
    v10 = ndimage.binary_erosion(v10, iterations=border, border_value=1)
    r_reg = subpixel_shift(r_raw, -dy, -dx)
    with torch.no_grad():
        d20 = degrade(torch.from_numpy(r_reg).float()[None], band_sigmas(mtf, TARGET_BANDS, 2), 2)[0].double().numpy()
    v20 = v10.reshape(128, 2, 128, 2).all(axis=(1, 3)) & rep_ok
    a_hat = lowpass_gain(x, d20, v20)
    xbar = np.array([x[b][v20].mean() for b in range(6)])
    rbar = np.array([r_reg[b][v10].mean() for b in range(6)])
    y = xbar[:, None, None] + (r_reg - rbar[:, None, None]) / a_hat[:, None, None]
    return dict(x=x, g=g, r_raw=r_raw, r_reg=r_reg, y=y, valid10=v10, valid20=v20, a_hat=a_hat)


def _job(args):
    root, tile, lr_rel, hr_rel, shift10, mtf = args
    torch.set_num_threads(1)
    with rasterio.open(root / lr_rel) as d:
        lr = d.read()
    with rasterio.open(root / hr_rel) as d:
        hr = d.read()
    try:
        return tile, prepare(lr, hr, shift10, mtf)
    except ValueError as e:
        return tile, str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sen2neon", default="data/sen2neon")
    ap.add_argument("--split", default="configs/p6_split.yaml")
    ap.add_argument("--out", default="data/p6_domain_r")
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args()
    root = Path(a.sen2neon)
    split = yaml.safe_load((ROOT / a.split).read_text())
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    e0 = pd.read_csv(ROOT / split["e0_tiles"]).set_index("tile")
    meta = pd.read_csv(root / "metadata.csv").set_index("id")
    roles = {"train": split["tiles"]["train_R"], "val": split["tiles"]["val_eval"], "test": split["tiles"]["test_eval"],
             "train_dt40": split["tiles"]["train_R_dt40"]}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    q16 = lambda v: np.clip(np.round(v / Q), 0, 65534).astype(np.uint16)  # noqa: E731
    report, gains = {}, []
    for role, tiles in roles.items():
        jobs = [(root, t, meta.loc[t, "lr"], meta.loc[t, "hr_2_5m_path"],
                 (float(e0.loc[t, "t_dy"]), float(e0.loc[t, "t_dx"])), mtf) for t in tiles]
        n_ok, skipped = 0, []
        with h5py.File(out / f"{role}.h5", "w") as h5, ProcessPoolExecutor(a.workers) as ex:
            for tile, o in ex.map(_job, jobs, chunksize=4):
                if isinstance(o, str):
                    skipped.append((tile, o))
                    continue
                grp = h5.create_group(tile)
                for k in ("x", "g", "r_raw", "r_reg", "y"):
                    grp.create_dataset(k, data=q16(o[k]), compression="lzf")
                grp.create_dataset("valid10", data=o["valid10"], compression="lzf")
                grp.create_dataset("valid20", data=o["valid20"], compression="lzf")
                grp.attrs["a_hat"] = o["a_hat"].tolist()
                grp.attrs["site"] = tile.split("_")[1]
                gains += [dict(role=role, tile=tile, band=b, a_hat=v) for b, v in zip(TARGET_BANDS, o["a_hat"])]
                n_ok += 1
        report[role] = dict(n_tiles=n_ok, n_sites=len({t.split("_")[1] for t in tiles}), skipped=skipped)
        print(role, n_ok, "tiles", len(skipped), "skipped", flush=True)
    g = pd.DataFrame(gains)
    g.to_csv(out / "gains.csv", index=False)
    report["a_hat_median"] = g.groupby("band").a_hat.median().round(4).to_dict()
    (out / "build_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report["a_hat_median"], indent=2))


if __name__ == "__main__":
    main()
