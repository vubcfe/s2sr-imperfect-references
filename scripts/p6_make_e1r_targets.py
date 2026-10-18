"""P6 E1-R: emulate the SEN2NEON radiometric harmonisation on domain S, where the true target y is known.

SEN2NEON (Donike et al., GRSL 2026, Eq. 2) fits, per band, rho_S2 = m * rho_synth + n by least squares with the
Sentinel-2 image linearly interpolated to the 2.5 m grid, and uses the fitted values as the reference. On domain S,
x_S plays the role of Sentinel-2 and the 2.5 m AVIRIS-derived reference z that of rho_synth. Per tile and band:
  x_up  = bilinear upsampling of x_S (20 m) to 2.5 m (x8)
  m, n  = OLS of x_up on z over valid pixels  (slope = R * sd(x_up) / sd(z): contrast compressed by about R)
  y_h   = m * y + n                            (y = A(z) at 10 m, so this equals A(m z + n))
Inversion proposed in the plan (§4.1): one frequency-independent gain a_hat per tile and band, estimated
symmetrically on the low-pass (sigma 60 m) between x_S and D(y_h) on the 20 m grid:
  y_inv = mean(x) + (y_h - mean(y_h)) / a_hat
Stored per tile (train and val roles, uint16 x 1e4): y_h, y_inv; attributes m, n, a_hat.

  python scripts/p6_make_e1r_targets.py --s data/p6_domain_s --sen2neon data/sen2neon --out data/p6_e1r
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402
import yaml  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.data.sen2neon import HR_NODATA, SCALE, TARGET_IDX  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402

Q = 1e-4
LOWPASS_M = 60.0


def harmonise(x20: np.ndarray, z: np.ndarray, ok: np.ndarray):
    """SEN2NEON Eq. 2 emulation. x20: (C, h, w) 20 m; z: (C, 8h, 8w) 2.5 m; ok: (8h, 8w). Returns m, n, R."""
    xu = F.interpolate(torch.from_numpy(x20).float()[None], scale_factor=8, mode="bilinear",
                       align_corners=False)[0].double().numpy()
    m, n, R = [], [], []
    for b in range(x20.shape[0]):
        xv, zv = xu[b][ok], z[b][ok]
        slope = np.cov(xv, zv)[0, 1] / zv.var()
        m.append(slope)
        n.append(xv.mean() - slope * zv.mean())
        R.append(np.corrcoef(xv, zv)[0, 1])
    return np.array(m), np.array(n), np.array(R)


def symmetric_lowpass_gain(x20: np.ndarray, r20: np.ndarray, valid20: np.ndarray, grid_m: float = 20.0):
    """Plan §4.1 gain: sd ratio of the low-pass (sigma 60 m) of r vs x, per band, on valid pixels."""
    s = LOWPASS_M / grid_m
    m = ndimage.binary_erosion(valid20, iterations=int(np.ceil(4 * s)), border_value=0)
    out = []
    for b in range(x20.shape[0]):
        lx, lr = ndimage.gaussian_filter(x20[b], s)[m], ndimage.gaussian_filter(r20[b], s)[m]
        out.append(lr.std() / max(lx.std(), 1e-12))
    return np.array(out)


def _job(args):
    s_path, tile, hr_path, mtf = args
    torch.set_num_threads(1)
    with h5py.File(s_path, "r") as f:
        y = f[tile]["y"][:] * Q
        x = f[tile]["x"][:] * Q
        v20 = f[tile]["valid20"][:]
    with rasterio.open(hr_path) as d:
        hr = d.read()
    z = hr[TARGET_IDX].astype(np.float64) / SCALE
    ok = np.all(hr[TARGET_IDX] != HR_NODATA, axis=0)
    m, n, R = harmonise(x, z, ok)
    yh = m[:, None, None] * y + n[:, None, None]
    with torch.no_grad():
        d20 = degrade(torch.from_numpy(yh).float()[None], band_sigmas(mtf, TARGET_BANDS, 2), 2)[0].double().numpy()
    a_hat = symmetric_lowpass_gain(x, d20, v20)
    vm = ndimage.zoom(v20, 2, order=0).astype(bool)
    xbar = np.array([x[b][v20].mean() for b in range(6)])
    yhbar = np.array([yh[b][vm].mean() for b in range(6)])
    yinv = xbar[:, None, None] + (yh - yhbar[:, None, None]) / a_hat[:, None, None]
    to16 = lambda a: np.clip(np.round(a / Q), 0, 65534).astype(np.uint16)  # noqa: E731
    return tile, to16(yh), to16(yinv), dict(m=m.tolist(), n=n.tolist(), R=R.tolist(), a_hat=a_hat.tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s", default="data/p6_domain_s")
    ap.add_argument("--sen2neon", default="data/sen2neon")
    ap.add_argument("--out", default="data/p6_e1r")
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    meta = pd.read_csv(Path(a.sen2neon) / "metadata.csv").set_index("id")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for role in ("train", "val"):
        s_path = Path(a.s) / f"{role}.h5"
        with h5py.File(s_path, "r") as f:
            tiles = sorted(f)[:a.limit] if a.limit else sorted(f)
        jobs = [(s_path, t, Path(a.sen2neon) / meta.loc[t, "hr_2_5m_path"], mtf) for t in tiles]
        name = "targets_S4_harmonized.h5" if role == "train" else "targets_val.h5"  # train file read by p6_e1.py
        with h5py.File(out / name, "w") as h5, ProcessPoolExecutor(a.workers) as ex:
            for tile, yh, yinv, p in ex.map(_job, jobs, chunksize=4):
                g = h5.create_group(tile)
                g.create_dataset("y_h", data=yh, compression="lzf")
                g.create_dataset("y_inv", data=yinv, compression="lzf")
                for k, v in p.items():
                    g.attrs[k] = v
                for b, band in enumerate(TARGET_BANDS):
                    rows.append(dict(role=role, tile=tile, band=band, m=p["m"][b], R=p["R"][b], a_hat=p["a_hat"][b]))
    df = pd.DataFrame(rows)
    df.to_csv(out / "harmonisation_params.csv", index=False)
    summ = df.groupby("band")[["m", "R", "a_hat"]].median().reindex(TARGET_BANDS)
    print(summ.round(3).to_string())
    (out / "summary.json").write_text(json.dumps(summ.round(4).to_dict(), indent=2))


if __name__ == "__main__":
    main()
