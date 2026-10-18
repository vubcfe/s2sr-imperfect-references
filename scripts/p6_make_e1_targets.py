"""P6 E1: corrupted training targets on domain S, one HDF5 per scenario (configs/p6_e1.yaml).

For every training tile of domain S the clean 10 m target y is corrupted as a cross-sensor reference would be
(global shift drawn from the E0 tile-shift distribution or fixed; optional smooth non-rigid field; optional local
content change taken from another tile), with per-tile parameters from a seeded generator. Then the same
registration as for real references is applied: one joint target-band shift estimated by phase correlation
between x and D(y_corrupted) on the 20 m grid (E0 procedure), and y_reg = subpixel_shift(y_corrupted, -shift).

Stored per tile: y_c (corrupted), y_reg (registered), params, estimated shift (uint16 x 1e4 for images).

  python scripts/p6_make_e1_targets.py --config configs/p6_e1.yaml --s data/p6_domain_s --out data/p6_e1
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.p6 import corrupt  # noqa: E402
from s2sr.registration import phase_corr, subpixel_shift  # noqa: E402

Q = 1e-4


def _seed(seed, *keys):
    return int(hashlib.sha256(":".join(map(str, (seed,) + keys)).encode()).hexdigest()[:8], 16)


def corrupt_tile(y, other, scenario, rng, pool):
    """y, other: (6, H, W) float reflectance. Returns corrupted target and the drawn parameters."""
    p = corrupt.draw_tile_params(scenario, rng, pool)
    yt = torch.from_numpy(y).float()[None]
    h, w = y.shape[-2:]
    dy = torch.full((1, h, w), p["dy"])
    dx = torch.full((1, h, w), p["dx"])
    if p["field_rms"] > 0:
        fy, fx = corrupt.smooth_field(1, h, w, p["field_rms"], scenario.get("nonrigid_corr_px", 16), rng)
        dy, dx = dy + fy, dx + fx
    yc = corrupt.warp(yt, dy, dx)
    if p["change"]:
        yc = corrupt.content_change(yc, torch.from_numpy(other).float()[None], scenario["change_area"], rng)
    return yc[0].numpy().astype(np.float64), p


def register(x20, yc, valid20, mtf):
    """E0 procedure: joint target-band shift between x and D(yc) on the 20 m grid; returns (dy, dx) in 10 m px."""
    with torch.no_grad():
        d = degrade(torch.from_numpy(yc).float()[None], band_sigmas(mtf, TARGET_BANDS, 2), 2)[0].numpy()
    ty, tx, pk = phase_corr(x20, d, valid20)
    if not np.isfinite(ty):
        return 0.0, 0.0, 0.0
    return 2 * ty, 2 * tx, pk


def _job(args):
    s_path, tile, other_tile, scenario, pool, seed, mtf = args
    torch.set_num_threads(1)
    with h5py.File(s_path, "r") as f:
        y = f[tile]["y"][:] * Q
        x = f[tile]["x"][:] * Q
        v20 = f[tile]["valid20"][:]
        other = f[other_tile]["y"][:] * Q
    rng = np.random.default_rng(_seed(seed, scenario["name"], tile))
    yc, p = corrupt_tile(y, other, scenario, rng, pool)
    dy, dx, pk = register(x, yc, v20, mtf)
    yreg = subpixel_shift(yc, -dy, -dx)
    to16 = lambda a: np.clip(np.round(a / Q), 0, 65534).astype(np.uint16)  # noqa: E731
    return tile, to16(yc), to16(yreg), dict(p, est_dy=dy, est_dx=dx, est_peak=pk, other=other_tile)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/p6_e1.yaml")
    ap.add_argument("--s", default="data/p6_domain_s")
    ap.add_argument("--out", default="data/p6_e1")
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    cfg = yaml.safe_load((ROOT / a.config).read_text())
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    e0 = pd.read_csv(ROOT / cfg["e0_tiles"])
    e0 = e0[(e0.status == "ok") & (e0.g_peak >= 0.1)]
    pool = e0[["t_dy", "t_dx"]].to_numpy()   # empirical joint target-band tile shifts (10 m px)
    s_path = Path(a.s) / "train.h5"
    with h5py.File(s_path, "r") as f:
        tiles = sorted(f)
    tiles = tiles[:a.limit] if a.limit else tiles
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for sc in cfg["scenarios"]:
        rng = np.random.default_rng(_seed(cfg["seed"], sc["name"], "other"))
        others = [tiles[(i + 1 + rng.integers(len(tiles) - 1)) % len(tiles)] for i in range(len(tiles))]
        jobs = [(s_path, t, o, sc, pool, cfg["seed"], mtf) for t, o in zip(tiles, others)]
        rows = []
        with h5py.File(out / f"targets_{sc['name']}.h5", "w") as h5, ProcessPoolExecutor(a.workers) as ex:
            h5.attrs["scenario"] = json.dumps(sc)
            for tile, yc, yreg, p in ex.map(_job, jobs, chunksize=4):
                g = h5.create_group(tile)
                g.create_dataset("y_c", data=yc, compression="lzf")
                g.create_dataset("y_reg", data=yreg, compression="lzf")
                g.attrs.update({k: v for k, v in p.items()})
                rows.append(dict(tile=tile, **p))
        df = pd.DataFrame(rows)
        df.to_csv(out / f"targets_{sc['name']}.csv", index=False)
        err = np.hypot(df.est_dy - df.dy, df.est_dx - df.dx)
        print(f"{sc['name']}: {len(df)} tiles; registration error (10 m px) median {err.median():.3f}, "
              f"q90 {err.quantile(0.9):.3f}; non-rigid {int((df.field_rms > 0).sum())}, change {int(df.change.sum())}",
              flush=True)


if __name__ == "__main__":
    main()
