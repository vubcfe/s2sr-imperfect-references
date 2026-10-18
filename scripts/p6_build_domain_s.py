"""P6 domain S: build self-aligned Sentinel-2-like pairs from the SEN2NEON AVIRIS-NG reference (plan §3).

For every clean tile of the train and validation sites (configs/p6_split.yaml; test sites are refused):
  y  (6, 256, 256)  target bands at 10 m   = A(z_target), MTF aggregation ratio 4 (as the real reference r)
  g  (4, 256, 256)  guide bands at 10 m    = A(z_guide), 10 m-band MTF, + sensor noise
  x  (6, 128, 128)  target bands at 20 m   = direct ratio-8 degradation with the 20 m-band MTF, per-band
                                             inter-band shift, + sensor noise
  valid10 / valid20 masks (reference nodata, eroded by the kernel reach)
All products are quantised to 1e-4 and stored as uint16 (reflectance x 10000) in one HDF5 per role.

The script also reports how far x is from D(y) (the LFC operator), since x is simulated directly from 2.5 m.

  python scripts/p6_build_domain_s.py --root data/sen2neon --out data/p6_domain_s
  python scripts/p6_build_domain_s.py --root data/sen2neon --out data/p6_domain_s_noise2 --noise-scale 2
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
import rasterio  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import GUIDE_BANDS, TARGET_BANDS  # noqa: E402
from s2sr.data.sen2neon import GUIDE_IDX, HR_NODATA, SCALE, TARGET_IDX, reference_valid  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.registration import subpixel_shift  # noqa: E402


def _site(tile: str) -> str:
    return tile.split("_")[1]


def tile_seed(seed: int, tile: str) -> int:
    """Per-tile seed independent of processing order (reproducible in parallel)."""
    return int(hashlib.sha256(f"{seed}:{tile}".encode()).hexdigest()[:8], 16)


def simulate(z: np.ndarray, ok: np.ndarray, cfg: dict, mtf: dict, rng: np.random.Generator) -> dict:
    """z: (12, 1024, 1024) reflectance at 2.5 m in SEN2NEON band order; ok: (1024, 1024) valid mask."""
    z = z.astype(np.float32).copy()
    for b in range(z.shape[0]):  # fill nodata before blurring; the masks exclude affected pixels
        if (~ok).any() and ok.any():
            z[b][~ok] = z[b][ok].mean()
    zt = torch.from_numpy(z)[None]
    s_t10, s_g10 = band_sigmas(mtf, TARGET_BANDS, ratio=4), band_sigmas(mtf, GUIDE_BANDS, ratio=4)
    s_t20 = band_sigmas(mtf, TARGET_BANDS, ratio=8)
    with torch.no_grad():
        y = degrade(zt[:, TARGET_IDX], s_t10, ratio=4)[0].numpy().astype(np.float64)
        g = degrade(zt[:, GUIDE_IDX], s_g10, ratio=4)[0].numpy().astype(np.float64)
        x = degrade(zt[:, TARGET_IDX], s_t20, ratio=8)[0].numpy().astype(np.float64)
    smax = cfg["interband_shift_max_px20"]
    shifts = rng.uniform(-smax, smax, size=(len(TARGET_BANDS), 2))
    x = np.stack([subpixel_shift(x[k], *shifts[k]) for k in range(len(TARGET_BANDS))])
    sig = cfg["noise_scale"]
    x += rng.normal(0, 1, x.shape) * np.array([cfg["noise_sigma"][b] for b in TARGET_BANDS])[:, None, None] * sig
    g += rng.normal(0, 1, g.shape) * np.array([cfg["noise_sigma"][b] for b in GUIDE_BANDS])[:, None, None] * sig
    q = cfg["quantisation"]
    to_u16 = lambda a: np.clip(np.round(a / q), 0, 65534).astype(np.uint16)  # noqa: E731
    v10 = reference_valid(ok, s_t10 + s_g10, 4)
    v20 = reference_valid(ok, s_t20, 8)
    return dict(y=to_u16(y), g=to_u16(g), x=to_u16(x), valid10=v10, valid20=v20, interband_shift=shifts)


def lfc_gap(y: np.ndarray, x: np.ndarray, v20: np.ndarray, mtf: dict) -> np.ndarray:
    """Per-band RMSE between D(y) (LFC operator, ratio 2) and the directly simulated x, on valid pixels."""
    with torch.no_grad():
        dy = degrade(torch.from_numpy(y.astype(np.float32))[None], band_sigmas(mtf, TARGET_BANDS, 2), 2)[0].numpy()
    d = (dy - x)[:, v20]
    return np.sqrt((d ** 2).mean(axis=1))


def _job(args):
    root, tile, lr_rel, hr_rel, cfg, mtf, seed = args
    torch.set_num_threads(1)
    with rasterio.open(root / hr_rel) as d:
        hr = d.read()
    ok = np.all(hr[TARGET_IDX + GUIDE_IDX] != HR_NODATA, axis=0)
    out = simulate(hr / SCALE, ok, cfg, mtf, np.random.default_rng(tile_seed(seed, tile)))
    q = cfg["quantisation"]
    out["lfc_gap"] = lfc_gap(out["y"] * q, out["x"] * q, out["valid20"], mtf)
    return tile, out


def build(root: Path, out_dir: Path, cfg: dict, mtf: dict, workers: int, limit: int | None = None) -> dict:
    split = yaml.safe_load((ROOT / cfg["split"]).read_text())
    allowed = {"train": set(split["sites"]["train"]), "val": set(split["sites"]["val"])}
    test_sites = set(split["sites"]["test"])
    tiles = {"train": split["tiles"]["train_S"]}
    # validation S tiles: every clean tile of the validation sites (S does not depend on |dt|)
    e0 = pd.read_csv(ROOT / split["e0_tiles"])
    f = split["filter"]
    clean = e0[(e0.status == "ok") & (e0.g_peak >= f["min_peak"]) & (e0.cloud_score <= f["max_cloud_score"])]
    tiles["val"] = sorted(map(str, clean[clean.site.isin(allowed["val"])].tile))
    meta = pd.read_csv(root / "metadata.csv").set_index("id")
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {}
    for role, names in tiles.items():
        names = names[:limit] if limit else names
        bad = [t for t in names if _site(t) in test_sites or _site(t) not in allowed[role]]
        if bad:  # leakage guard (plan §3, R3-2.4)
            raise ValueError(f"{role}: tiles from disallowed sites: {bad[:5]}")
        jobs = [(root, t, meta.loc[t, "lr"], meta.loc[t, "hr_2_5m_path"], cfg, mtf, cfg["seed"]) for t in names]
        gaps = []
        with h5py.File(out_dir / f"{role}.h5", "w") as h5, ProcessPoolExecutor(workers) as pool:
            h5.attrs["config"] = json.dumps(cfg)
            h5.attrs["target_bands"] = json.dumps(TARGET_BANDS)
            h5.attrs["guide_bands"] = json.dumps(GUIDE_BANDS)
            for i, (tile, o) in enumerate(pool.map(_job, jobs, chunksize=4)):
                grp = h5.create_group(tile)
                for k in ("y", "g", "x"):
                    grp.create_dataset(k, data=o[k], compression="lzf")
                for k in ("valid10", "valid20"):
                    grp.create_dataset(k, data=o[k], compression="lzf")
                grp.attrs["interband_shift_px20"] = o["interband_shift"].tolist()
                grp.attrs["site"] = _site(tile)
                gaps.append(o["lfc_gap"])
                if (i + 1) % 100 == 0:
                    print(f"{role}: {i + 1}/{len(jobs)}", flush=True)
        gaps = np.array(gaps)
        report[role] = dict(n_tiles=len(names), n_sites=len({_site(t) for t in names}),
                            lfc_gap_rmse_median=dict(zip(TARGET_BANDS, np.median(gaps, axis=0).round(5).tolist())))
    (out_dir / "build_report.json").write_text(json.dumps(report, indent=2))
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="data/sen2neon")
    ap.add_argument("--out", default="data/p6_domain_s")
    ap.add_argument("--config", default="configs/p6_domain_s.yaml")
    ap.add_argument("--noise-scale", type=float, default=None)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    cfg = yaml.safe_load((ROOT / a.config).read_text())
    if a.noise_scale is not None:
        cfg["noise_scale"] = a.noise_scale
    mtf = yaml.safe_load((ROOT / cfg["mtf_config"]).read_text())["mtf_nyquist"]
    print(json.dumps(build(Path(a.root), Path(a.out), cfg, mtf, a.workers, a.limit), indent=2))


if __name__ == "__main__":
    main()
