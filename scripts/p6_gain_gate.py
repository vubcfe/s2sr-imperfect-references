"""P6 pre-E2 gate: does the E0 gain estimator recover a known contrast gain when the Sentinel-2 side carries
sensor noise, aliasing and inter-band misregistration? (plan §8, review rounds 2-3)

On domain-S tiles the truth is known: x is simulated from the AVIRIS reference exactly as in domain S (MTF, ratio 8,
inter-band shift, Gaussian noise at `noise_scale` x the configured S2 noise), and the 20 m reference r is the
noise-free MTF degradation of the same scene. A known uniform contrast gain g is applied to r about its tile mean
(what a per-tile linear harmonisation would do), then E0's `_radiometry` estimates the low-pass gain (a_sd) and the
high-pass gain (hp_gain_sd, hp_gain_ols). The ratio estimate / g shows the bias caused by the S2 side alone.

  python scripts/p6_gain_gate.py --root data/sen2neon --every 4 --out results/p6_gain_gate
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


E0 = _load("e0", "p6_e0_reference_errors.py")
DS = _load("dom_s", "p6_build_domain_s.py")
from s2sr.bands import TARGET_BANDS  # noqa: E402
from s2sr.data.sen2neon import GUIDE_IDX, HR_NODATA, SCALE, TARGET_IDX, reference_valid  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402

GAINS = (1.0, 0.75, 0.5)
NOISE = (0.0, 1.0, 2.0)


def _job(args):
    root, tile, hr_rel, cfg, mtf = args
    torch.set_num_threads(1)
    with rasterio.open(root / hr_rel) as d:
        hr = d.read()
    ok = np.all(hr[TARGET_IDX + GUIDE_IDX] != HR_NODATA, axis=0)
    z = hr / SCALE
    s_t20 = band_sigmas(mtf, TARGET_BANDS, ratio=8)
    zf = z.astype(np.float32).copy()
    for b in range(zf.shape[0]):
        if (~ok).any() and ok.any():
            zf[b][~ok] = zf[b][ok].mean()
    with torch.no_grad():
        r = degrade(torch.from_numpy(zf[TARGET_IDX])[None], s_t20, ratio=8)[0].numpy().astype(np.float64)
    v20 = reference_valid(ok, s_t20, 8)
    m = np.ones_like(v20) & v20
    ms = E0._stat_mask(m, 20)
    rows = []
    if ms.sum() < 200:
        return rows
    for ns in NOISE:
        o = DS.simulate(z, ok, {**cfg, "noise_scale": ns}, mtf,
                        np.random.default_rng(DS.tile_seed(cfg["seed"], tile)))
        x = o["x"] * cfg["quantisation"]
        for g in GAINS:
            mean = r[:, ms].mean(axis=1)[:, None, None]
            rg = mean + g * (r - mean)
            for k, band in enumerate(TARGET_BANDS):
                rad = E0._radiometry(x[k], rg[k], ms, 20, m)
                rows.append(dict(tile=tile, noise_scale=ns, gain=g, band=band, a_sd=rad["a_sd"],
                                 hp_gain_sd=rad["hp_gain_sd"], hp_gain_ols=rad["hp_gain_ols"], hp_corr=rad["hp_corr"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="data/sen2neon")
    ap.add_argument("--every", type=int, default=4)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default="results/p6_gain_gate")
    a = ap.parse_args()
    cfg = yaml.safe_load((ROOT / "configs" / "p6_domain_s.yaml").read_text())
    mtf = yaml.safe_load((ROOT / cfg["mtf_config"]).read_text())["mtf_nyquist"]
    split = yaml.safe_load((ROOT / cfg["split"]).read_text())
    tiles = split["tiles"]["train_S"][::a.every]
    meta = pd.read_csv(Path(a.root) / "metadata.csv").set_index("id")
    jobs = [(Path(a.root), t, meta.loc[t, "hr_2_5m_path"], cfg, mtf) for t in tiles]
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for r in ex.map(_job, jobs, chunksize=2):
            rows += r
    df = pd.DataFrame(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "gate.csv", index=False)
    s = df.groupby(["noise_scale", "gain"])[["a_sd", "hp_gain_sd", "hp_gain_ols", "hp_corr"]].median()
    for c in ("a_sd", "hp_gain_sd", "hp_gain_ols"):
        s[c + "/gain"] = s[c] / s.index.get_level_values("gain")
    per_band = df[df.gain == 1.0].groupby(["noise_scale", "band"]).hp_gain_sd.median().unstack()
    text = (f"# Gain-estimator gate ({df.tile.nunique()} domain-S tiles, every {a.every}th train tile)\n\n"
            "Medians over tiles and bands. Truth: low-pass and high-pass gain = `gain`.\n\n"
            + s.round(4).to_markdown() + "\n\nhp_gain_sd at true gain 1, per band:\n\n"
            + per_band.reindex(columns=TARGET_BANDS).round(4).to_markdown() + "\n")
    (out / "summary.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
