"""P6 experiment E0: characterise geometric and radiometric errors of the SEN2NEON reference.

For every tile the AVIRIS-NG reference (2.5 m) is degraded with the Sentinel-2 MTF to the grid of
the observation it is compared with, then:

  geometry   global sub-pixel shift on the 10 m guide bands (B2 B3 B4 B8, most texture), a joint
             shift on the six 20 m target bands, a shift per target band, and block-wise shifts on
             both band groups (640 m blocks) to measure what one translation or an affine field leaves.
  radiometry after removing the global shift (and, for target bands, also after each band's own
             shift): gain/offset r = a x + c from low-pass images (OLS, symmetric, and cross-fitted
             on tile halves), RMSE before/after, high-pass gain and correlation, and the full-
             resolution correlation R and standard-deviation ratio (to test whether the reference
             was harmonised by a per-tile regression, which shrinks its contrast to about R).

All shifts are in 10 m pixels, sign convention of s2sr.registration (r is x shifted by d). Filter
scales are in metres, so 10 m and 20 m bands are comparable. Alignment uses a non-periodic cubic-
spline shift (an FFT shift would ring across the whole tile from the opposite edge).

Round-1 review fixes (review/round1_*.md, review/response_round1.md) are marked "R1-<reviewer>.<n>".

  python scripts/p6_e0_reference_errors.py --root data/sen2neon --out results/p6_e0
  python scripts/p6_e0_reference_errors.py --root data/sen2neon --out results/p6_e0_mtf+0.05 --mtf-delta 0.05 --every 5
  python scripts/p6_e0_reference_errors.py --root data/sen2neon --out results/p6_e0 --summary-only
"""
from __future__ import annotations

import os

# one BLAS/FFT thread per worker process: with the default, 24 workers drove the load of the shared
# 32-core server to ~390 (must be set before numpy/scipy/torch are imported)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import rasterio
import torch
import yaml
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import GUIDE_BANDS, TARGET_BANDS  # noqa: E402
from s2sr.data.sen2neon import (GUIDE_IDX, HR_NODATA, LR_NODATA, SCALE, TARGET_IDX,  # noqa: E402
                                reference_valid)
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.registration import (block_shifts, common_noise_decomposition, linear_fit,  # noqa: E402
                               match_blocks, phase_corr, shift_field_residual, subpixel_shift)

DT_BINS = [-0.5, 3, 7, 15, 40]
DT_LABELS = ["0-3", "4-7", "8-15", "16-40"]
LOWPASS_M = 60.0       # Gaussian sigma (m) of the low-pass used for gain/offset   (R1-1.M1, R1-3.6)
HIGHPASS_M = 20.0      # Gaussian sigma (m) removed to form the high-pass detail
BLOCK_M = 640.0        # block size (m) for block-wise shifts on both band groups (R1-1.C5)
MAX_RESAMPLED = 0.01   # share of non-replicated 2x2 blocks above which a tile counts as resampled
H0_THRESHOLD = 0.25    # px (10 m) residual after one translation


def _read(path: Path) -> np.ndarray:
    with rasterio.open(path) as d:
        return d.read()


def _simulate(hr: np.ndarray, ok: np.ndarray, sigmas: list[float], ratio: int):
    """Degrade 2.5 m reference bands to a coarse grid; nodata is filled before blurring and the
    coarse mask is eroded by the kernel reach so no output pixel sees filled values."""
    hr = hr.astype(np.float32) / SCALE
    for b in range(hr.shape[0]):
        if (~ok).any() and ok.any():
            hr[b][~ok] = hr[b][ok].mean()
    with torch.no_grad():
        sim = degrade(torch.from_numpy(hr)[None], sigmas, ratio=ratio)[0].numpy().astype(np.float64)
    return sim, reference_valid(ok, sigmas, ratio)


def native_20m(stack: np.ndarray):
    """(C, H, W) 2x2-replicated stack -> native (C, H/2, W/2) values, validity mask and the share of
    2x2 blocks that are not exact replicas. A few such blocks occur at tile edges and are masked;
    a large share means the 20 m bands were resampled and have no native values."""
    c, h, w = stack.shape
    blk = stack.reshape(c, h // 2, 2, w // 2, 2)
    bad = (blk.max(axis=(2, 4)) != blk.min(axis=(2, 4))).any(axis=0)
    return stack[:, ::2, ::2], ~bad, float(bad.mean())


def _smooth(img, sigma, valid):
    """Gaussian smoothing normalised by the valid mask, so invalid pixels (e.g. LR nodata 65535,
    i.e. 6.55 in reflectance) contribute nothing to their valid neighbours (R2-3.N1)."""
    w = valid.astype(np.float64)
    num = ndimage.gaussian_filter(np.where(valid, img, 0.0), sigma)
    den = ndimage.gaussian_filter(w, sigma)
    return np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)


def _lp(img, grid_m, valid):
    return _smooth(img, LOWPASS_M / grid_m, valid)


def _hp(img, grid_m, valid):
    return img - _smooth(img, HIGHPASS_M / grid_m, valid)


def _rmse(a, b):
    return float(np.sqrt(((a - b) ** 2).mean()))


def _radiometry(x: np.ndarray, r_al: np.ndarray, m: np.ndarray, grid_m: int, valid: np.ndarray) -> dict:
    """Compare one band of the observation x with the aligned reference r_al on mask m; `valid` (a
    superset of m) is the mask used by the normalised filters."""
    xk, rk = x[m], r_al[m]
    lxf, lrf = _lp(x, grid_m, valid), _lp(r_al, grid_m, valid)
    lx, lr_ = lxf[m], lrf[m]
    ols = linear_fit(lx, lr_)
    # OLS gain is attenuated when x has variation not shared with r; the standard-deviation ratio
    # (geometric-mean regression) is symmetric. Report residuals for both (R1-3.3).
    a_sd = float(lr_.std() / max(lx.std(), 1e-12))
    c_sd = float(lr_.mean() - a_sd * lx.mean())
    # cross-fitted OLS residual: fit on one half of the tile, evaluate on the other (R1-1.C3)
    cols = np.broadcast_to(np.arange(m.shape[1])[None, :], m.shape)
    half = cols < m.shape[1] // 2
    se, n = 0.0, 0
    for fit_half in (half, ~half):
        mf, me = m & fit_half, m & ~fit_half
        if mf.sum() < 50 or me.sum() < 50:
            se, n = np.nan, 1
            break
        f = linear_fit(lxf[mf], lrf[mf])
        se += float(((r_al[me] - (f["a"] * x[me] + f["c"])) ** 2).sum())
        n += int(me.sum())
    hx, hr = _hp(x, grid_m, valid)[m], _hp(r_al, grid_m, valid)[m]
    return dict(n_px=int(m.sum()), mean_x=float(xk.mean()), mean_r=float(rk.mean()),
                rmse_al=_rmse(rk, xk),
                rmse_lin_ols=_rmse(rk, ols["a"] * xk + ols["c"]),
                rmse_lin_sym=_rmse(rk, a_sd * xk + c_sd),
                rmse_lin_cv=float(np.sqrt(se / n)),
                a_ols=ols["a"], c_ols=ols["c"], a_sd=a_sd, c_sd=c_sd, r2_lowpass=ols["r2"],
                # high-pass gain and agreement (R1-1.C2): is detail amplitude in r lower than in x?
                hp_gain_sd=float(hr.std() / max(hx.std(), 1e-12)),
                hp_gain_ols=float(np.cov(hx, hr)[0, 1] / max(hx.var(), 1e-24)),
                hp_corr=float(np.corrcoef(hx, hr)[0, 1]),
                # harmonisation check: a per-tile regression of x on r would give sd ratio ~ R
                R_full=float(np.corrcoef(xk, rk)[0, 1]),
                sd_ratio_full=float(rk.std() / max(xk.std(), 1e-12)))


def _stat_mask(m: np.ndarray, grid_m: int) -> np.ndarray:
    """Pixels whose whole low-pass footprint lies inside m: statistics use only these, while the normalised
    filters use all of m, so no retained pixel sees a truncated kernel (R3-3.NEW-1). scipy truncates the
    Gaussian at 4 sigma, hence a 4-sigma erosion."""
    return ndimage.binary_erosion(m, iterations=math.ceil(4 * LOWPASS_M / grid_m), border_value=0)


def _band_rows(tile: str, bands, x, r, valid, d10, grid_m: int, own_shifts=None, min_peak=0.1):
    """Per-band radiometry after the tile shift d10 (10 m px) and, if given, after each band's own shift.

    Global and own-shift metrics use separate masks, so an unreliable own shift (peak < min_peak) cannot
    shrink the mask of the global metrics (R3-3.NEW-4)."""
    f = 10.0 / grid_m

    def masks(shifts):
        # pixels within the shift reach of invalid pixels hold mirrored or filled values after the shift
        m = ndimage.binary_erosion(valid, iterations=math.ceil(max(shifts) * f) + 3, border_value=0)
        return m, _stat_mask(m, grid_m)

    m, ms = masks([abs(v) for v in d10])
    if ms.sum() < 200:
        return []
    r_al = subpixel_shift(r, -d10[0] * f, -d10[1] * f)
    rows = []
    for k, band in enumerate(bands):
        row = dict(tile=tile, band=band, grid_m=grid_m,
                   rmse_raw=_rmse(r[k][ms], x[k][ms]),
                   hp_corr_raw=float(np.corrcoef(_hp(x[k], grid_m, m)[ms], _hp(r[k], grid_m, m)[ms])[0, 1]))
        row.update(_radiometry(x[k], r_al[k], ms, grid_m, m))
        if own_shifts is not None:
            dy, dx, pk = own_shifts[k]
            row.update(dy=dy, dx=dx, peak=pk)
            if pk >= min_peak and np.isfinite(dy):  # alignment with the band's own shift
                mo, mso = masks([abs(v) for v in d10] + [abs(dy), abs(dx)])
                if mso.sum() >= 200:
                    rb = subpixel_shift(r[k], -dy * f, -dx * f)
                    ref = _radiometry(x[k], r_al[k], mso, grid_m, mo)  # paired: same pixels as `own`
                    own = _radiometry(x[k], rb, mso, grid_m, mo)
                    row.update(rmse_al_pair=ref["rmse_al"], hp_corr_pair=ref["hp_corr"],
                               rmse_al_own=own["rmse_al"], rmse_lin_ols_own=own["rmse_lin_ols"],
                               hp_corr_own=own["hp_corr"])
        rows.append(row)
    return rows


def _rms(d: np.ndarray) -> float:
    return float(np.sqrt((d ** 2).sum(axis=1).mean())) if len(d) else float("nan")


def _block_rms_split(blocks: list[dict], shape, block: int, scale: float) -> tuple[float, float]:
    """rms of block shifts (around 0) for edge blocks and interior blocks separately (R3-3.NEW-3)."""
    h, w = shape
    edge, inner = [], []
    for b in blocks:
        e = min(b["cy"], b["cx"], h - b["cy"], w - b["cx"]) <= block / 2 + 1e-6
        (edge if e else inner).append([b["dy"] * scale, b["dx"] * scale])
    return _rms(np.array(edge).reshape(-1, 2)), _rms(np.array(inner).reshape(-1, 2))


def _residual_after_tile_shift(xg, rg, v10, xt, rt, v20, dg, dt, upsample, min_peak) -> dict:
    """Block shifts left after the shifts the correction will use: guide bands after the guide tile shift
    dg, target bands after the joint target-band shift dt (R3-1, R3-2); target bands also after dg as a
    sensitivity. rms around zero, in 10 m px.

    Target blocks are measured on two independent band groups (red edge B5-B8A, SWIR B11-B12) so that
    estimation noise can be separated from a field shared by both groups. H0' thresholds the unbiased
    `common2` (s2sr.registration.common_noise_decomposition)."""
    bt_px, bg_px = int(BLOCK_M / 20), int(BLOCK_M / 10)

    def aligned(r, v, d, grid_f):
        b = math.ceil(max(abs(d[0]), abs(d[1])) * grid_f) + 3
        # erode only around interior nodata (border_value=1): eroding from the tile edge would drop every
        # edge block below the 95% validity threshold; edges are already down-weighted by the window
        return subpixel_shift(r, -d[0] * grid_f, -d[1] * grid_f), ndimage.binary_erosion(v, iterations=b,
                                                                                            border_value=1)
    keep = lambda bl: [b for b in bl if b["peak"] >= min_peak]  # noqa: E731
    rga, vg = aligned(rg, v10, dg, 1.0)
    bg = keep(block_shifts(xg, rga, vg, block=bg_px, upsample=upsample))
    out = dict(n_blocks_gs=len(bg), rms_tileshift_g=_rms(match_blocks(bg)[0]))
    out["rms_edge_g"], out["rms_inner_g"] = _block_rms_split(bg, xg.shape[-2:], bg_px, 1.0)
    for tag, d in (("", dt), ("_gshift", dg)):
        rta, vt = aligned(rt, v20, d, 0.5)
        bt = keep(block_shifts(xt, rta, vt, block=bt_px, upsample=upsample))
        ba = keep(block_shifts(xt[:4], rta[:4], vt, block=bt_px, upsample=upsample))
        bb = keep(block_shifts(xt[4:], rta[4:], vt, block=bt_px, upsample=upsample))
        out[f"n_blocks_ts{tag}"] = len(bt)
        out[f"rms_tileshift_t{tag}"] = _rms(2.0 * match_blocks(bt)[0])
        if tag == "":
            out["rms_group_rededge"] = _rms(2.0 * match_blocks(ba)[0])
            out["rms_group_swir"] = _rms(2.0 * match_blocks(bb)[0])
            out["rms_edge_t"], out["rms_inner_t"] = _block_rms_split(bt, xt.shape[-2:], bt_px, 2.0)
        fg, fa, fb = match_blocks(bg, ba, bb, scale=2.0)
        if len(fa) >= 6:
            ab = common_noise_decomposition(fa, fb, centre="none")
            tg = common_noise_decomposition(fg, (fa + fb) / 2, centre="none")
            out.update({f"n_blocks_dec{tag}": ab["n"], f"common2_t{tag}": ab["common2"],
                        f"common_t{tag}": ab["common"], f"noise_t{tag}": ab["noise"],
                        f"common2_tg{tag}": tg["common2"], f"common_tg{tag}": tg["common"]})
        else:
            out.update({f"n_blocks_dec{tag}": len(fa), f"common2_t{tag}": np.nan, f"common_t{tag}": np.nan,
                        f"noise_t{tag}": np.nan, f"common2_tg{tag}": np.nan, f"common_tg{tag}": np.nan})
    return out


def analyse_tile(root: Path, row, sig_g, sig_t, upsample: int, min_peak: float = 0.1):
    lr = _read(root / row.lr)
    hr = _read(root / row.hr_2_5m_path)
    lr_ok = np.all(lr != LR_NODATA, axis=0)
    tile = dict(tile=row.id, site=row.id.split("_")[1], flight=row.id.split("__")[0],
                year=int(row.id[:4]), dt_days=abs(float(row.temporal_difference_days)),
                land_cover=row.land_cover_superclass,
                cloud_score=float(getattr(row, "cloud_score_plus_cdf", np.nan)))

    # 20 m target bands; tiles whose 20 m bands were resampled have no native values (R1-1.m2)
    xt, rep_ok, bad_frac = native_20m(lr[TARGET_IDX])
    tile["nonreplicated_frac"] = bad_frac
    if bad_frac > MAX_RESAMPLED:
        tile["status"] = "skipped: 20 m bands resampled"
        return tile, []
    xt = xt.astype(np.float64) / SCALE
    rt, ok_t = _simulate(hr[TARGET_IDX], np.all(hr[TARGET_IDX] != HR_NODATA, axis=0), sig_t, 8)
    v20 = lr_ok.reshape(lr_ok.shape[0] // 2, 2, -1, 2).all(axis=(1, 3)) & ok_t & rep_ok

    # 10 m guide bands
    xg = lr[GUIDE_IDX].astype(np.float64) / SCALE
    rg, ok_g = _simulate(hr[GUIDE_IDX], np.all(hr[GUIDE_IDX] != HR_NODATA, axis=0), sig_g, 4)
    v10 = lr_ok & ok_g
    tile["valid_frac_10m"] = float(v10.mean())
    if v10.mean() < 0.5:
        tile["status"] = "skipped: <50% valid"
        return tile, []

    gy, gx, gp = phase_corr(xg, rg, v10, upsample)
    if not (np.isfinite(gy) and np.isfinite(gx)):  # R3-1: flat or fully invalid guide bands
        tile["status"] = "skipped: no global shift"
        return tile, []
    ty, tx, tp = phase_corr(xt, rt, v20, upsample)
    if not (np.isfinite(ty) and np.isfinite(tx)):  # fall back to the guide shift (20 m px)
        ty, tx, tp = gy / 2, gx / 2, 0.0
    res_g = shift_field_residual(block_shifts(xg, rg, v10, block=int(BLOCK_M / 10), upsample=upsample),
                                 min_peak=min_peak)
    res_t = shift_field_residual(block_shifts(xt, rt, v20, block=int(BLOCK_M / 20), upsample=upsample),
                                 min_peak=min_peak)
    tile.update(status="ok", g_dy=gy, g_dx=gx, g_mag=math.hypot(gy, gx), g_peak=gp,
                t_dy=2 * ty, t_dx=2 * tx, t_mag=2 * math.hypot(ty, tx), t_peak=tp,
                gt_diff=math.hypot(gy - 2 * ty, gx - 2 * tx),
                n_blocks_g=res_g["n_blocks"], rms_translation_g=res_g["rms_translation"],
                rms_affine_g=res_g["rms_affine"],
                # target blocks are in 20 m px: convert to 10 m px
                n_blocks_t=res_t["n_blocks"], rms_translation_t=2 * res_t["rms_translation"],
                rms_affine_t=2 * res_t["rms_affine"])
    tile.update(_residual_after_tile_shift(xg, rg, v10, xt, rt, v20, (gy, gx), (2 * ty, 2 * tx), upsample,
                                           min_peak))
    own = []
    for k in range(len(TARGET_BANDS)):
        dy, dx, pk = phase_corr(xt[k], rt[k], v20, upsample)
        own.append((2 * dy, 2 * dx, pk))
    # each band group is aligned with the shift the correction will use for it (§4.1 of the plan)
    bands = (_band_rows(row.id, GUIDE_BANDS, xg, rg, v10, (gy, gx), 10)
             + _band_rows(row.id, TARGET_BANDS, xt, rt, v20, (2 * ty, 2 * tx), 20, own, min_peak))
    return tile, bands


def _init_worker():
    torch.set_num_threads(1)


def _job(args):
    row = args[1]
    try:
        return analyse_tile(*args)
    except Exception as e:  # keep the run going; the failure is recorded per tile
        return dict(tile=row.id, site=row.id.split("_")[1], status=f"error: {type(e).__name__}: {e}"), []


# ------------------------------------------------------------------------------------- summary
def _cluster_ci(df: pd.DataFrame, stat, n_boot: int = 2000, seed: int = 0):
    """Percentile CI of stat(df) under resampling of sites with replacement."""
    rng = np.random.default_rng(seed)
    groups = [g for _, g in df.groupby("site")]
    vals = [stat(pd.concat([groups[i] for i in rng.integers(0, len(groups), len(groups))]))
            for _ in range(n_boot)]
    return np.nanpercentile(vals, [2.5, 97.5])


def _nested_ss(df: pd.DataFrame, col: str) -> dict:
    """Sum-of-squares decomposition of col into site, flight-within-site and tile-within-flight."""
    y, gm = df[col], df[col].mean()
    site_m = df.groupby("site")[col].transform("mean")
    flight_m = df.groupby("flight")[col].transform("mean")
    ss = dict(site=((site_m - gm) ** 2).sum(), flight_in_site=((flight_m - site_m) ** 2).sum(),
              tile_in_flight=((y - flight_m) ** 2).sum())
    tot = sum(ss.values())
    return {k: v / tot for k, v in ss.items()}


def _h0_lines(ok: pd.DataFrame, col: str, label: str, n_col: str | None = None,
              thr: float = H0_THRESHOLD) -> list[str]:
    """Share of tiles with col <= thr: complete-case with a site-bootstrap CI, and the same
    over all reliable tiles counting a missing value as failing (R2-1: CI over all sites)."""
    fin = ok[np.isfinite(ok[col])]
    rate = (fin[col] <= thr).mean()
    lo, hi = _cluster_ci(fin, lambda d: (d[col] <= thr).mean())
    fail_missing = lambda d: (d[col] <= thr).sum() / len(d)  # noqa: E731  NaN <= x is False
    worst = fail_missing(ok)
    wlo, whi = _cluster_ci(ok, fail_missing)
    nb = f"; blocks per tile (median, same tiles): {fin[n_col].median():.0f}" if n_col else ""
    return [f"- {label}: available for {len(fin)} of {len(ok)} tiles{nb}; median {fin[col].median():.3f}, "
            f"q90 {fin[col].quantile(0.9):.3f}; share <= {thr:g}: **{rate:.1%}** "
            f"(site-bootstrap CI {lo:.1%}-{hi:.1%}); missing counted as failing: {worst:.1%} "
            f"(CI {wlo:.1%}-{whi:.1%})"]


def summarise(tiles: pd.DataFrame, bands: pd.DataFrame, out: Path, min_peak: float) -> str:
    ok = tiles[(tiles.status == "ok") & (tiles.g_peak >= min_peak)].copy()
    ok["dt_bin"] = pd.cut(ok.dt_days, DT_BINS, labels=DT_LABELS)
    q = [0.1, 0.25, 0.5, 0.75, 0.9]
    qcols = [f"q{int(100 * c)}" for c in q]
    lines = [f"# E0 summary: {len(ok)} reliable tiles of {len(tiles)} "
             f"(status ok: {(tiles.status == 'ok').sum()}, of which g_peak >= {min_peak}: {len(ok)}), "
             f"{ok.site.nunique()} sites, {ok.flight.nunique()} flights among reliable tiles "
             f"({tiles.flight.nunique() if 'flight' in tiles else 'n/a'} in all tiles)", "",
             "Filter scales: low-pass sigma %.0f m, high-pass removes sigma %.0f m; blocks %.0f m." %
             (LOWPASS_M, HIGHPASS_M, BLOCK_M), "",
             "## Tile status", "", tiles.status.value_counts().to_markdown(), ""]

    geo = ok[["g_mag", "t_mag", "gt_diff", "rms_translation_g", "rms_affine_g",
              "rms_translation_t", "rms_affine_t"]].quantile(q).T
    geo.columns = qcols
    lines += ["## Geometry (10 m pixels)", "",
              "g/t: guide (10 m) / target (20 m) bands; gt_diff: difference of their global shifts;",
              "rms_translation / rms_affine: block shifts left after one translation / an affine field.",
              "", geo.round(3).to_markdown(), ""]
    for thr in (0.25, 0.5, 1.0):
        lines.append(f"- tiles with |global shift| > {thr} px: {(ok.g_mag > thr).mean():.1%}")
    lines += ["", "### H0 (one translation per tile leaves <= %.2f px)" % H0_THRESHOLD, "",
              "Pre-specified measure (rms of block shifts around their median):", ""]
    lines += _h0_lines(ok, "rms_translation_g", "guide bands (10 m)", "n_blocks_g")
    lines += _h0_lines(ok, "rms_translation_t", "target bands (20 m)", "n_blocks_t")
    H2 = H0_THRESHOLD ** 2
    lines += ["", "After the shifts the correction uses (guide bands: guide tile shift; target bands: joint "
              "target-band shift); rms of block shifts around 0:", ""]
    lines += _h0_lines(ok, "rms_tileshift_g", "guide bands", "n_blocks_gs")
    lines += _h0_lines(ok, "rms_tileshift_t", "target bands (all six, raw; includes estimation noise)",
                       "n_blocks_ts")
    lines += _h0_lines(ok, "rms_group_rededge", "target group red edge B5-B8A (raw)")
    lines += _h0_lines(ok, "rms_group_swir", "target group SWIR B11-B12 (raw)")
    lines += ["", "H0' (post hoc): field common to the two target groups, unbiased squared estimate common2 = "
              "mean(a . b) thresholded at 0.25^2 (the clipped sqrt `common` is biased up by ~0.04 px at 14 blocks):", ""]
    lines += _h0_lines(ok, "common2_t", "H0' common target field (common2 <= 0.0625)", "n_blocks_dec", thr=H2)
    lines += _h0_lines(ok, "common2_tg", "field common to target and guide (common2 <= 0.0625)", "n_blocks_dec",
                       thr=H2)
    lines += _h0_lines(ok, "noise_t", "estimation noise of one target group (rms)", "n_blocks_dec")
    lines += ["", "Sensitivity: target bands aligned with the guide shift instead:", ""]
    lines += _h0_lines(ok, "rms_tileshift_t_gshift", "target bands (raw)")
    lines += _h0_lines(ok, "common2_t_gshift", "H0' measure", thr=H2)
    top = ok.cloud_score > ok.cloud_score.quantile(0.75)
    lines += ["", f"Sensitivity, excluding the top cloud-score quartile ({top.sum()} tiles):", ""]
    lines += _h0_lines(ok[~top], "rms_tileshift_g", "guide bands")
    lines += _h0_lines(ok[~top], "common2_t", "H0' measure", thr=H2)
    geo_dt = ok.groupby("dt_bin", observed=False).agg(
        n=("tile", "size"), rms_tileshift_g=("rms_tileshift_g", "median"),
        rms_tileshift_t=("rms_tileshift_t", "median"), common_t=("common_t", "median"),
        H0prime_share=("common2_t", lambda v: (v.dropna() <= H2).mean()),
        common_tg=("common_tg", "median"))
    lines += ["", "Geometry by |dt| (days): if residual fields were surface change they would grow with |dt| "
              "(R3-1):", "", geo_dt.round(3).to_markdown(), ""]
    ei = ok[["rms_edge_g", "rms_inner_g", "rms_edge_t", "rms_inner_t"]].median()
    lines += ["Edge vs interior blocks after the tile shift (medians, 10 m px; reflect-padded tile edges, "
              "R3-3.NEW-3): " + ", ".join(f"{k} {v:.3f}" for k, v in ei.items()), ""]
    pk = ok[["g_peak", "t_peak"]].quantile(q).T
    pk.columns = qcols
    lines += ["", "Phase-correlation peaks:", "", pk.round(3).to_markdown(), ""]
    ns = pd.DataFrame({c: _nested_ss(ok, c) for c in ["g_dy", "g_dx"]}).T
    lines += ["### Where shifts vary (share of sum of squares; R1-1.M2)", "", ns.round(3).to_markdown(), ""]

    clean = ok[ok.cloud_score <= ok.cloud_score.quantile(0.75)]
    b_all = bands[bands.tile.isin(ok.tile)].merge(ok[["tile", "dt_bin", "land_cover"]], on="tile")
    b = b_all[b_all.tile.isin(clean.tile)]  # R3-1: main radiometric tables on the cloud-clean subset
    order = GUIDE_BANDS + TARGET_BANDS
    cols = ["rmse_raw", "rmse_al", "rmse_lin_ols", "rmse_lin_sym", "rmse_lin_cv", "a_ols", "a_sd",
            "r2_lowpass", "hp_gain_sd", "hp_gain_ols", "hp_corr_raw", "hp_corr"]
    per_band = b.groupby("band")[cols].median().reindex(order)
    lines += [f"## Radiometry per band, cloud-clean subset ({len(clean)} tiles, cloud score <= q75; medians; "
              "reflectance units)", "",
              "Statistics use only pixels whose whole 4-sigma low-pass footprint is valid. hp_gain is an "
              "*apparent* gain: S2 noise and aliasing in x lower it; not yet corrected (pre-E2 gate).", "",
              "rmse_al: after global shift; rmse_lin_*: after shift + gain/offset fitted on the low-pass "
              "(ols, symmetric, or cross-fitted on the other tile half); hp_gain: std / OLS ratio of "
              "high-pass detail r vs x.", "", per_band.round(4).to_markdown(), ""]
    agg = pd.DataFrame({
        "align_pct (ratio of medians)": 100 * (1 - per_band.rmse_al / per_band.rmse_raw),
        "lin_cv_pct (ratio of medians)": 100 * (1 - per_band.rmse_lin_cv / per_band.rmse_al)})
    lines += ["RMSE reduction computed from the medians above (R1-3.7):", "", agg.round(1).to_markdown(), ""]
    all_t = b_all.groupby("band")[["rmse_lin_cv", "a_sd", "hp_gain_sd", "hp_corr"]].median().reindex(order)
    lines += [f"All reliable tiles incl. cloudy ones ({len(ok)}; sensitivity):", "", all_t.round(4).to_markdown(), ""]

    harm = b.groupby("band").apply(lambda g: pd.Series({
        "median R": g.R_full.median(), "median sd ratio": g.sd_ratio_full.median(),
        "median hp_gain_sd": g.hp_gain_sd.median(),
        "corr(sd ratio, R) over tiles": np.corrcoef(g.sd_ratio_full, g.R_full)[0, 1],
        "median |mean_r - mean_x|": (g.mean_r - g.mean_x).abs().median()}), include_groups=False).reindex(order)
    lines += ["## Harmonisation check", "",
              "A per-tile regression of S2 on the reference would match tile means exactly and shrink the "
              "reference's standard deviation to about R times that of S2 at the fitted scale.", "",
              harm.round(4).to_markdown(), ""]

    tb = b[(b.grid_m == 20)]
    if "rmse_al_own" in tb:
        tb = tb[tb.rmse_al_own.notna()].merge(ok[["tile", "gt_diff"]], on="tile")  # paired rows only (R2-3.N2)
        # baseline = joint target-band shift, on exactly the same pixels as the own-shift metrics
        tb = tb.assign(d_rmse=tb.rmse_al_pair - tb.rmse_al_own, d_hpcorr=tb.hp_corr_own - tb.hp_corr_pair,
                       mag=np.hypot(tb.dy, tb.dx))
        own = tb.groupby("band").agg(
            n=("tile", "size"), own_shift_q50=("mag", "median"),
            own_shift_q90=("mag", lambda v: v.quantile(0.9)),
            d_rmse_median=("d_rmse", "median"), d_rmse_q90=("d_rmse", lambda v: v.quantile(0.9)),
            share_rmse_gain_gt_0_001=("d_rmse", lambda v: (v > 0.001).mean()),
            d_hpcorr_median=("d_hpcorr", "median"), d_hpcorr_q90=("d_hpcorr", lambda v: v.quantile(0.9)),
            share_hpcorr_gain_gt_0_01=("d_hpcorr", lambda v: (v > 0.01).mean())).reindex(TARGET_BANDS)
        lines += ["## Target bands: paired gain of each band's own shift over the joint target-band shift", "",
                  "d_rmse = rmse(tile shift) - rmse(own shift), d_hpcorr = hp_corr(own) - hp_corr(tile); "
                  "positive = own shift better. Rows: tiles where the band's own peak >= min_peak.", "",
                  own.round(4).to_markdown(), ""]
        q4 = tb.groupby(pd.qcut(tb.gt_diff, 4, duplicates="drop"), observed=True).agg(
            n=("tile", "size"), d_rmse_median=("d_rmse", "median"), d_hpcorr_median=("d_hpcorr", "median"),
            share_hpcorr_gain_gt_0_01=("d_hpcorr", lambda v: (v > 0.01).mean()))
        lines += ["By quartile of the guide-target shift difference gt_diff (10 m px):", "",
                  q4.round(4).to_markdown(), ""]

    by_dt = b.pivot_table(index="band", columns="dt_bin", values="rmse_lin_cv", aggfunc="median",
                          observed=False).reindex(order)
    n_dt = ok.dt_bin.value_counts().reindex(DT_LABELS)
    lines += ["## Cross-fitted residual RMSE by |dt| (days)", "",
              "tiles per bin: " + ", ".join(f"{k}: {v}" for k, v in n_dt.items()), "",
              by_dt.round(4).to_markdown(), ""]
    by_lc = b.pivot_table(index="band", columns="land_cover", values="rmse_lin_cv", aggfunc="median",
                          observed=False).reindex(order)
    n_lc = ok.land_cover.value_counts()
    lines += ["## Same, by land-cover superclass (medians; small classes are indicative only)", "",
              "tiles per class: " + ", ".join(f"{k}: {v}" for k, v in n_lc.items()), "",
              by_lc.round(4).to_markdown(), ""]
    gain_lc = b[b.grid_m == 20].pivot_table(index="band", columns="land_cover", values="hp_gain_sd",
                                            aggfunc="median", observed=False).reindex(TARGET_BANDS)
    lines += ["## High-pass gain (sd ratio) of target bands by land cover", "", gain_lc.round(3).to_markdown(), ""]
    cs = ok.groupby(pd.qcut(ok.cloud_score, 4, duplicates="drop"), observed=True).agg(
        n=("tile", "size"), g_peak=("g_peak", "median"), rms_tileshift_g=("rms_tileshift_g", "median"))
    bt20 = b[b.grid_m == 20].merge(ok[["tile", "cloud_score"]], on="tile")
    rad = bt20.groupby(pd.qcut(bt20.cloud_score, 4, duplicates="drop"), observed=True)[
        ["rmse_lin_cv", "a_sd", "hp_gain_sd", "hp_corr"]].median()
    cs = cs.join(rad)
    lines += ["## By tile cloud score (cloud_score_plus_cdf quartiles; SEN2NEON has no per-pixel cloud "
              "mask). Radiometry columns: medians over the six target bands and tiles.", "",
              cs.round(4).to_markdown(), ""]
    text = "\n".join(lines)
    (out / "summary.md").write_text(text)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="data/sen2neon", help="folder with metadata.csv and the tiles")
    ap.add_argument("--out", default="results/p6_e0")
    ap.add_argument("--max-dt", type=float, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--every", type=int, default=1, help="use every n-th tile (sensitivity runs)")
    ap.add_argument("--mtf-delta", type=float, default=0.0, help="add to every MTF value (sensitivity)")
    ap.add_argument("--upsample", type=int, default=100,
                    help="phase-correlation refinement: resolution is 1/upsample px of each grid; 20 "
                         "quantised 20 m-band shifts to 0.1 px of 10 m (found by the end-to-end test)")
    ap.add_argument("--min-peak", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    if not a.summary_only:
        meta = pd.read_csv(root / "metadata.csv")
        meta = meta[[(root / p).exists() and (root / q).exists()
                     for p, q in zip(meta.lr, meta.hr_2_5m_path)]]
        if a.max_dt is not None:
            meta = meta[meta.temporal_difference_days.abs() <= a.max_dt]
        meta = meta.iloc[::a.every]
        if a.limit:
            meta = meta.head(a.limit)
        mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
        mtf = {k: v + a.mtf_delta for k, v in mtf.items()}
        sig_g = band_sigmas(mtf, GUIDE_BANDS, ratio=4)
        sig_t = band_sigmas(mtf, TARGET_BANDS, ratio=8)
        jobs = [(root, SimpleNamespace(**r), sig_g, sig_t, a.upsample, a.min_peak)
                for r in meta.to_dict("records")]
        tiles, bands = [], []
        with ProcessPoolExecutor(a.workers, initializer=_init_worker) as pool:
            for i, (t, bs) in enumerate(pool.map(_job, jobs, chunksize=4)):
                tiles.append(t)
                bands += bs
                if t.get("status") == "ok":
                    print(f"[{i + 1}/{len(meta)}] {t['tile']}: shift {t['g_dy']:+.2f},{t['g_dx']:+.2f} px "
                          f"(peak {t['g_peak']:.2f}), residual {t['rms_translation_g']:.2f} px", flush=True)
                else:
                    print(f"[{i + 1}/{len(meta)}] {t['tile']}: {t.get('status')}", flush=True)
        pd.DataFrame(tiles).to_csv(out / "tiles.csv", index=False)
        pd.DataFrame(bands).to_csv(out / "bands.csv", index=False)
    print(summarise(pd.read_csv(out / "tiles.csv"), pd.read_csv(out / "bands.csv"), out, a.min_peak))


if __name__ == "__main__":
    main()
