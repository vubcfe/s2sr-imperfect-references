"""P6 endpoints (plan §6, §8): high-pass PSNR, high-pass correlation and edge-position error.

High-pass: HP(z) = z - G_sigma * z with sigma = 20 m (2 px on the 10 m grid), Gaussian smoothing normalised by
the valid mask; statistics only on pixels whose whole 4-sigma footprint is valid (as in E0). All functions take
(C, H, W) numpy arrays of one tile and return one value per band.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

HP_SIGMA_M = 20.0
GRID_M = 10.0


def _smooth(img, sigma, valid):
    num = ndimage.gaussian_filter(np.where(valid, img, 0.0), sigma)
    den = ndimage.gaussian_filter(valid.astype(np.float64), sigma)
    return np.where(den > 1e-6, num / np.maximum(den, 1e-6), 0.0)


def stat_mask(valid: np.ndarray, sigma_px: float = HP_SIGMA_M / GRID_M, border: int = 8) -> np.ndarray:
    m = ndimage.binary_erosion(valid, iterations=int(np.ceil(4 * sigma_px)), border_value=0)
    if border:
        m[:border], m[-border:], m[:, :border], m[:, -border:] = False, False, False, False
    return m


def highpass(z: np.ndarray, valid: np.ndarray, sigma_px: float = HP_SIGMA_M / GRID_M) -> np.ndarray:
    z = np.asarray(z, dtype=np.float64)
    return np.stack([z[k] - _smooth(z[k], sigma_px, valid) for k in range(z.shape[0])])


def hp_psnr(pred: np.ndarray, ref: np.ndarray, valid: np.ndarray, peak: float = 1.0) -> np.ndarray:
    m = stat_mask(valid)
    hp, hr = highpass(pred, valid), highpass(ref, valid)
    mse = ((hp - hr) ** 2)[:, m].mean(axis=1)
    return 10 * np.log10(peak ** 2 / np.maximum(mse, 1e-20))


def hp_corr(pred: np.ndarray, ref: np.ndarray, valid: np.ndarray) -> np.ndarray:
    m = stat_mask(valid)
    hp, hr = highpass(pred, valid)[:, m], highpass(ref, valid)[:, m]
    hp = hp - hp.mean(axis=1, keepdims=True)
    hr = hr - hr.mean(axis=1, keepdims=True)
    return (hp * hr).sum(1) / np.maximum(np.sqrt((hp ** 2).sum(1) * (hr ** 2).sum(1)), 1e-20)


def _edges(band: np.ndarray, m: np.ndarray, q: float) -> np.ndarray:
    """Edge map: local maxima of the Sobel gradient magnitude along x or y, above the q-quantile."""
    gy, gx = ndimage.sobel(band, 0), ndimage.sobel(band, 1)
    mag = np.hypot(gx, gy)
    thr = np.quantile(mag[m], q) if m.any() else np.inf
    nms = (mag >= ndimage.maximum_filter(mag, size=3)) & (mag > thr)
    return nms & m


def edge_error(pred: np.ndarray, ref: np.ndarray, valid: np.ndarray, q: float = 0.9,
               max_dist: float = 3.0) -> np.ndarray:
    """Symmetric mean distance (px) between edges of pred and edges of ref, distances capped at max_dist so
    that unmatched edges count as a fixed penalty rather than dominating. Lower is better."""
    m = stat_mask(valid)
    out = []
    for k in range(pred.shape[0]):
        ep, er = _edges(pred[k], m, q), _edges(ref[k], m, q)
        if not ep.any() or not er.any():
            out.append(np.nan)
            continue
        d_to_r = ndimage.distance_transform_edt(~er)
        d_to_p = ndimage.distance_transform_edt(~ep)
        out.append(0.5 * (np.minimum(d_to_r[ep], max_dist).mean() + np.minimum(d_to_p[er], max_dist).mean()))
    return np.array(out)


def hp_gain(pred: np.ndarray, ref: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Amplitude of the predicted detail relative to the reference: sd(HP(pred)) / sd(HP(ref)) per band
    (1 = correct contrast, < 1 = too smooth). HP-corr is blind to this, so E1-R reports both."""
    m = stat_mask(valid)
    hp, hr = highpass(pred, valid)[:, m], highpass(ref, valid)[:, m]
    return hp.std(1) / np.maximum(hr.std(1), 1e-20)
