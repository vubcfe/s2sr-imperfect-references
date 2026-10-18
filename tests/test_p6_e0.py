"""End-to-end and component tests for P6 experiment E0 (added after review round 1, R1-3.5)."""
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import rasterio
import torch

from s2sr.data.sen2neon import GUIDE_IDX, TARGET_IDX
from s2sr.degrade import band_sigmas, degrade
from s2sr.registration import _upsampled_dft, fourier_shift, phase_corr, subpixel_shift

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("e0", ROOT / "scripts" / "p6_e0_reference_errors.py")
E0 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E0)

MTF = {b: v for b, v in zip(["B02", "B03", "B04", "B08", "B05", "B06", "B07", "B8A", "B11", "B12"],
                            [0.24] * 4 + [0.29] * 4 + [0.195] * 2)}


def _scene(n, c, seed):
    """Periodic scene-like texture (1/f^2 spectrum) on an n x n grid, c bands."""
    rng = np.random.default_rng(seed)
    fy, fx = np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :]
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((c, n, n))) / (1 + (fy ** 2 + fx ** 2) / 0.01 ** 2)))
    return t / t.std(axis=(1, 2), keepdims=True)


def test_upsampled_dft_matches_direct_inverse_dft():
    rng = np.random.default_rng(0)
    X = np.fft.fft2(rng.standard_normal((16, 12)))
    up, size, off = 4, 6, (3.0, 2.0)
    got = _upsampled_dft(np.conj(X), size, up, off).conj()
    ky, kx = np.fft.fftfreq(16) * 16, np.fft.fftfreq(12) * 12
    for i in range(size):
        for j in range(size):
            uy, ux = (i - off[0]) / up, (j - off[1]) / up
            ref = (X * np.exp(2j * np.pi * (ky[:, None] * uy / 16 + kx[None, :] * ux / 12))).sum()
            assert abs(got[i, j] - ref) < 1e-8 * abs(X).sum()


@pytest.mark.parametrize("ratio", [2, 4, 8])
def test_degrade_keeps_grids_aligned(ratio):
    """A blob centred on a coarse pixel centre must stay centred after degradation."""
    n = 32 * ratio
    yy, xx = np.mgrid[0:n, 0:n] + 0.5          # fine pixel centres, in fine pixels
    cy = cx = (10 + 0.5) * ratio               # centre of coarse pixel (10, 10)
    blob = np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * (2.0 * ratio) ** 2))
    out = degrade(torch.from_numpy(blob[None, None]).float(), [1.2 * ratio / 2], ratio=ratio)[0, 0].numpy()
    w = out / out.sum()
    iy, ix = np.mgrid[0:32, 0:32]
    assert abs((w * iy).sum() - 10) < 0.01 and abs((w * ix).sum() - 10) < 0.01


def test_subpixel_shift_has_no_wraparound_from_opposite_edge():
    img = np.zeros((64, 64))
    img[:, -4:] = 1.0  # strong edge at the right border
    out = subpixel_shift(img, 0.0, 0.4)
    # far from the right border nothing may change (an FFT shift rings across the whole tile)
    assert np.abs(out[:, 4:40]).max() < 1e-6
    assert np.abs(fourier_shift(img, 0.0, 0.4)[:, 4:40]).max() > 1e-3


def test_subpixel_shift_fractional_accuracy():
    tex = _scene(256, 1, 5)[0]
    d = (0.37, -0.61)
    a = subpixel_shift(tex, *d)[20:-20, 20:-20]
    b = fourier_shift(tex, *d)[20:-20, 20:-20]  # exact for this periodic band-limited image
    err = np.abs(a - b) / tex.std()
    assert err.mean() < 0.015 and err.max() < 0.06


def test_native_20m_masks_isolated_blocks_and_flags_resampled():
    rng = np.random.default_rng(0)
    nat = rng.integers(100, 5000, (6, 16, 16))
    rep = nat.repeat(2, axis=1).repeat(2, axis=2)
    x, ok, frac = E0.native_20m(rep)
    assert np.array_equal(x, nat) and ok.all() and frac == 0
    rep[0, 0, 1] += 1                           # one broken block (edge artefact)
    _, ok, frac = E0.native_20m(rep)
    assert not ok[0, 0] and ok.sum() == 16 * 16 - 1 and frac < E0.MAX_RESAMPLED
    resampled = rng.integers(100, 5000, (6, 32, 32))
    assert E0.native_20m(resampled)[2] > 0.9


def _write(path, arr):
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[1], width=arr.shape[2],
                       count=arr.shape[0], dtype="uint16") as d:
        d.write(arr)


def _synthetic_tile(tmp_path, warp_amp10: float = 0.0, d10=(0.40, -0.30), gain_t: float = 0.6):
    """Synthetic SEN2NEON tile. The reference is the scene shifted by d10 (10 m px) plus, optionally,
    a smooth non-rigid field dy(x) = warp_amp10 * sin(2 pi x / width), with target-band gain gain_t."""
    from scipy import ndimage
    n_hr, n_lr = 1024, 256
    big = _scene(2 * n_hr, 12, 1)
    base = 0.15 + 0.02 * big[:, :n_hr, :n_hr]
    shifted = 0.15 + 0.02 * fourier_shift(big, 4 * d10[0], 4 * d10[1])[:, :n_hr, :n_hr]
    if warp_amp10:
        yy, xx = np.mgrid[0:n_hr, 0:n_hr].astype(float)
        dy = 4 * warp_amp10 * np.sin(2 * np.pi * xx / n_hr)
        shifted = np.stack([ndimage.map_coordinates(b, [yy - dy, xx], order=3, mode="mirror") for b in shifted])
    gain = np.ones(12)
    gain[TARGET_IDX] = gain_t
    ref = 0.15 + gain[:, None, None] * (shifted - 0.15)
    s_g, s_t = band_sigmas(MTF, E0.GUIDE_BANDS, 4), band_sigmas(MTF, E0.TARGET_BANDS, 8)
    lr = np.zeros((12, n_lr, n_lr))
    lr[GUIDE_IDX] = degrade(torch.from_numpy(base[GUIDE_IDX]).float()[None], s_g, 4)[0].numpy()
    lr[TARGET_IDX] = degrade(torch.from_numpy(base[TARGET_IDX]).float()[None], s_t, 8)[0].numpy() \
        .repeat(2, axis=1).repeat(2, axis=2)
    _write(tmp_path / "lr/t.tif", np.round(np.clip(lr, 1e-3, 1.0) * 1e4).astype(np.uint16))
    _write(tmp_path / "hr/t.tif", np.round(np.clip(ref, 1e-3, 1.0) * 1e4).astype(np.uint16))
    row = SimpleNamespace(id="2020_TEST_1__0_0", lr="lr/t.tif", hr_2_5m_path="hr/t.tif",
                          temporal_difference_days=2, land_cover_superclass="Rural", cloud_score_plus_cdf=0.9)
    return E0.analyse_tile(tmp_path, row, s_g, s_t, upsample=100)


def test_e0_decomposition_recovers_known_nonrigid_field(tmp_path):
    """R2-1/R2-2: validate the H0' measure on an image with known local distortion. A sinusoidal field of
    amplitude 0.4 px has rms 0.28 px; averaging over 640 m blocks (1/4 period) lowers it to about 0.25."""
    tile, _ = _synthetic_tile(tmp_path, warp_amp10=0.4)
    assert 0.18 < tile["common_t"] < 0.32
    assert 0.18 < tile["common_tg"] < 0.32
    assert 0.18 < tile["rms_tileshift_g"] < 0.32
    assert tile["noise_t"] < 0.1


def test_e0_end_to_end_recovers_known_shift_gain_and_offset(tmp_path):
    """Synthetic tile: S2 is the MTF-degraded scene; the 2.5 m reference is the same scene shifted by
    d (10 m px), with gain `a` and an offset that keeps the tile mean. E0 must recover all three."""
    n_hr, n_lr = 1024, 256                      # real SEN2NEON tile size
    big = _scene(2 * n_hr, 12, 1)
    d10 = (0.40, -0.30)
    gain = np.ones(12)                          # guide bands keep their contrast,
    gain[TARGET_IDX] = 0.6                      # target bands have 0.6x contrast
    base = 0.15 + 0.02 * big[:, :n_hr, :n_hr]
    shifted = 0.15 + 0.02 * fourier_shift(big, 4 * d10[0], 4 * d10[1])[:, :n_hr, :n_hr]
    assert base.min() > 0 and shifted.min() > 0
    ref = 0.15 + gain[:, None, None] * (shifted - 0.15)

    def to_s2(z, idx, ratio, mtfs):
        t = torch.from_numpy(z[idx]).float()[None]
        return degrade(t, mtfs, ratio=ratio)[0].numpy()

    s_g, s_t = band_sigmas(MTF, E0.GUIDE_BANDS, 4), band_sigmas(MTF, E0.TARGET_BANDS, 8)
    lr = np.zeros((12, n_lr, n_lr))
    lr[GUIDE_IDX] = to_s2(base, GUIDE_IDX, 4, s_g)
    lr[TARGET_IDX] = to_s2(base, TARGET_IDX, 8, s_t).repeat(2, axis=1).repeat(2, axis=2)
    # clip before the uint16 cast: negative values would wrap around to ~65000
    lr_u16 = np.round(np.clip(lr, 1e-3, 1.0) * 1e4).astype(np.uint16)
    hr_u16 = np.round(np.clip(ref, 1e-3, 1.0) * 1e4).astype(np.uint16)
    hr_u16[:, :40, :40] = 0                     # some reference nodata
    lr_u16[:, -20:, -20:] = 65535               # and some observation nodata (R2-3.N1)
    _write(tmp_path / "lr/t.tif", lr_u16)
    _write(tmp_path / "hr/t.tif", hr_u16)
    row = SimpleNamespace(id="2020_TEST_1__0_0", lr="lr/t.tif", hr_2_5m_path="hr/t.tif",
                          temporal_difference_days=2, land_cover_superclass="Rural", cloud_score_plus_cdf=0.9)

    tile, bands = E0.analyse_tile(tmp_path, row, s_g, s_t, upsample=100)
    assert tile["status"] == "ok"
    assert abs(tile["g_dy"] - d10[0]) < 0.05 and abs(tile["g_dx"] - d10[1]) < 0.05
    assert abs(tile["t_dy"] - d10[0]) < 0.05 and abs(tile["t_dx"] - d10[1]) < 0.05
    assert tile["rms_translation_g"] < 0.05 and tile["rms_translation_t"] < 0.1
    # after the guide tile shift nothing should be left; the decomposition sees no common field
    assert tile["rms_tileshift_g"] < 0.05 and tile["rms_tileshift_t"] < 0.1
    assert tile["n_blocks_gs"] >= 12 and tile["n_blocks_dec"] >= 10   # edge blocks are kept
    assert tile["common_t"] < 0.05
    b = pd.DataFrame(bands).set_index("band")
    for band in E0.TARGET_BANDS:
        assert abs(b.loc[band, "a_sd"] - 0.6) < 0.05
        assert abs(b.loc[band, "hp_gain_sd"] - 0.6) < 0.08
        expected_mean_r = 0.15 + 0.6 * (b.loc[band, "mean_x"] - 0.15)  # gain applied about 0.15
        assert abs(b.loc[band, "mean_r"] - expected_mean_r) < 1e-3
        assert b.loc[band, "rmse_lin_cv"] < 0.3 * b.loc[band, "rmse_al"]
        assert abs(b.loc[band, "dy"] - d10[0]) < 0.05
    for band in E0.GUIDE_BANDS:
        assert abs(b.loc[band, "a_sd"] - 1.0) < 0.05
        assert b.loc[band, "hp_corr"] > 0.95 > b.loc[band, "hp_corr_raw"]


def test_masked_smoothing_ignores_nodata_values():
    img = np.full((64, 64), 0.2)
    valid = np.ones((64, 64), bool)
    img[30:34, 30:34], valid[30:34, 30:34] = 6.5535, False   # LR nodata scaled to reflectance
    lp = E0._lp(img, 10, valid)
    hp = E0._hp(img, 10, valid)
    assert np.allclose(lp[valid], 0.2, atol=1e-9) and np.allclose(hp[valid], 0.0, atol=1e-9)


def test_stat_mask_excludes_truncated_low_pass():
    """Retained pixels must have the same low-pass as without any mask edge (R3-3.NEW-1)."""
    yy, xx = np.mgrid[0:128, 0:128]
    img = 0.1 + 0.001 * xx + 0.02 * np.sin(xx / 5.0) * np.cos(yy / 7.0)
    m = np.zeros((128, 128), bool)
    m[20:110, 20:110] = True
    ms = E0._stat_mask(m, 10)
    full = E0._lp(img, 10, np.ones_like(m))
    assert ms.sum() > 0
    assert np.abs(E0._lp(img, 10, m)[ms] - full[ms]).max() < 1e-6 * img.std()
    edge = m & ~ms                                   # without the extra erosion the kernel is truncated
    assert np.abs(E0._lp(img, 10, m)[edge] - full[edge]).max() > 1e-3 * img.std()


def test_bogus_own_shift_does_not_change_global_metrics():
    rng = np.random.default_rng(0)
    x = 0.1 + 0.01 * rng.standard_normal((6, 128, 128))
    valid = np.ones((128, 128), bool)
    good = [(0.1, 0.0, 0.8)] * 6
    bad = good[:5] + [(25.0, -25.0, 0.05)]          # huge shift with a low peak
    a = pd.DataFrame(E0._band_rows("t", E0.TARGET_BANDS, x, x, valid, (0.0, 0.0), 20, good)).set_index("band")
    b = pd.DataFrame(E0._band_rows("t", E0.TARGET_BANDS, x, x, valid, (0.0, 0.0), 20, bad)).set_index("band")
    assert len(b) == 6
    assert np.allclose(a["rmse_al"], b["rmse_al"]) and np.allclose(a["hp_corr"], b["hp_corr"])
    assert np.isnan(b.loc["B12"].get("rmse_al_own", np.nan))
