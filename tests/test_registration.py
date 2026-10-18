import numpy as np
import pytest

from s2sr.registration import (block_shifts, common_noise_decomposition, fourier_shift, linear_fit,
                               match_blocks, phase_corr, shift_field_residual)


def _texture(h=128, w=128, c=1, seed=0):
    """Scene-like texture (1/f^2 spectrum), cut from a larger periodic image so tile edges are not
    periodic. Returns a function mapping a shift (dy, dx) to the cropped, shifted tile."""
    rng = np.random.default_rng(seed)
    n = 4 * max(h, w)
    fy, fx = np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :]
    big = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((c, n, n)))
                               / (1 + (fy ** 2 + fx ** 2) / 0.05 ** 2)))
    return lambda dy=0.0, dx=0.0: fourier_shift(big, dy, dx)[:, n // 3:n // 3 + h, n // 3:n // 3 + w]


@pytest.mark.parametrize("d", [(0.0, 0.0), (0.3, -0.7), (-1.25, 0.5), (2.6, 1.1)])
def test_phase_corr_recovers_subpixel_shift_and_sign(d):
    tex = _texture(c=3)
    ref, mov = tex(), tex(*d)
    dy, dx, peak = phase_corr(ref, mov, upsample=50)
    assert abs(dy - d[0]) < 0.05 and abs(dx - d[1]) < 0.05
    assert peak > 0.3


def test_phase_corr_unbiased_on_smooth_coarse_tiles():
    """Aliased near-Nyquist frequencies bias smooth (MTF-blurred, decimated) tiles towards zero."""
    import torch
    from s2sr.degrade import degrade
    rng = np.random.default_rng(4)
    n = 1024
    fy, fx = np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :]
    big = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((2, n, n))) / (1 + (fy ** 2 + fx ** 2) / 0.01 ** 2)))
    d = (1.6, -1.2)  # in 2.5 m px, i.e. (0.2, -0.15) px of a 20 m grid
    to20 = lambda z: degrade(torch.from_numpy(z[:, :512, :512]).float()[None], [4.0, 4.0], 8)[0].numpy()
    ref, mov = to20(big), to20(fourier_shift(big, *d))
    dy0, dx0, _ = phase_corr(ref, mov, fmax=0.75)           # all frequencies, incl. aliased ones
    dy, dx, _ = phase_corr(ref, mov)
    assert abs(dy0 - 0.2) > 0.02                             # the bias is real ...
    assert abs(dy - 0.2) < 0.01 and abs(dx + 0.15) < 0.01    # ... and the frequency limit removes it


def test_fourier_shift_convention():
    img = np.zeros((32, 32))
    img[10, 12] = 1.0
    out = fourier_shift(img, 3, -2)
    assert np.unravel_index(np.argmax(out), out.shape) == (13, 10)


def test_phase_corr_with_nodata_mask():
    tex = _texture(c=2, seed=1)
    ref, mov = tex(), tex(0.5, -0.4)
    valid = np.ones(ref.shape[-2:], bool)
    valid[:30, :40] = False
    mov[:, ~valid] = 0.0
    dy, dx, _ = phase_corr(ref, mov, valid, upsample=50)
    assert abs(dy - 0.5) < 0.1 and abs(dx + 0.4) < 0.1


def test_block_shifts_detect_translation_vs_affine_field():
    tex = _texture(256, 256, c=2, seed=2)
    ref, mov = tex(), tex(0.6, -0.3)
    valid = np.ones((256, 256), bool)
    bl = block_shifts(ref, mov, valid, block=64, min_std=0.0)
    res = shift_field_residual(bl, min_peak=0.0)
    assert res["n_blocks"] == 16
    assert abs(res["med_dy"] - 0.6) < 0.05 and abs(res["med_dx"] + 0.3) < 0.05
    assert res["rms_translation"] < 0.05


def test_block_shifts_keep_blocks_when_only_some_bands_are_flat():
    tex = _texture(128, 128, c=2, seed=3)
    ref, mov = tex(), tex(0.4, 0.2)
    ref[0] *= 1e-3  # one nearly flat band, as visible bands over dense forest
    mov[0] *= 1e-3
    ref[1] *= 0.05 / ref[1].std()
    mov[1] *= 0.05 / mov[1].std()
    bl = block_shifts(ref, mov, np.ones((128, 128), bool), block=64, min_std=0.005)
    assert len(bl) == 4


def test_linear_fit_recovers_gain_offset():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 0.4, 5000)
    r = 1.1 * x - 0.01
    f = linear_fit(x, r)
    assert abs(f["a"] - 1.1) < 1e-6 and abs(f["c"] + 0.01) < 1e-6
    assert f["rmse_lin"] < 1e-9 < f["rmse_raw"]


def test_phase_corr_returns_nan_for_flat_or_invalid_input():
    flat = np.ones((2, 64, 64))
    dy, dx, pk = phase_corr(flat, flat)
    assert np.isnan(dy) and np.isnan(dx) and pk == 0.0
    tex = _texture(64, 64, c=2)()
    dy, dx, pk = phase_corr(tex, tex, np.zeros((64, 64), bool))
    assert np.isnan(dy) and pk == 0.0


@pytest.mark.parametrize("true_rms,noise_rms", [(0.0, 0.15), (0.2, 0.15), (0.3, 0.05)])
def test_common_noise_decomposition_recovers_known_field_and_noise(true_rms, noise_rms):
    rng = np.random.default_rng(1)
    n = 4000
    s = rng.normal(0, true_rms / np.sqrt(2), (n, 2)) + np.array([0.4, -0.2])  # field + translation
    a = s + rng.normal(0, noise_rms / np.sqrt(2), (n, 2))
    b = s + rng.normal(0, noise_rms / np.sqrt(2), (n, 2))
    d = common_noise_decomposition(a, b)
    assert abs(d["common"] - true_rms) < 0.02
    assert abs(d["noise"] - noise_rms) < 0.02


def test_match_blocks_pairs_grids_and_scales_units():
    g = [dict(cy=32, cx=32, dy=0.1, dx=0.0), dict(cy=96, cx=32, dy=0.2, dx=0.0)]
    t = [dict(cy=16, cx=16, dy=0.05, dx=0.0), dict(cy=48, cx=48, dy=0.3, dx=0.0)]  # 20 m grid
    fg, ft = match_blocks(g, t, scale=2.0)
    assert fg.shape == (1, 2) and np.allclose(fg, [[0.1, 0.0]]) and np.allclose(ft, [[0.1, 0.0]])


def test_common_noise_decomposition_small_n_bias():
    """At E0's block count (n = 14) and noise (0.225 px) the clipped sqrt is biased upwards under a zero
    field, while common2 = mean(a . b) stays unbiased; H0' therefore thresholds common2."""
    rng = np.random.default_rng(7)
    c2, c = [], []
    for _ in range(4000):
        a = rng.normal(0, 0.225 / np.sqrt(2), (14, 2))
        b = rng.normal(0, 0.225 / np.sqrt(2), (14, 2))
        d = common_noise_decomposition(a, b, centre="none")
        c2.append(d["common2"])
        c.append(d["common"])
    assert abs(np.mean(c2)) < 0.002            # unbiased
    assert np.mean(c) > 0.03                   # the clipped sqrt is not (about 0.04 px under a zero field)
    assert np.mean(np.array(c2) <= 0.25 ** 2) > 0.99
