"""Tests for the E1-R harmonisation emulation (scripts/p6_make_e1r_targets.py)."""
import importlib.util
from pathlib import Path

import numpy as np
import torch

from s2sr.degrade import degrade

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("e1r", ROOT / "scripts" / "p6_make_e1r_targets.py")
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)


def _z(n=512, seed=0):
    rng = np.random.default_rng(seed)
    fy, fx = np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :]
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((6, n, n))) / (1 + (fy ** 2 + fx ** 2) / 0.01 ** 2)))
    return 0.2 + 0.03 * t / t.std()


def test_ols_harmonisation_compresses_contrast_by_r_and_inversion_restores_it():
    z = _z()
    x20 = degrade(torch.from_numpy(z).float()[None], [4.0] * 6, 8)[0].double().numpy()   # smooth "S2"
    ok = np.ones(z.shape[-2:], bool)
    m, n, R = E.harmonise(x20, z, ok)
    assert np.all(m < 1.0) and np.all(R < 1.0)
    # OLS slope = R * sd(x_up) / sd(z): the harmonised reference has sd = R * sd(x_up)
    zh = m[:, None, None] * z + n[:, None, None]
    assert np.all(np.abs(zh.mean(axis=(1, 2)) - z.mean(axis=(1, 2))) < 2e-3)
    # inversion with the low-pass gain between x and D(z_h) restores the low-pass contrast of x
    d20 = degrade(torch.from_numpy(zh).float()[None], [4.0] * 6, 8)[0].double().numpy()
    a_hat = E.symmetric_lowpass_gain(x20, d20, np.ones(x20.shape[-2:], bool))
    assert np.allclose(a_hat, m, rtol=0.05)
    zinv = zh.mean() + (zh - zh.mean()) / a_hat[:, None, None]
    assert np.allclose(zinv.std(axis=(1, 2)), z.std(axis=(1, 2)), rtol=0.05)
