"""Tests for P6 losses, metrics and target corruptions (src/s2sr/p6)."""
import numpy as np
import torch

from s2sr.p6 import corrupt, losses, metrics
from s2sr.registration import fourier_shift

torch.manual_seed(0)


def _tex(n=2, c=6, h=64, w=64, seed=0):
    rng = np.random.default_rng(seed)
    fy, fx = np.fft.fftfreq(h)[:, None], np.fft.fftfreq(w)[None, :]
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, c, h, w))) / (1 + (fy ** 2 + fx ** 2) / 0.02 ** 2)))
    return torch.from_numpy(0.2 + 0.03 * t / t.std()).float()


def test_warp_matches_module_shift_convention():
    x = _tex(1, 1, 64, 64)
    d = (0.6, -0.4)
    ref = torch.from_numpy(fourier_shift(x[0, 0].numpy().astype(float), *d)).float()
    out = corrupt.warp(x, torch.tensor([[[d[0]]]]), torch.tensor([[[d[1]]]]))[0, 0]
    wrong = corrupt.warp(x, torch.tensor([[[-d[0]]]]), torch.tensor([[[-d[1]]]]))[0, 0]
    err = (out[8:-8, 8:-8] - ref[8:-8, 8:-8]).abs().mean()
    assert err < 0.03 * x.std()                                       # bicubic interpolation error only
    assert (wrong[8:-8, 8:-8] - ref[8:-8, 8:-8]).abs().mean() > 5 * err  # the opposite sign is clearly worse


def test_shift_l1_tolerates_shifts_within_radius_only():
    y = _tex()
    shifted = corrupt.warp(y, torch.full((2, 1, 1), 0.5), torch.full((2, 1, 1), -0.25))
    m = torch.ones(2, 64, 64)
    assert losses.shift_l1(y, shifted, m, radius=0.5) < 0.2 * losses.l1_masked(y, shifted, m)
    far = corrupt.warp(y, torch.full((2, 1, 1), 2.0), torch.zeros(2, 1, 1))
    assert losses.shift_l1(y, far, m, radius=0.5) > 0.5 * losses.l1_masked(y, far, m)


def test_reg_hp_is_gain_invariant():
    y = _tex()
    m = torch.ones(2, 64, 64)
    low_contrast = 0.2 + 0.5 * (y - 0.2)
    assert losses.reg_hp(y, low_contrast, m, ref_for_scale=y) < 0.02
    assert losses.l1_masked(y, low_contrast, m) > 1e-3


def test_cobi_zero_for_identical_and_more_tolerant_than_l1():
    y = _tex()
    g = torch.Generator().manual_seed(0)
    assert losses.cobi_patch(y, y, generator=g) < 1e-5
    sh = corrupt.warp(y, torch.full((2, 1, 1), 1.0), torch.zeros(2, 1, 1))
    other = _tex(seed=5)
    g = torch.Generator().manual_seed(0)
    c_sh = losses.cobi_patch(y, sh, generator=g)
    g = torch.Generator().manual_seed(0)
    c_ot = losses.cobi_patch(y, other, generator=g)
    assert c_sh < c_ot


def test_fdl_zero_for_identical_and_positive_otherwise():
    f = [torch.randn(2, 8, 16, 16), torch.randn(2, 8, 8, 8)]
    g = torch.Generator().manual_seed(0)
    assert losses.fdl(f, f, generator=g) < 1e-6
    g = torch.Generator().manual_seed(0)
    assert losses.fdl(f, [x + torch.randn_like(x) for x in f], generator=g) > 0.01


def test_vid_estimator_learns_teacher_features():
    est = losses.VIDEstimator(8, 8)
    fs = torch.randn(4, 8, 16, 16)
    ft = 2 * fs + 0.1
    opt = torch.optim.Adam(est.parameters(), lr=1e-2)
    first = None
    for _ in range(200):
        loss = est(fs, ft)
        first = loss.item() if first is None else first
        opt.zero_grad()
        loss.backward()
        opt.step()
    assert loss.item() < first - 1.0


def _proj(seed=0, n=3, d=16, h=32, w=32):
    g = torch.Generator().manual_seed(seed)
    return torch.nn.functional.normalize(torch.randn(n, d, h, w, generator=g), dim=1)


def test_pcrd_shapes_and_ordering():
    v = _proj()
    m = torch.ones(3, 32, 32)
    g = torch.Generator().manual_seed(1)
    matched = losses.pcrd(v, v, m, mode="global", generator=g)
    g = torch.Generator().manual_seed(1)
    random = losses.pcrd(v, _proj(seed=9), m, mode="global", generator=g)
    assert torch.isfinite(matched) and matched < random


def test_pcrd_global_tolerates_one_pixel_teacher_shift():
    v = _proj()
    vt = torch.roll(v, shifts=(1, 0), dims=(2, 3))   # teacher content one feature pixel down
    m = torch.ones(3, 32, 32)
    g = torch.Generator().manual_seed(2)
    aligned = losses.pcrd(v, v, m, mode="global", generator=g)
    g = torch.Generator().manual_seed(2)
    shifted = losses.pcrd(v, vt, m, mode="global", generator=g)
    assert abs(shifted.item() - aligned.item()) < 0.05


def test_metrics_identity_and_shift():
    y = _tex(1, 2, 128, 128)[0].numpy().astype(float)
    v = np.ones((128, 128), bool)
    assert np.allclose(metrics.hp_corr(y, y, v), 1.0)
    assert np.allclose(metrics.edge_error(y, y, v), 0.0)
    sh = np.stack([fourier_shift(b, 1.0, 0.0) for b in y])
    assert np.all(metrics.edge_error(sh, y, v) > 0.5)
    assert np.all(metrics.hp_corr(sh, y, v) < 0.9)
    assert np.all(metrics.hp_psnr(y, y, v) > 100)


def test_corruption_params_and_field():
    rng = np.random.default_rng(0)
    pool = np.array([[0.3, -0.1], [0.5, 0.2]])
    p = corrupt.draw_tile_params(dict(global_shift="e0", nonrigid_frac=1.0, nonrigid_rms=0.3), rng, pool)
    assert (p["dy"], p["dx"]) in {(0.3, -0.1), (0.5, 0.2)} and p["field_rms"] == 0.3
    dy, dx = corrupt.smooth_field(2, 64, 64, rms=0.3, corr_px=8, rng=rng)
    assert abs(float(torch.sqrt((dy ** 2 + dx ** 2).mean())) - 0.3) < 1e-4


def test_reg_hp_gradient_does_not_shrink_detail():
    """Scale invariance in the prediction: d loss / d c at c = 1 for loss(c * pred) is zero, so the Reg-HP loss
    exerts no force on the detail amplitude (the v1 detached scale shrank it)."""
    g = torch.Generator().manual_seed(0)
    y = torch.rand(2, 6, 64, 64, generator=g)
    pred = (0.6 * y + 0.4 * torch.rand(2, 6, 64, 64, generator=g)).requires_grad_(True)
    loss = losses.reg_hp(pred, 0.7 * y)
    (grad,) = torch.autograd.grad(loss, pred)
    assert abs(float((grad * pred).sum())) < 1e-4 * float(grad.abs().sum() * pred.abs().mean())
