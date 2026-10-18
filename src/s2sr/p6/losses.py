"""Training losses of P6 (plan §4.2–4.4, notes/baselines_fulltext.md). Tensors are (N, C, H, W); masks (N, H, W)
are 1 where the target is valid. Shift sign convention of s2sr.registration (content moves by +d).
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .corrupt import warp


# ----------------------------------------------------------------------------------------------- helpers
def masked_mean(err: torch.Tensor, mask: torch.Tensor | None) -> torch.Tensor:
    """Mean of err (N, C, H, W) over valid pixels, per sample -> (N,)."""
    if mask is None:
        return err.mean(dim=(1, 2, 3))
    m = mask[:, None].to(err.dtype)
    return (err * m).sum(dim=(1, 2, 3)) / (m.sum(dim=(1, 2, 3)) * err.shape[1]).clamp_min(1.0)


def gaussian_blur(x: torch.Tensor, sigma: float) -> torch.Tensor:
    r = max(1, math.ceil(4 * sigma))
    t = torch.arange(-r, r + 1, device=x.device, dtype=x.dtype)
    k = torch.exp(-0.5 * (t / sigma) ** 2)
    k = (k / k.sum()).view(1, 1, 1, -1)
    c = x.shape[1]
    x = F.pad(x, (r, r, r, r), mode="reflect")
    x = F.conv2d(x, k.expand(c, 1, 1, -1), groups=c)
    return F.conv2d(x, k.transpose(2, 3).expand(c, 1, -1, 1), groups=c)


def highpass(x: torch.Tensor, sigma_px: float = 2.0) -> torch.Tensor:
    """HP(z) = z - G_sigma * z; sigma 2 px = 20 m on the 10 m grid (same definition as the endpoints)."""
    return x - gaussian_blur(x, sigma_px)


# ----------------------------------------------------------------------------------------------- pixel losses
def l1_masked(pred, target, mask=None) -> torch.Tensor:
    return masked_mean((pred - target).abs(), mask).mean()


def shift_l1(pred, target, mask=None, radius: float = 0.5, step: float = 0.25) -> torch.Tensor:
    """min over shifts s, |s|_inf <= radius on a `step` grid, of L1(S_s(pred), target), per sample (HighRes-net /
    cPSNR style). The border affected by the shift is excluded through the mask."""
    vals = torch.arange(-radius, radius + 1e-6, step, device=pred.device, dtype=pred.dtype)
    if mask is not None:
        b = math.ceil(radius) + 2
        mask = mask.clone()
        mask[:, :b], mask[:, -b:], mask[:, :, :b], mask[:, :, -b:] = 0, 0, 0, 0
    best = None
    for sy in vals:
        for sx in vals:
            e = masked_mean((warp(pred, sy.expand(pred.shape[0])[:, None, None],
                                  sx.expand(pred.shape[0])[:, None, None]) - target).abs(), mask)
            best = e if best is None else torch.minimum(best, e)
    return best.mean()


def reg_hp(pred, target, mask=None, ref_for_scale=None, sigma_px: float = 2.0) -> torch.Tensor:
    """Gain-invariant high-pass L1 (Reg-HP): each high-pass divided by its own per-sample standard deviation, so a
    global contrast error of the reference does not set the detail amplitude (plan §4.2); the amplitude is left to
    the Wald loss and the LFC. The prediction's scale keeps its gradient (the loss is then exactly scale-invariant
    in the prediction); with a detached scale (v1, used before 15/10/2026) every step shrinks the detail by the
    L1 regression-to-the-mean factor and the amplitude collapses when the correlation is low. `ref_for_scale`
    (detached) is kept for tests only."""
    hp, ht = highpass(pred, sigma_px), highpass(target, sigma_px)
    ref = hp if ref_for_scale is None else highpass(ref_for_scale, sigma_px).detach()
    s = ref.flatten(2).std(dim=2).clamp_min(1e-4)[:, :, None, None]
    st = ht.flatten(2).std(dim=2).clamp_min(1e-4)[:, :, None, None]
    return masked_mean((hp / s - ht / st).abs(), mask).mean()


# ----------------------------------------------------------------------------------------------- CoBi
def cobi_patch(pred, target, mask=None, patch: int = 5, w_s: float = 0.05, n_samples: int = 1024,
               generator: torch.Generator | None = None) -> torch.Tensor:
    """Contextual bilateral loss (Zhang et al., CVPR 2019, Eq. 2) on n x n patches of all bands (no VGG: 6-band
    input). CoBi = mean_i min_j [ cosine_dist(p_i, q_j) + w_s * ||pos_i - pos_j||_2 (pixels) ], on a random
    subset of positions (the same subset for prediction and target)."""
    n, c, h, w = pred.shape
    pp = F.unfold(pred, patch, padding=patch // 2)       # (N, C*patch^2, H*W)
    pt = F.unfold(target, patch, padding=patch // 2)
    L = h * w
    if mask is not None:
        valid = mask.flatten(1) > 0
    else:
        valid = torch.ones(n, L, dtype=torch.bool, device=pred.device)
    ys = torch.arange(h, device=pred.device).repeat_interleave(w).to(pred.dtype)
    xs = torch.arange(w, device=pred.device).repeat(h).to(pred.dtype)
    losses = []
    for i in range(n):
        idx = valid[i].nonzero().squeeze(1)
        if idx.numel() == 0:
            continue
        sel = idx[torch.randperm(idx.numel(), generator=generator, device="cpu").to(idx.device)[:n_samples]]
        a = F.normalize(pp[i, :, sel].T, dim=1)
        b = F.normalize(pt[i, :, sel].T, dim=1)
        d = 1 - a @ b.T
        pos = torch.stack([ys[sel], xs[sel]], 1)
        d = d + w_s * torch.cdist(pos, pos)
        losses.append(d.min(dim=1).values.mean())
    return torch.stack(losses).mean()


# ----------------------------------------------------------------------------------------------- FDL
def sliced_wasserstein(a: torch.Tensor, b: torch.Tensor, n_proj: int = 256,
                       generator: torch.Generator | None = None) -> torch.Tensor:
    """a, b: (N, S, D) sets of S samples in D dims. Random unit projections, sort, mean |difference|."""
    d = a.shape[-1]
    proj = torch.randn(d, n_proj, generator=generator).to(a.device, a.dtype)
    proj = proj / proj.norm(dim=0, keepdim=True)
    pa, _ = torch.sort(a @ proj, dim=1)
    pb, _ = torch.sort(b @ proj, dim=1)
    return (pa - pb).abs().mean()


def fdl(feats_pred: list[torch.Tensor], feats_target: list[torch.Tensor], phase_weight: float = 1.0,
        n_proj: int = 256, generator: torch.Generator | None = None) -> torch.Tensor:
    """Frequency Distribution Loss (Ni et al., CVPR 2024): sliced Wasserstein distance between the
    distributions of DFT amplitude and phase of features; samples = frequencies, dims = channels."""
    total = 0.0
    for fp, ft in zip(feats_pred, feats_target):
        Fp, Ft = torch.fft.fft2(fp.float()), torch.fft.fft2(ft.float())
        to_set = lambda z: z.flatten(2).transpose(1, 2)  # noqa: E731  (N, H*W, C)
        total = total + sliced_wasserstein(to_set(Fp.abs()), to_set(Ft.abs()), n_proj, generator)
        total = total + phase_weight * sliced_wasserstein(to_set(Fp.angle()), to_set(Ft.angle()), n_proj, generator)
    return total / len(feats_pred)


# ----------------------------------------------------------------------------------------------- distillation
def feature_highpass(f: torch.Tensor, k: int = 5) -> torch.Tensor:
    """KD-HP: feature minus its k x k average (plan §4.4)."""
    return f - F.avg_pool2d(f, k, stride=1, padding=k // 2, count_include_pad=False)


class VIDEstimator(nn.Module):
    """PISR/VID estimator: location and scale maps from student features (two 1x1-PReLU-1x1 branches,
    softplus on the scale)."""

    def __init__(self, c_student: int, c_teacher: int, eps: float = 1e-6):
        super().__init__()
        self.mu = nn.Sequential(nn.Conv2d(c_student, c_teacher, 1), nn.PReLU(), nn.Conv2d(c_teacher, c_teacher, 1))
        self.scale = nn.Sequential(nn.Conv2d(c_student, c_teacher, 1), nn.PReLU(), nn.Conv2d(c_teacher, c_teacher, 1))
        self.eps = eps

    def forward(self, f_student, f_teacher, mask=None):
        mu = self.mu(f_student)
        b = F.softplus(self.scale(f_student)) + self.eps
        nll = torch.log(b) + (f_teacher.detach() - mu).abs() / b
        return masked_mean(nll, mask).mean()


class ProjectionHead(nn.Module):
    def __init__(self, c_in: int, c_out: int = 128):
        super().__init__()
        self.net = nn.Sequential(nn.Conv2d(c_in, c_in, 1), nn.ReLU(inplace=True), nn.Conv2d(c_in, c_out, 1))

    def forward(self, f):
        return F.normalize(self.net(f), dim=1)


def _sample_anchors(mask, n_anchor, border, generator):
    """n_anchor valid positions per sample (with replacement when fewer are valid) -> (N, A) rows, cols."""
    n, h, w = mask.shape
    m = mask.clone().bool()
    m[:, :border], m[:, -border:], m[:, :, :border], m[:, :, -border:] = False, False, False, False
    ys, xs = [], []
    for i in range(n):
        idx = m[i].flatten().nonzero().squeeze(1)
        if idx.numel() == 0:
            idx = torch.tensor([(h // 2) * w + w // 2], device=mask.device)
        sel = idx[torch.randint(0, idx.numel(), (n_anchor,), generator=generator).to(idx.device)]
        ys.append(sel // w)
        xs.append(sel % w)
    return torch.stack(ys), torch.stack(xs)


def pcrd(v_s: torch.Tensor, v_t: torch.Tensor, mask: torch.Tensor, mode: str = "global", radius: int = 1,
         tau: float = 0.1, n_anchor: int = 256, n_neg: int = 256, min_neg_dist: int = 8,
         false_neg_thr: float = 0.9, generator: torch.Generator | None = None) -> torch.Tensor:
    """Patch contrastive distillation (plan §4.4), vectorised over the batch. v_s, v_t: L2-normalised
    projections (N, D, H, W).

    mode "local" (pCRD-L): positives = all teacher positions within `radius` (inf-norm) of the anchor; the
      numerator sums over them (each anchor may pick its own shift).
    mode "global" (pCRD-G): one shift s* per tile, the one maximising the mean anchor-positive similarity over
      the (2 radius + 1)^2 candidates; a single positive at p + s*.
    Negatives: n_neg positions of the same tile at inf-distance >= min_neg_dist plus n_neg positions from the
    other tiles of the batch; negatives whose teacher vector is more similar than `false_neg_thr` to the
    positive's teacher vector are dropped (false negatives in homogeneous areas)."""
    n, d, h, w = v_s.shape
    dev = v_s.device
    ay, ax = _sample_anchors(mask, n_anchor, radius + 1, generator)            # (N, A)
    ay, ax = ay.to(dev), ax.to(dev)
    bi = torch.arange(n, device=dev)[:, None]
    S = v_s.permute(0, 2, 3, 1)                                               # (N, H, W, D)
    Tt = v_t.permute(0, 2, 3, 1)
    s_vec = S[bi, ay, ax]                                                     # (N, A, D)
    offs = [(oy, ox) for oy in range(-radius, radius + 1) for ox in range(-radius, radius + 1)]
    pos = torch.stack([Tt[bi, ay + oy, ax + ox] for oy, ox in offs], 2)       # (N, A, K, D)
    sim_pos = (pos * s_vec[:, :, None]).sum(-1) / tau                         # (N, A, K)
    if mode == "global":
        k = sim_pos.mean(1).argmax(1)                                         # (N,)
        sim_pos = sim_pos.gather(2, k[:, None, None].expand(n, n_anchor, 1))
        pos_ref = pos[torch.arange(n, device=dev), :, k]                      # (N, A, D)
    else:
        pos_ref = pos[:, :, len(offs) // 2]
    r = lambda hi, size: torch.randint(0, hi, size, generator=generator).to(dev)  # noqa: E731
    cy, cx = r(h, (n, n_anchor, n_neg)), r(w, (n, n_anchor, n_neg))
    far = torch.maximum((cy - ay[..., None]).abs(), (cx - ax[..., None]).abs()) >= min_neg_dist
    neg = Tt[bi[..., None], cy, cx]                                           # (N, A, n_neg, D)
    if n > 1:
        j = r(n - 1, (n, n_anchor, n_neg))
        j = j + (j >= bi[..., None]).long()
        neg = torch.cat([neg, Tt[j, r(h, j.shape), r(w, j.shape)]], 2)
        far = torch.cat([far, torch.ones_like(far)], 2)
    keep = far & ((neg * pos_ref[:, :, None]).sum(-1) < false_neg_thr)
    sim_neg = ((neg * s_vec[:, :, None]).sum(-1) / tau).masked_fill(~keep, float("-inf"))
    num = torch.logsumexp(sim_pos, 2)
    den = torch.logsumexp(torch.cat([sim_pos, sim_neg], 2), 2)
    return (den - num).mean()
