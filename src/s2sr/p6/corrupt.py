"""Controlled target corruptions for experiment E1 (plan §7): the training target of domain S is corrupted the
way a real cross-sensor reference is, while evaluation uses the clean target.

All operations work on torch tensors (N, C, H, W) with shifts in pixels of that grid, sign convention of
s2sr.registration (content moves by +d). Per-tile parameters are drawn once from a seeded generator so that a
tile is corrupted identically in every epoch and for every method.
"""
from __future__ import annotations

import math

import numpy as np
import torch
import torch.nn.functional as F


def _base_grid(n, h, w, device, dtype):
    ys = (torch.arange(h, device=device, dtype=dtype) + 0.5) / h * 2 - 1
    xs = (torch.arange(w, device=device, dtype=dtype) + 0.5) / w * 2 - 1
    gy, gx = torch.meshgrid(ys, xs, indexing="ij")
    return torch.stack([gx, gy], dim=-1).expand(n, h, w, 2)


def warp(img: torch.Tensor, dy: torch.Tensor, dx: torch.Tensor) -> torch.Tensor:
    """Move content by the displacement field (dy, dx) in pixels: out(p) = img(p - d(p)).
    dy, dx: (N, H, W) or broadcastable (N, 1, 1). Bicubic sampling, border padding."""
    n, _, h, w = img.shape
    grid = _base_grid(n, h, w, img.device, img.dtype).clone()
    grid[..., 0] = grid[..., 0] - dx * 2 / w
    grid[..., 1] = grid[..., 1] - dy * 2 / h
    return F.grid_sample(img, grid, mode="bicubic", padding_mode="border", align_corners=False)


def smooth_field(n: int, h: int, w: int, rms: float, corr_px: float, rng: np.random.Generator):
    """Random smooth displacement field (dy, dx), each (n, h, w), with total rms `rms` px and spatial
    correlation length about corr_px (Gaussian-filtered white noise, renormalised)."""
    from scipy import ndimage
    out = []
    for _ in range(n):
        f = np.stack([ndimage.gaussian_filter(rng.standard_normal((h, w)), corr_px, mode="wrap") for _ in range(2)])
        f -= f.mean(axis=(1, 2), keepdims=True)
        f *= rms / max(np.sqrt((f ** 2).sum(0).mean()), 1e-12)
        out.append(f)
    out = torch.from_numpy(np.stack(out)).float()
    return out[:, 0], out[:, 1]


def content_change(img: torch.Tensor, other: torch.Tensor, frac: float, rng: np.random.Generator,
                   soft_px: float = 4.0) -> torch.Tensor:
    """Replace a random elliptical region covering about `frac` of each image by `other` (another tile's
    content, same bands), with a soft edge: mimics surface change between acquisitions."""
    from scipy import ndimage
    n, _, h, w = img.shape
    masks = []
    for _ in range(n):
        cy, cx = rng.uniform(0.2, 0.8) * h, rng.uniform(0.2, 0.8) * w
        a = math.sqrt(frac * h * w / math.pi) * rng.uniform(0.7, 1.3)
        b = frac * h * w / (math.pi * a)
        yy, xx = np.mgrid[0:h, 0:w]
        mk = (((yy - cy) / a) ** 2 + ((xx - cx) / b) ** 2 <= 1).astype(np.float64)
        masks.append(ndimage.gaussian_filter(mk, soft_px))
    m = torch.from_numpy(np.stack(masks)).float().to(img.device)[:, None]
    return img * (1 - m) + other * m


def draw_tile_params(scenario: dict, rng: np.random.Generator, shift_pool: np.ndarray | None = None) -> dict:
    """Per-tile corruption parameters for a scenario (see configs/p6_e1.yaml)."""
    p = dict(dy=0.0, dx=0.0, field_rms=0.0, change=False)
    g = scenario.get("global_shift")
    if g == "e0":            # draw (dy, dx) from the empirical E0 distribution of tile shifts
        dy, dx = shift_pool[rng.integers(len(shift_pool))]
        p.update(dy=float(dy), dx=float(dx))
    elif g is not None:      # fixed magnitude, random direction
        ang = rng.uniform(0, 2 * math.pi)
        p.update(dy=float(g) * math.sin(ang), dx=float(g) * math.cos(ang))
    if rng.uniform() < scenario.get("nonrigid_frac", 0.0):
        p["field_rms"] = float(scenario["nonrigid_rms"])
    if rng.uniform() < scenario.get("change_frac", 0.0):
        p["change"] = True
    return p
