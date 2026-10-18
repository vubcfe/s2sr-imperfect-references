"""P6 datasets on domain S (plan §3, §7 E1).

Training items: random aligned crops (10 m crop `crop`, 20 m crop `crop/2`) of a domain-S tile with the training
target chosen by the method: clean y, corrupted y_c or registered y_reg (data/p6_e1/targets_<scenario>.h5), plus
random flips / 90-degree rotations applied identically to every product. Evaluation items: full tiles with the
clean target.
"""
from __future__ import annotations

import h5py
import numpy as np
import torch
from torch.utils.data import Dataset

Q = 1e-4


class _H5:
    """Lazily opened HDF5 handle (one per worker process)."""

    def __init__(self, path):
        self.path, self._f = str(path), None

    @property
    def f(self):
        if self._f is None:
            self._f = h5py.File(self.path, "r")
        return self._f


class SPatches(Dataset):
    def __init__(self, s_h5, tiles, targets_h5=None, target="y", crop: int = 128, augment: bool = True,
                 length: int = 10 ** 7, seed: int = 0, offset: int = 0):
        self.s, self.t = _H5(s_h5), (_H5(targets_h5) if targets_h5 else None)
        self.tiles, self.target, self.crop, self.augment = list(tiles), target, crop, augment
        self.length, self.seed, self.offset = length, seed, offset  # offset: items already used (resume)
        if target != "y" and targets_h5 is None:
            raise ValueError(f"target {target} needs a targets file")

    def __len__(self):
        return self.length

    def __getitem__(self, i):
        rng = np.random.default_rng((self.seed, i + self.offset))
        tile = self.tiles[rng.integers(len(self.tiles))]
        g = self.s.f[tile]
        h = g["y"].shape[-1]
        c = self.crop
        oy, ox = 2 * rng.integers(0, (h - c) // 2 + 1), 2 * rng.integers(0, (h - c) // 2 + 1)
        sl10, sl20 = np.s_[oy:oy + c, ox:ox + c], np.s_[oy // 2:(oy + c) // 2, ox // 2:(ox + c) // 2]
        x = g["x"][(slice(None),) + sl20] * Q
        gd = g["g"][(slice(None),) + sl10] * Q
        y = g["y"][(slice(None),) + sl10] * Q
        yt = y if self.target == "y" else self.t.f[tile][self.target][(slice(None),) + sl10] * Q
        v = g["valid10"][sl10]
        if self.augment:
            k, flip = int(rng.integers(4)), bool(rng.integers(2))
            aug = lambda a: np.rot90(a[..., ::-1] if flip else a, k, axes=(-2, -1))  # noqa: E731
            x, gd, y, yt, v = map(aug, (x, gd, y, yt, v))
        t = lambda a: torch.from_numpy(np.ascontiguousarray(a)).float()  # noqa: E731
        return dict(x=t(x), g=t(gd), y=t(y), target=t(yt), valid=t(v))


class STiles(Dataset):
    """Full tiles with the clean target, for evaluation."""

    def __init__(self, s_h5, tiles, ref_key: str = "y"):
        self.s, self.tiles, self.ref_key = _H5(s_h5), list(tiles), ref_key

    def __len__(self):
        return len(self.tiles)

    def __getitem__(self, i):
        g = self.s.f[self.tiles[i]]
        t = lambda a: torch.from_numpy(a[:] * Q).float()  # noqa: E731
        return dict(x=t(g["x"]), g=t(g["g"]), y=t(g[self.ref_key]), valid=torch.from_numpy(g["valid10"][:]).float(),
                    tile=self.tiles[i])
