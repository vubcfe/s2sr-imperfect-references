"""Model registry: build any benchmark method from its configs/methods.yaml entry."""
from __future__ import annotations

import torch.nn as nn

from ..bands import TARGET_BANDS
from ..degrade import band_sigmas, upsample
from .common import LowFrequencyConstraint, count_params
from .dsen2 import DSen2
from .edsr_rcan import EDSR, RCAN


class Bicubic(nn.Module):
    trainable = False

    def forward(self, lr, guide):
        return upsample(lr)


class ATPRK(nn.Module):
    """Placeholder. ATPRK (Wang et al., 2016) is not re-implemented here.

    TODO: run it with an existing implementation (the authors' MATLAB code, or the version in
    the Sentinel-2 sharpening benchmark of Ciotola et al., 2025, if included) and import its
    outputs with scripts/import_external.py. Do not report ATPRK numbers from anything else.
    """
    trainable = False

    def forward(self, lr, guide):
        raise NotImplementedError("ATPRK must be run with an external implementation")


_ARCH = {"dsen2": DSen2, "edsr": EDSR, "rcan": RCAN}


def build_model(spec: dict, mtf: dict | None = None) -> nn.Module:
    """spec: {"arch": ..., "params": {...}, "lfc": bool}."""
    arch = spec["arch"]
    if arch == "bicubic":
        return Bicubic()
    if arch == "atprk":
        return ATPRK()
    model = _ARCH[arch](**spec.get("params", {}))
    if spec.get("lfc"):
        if mtf is None:
            raise ValueError("LFC needs the MTF table")
        model = LowFrequencyConstraint(model, band_sigmas(mtf, TARGET_BANDS))
    model.trainable = True
    return model


__all__ = ["build_model", "count_params", "Bicubic", "ATPRK", "DSen2", "EDSR", "RCAN"]
