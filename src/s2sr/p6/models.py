"""P6 networks: the common student (guided EDSR + LFC, the P1 `lfc_guided` model, ~1.8 M parameters), the PISR-style
teacher (encoder + decoder with the student's architecture) and feature taps for distillation (plan §4.3–4.4).
"""
from __future__ import annotations

import torch
import torch.nn as nn

from ..bands import N_TARGET
from ..degrade import band_sigmas, degrade
from ..models import build_model
from ..bands import TARGET_BANDS

STUDENT_SPEC = {"arch": "edsr", "params": {"nf": 80, "n_blocks": 16, "guided": True}, "lfc": True}


def build_student(mtf: dict) -> nn.Module:
    return build_model(STUDENT_SPEC, mtf)


def _body(model: nn.Module) -> nn.Sequential:
    inner = model.model if hasattr(model, "model") else model  # unwrap the LFC wrapper
    return inner.body


class FeatureTap:
    """Collects the outputs of the last `n` residual blocks of a student-architecture model during forward."""

    def __init__(self, model: nn.Module, n: int = 2):
        self.feats: list[torch.Tensor] = []
        blocks = list(_body(model))[-n:]
        self.handles = [b.register_forward_hook(self._hook) for b in blocks]

    def _hook(self, _m, _i, out):
        self.feats.append(out)

    def pop(self) -> list[torch.Tensor]:
        f, self.feats = self.feats, []
        return f

    def remove(self):
        for h in self.handles:
            h.remove()


class TeacherEncoder(nn.Module):
    """PISR encoder: 4 conv blocks (kernels 5, 5, 3, 3; second strided by 2) with PReLU; maps the 10 m reference
    (6 bands) to a 20 m 'compact feature' c with the same shape as the observation x."""

    def __init__(self, nf: int = 64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(N_TARGET, nf, 5, padding=2), nn.PReLU(nf),
            nn.Conv2d(nf, nf, 5, stride=2, padding=2), nn.PReLU(nf),
            nn.Conv2d(nf, nf, 3, padding=1), nn.PReLU(nf),
            nn.Conv2d(nf, N_TARGET, 3, padding=1))

    def forward(self, y):
        return self.net(y)


class Teacher(nn.Module):
    """y (privileged reference, 10 m) -> c = E(y) (20 m) -> decoder(c, guide) -> reconstruction of y.
    Loss: L1(recon, y) + lambda_T * L1(c, D(y)) (imitation); the decoder has the student's architecture."""

    def __init__(self, mtf: dict):
        super().__init__()
        self.encoder = TeacherEncoder()
        self.decoder = build_student(mtf)
        self.sigmas = band_sigmas(mtf, TARGET_BANDS)

    def forward(self, y, guide):
        c = self.encoder(y)
        return c, self.decoder(c, guide)

    def imitation_target(self, y):
        return degrade(y, self.sigmas)
