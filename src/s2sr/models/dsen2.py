"""DSen2 (Lanaras et al., 2018) and a re-implementation of DSen2+CA+HPF (Nguyen-Vi et al., 2025).

DSen2: concatenate the guide bands with the bicubic-upsampled 20 m bands, one conv + ReLU,
a stack of residual blocks (scaling 0.1), one conv to the output bands, plus a global skip
from the upsampled input. The paper uses 128 features and 6 blocks.

DSen2+CA+HPF: re-implemented from the paper's description (channel attention in every
residual block and a high-pass enhancement of the guide). TODO: check against the authors'
code or the full text before reporting it under their name.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from ..bands import N_GUIDE, N_TARGET
from ..degrade import upsample
from .common import ChannelAttention, ResBlock, conv3


class DSen2(nn.Module):
    def __init__(self, nf: int = 128, n_blocks: int = 6, channel_attention: bool = False,
                 high_pass: bool = False):
        super().__init__()
        self.high_pass = high_pass
        cin = N_GUIDE + N_TARGET + (N_GUIDE if high_pass else 0)
        self.head = nn.Sequential(conv3(cin, nf), nn.ReLU(inplace=True))
        self.body = nn.Sequential(*[
            ResBlock(nf, 0.1, ChannelAttention(nf) if channel_attention else None)
            for _ in range(n_blocks)])
        self.tail = conv3(nf, N_TARGET)

    @staticmethod
    def _high_pass(g):
        # guide minus its 3x3 box-filtered version: emphasises edges of the 10 m bands
        return g - F.avg_pool2d(F.pad(g, (1, 1, 1, 1), mode="reflect"), 3, stride=1)

    def forward(self, lr, guide):
        up = upsample(lr)
        feats = [guide, up] + ([self._high_pass(guide)] if self.high_pass else [])
        x = self.head(torch.cat(feats, dim=1))
        return up + self.tail(self.body(x))
