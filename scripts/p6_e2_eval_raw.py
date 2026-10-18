"""Score every finished E2 run (and the E2 post-hoc runs) against the reference exactly as released (r_raw: no
registration, no gain inversion), as a sensitivity analysis free of any correction made by the authors.
Writes test_r_raw.csv next to test_y.csv in each run directory; skips runs that already have it.

  python scripts/p6_e2_eval_raw.py --gpu 0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.p6 import train as T  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", type=int, default=0)
    a = ap.parse_args()
    dev = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    test_h5 = ROOT / "data/p6_domain_r/test.h5"
    with h5py.File(test_h5, "r") as f:
        tiles = sorted(f)
    for runs in ("runs/p6_e2", "runs/p6_e2_posthoc"):
        for d in sorted((ROOT / runs).glob("*/seed*")):
            if d.parts[-2].startswith("_") or not (d / "student.pt").exists() or not (d / "test_y.csv").exists():
                continue
            if (d / "test_r_raw.csv").exists():
                continue
            T.evaluate(d / "student.pt", {}, mtf, tiles, d / "test_r_raw.csv", dev, eval_h5=test_h5, ref_key="r_raw")
            print("scored", d.relative_to(ROOT), flush=True)


if __name__ == "__main__":
    main()
