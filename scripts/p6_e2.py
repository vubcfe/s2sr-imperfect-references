"""P6 experiment E2 driver (configs/p6_e2.yaml): train on domain R (real Sentinel-2 + corrected SEN2NEON
reference) or domain S, evaluate on the test sites.

  python scripts/p6_e2.py run --method shift_l1 --seed 0 --gpu 0
  python scripts/p6_e2.py teacher --seed 0 --gpu 0

Outputs runs/p6_e2/<method>/seed<k>/: student.pt, test_y.csv (corrected reference), test_r_reg.csv
(contrast as delivered), wald_test.csv (Wald self-consistency), run.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.bands import GUIDE_BANDS, TARGET_BANDS  # noqa: E402
from s2sr.degrade import band_sigmas, degrade  # noqa: E402
from s2sr.p6 import train as T  # noqa: E402
from s2sr.p6.data import STiles  # noqa: E402
from s2sr.p6.models import build_student  # noqa: E402


def load_cfg(path="configs/p6_e2.yaml"):
    return yaml.safe_load((ROOT / path).read_text())


def _tiles(h5):
    with h5py.File(h5, "r") as f:
        return sorted(f)


def paths(cfg, domain, seed):
    d = ROOT / (cfg["domain_r"] if domain == "R" else cfg["domain_s"])
    runs = ROOT / cfg["runs_dir"]
    tr = d / "train.h5"
    return dict(s_train=tr, s_val=d / "val.h5", train_tiles=_tiles(tr), targets=tr,
                teacher=runs / "_teacher" / f"seed{seed}" / "teacher.pt")


@torch.no_grad()
def wald_consistency(ckpt, mtf, test_h5, tiles, device):
    """PSNR (peak 1) of f(D(x), D(g)) vs x on valid 20 m pixels, per tile (H1 non-inferiority metric)."""
    model = build_student(mtf).to(device).eval()
    model.load_state_dict(torch.load(ckpt, map_location=device))
    st, sg = band_sigmas(mtf, TARGET_BANDS), band_sigmas(mtf, GUIDE_BANDS)
    rows = []
    for b in torch.utils.data.DataLoader(STiles(test_h5, tiles), batch_size=4):
        x, g = b["x"].to(device), b["g"].to(device)
        pred = model(degrade(x, st), degrade(g, sg))
        v20 = (b["valid"].reshape(-1, 128, 2, 128, 2).amin(dim=(2, 4)) > 0).to(device)
        for i, tile in enumerate(b["tile"]):
            e = ((pred[i] - x[i]) ** 2)[:, v20[i]].mean()
            rows.append(dict(tile=tile, site=tile.split("_")[1], wald_psnr=float(10 * torch.log10(1 / e))))
    return pd.DataFrame(rows)


def evaluate_val(ck, cfg, mtf, out, device):
    val_h5 = ROOT / cfg["domain_r"] / "val.h5"
    tiles = _tiles(val_h5)
    T.evaluate(ck, {}, mtf, tiles, out / "val_y.csv", device, eval_h5=val_h5, ref_key="y")
    wald_consistency(ck, mtf, val_h5, tiles, device).to_csv(out / "wald_val.csv", index=False)


def run(cfg, name, seed, device):
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    spec = cfg["methods"][name]
    out = ROOT / cfg["runs_dir"] / name / f"seed{seed}"
    out.mkdir(parents=True, exist_ok=True)
    if (out / "run.json").exists():
        return out
    t0 = time.time()
    n = cfg["training"]["iterations"]
    if "stages" in spec:  # B-sim+FT: stage 1 on domain S, stage 2 on domain R from the stage-1 weights
        prev = None
        for k, st in enumerate(spec["stages"]):
            st = dict(st)
            if prev is not None:
                st["init_ckpt"] = str(prev)
            prev = T.train_student(st, cfg, paths(cfg, st["domain"], seed), mtf, out / f"stage{k}", seed, device,
                                   int(n * st["fraction"]))
        ck = out / "student.pt"
        ck.write_bytes(Path(prev).read_bytes())
    else:
        p = paths(cfg, spec["domain"], seed)
        if spec.get("teacher"):
            T.train_teacher(spec["target"], cfg, p, mtf, p["teacher"].parent, seed, device)
        ck = T.train_student(dict(spec), cfg, p, mtf, out, seed, device)
    if cfg.get("eval_val"):  # post-hoc selection on the validation sites (never on test)
        evaluate_val(ck, cfg, mtf, out, device)
    if cfg.get("val_only"):
        T.save_json(out / "run.json", dict(method=name, spec=spec, seed=seed, val_only=True,
                                           minutes=round((time.time() - t0) / 60, 1)))
        return out
    test_h5 = ROOT / cfg["domain_r"] / "test.h5"
    tiles = _tiles(test_h5)
    T.evaluate(ck, {}, mtf, tiles, out / "test_y.csv", device, eval_h5=test_h5, ref_key="y")
    T.evaluate(ck, {}, mtf, tiles, out / "test_r_reg.csv", device, eval_h5=test_h5, ref_key="r_reg")
    wald_consistency(ck, mtf, test_h5, tiles, device).to_csv(out / "wald_test.csv", index=False)
    T.save_json(out / "run.json", dict(method=name, spec=spec, seed=seed, minutes=round((time.time() - t0) / 60, 1),
                                       frozen=cfg.get("frozen")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "teacher", "evalval"])
    ap.add_argument("--ckpt_dir", help="evalval: run directory with student.pt (writes val_y.csv, wald_val.csv)")
    ap.add_argument("--config", default="configs/p6_e2.yaml")
    ap.add_argument("--method", default="reg_l1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gpu", type=int, default=0)
    a = ap.parse_args()
    cfg = load_cfg(a.config)
    device = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    if a.cmd == "evalval":
        mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
        d = ROOT / a.ckpt_dir
        evaluate_val(d / "student.pt", cfg, mtf, d, device)
    elif a.cmd == "teacher":
        mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
        p = paths(cfg, "R", a.seed)
        T.train_teacher("y", cfg, p, mtf, p["teacher"].parent, a.seed, device)
    else:
        run(cfg, a.method, a.seed, device)


if __name__ == "__main__":
    main()
