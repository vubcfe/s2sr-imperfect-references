"""Smoke test of the E2 driver (scripts/p6_e2.py) on tiny synthetic domains R and S (CPU)."""
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("e2", ROOT / "scripts" / "p6_e2.py")
E2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E2)


def _domain(d, roles, keys):
    rng = np.random.default_rng(0)
    d.mkdir()
    for role, names in roles.items():
        with h5py.File(d / f"{role}.h5", "w") as f:
            for t in names:
                g = f.create_group(t)
                g["x"] = rng.integers(1000, 3000, (6, 128, 128)).astype(np.uint16)
                g["g"] = rng.integers(1000, 3000, (4, 256, 256)).astype(np.uint16)
                for k in keys:
                    g[k] = rng.integers(1000, 3000, (6, 256, 256)).astype(np.uint16)
                g["valid10"] = np.ones((256, 256), bool)
                g["valid20"] = np.ones((128, 128), bool)


@pytest.mark.parametrize("method", ["wald", "naive_l1", "reg_l1_nogain", "reg_l1", "shift_l1", "reg_hp", "b_sim",
                                    "b_sim_ft", "p6"])
def test_e2_methods_train_and_evaluate(tmp_path, method):
    _domain(tmp_path / "r", {"train": ["2020_KONZ_1__0_0", "2020_RMNP_1__0_0"], "val": ["2020_SRER_1__0_0"],
                             "test": ["2020_OAES_1__0_0", "2020_JORN_1__0_0"]}, ["y", "r_raw", "r_reg"])
    _domain(tmp_path / "s", {"train": ["2020_KONZ_1__0_1"], "val": ["2020_SRER_1__0_1"]}, ["y"])
    cfg = yaml.safe_load((ROOT / "configs" / "p6_e2.yaml").read_text())
    cfg.update(domain_r=str(tmp_path / "r"), domain_s=str(tmp_path / "s"), runs_dir=str(tmp_path / "runs"))
    cfg["training"].update(iterations=2, teacher_iterations=2, batch=2, crop=64, workers=0, log_every=1, warmup=1)
    out = E2.run(cfg, method, 0, torch.device("cpu"))
    t = pd.read_csv(out / "test_y.csv")
    assert set(t.site) == {"OAES", "JORN"} and {"hp_corr", "hp_gain"} <= set(t.columns)
    assert len(pd.read_csv(out / "test_r_reg.csv")) == len(t)
    assert pd.read_csv(out / "wald_test.csv").wald_psnr.notna().all()
