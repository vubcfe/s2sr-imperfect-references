"""Smoke test of the P6 training code on a tiny synthetic domain S (CPU, a few steps per method)."""
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("e1", ROOT / "scripts" / "p6_e1.py")
E1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E1)


def _fake(tmp):
    rng = np.random.default_rng(0)
    s = tmp / "s"
    s.mkdir()
    tiles = {"train": [f"2020_TRA{i}_1__0_{i}" for i in range(3)],
             "val": ["2020_CLBJ_1__0_0", "2020_NOGP_1__0_0"]}
    for role, names in tiles.items():
        with h5py.File(s / f"{role}.h5", "w") as f:
            for t in names:
                g = f.create_group(t)
                g["y"] = rng.integers(1000, 3000, (6, 128, 128)).astype(np.uint16)
                g["g"] = rng.integers(1000, 3000, (4, 128, 128)).astype(np.uint16)
                g["x"] = rng.integers(1000, 3000, (6, 64, 64)).astype(np.uint16)
                g["valid10"] = np.ones((128, 128), bool)
                g["valid20"] = np.ones((64, 64), bool)
    tg = tmp / "targets"
    tg.mkdir()
    with h5py.File(tg / "targets_S2_realistic.h5", "w") as f:
        for t in tiles["train"]:
            g = f.create_group(t)
            g["y_c"] = rng.integers(1000, 3000, (6, 128, 128)).astype(np.uint16)
            g["y_reg"] = rng.integers(1000, 3000, (6, 128, 128)).astype(np.uint16)
    return s, tg


@pytest.mark.parametrize("method", ["wald", "naive_l1", "reg_l1", "shift_l1", "cobi", "fdl", "pisr", "p6_local",
                                    "p6_global"])
def test_every_method_trains_and_evaluates(tmp_path, method):
    s, tg = _fake(tmp_path)
    cfg = yaml.safe_load((ROOT / "configs" / "p6_e1.yaml").read_text())
    cfg.update(domain_s=str(s), targets_dir=str(tg), runs_dir=str(tmp_path / "runs"))
    cfg["training"].update(batch=2, crop=64, workers=0, log_every=1, warmup=1)
    out = E1.run(cfg, "S2_realistic", method, 0, torch.device("cpu"), iterations=2)
    ev = pd.read_csv(out / "eval.csv")
    assert set(ev.site) == {"NOGP"} and len(ev) == 6 and ev.hp_corr.notna().all()
    assert pd.read_csv(out / "tune.csv").site.unique().tolist() == ["CLBJ"]


def test_interrupted_run_resumes_to_identical_weights(tmp_path, monkeypatch):
    from s2sr.p6 import losses, train as T
    s, tg = _fake(tmp_path)
    cfg = yaml.safe_load((ROOT / "configs" / "p6_e1.yaml").read_text())
    cfg.update(domain_s=str(s), targets_dir=str(tg), runs_dir=str(tmp_path / "runs"))
    cfg["training"].update(batch=2, crop=64, workers=0, log_every=1, warmup=1, resume_every=2)
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    p = E1.paths_for(cfg, "S2_realistic", 0)
    method = dict(cfg["methods"]["reg_l1"])
    full = T.train_student(method, cfg, p, mtf, tmp_path / "full", 0, torch.device("cpu"), iterations=4)
    calls = {"n": 0}
    orig = losses.l1_masked

    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] == 5:  # two l1 calls per step: interrupt at the start of step 3
            raise KeyboardInterrupt
        return orig(*a, **k)

    monkeypatch.setattr(T.L, "l1_masked", flaky)
    with pytest.raises(KeyboardInterrupt):
        T.train_student(method, cfg, p, mtf, tmp_path / "cut", 0, torch.device("cpu"), iterations=4)
    assert (tmp_path / "cut" / "resume.pt").exists()
    monkeypatch.setattr(T.L, "l1_masked", orig)
    cut = T.train_student(method, cfg, p, mtf, tmp_path / "cut", 0, torch.device("cpu"), iterations=4)
    a, b = torch.load(full), torch.load(cut)
    assert all(torch.allclose(a[k], b[k], atol=1e-6) for k in a)
    assert not (tmp_path / "cut" / "resume.pt").exists()


def test_corrupt_resume_file_is_moved_aside(tmp_path):
    from s2sr.p6 import train as T
    p = tmp_path / "resume.pt"
    p.write_bytes(b"")
    assert T._load_resume(p, "cpu") is None
    assert not p.exists() and (tmp_path / "resume.corrupt").exists()
    T._save_resume(p, step=3)
    assert T._load_resume(p, "cpu")["step"] == 3
