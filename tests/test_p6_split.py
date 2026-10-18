"""Tests for the P6 site split (scripts/p6_make_split.py)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("split", ROOT / "scripts" / "p6_make_split.py")
SPLIT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SPLIT)


def _fake_tiles(seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i, site in enumerate(sorted(SPLIT.DOMAIN)[:30]):
        for k in range(int(rng.integers(20, 60))):
            rows.append(dict(tile=f"2020_{site}_1__{k}_0", site=site, status="ok", g_peak=0.6,
                             cloud_score=float(rng.uniform(0.86, 0.94)), dt_days=float(rng.integers(0, 41)),
                             land_cover=str(rng.choice(["Rural", "Forest", "Water"]))))
    return pd.DataFrame(rows)


def test_choose_is_deterministic_and_respects_rules(monkeypatch):
    monkeypatch.setattr(SPLIT, "RANDOM_DRAWS", 500)
    ok, _ = SPLIT.clean_tiles(_fake_tiles())
    (t1, v1), k1, big, eligible = SPLIT.choose(ok)
    (t2, v2), k2, _, _ = SPLIT.choose(ok)
    assert (t1, v1, k1) == (t2, v2, k2)
    assert not set(t1) & set(v1)
    assert set(t1) | set(v1) <= set(eligible)
    assert len(big & set(t1)) <= SPLIT.MAX_BIG_PER_ROLE and len(big & set(v1)) <= SPLIT.MAX_BIG_PER_ROLE


@pytest.mark.skipif(not (ROOT / "configs" / "p6_split.yaml").exists(), reason="split not generated")
def test_generated_split_has_no_leakage():
    cfg = yaml.safe_load((ROOT / "configs" / "p6_split.yaml").read_text())
    sites = cfg["sites"]
    assert not set(sites["test"]) & set(sites["val"])
    assert not (set(sites["test"]) | set(sites["val"])) & set(sites["train"])
    site_of = lambda t: t.split("_")[1]  # noqa: E731
    for key in ("train_R", "train_R_dt7", "train_R_dt40", "train_S"):
        assert {site_of(t) for t in cfg["tiles"][key]} <= set(sites["train"])
    assert {site_of(t) for t in cfg["tiles"]["test_eval"]} <= set(sites["test"])
    assert {site_of(t) for t in cfg["tiles"]["val_eval"]} <= set(sites["val"])
