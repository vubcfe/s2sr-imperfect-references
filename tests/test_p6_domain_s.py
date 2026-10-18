"""Tests for P6 domain S (scripts/p6_build_domain_s.py)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("dom_s", ROOT / "scripts" / "p6_build_domain_s.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

CFG = yaml.safe_load((ROOT / "configs" / "p6_domain_s.yaml").read_text())
MTF = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]


def _scene(seed=0, flat=False):
    rng = np.random.default_rng(seed)
    if flat:
        return np.full((12, 512, 512), 0.2)
    n = 512
    fy, fx = np.fft.fftfreq(n)[:, None], np.fft.fftfreq(n)[None, :]
    t = np.real(np.fft.ifft2(np.fft.fft2(rng.standard_normal((12, n, n))) / (1 + (fy ** 2 + fx ** 2) / 0.01 ** 2)))
    return 0.2 + 0.03 * t / t.std()


def test_shapes_types_and_masks():
    z = _scene()
    ok = np.ones((512, 512), bool)
    ok[:20, :20] = False
    o = S.simulate(z, ok, CFG, MTF, np.random.default_rng(0))
    assert o["y"].shape == (6, 128, 128) and o["g"].shape == (4, 128, 128) and o["x"].shape == (6, 64, 64)
    assert o["y"].dtype == np.uint16 and o["valid10"].shape == (128, 128) and o["valid20"].shape == (64, 64)
    assert not o["valid10"][0, 0] and not o["valid20"][0, 0] and o["valid10"][-1, -1]
    assert np.all(np.abs(o["interband_shift"]) <= CFG["interband_shift_max_px20"])


def test_noise_level_matches_config():
    o = S.simulate(_scene(flat=True), np.ones((512, 512), bool), CFG, MTF, np.random.default_rng(1))
    q = CFG["quantisation"]
    for k, b in enumerate(S.TARGET_BANDS):
        assert abs((o["x"][k] * q).std() - CFG["noise_sigma"][b]) < 0.15 * CFG["noise_sigma"][b] + q
    assert np.all(o["y"] == round(0.2 / q))  # the target is noise-free


def test_per_tile_seed_is_order_independent_and_reproducible():
    a = S.tile_seed(2026, "2019_KONZ_3__1_2")
    assert a == S.tile_seed(2026, "2019_KONZ_3__1_2") != S.tile_seed(2026, "2019_KONZ_3__1_3")
    z, ok = _scene(), np.ones((512, 512), bool)
    o1 = S.simulate(z, ok, CFG, MTF, np.random.default_rng(a))
    o2 = S.simulate(z, ok, CFG, MTF, np.random.default_rng(a))
    assert np.array_equal(o1["x"], o2["x"]) and np.array_equal(o1["g"], o2["g"])


def test_x_is_close_to_lfc_operator_of_y():
    """x is simulated directly from 2.5 m; it must stay close to D(y) (LFC) but need not be identical."""
    o = S.simulate(_scene(), np.ones((512, 512), bool), {**CFG, "noise_scale": 0.0, "interband_shift_max_px20": 0.0},
                   MTF, np.random.default_rng(0))
    q = CFG["quantisation"]
    gap = S.lfc_gap(o["y"] * q, o["x"] * q, o["valid20"], MTF)
    assert np.all(gap < 0.1 * 0.03)  # well below the scene contrast


def test_build_refuses_test_sites(tmp_path, monkeypatch):
    split = yaml.safe_load((ROOT / "configs" / "p6_split.yaml").read_text())
    bad = split["tiles"]["test_eval"][:1]
    fake = dict(split, tiles=dict(split["tiles"], train_S=bad))
    p = tmp_path / "split.yaml"
    p.write_text(yaml.safe_dump(fake))
    monkeypatch.setattr(S, "ROOT", Path("/"))
    cfg = dict(CFG, split=str(p), e0_tiles=None)
    fake["e0_tiles"] = str(ROOT / split["e0_tiles"])
    p.write_text(yaml.safe_dump(fake))
    (tmp_path / "metadata.csv").write_text("id,lr,hr_2_5m_path\n")
    with pytest.raises(ValueError, match="disallowed sites"):
        S.build(tmp_path, tmp_path / "out", cfg, MTF, workers=1)
