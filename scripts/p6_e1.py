"""P6 experiment E1 driver (configs/p6_e1.yaml).

  python scripts/p6_e1.py targets                      # corrupted/registered targets for every scenario
  python scripts/p6_e1.py run --scenario S2_realistic --method p6_global --seed 0 [--iterations N] [--gpu 0]
  python scripts/p6_e1.py smoke --gpu 0                # short runs of every method, reports seconds per step

`run` trains (teacher first when needed, cached per scenario/seed), then evaluates on the `evaluate` sites and on
the `tune` sites of domain-S validation, writing runs/p6_e1/<scenario>/<method>/seed<k>/{eval,tune}.csv.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import h5py
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s2sr.p6 import train as T  # noqa: E402


def load_cfg(path):
    return yaml.safe_load((ROOT / path).read_text())


def paths_for(cfg, scenario, seed, method_name=None):
    s = ROOT / cfg["domain_s"]
    with h5py.File(s / "train.h5", "r") as f:
        train_tiles = sorted(f)
    with h5py.File(s / "val.h5", "r") as f:
        val_tiles = sorted(f)
    runs = ROOT / cfg["runs_dir"] / scenario
    return dict(s_train=s / "train.h5", s_val=s / "val.h5", train_tiles=train_tiles,
                targets=ROOT / cfg["targets_dir"] / f"targets_{scenario}.h5",
                teacher=runs / "_teacher" / f"seed{seed}" / "teacher.pt",
                wald_model=runs / "wald" / f"seed{seed}" / "student.pt",
                tune_tiles=[t for t in val_tiles if t.split("_")[1] in cfg["val_sites"]["tune"]],
                eval_tiles=[t for t in val_tiles if t.split("_")[1] in cfg["val_sites"]["evaluate"]])


def run(cfg, scenario, method_name, seed, device, iterations=None, out=None, weights=None):
    mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
    method = dict(cfg["methods"][method_name])
    if weights:
        method["weights"] = {**method.get("weights", {}), **weights}
    p = paths_for(cfg, scenario, seed)
    out = out or ROOT / cfg["runs_dir"] / scenario / method_name / f"seed{seed}"
    out.mkdir(parents=True, exist_ok=True)
    if (out / "run.json").exists():
        return out  # finished earlier (queue restart)
    t0 = time.time()
    if method.get("teacher"):
        T.train_teacher(method["target"], cfg, p, mtf, p["teacher"].parent, seed, device, iterations)
    if "fdl" in method.get("weights", {}) and not p["wald_model"].exists():
        run(cfg, scenario, "wald", seed, device, iterations)
    ck = T.train_student(method, cfg, p, mtf, out, seed, device, iterations)
    T.evaluate(ck, p, mtf, p["eval_tiles"], out / "eval.csv", device)
    T.evaluate(ck, p, mtf, p["tune_tiles"], out / "tune.csv", device)
    T.save_json(out / "run.json", dict(scenario=scenario, method=method_name, spec=method, seed=seed,
                                       iterations=iterations or cfg["training"]["iterations"],
                                       minutes=round((time.time() - t0) / 60, 1), frozen=cfg.get("frozen")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["targets", "run", "smoke", "curve", "teacher"])
    ap.add_argument("--config", default="configs/p6_e1.yaml")
    ap.add_argument("--scenario", default="S2_realistic")
    ap.add_argument("--method", default="reg_l1")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--iterations", type=int, default=None)
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--weights", default=None, help='JSON dict overriding auxiliary weights, e.g. {"kd": 1e-3}')
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cfg = load_cfg(a.config)
    device = torch.device(f"cuda:{a.gpu}" if torch.cuda.is_available() else "cpu")
    if a.cmd == "targets":
        subprocess.run([sys.executable, str(ROOT / "scripts" / "p6_make_e1_targets.py"), "--config", a.config,
                        "--s", cfg["domain_s"], "--out", cfg["targets_dir"]], check=True, cwd=ROOT)
    elif a.cmd == "run":
        run(cfg, a.scenario, a.method, a.seed, device, a.iterations,
            Path(a.out) if a.out else None, json.loads(a.weights) if a.weights else None)
    elif a.cmd == "teacher":  # train (or resume) the scenario's teacher only
        mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
        p = paths_for(cfg, a.scenario, a.seed)
        T.train_teacher(cfg["methods"]["p6_global"]["target"], cfg, p, mtf, p["teacher"].parent, a.seed, device)
    elif a.cmd == "curve":  # learning curve to fix the step budget before freezing the config
        steps = [5000, 10000, 20000, 30000]
        cfg["training"]["save_steps"] = steps
        mtf = yaml.safe_load((ROOT / "configs" / "mtf.yaml").read_text())["mtf_nyquist"]
        out = ROOT / "runs" / "p6_e1_curve" / a.scenario / a.method
        run(cfg, a.scenario, a.method, a.seed, device, iterations=max(steps), out=out)
        p = paths_for(cfg, a.scenario, a.seed)
        import pandas as pd
        rows = []
        for st in steps:
            ck = out / (f"student_{st}.pt")
            df = T.evaluate(ck, p, mtf, p["tune_tiles"], out / f"tune_{st}.csv", device)
            rows.append(dict(step=st, hp_corr=df.hp_corr.mean(), hp_psnr=df.hp_psnr.mean(),
                             edge_error=df.edge_error.mean()))
        pd.DataFrame(rows).to_csv(out / "curve.csv", index=False)
        print(pd.DataFrame(rows).to_string(index=False), flush=True)
    elif a.cmd == "smoke":
        res = {}
        for m in cfg["methods"]:
            t0 = time.time()
            out = ROOT / "runs" / "p6_e1_smoke" / m
            run(cfg, "S2_realistic", m, 99, device, iterations=200, out=out)
            res[m] = round(time.time() - t0, 1)
            print(m, res[m], "s for 200 steps + eval", flush=True)
        (ROOT / "runs" / "p6_e1_smoke" / "timing.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
