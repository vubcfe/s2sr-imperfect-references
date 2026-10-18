"""P6 E1 analysis: apply the pre-specified H2 rule (configs/p6_e1.yaml) to runs/p6_e1/<scenario>/<method>/seed0.

Per tile, endpoints are averaged over the six bands; per method x scenario we report means and tile-bootstrap 95%
CIs (2,000 resamples, same resamples for every method, so differences are paired).
  robust(method, scenario): HP-corr drop vs the method's own S0_clean value < 0.01 and CI upper bound < 0.02
  naive degrades: drop > 0.01 with CI excluding 0
  E2 loss: highest HP-corr on S2_realistic among losses whose edge-position error does not exceed Reg-L1's by
    more than 0.05 px; paired bootstrap tests of that loss against each other loss, Holm-adjusted.
  pCRD-L enters E2 only if it beats Shift-L1 (paired, Holm-adjusted p < 0.05 and positive difference).

  python scripts/p6_e1_analyze.py --runs runs/p6_e1 --split eval          # final (evaluate sites)
  python scripts/p6_e1_analyze.py --runs runs/p6_e1_tune --split tune     # weight tuning (tune sites)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_BOOT, SEED = 2000, 2026
LOSSES = ["naive_l1", "reg_l1", "shift_l1", "cobi", "fdl", "pisr", "p6_local", "p6_global"]


def load(runs: Path, split: str, seed: str = "seed0") -> pd.DataFrame:
    rows = []
    for f in runs.glob(f"*/*/{seed}/{split}.csv"):
        d = pd.read_csv(f)
        d = d.groupby("tile", as_index=False)[["hp_corr", "hp_psnr", "edge_error"]].mean()
        d["scenario"], d["method"], d["seed"] = f.parts[-4], f.parts[-3], f.parts[-2]
        rows.append(d)
    return pd.concat(rows) if rows else pd.DataFrame()


def boot_idx(n, rng):
    return rng.integers(0, n, (N_BOOT, n))


def holm(p: dict) -> dict:
    keys = sorted(p, key=p.get)
    out, m, run = {}, len(keys), 0.0
    for i, k in enumerate(keys):
        run = max(run, min(1.0, (m - i) * p[k]))
        out[k] = run
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs/p6_e1")
    ap.add_argument("--split", default="eval")
    a = ap.parse_args()
    df = load(Path(a.runs), a.split)
    if df.empty:
        raise SystemExit("no results")
    rng = np.random.default_rng(SEED)
    wide = {m: df.pivot_table(index="tile", columns=["scenario", "method"], values=m) for m in
            ("hp_corr", "hp_psnr", "edge_error")}
    tiles = wide["hp_corr"].index
    idx = boot_idx(len(tiles), rng)
    lines = [f"# E1 results ({a.split} sites, {len(tiles)} tiles)", ""]
    summ = []
    for (sc, me) in wide["hp_corr"].columns:
        v = wide["hp_corr"][(sc, me)].to_numpy()
        bs = np.nanmean(v[idx], axis=1)
        summ.append(dict(scenario=sc, method=me, hp_corr=np.nanmean(v), ci_lo=np.percentile(bs, 2.5),
                         ci_hi=np.percentile(bs, 97.5), hp_psnr=np.nanmean(wide["hp_psnr"][(sc, me)]),
                         edge_error=np.nanmean(wide["edge_error"][(sc, me)])))
    summ = pd.DataFrame(summ)
    lines += ["## Endpoints (mean over tiles; HP-corr with tile-bootstrap 95% CI)", "",
              summ.round(4).to_markdown(index=False), ""]

    # robustness relative to the method's own clean (S0) result
    rob = []
    for (sc, me) in wide["hp_corr"].columns:
        if sc == "S0_clean" or ("S0_clean", me) not in wide["hp_corr"].columns:
            continue
        d = (wide["hp_corr"][("S0_clean", me)] - wide["hp_corr"][(sc, me)]).to_numpy()
        bs = np.nanmean(d[idx], axis=1)
        drop, hi, lo = np.nanmean(d), np.percentile(bs, 97.5), np.percentile(bs, 2.5)
        rob.append(dict(scenario=sc, method=me, drop=drop, ci_lo=lo, ci_hi=hi,
                        robust=bool(drop < 0.01 and hi < 0.02), degrades=bool(drop > 0.01 and lo > 0)))
    if rob:
        lines += ["## Drop of HP-corr vs the method's own S0_clean (H2)", "",
                  pd.DataFrame(rob).round(4).to_markdown(index=False), ""]

    # E2 loss selection on S2_realistic
    sc = "S2_realistic"
    cand = [m for m in LOSSES if (sc, m) in wide["hp_corr"].columns]
    if cand and (sc, "reg_l1") in wide["edge_error"].columns:
        ee_reg = np.nanmean(wide["edge_error"][(sc, "reg_l1")])
        ok = [m for m in cand if np.nanmean(wide["edge_error"][(sc, m)]) <= ee_reg + 0.05]
        best = max(ok, key=lambda m: np.nanmean(wide["hp_corr"][(sc, m)]))
        pv, diff = {}, {}
        for m in cand:
            if m == best:
                continue
            d = (wide["hp_corr"][(sc, best)] - wide["hp_corr"][(sc, m)]).to_numpy()
            bs = np.nanmean(d[idx], axis=1)
            diff[m] = np.nanmean(d)
            pv[m] = min(1.0, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))
        ph = holm(pv)
        sel = pd.DataFrame([dict(vs=m, diff=diff[m], p=pv[m], p_holm=ph[m]) for m in pv])
        lines += [f"## E2 loss selection on {sc}", "",
                  f"Edge-error constraint: <= Reg-L1 + 0.05 px = {ee_reg + 0.05:.3f}; eligible: {ok}",
                  f"**Selected: {best}**", "", sel.round(4).to_markdown(index=False), ""]
        if ("S2_realistic", "p6_local") in wide["hp_corr"].columns and (sc, "shift_l1") in wide["hp_corr"].columns:
            d = (wide["hp_corr"][(sc, "p6_local")] - wide["hp_corr"][(sc, "shift_l1")]).to_numpy()
            bs = np.nanmean(d[idx], axis=1)
            p = min(1.0, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))
            lines.append(f"pCRD-L vs Shift-L1: diff {np.nanmean(d):+.4f}, p {p:.4f} -> pCRD-L "
                         f"{'enters' if (np.nanmean(d) > 0 and p * len(pv) < 0.05) else 'does not enter'} E2")
    # seed replication (methods with several seeds): mean, sd and paired per-seed differences
    allseed = load(Path(a.runs), a.split, seed="seed*")
    per = allseed.groupby(["scenario", "method", "seed"]).hp_corr.mean().reset_index()
    multi = per.groupby(["scenario", "method"]).filter(lambda g: len(g) > 1)
    if not multi.empty:
        rep = multi.groupby(["scenario", "method"]).hp_corr.agg(["mean", "std", "count"]).reset_index()
        lines += ["## Seed replication (mean HP-corr per seed; mean, sd over seeds)", "",
                  rep.round(5).to_markdown(index=False), ""]
        w = multi.pivot_table(index=["scenario", "seed"], columns="method", values="hp_corr")
        pairs = [("shift_l1", "reg_l1"), ("reg_l1", "naive_l1"), ("p6_global", "reg_l1"), ("shift_l1", "naive_l1")]
        rows = []
        for sc, g in w.groupby(level=0):
            for a_, b_ in pairs:
                if a_ in g and b_ in g:
                    dd = (g[a_] - g[b_]).dropna()
                    rows.append(dict(scenario=sc, diff=f"{a_} - {b_}", n_seeds=len(dd), mean=dd.mean(),
                                     sd=dd.std(), min=dd.min(), max=dd.max()))
        lines += ["Paired differences across seeds (same seed, same scenario):", "",
                  pd.DataFrame(rows).round(5).to_markdown(index=False), ""]
    text = "\n".join(lines)
    out = Path(a.runs) / f"summary_{a.split}.md"
    out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
