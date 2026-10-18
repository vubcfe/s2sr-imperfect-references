"""P6: fixed, seeded site split of SEN2NEON after E0 filtering (pre-E1 gate, step 2).

Tile filter (from E0 v4, results/p6_e0/tiles.csv): status ok, phase-correlation peak >= 0.1 (excludes tiles with
resampled 20 m bands and < 50% valid), cloud score <= the 75th percentile of reliable tiles.

Roles:
  test  sealed; evaluated on its |dt| <= 7 day tiles only (where E0 shows translation-only geometry)
  val   checkpoint/hyper-parameter selection; |dt| <= 7 day tiles
  train real-reference (R) training on |dt| <= 15 days (primary; 7 and 40 are sensitivities); domain S from all
        clean tiles of train sites
Only sites with >= MIN_TEST_TILES tiles at |dt| <= 7 days can be test or val sites (18 sites after filtering).

Selection: among RANDOM_DRAWS random assignments (seed fixed), keep the one that maximises the number of distinct
NEON domains in test, then minimises the land-cover distance between test and all eligible tiles, subject to at
most MAX_BIG_PER_ROLE of the BIG_SITES in test and in val, and a test share of eligible |dt| <= 7 tiles within
TEST_SHARE. The domain table is used only for stratification.

  python scripts/p6_make_split.py --e0 results/p6_e0/tiles.csv --out configs/p6_split.yaml
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SEED = 2026
RANDOM_DRAWS = 20000
N_TEST, N_VAL = 9, 4
MIN_TEST_TILES = 5
MAX_BIG_PER_ROLE = 2
TEST_SHARE = (0.35, 0.50)
DT_EVAL, DT_TRAIN = 7, 15

# NEON site -> ecoclimatic domain (stratification only)
DOMAIN = {
    "HARV": "D01", "BART": "D01", "SCBI": "D02", "SERC": "D02", "BLAN": "D02", "OSBS": "D03", "DSNY": "D03",
    "JERC": "D03", "GUAN": "D04", "LAJA": "D04", "CUPE": "D04", "GUIL": "D04", "UNDE": "D05", "STEI": "D05",
    "CHEQ": "D05", "KONZ": "D06", "UKFS": "D06", "MCDI": "D06", "MLBS": "D07", "GRSM": "D07", "TALL": "D08",
    "DELA": "D08", "LENO": "D08", "NOGP": "D09", "CPER": "D10", "STER": "D10", "RMNP": "D10", "ARIK": "D10",
    "CLBJ": "D11", "OAES": "D11", "PRIN": "D11", "YELL": "D12", "NIWO": "D13", "MOAB": "D13", "WLOU": "D13",
    "SRER": "D14", "JORN": "D14", "ONAQ": "D15", "REDB": "D15", "ABBY": "D16", "WREF": "D16", "MCRA": "D16",
    "SJER": "D17", "SOAP": "D17", "TEAK": "D17", "TOOL": "D18", "BONA": "D19", "DEJU": "D19", "HEAL": "D19",
}


def clean_tiles(tiles: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    rel = tiles[(tiles.status == "ok") & (tiles.g_peak >= 0.1)]
    q75 = float(rel.cloud_score.quantile(0.75))
    ok = rel[rel.cloud_score <= q75].copy()
    ok["land_cover"] = ok.land_cover.fillna("Unknown")
    ok["domain"] = ok.site.map(DOMAIN).fillna("unknown")
    return ok, q75


def choose(ok: pd.DataFrame, seed: int = SEED):
    per = ok.groupby("site").agg(n7=("dt_days", lambda v: int((v <= DT_EVAL).sum())))
    eligible = sorted(map(str, per.index[per.n7 >= MIN_TEST_TILES]))
    big = set(map(str, per.loc[eligible].n7.sort_values(ascending=False).index[:5]))
    ev = ok[ok.dt_days <= DT_EVAL]
    target_lc = ev[ev.site.isin(eligible)].land_cover.value_counts(normalize=True)
    total7 = per.loc[eligible].n7.sum()
    rng = np.random.default_rng(seed)
    best, best_key = None, None
    for _ in range(RANDOM_DRAWS):
        perm = rng.permutation(eligible)
        test, val = sorted(map(str, perm[:N_TEST])), sorted(map(str, perm[N_TEST:N_TEST + N_VAL]))
        if len(big & set(test)) > MAX_BIG_PER_ROLE or len(big & set(val)) > MAX_BIG_PER_ROLE:
            continue
        share = per.loc[test].n7.sum() / total7
        if not TEST_SHARE[0] <= share <= TEST_SHARE[1]:
            continue
        lc = ev[ev.site.isin(test)].land_cover.value_counts(normalize=True).reindex(target_lc.index, fill_value=0)
        key = (-len({DOMAIN.get(s) for s in test}), float((lc - target_lc).abs().sum()), tuple(test))
        if best_key is None or key < best_key:
            best, best_key = (test, val), key
    if best is None:
        raise ValueError(f"no assignment satisfies the rules ({len(eligible)} eligible sites)")
    return best, best_key, big, eligible


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e0", default="results/p6_e0/tiles.csv")
    ap.add_argument("--out", default="configs/p6_split.yaml")
    a = ap.parse_args()
    ok, q75 = clean_tiles(pd.read_csv(a.e0))
    (test, val), key, big, eligible = choose(ok)
    role = {str(s): "train" for s in ok.site.unique()}
    role.update({s: "test" for s in test})
    role.update({s: "val" for s in val})
    ok["role"] = ok.site.map(role)

    def tiles_for(r, dt):
        return sorted(map(str, ok[(ok.role == r) & (ok.dt_days <= dt)].tile))

    table = ok.groupby(["role", "site"]).agg(
        domain=("domain", "first"), n_clean=("tile", "size"),
        n_dt7=("dt_days", lambda v: int((v <= 7).sum())), n_dt15=("dt_days", lambda v: int((v <= 15).sum())),
        main_land_cover=("land_cover", lambda v: v.value_counts().index[0])).reset_index()
    cfg = dict(
        created="by scripts/p6_make_split.py", seed=SEED, e0_tiles=a.e0,
        filter=dict(status="ok", min_peak=0.1, max_cloud_score=round(q75, 4)),
        rules=dict(n_test=N_TEST, n_val=N_VAL, min_dt7_tiles_for_eval_sites=MIN_TEST_TILES,
                   big_sites=sorted(big), max_big_per_role=MAX_BIG_PER_ROLE, test_share=list(TEST_SHARE),
                   eval_dt_days=DT_EVAL, train_dt_days=DT_TRAIN),
        selection=dict(distinct_test_domains=-key[0], test_landcover_l1=round(key[1], 4),
                       eligible_sites=eligible),
        sites=dict(test=list(test), val=list(val), train=sorted(s for s, r in role.items() if r == "train")),
        tiles=dict(test_eval=tiles_for("test", DT_EVAL), val_eval=tiles_for("val", DT_EVAL),
                   train_R=tiles_for("train", DT_TRAIN), train_R_dt7=tiles_for("train", 7),
                   train_R_dt40=tiles_for("train", 40), train_S=tiles_for("train", 40)),
    )
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False, width=120))
    table.to_csv(out.with_suffix(".sites.csv"), index=False)
    counts = {k: len(v) for k, v in cfg["tiles"].items()}
    print(yaml.safe_dump(dict(sites=cfg["sites"], counts=counts, selection=cfg["selection"]), sort_keys=False))
    lc = ok[ok.dt_days <= DT_EVAL].groupby("role").land_cover.value_counts().unstack(fill_value=0)
    print(lc.to_string())
    print("test domains:", Counter(DOMAIN.get(s) for s in test))


if __name__ == "__main__":
    main()
