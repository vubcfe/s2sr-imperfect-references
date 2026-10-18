"""Figures for the P6 paper (PFG). Every number is read from the result files; nothing is typed in.

  python paper/make_figures.py   # from the repository root

Writes fig*.pdf (vector) and fig*.png (preview) next to this script. Widths follow PFG: 62 mm (one column) or
130 mm (two columns).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent  # repository root
E0 = BENCH / "results" / "p6_e0"
MM = 1 / 25.4
W1, W2 = 62 * MM, 130 * MM
# validated categorical slots 1-3 (dataviz reference palette, light), context gray, ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRAY, INK, MUTED = "#b5b3ad", "#0b0b0b", "#52514e"
TARGET = ["B05", "B06", "B07", "B8A", "B11", "B12"]

plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "grid.color": "#e6e5e1", "grid.linewidth": 0.5,
    "pdf.fonttype": 42, "font.family": "DejaVu Sans",
})


def save(fig, name):
    (HERE / "figures").mkdir(exist_ok=True)
    fig.savefig(HERE / "figures" / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(HERE / "figures" / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_geometry():
    t = pd.read_csv(E0 / "tiles.csv")
    t = t[t.status == "ok"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 48 * MM), gridspec_kw={"wspace": 0.45})
    bins = np.linspace(0, 1.5, 31)
    for col, c, lab in (("g_mag", BLUE, "guide bands (10 m)"), ("t_mag", ORANGE, "target bands (20 m)")):
        v = t[col].clip(upper=1.5)
        a.hist(v, bins=bins, histtype="step", color=c, linewidth=1.4, label=f"{lab}\nmedian {t[col].median():.2f} px")
    a.set_xlabel("tile shift (10 m px; last bin ≥ 1.5)")
    a.set_ylabel("tiles")
    a.set_ylim(0, None)
    a.legend(frameon=False, loc="upper right", handlelength=1.2)
    a.set_title("(a) per-tile misregistration", loc="left")
    edges = [0, 7, 15, 25, 40]
    t["dt_bin"] = pd.cut(t.dt_days.abs(), edges, include_lowest=True)
    g = t.groupby("dt_bin", observed=True)
    for col, c, lab in (("rms_tileshift_g", BLUE, "guide"), ("rms_tileshift_t", ORANGE, "target")):
        q = g[col].quantile([0.25, 0.5, 0.75]).unstack()
        x = np.arange(len(q)) + (-0.1 if lab == "guide" else 0.1)
        b.errorbar(x, q[0.5], yerr=[q[0.5] - q[0.25], q[0.75] - q[0.5]], fmt="o", ms=3.5, color=c, lw=1, capsize=0,
                   label=lab)
    b.set_xticks(range(len(edges) - 1), [f"{lo}–{hi}" for lo, hi in zip(edges[:-1], edges[1:])])
    b.set_xlabel("|Δt| S2 vs airborne (days)")
    b.set_ylabel("residual after translation (px)")
    b.set_title("(b) residual grows with |Δt| (median, IQR)", loc="left")
    b.legend(frameon=False, loc="upper left")
    b.grid(axis="y")
    save(fig, "fig2_geometry")


def fig_harmonisation():
    h = pd.read_csv(E0 / "harmonisation_regression_check.csv")
    h = h[h.band.isin(TARGET)]
    fig, axs = plt.subplots(2, 3, figsize=(W2, 84 * MM), sharex=True, sharey=True,
                            gridspec_kw={"wspace": 0.12, "hspace": 0.12})
    for ax, band in zip(axs.flat, TARGET):
        d = h[h.band == band]
        ax.plot([0, 1.05], [0, 1.05], color=MUTED, lw=0.7, ls="--", zorder=1)
        ax.scatter(d.R25, d.sd_ratio25, s=2, color=BLUE, alpha=0.25, linewidths=0, zorder=2)
        r = np.corrcoef(d.R25, d.sd_ratio25)[0, 1]
        ax.text(0.04, 0.96, band, transform=ax.transAxes, va="top", fontsize=7, fontweight="bold", color=INK)
        ax.text(0.96, 0.04, f"r = {r:.2f}\nmedian gap {d.gap.median():+.3f}", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6, color=INK)
        ax.set_xlim(0.2, 1.02)
        ax.set_ylim(0.2, 1.02)
        ax.grid(True)
    fig.supxlabel("R: correlation of Sentinel-2 (interpolated to 2.5 m) with the AVIRIS simulation", fontsize=7,
                  y=0.0)
    fig.supylabel("sd(SEN2NEON reference) / sd(Sentinel-2)", fontsize=7, x=0.04)
    n = h.groupby("band").size().min()
    fig.suptitle(f"Each point is one tile (n = {n} per band); dashed line: sd ratio = R", fontsize=6.5, x=0.12,
                 ha="left", y=0.99, color=MUTED)
    save(fig, "fig3_harmonisation")


def fig_rank_flip():
    p = BENCH / "runs" / "p6_e1r" / "analysis" / "scores_by_reference.csv"
    s = pd.read_csv(p)
    # band mean per tile, tile mean per model and seed, seed mean per model (as runs/p6_e1r/analysis/tau_ci.md)
    m = (s.groupby(["model", "seed", "ref", "tile"]).hp_psnr.mean().groupby(["model", "seed", "ref"]).mean()
         .groupby(["model", "ref"]).mean().unstack("ref"))
    ry, rh = m["y"].rank(ascending=False), m["y_h"].rank(ascending=False)
    fig, ax = plt.subplots(figsize=(W1 * 1.25, 80 * MM))
    hi = {"E1R/reg_l1_h": (ORANGE, "trained on harmonised ref."), "S0_clean/fdl": (AQUA, "FDL"),
          "S0_clean/shift_l1": (BLUE, "Shift-L1")}
    for k in m.index:
        if k in hi:
            continue
        ax.plot([0, 1], [ry[k], rh[k]], color=GRAY, lw=0.8, marker="o", ms=2.5)
    for k, (c, lab) in hi.items():
        ax.plot([0, 1], [ry[k], rh[k]], color=c, lw=1.6, marker="o", ms=3.5, zorder=3)
        ax.text(1.04, rh[k], lab, va="center", fontsize=6.3, color=INK)
    ax.set_xticks([0, 1], ["scored vs\ntrue y", "scored vs\nharmonised y_h"])
    ax.set_ylim(len(m) + 0.5, 0.5)
    ax.set_yticks([1, 5, 10, 15, 20])
    ax.set_ylabel("rank by HP-PSNR (1 = best)")
    ax.set_xlim(-0.15, 1.75)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    from scipy.stats import kendalltau
    tau = kendalltau(m["y"], m["y_h"])[0]
    ax.set_title(f"Same {len(m)} models, two references\n(Kendall τ = {tau:.2f})", loc="left")
    print("rank flip:", {k: (int(ry[k]), int(rh[k])) for k in hi}, "tau", round(tau, 3))
    save(fig, "fig4_rank_flip")


E2_LABELS = {"wald": "Wald", "naive_l1": "Naive-L1", "reg_l1_nogain": "Reg-L1, no gain", "b_sim": "B-sim",
             "b_sim_ft": "B-sim+FT", "shift_l1": "Shift-L1", "p6": "P6 (distillation)", "reg_l1": "Reg-L1",
             "reg_hp": "Reg-HP (λW = 1)", "reg_hp_w30": "Reg-HP (λW = 30)*"}


def fig_e2():
    import sys
    sys.path.insert(0, str(BENCH / "scripts"))
    import p6_e2_analyze as A
    df = pd.concat([A.load(BENCH / "runs" / "p6_e2"), A.load(BENCH / "runs" / "p6_e2_posthoc")])
    df = df[df.method.isin(E2_LABELS)]
    per_seed = df.groupby(["ref", "method", "seed"])[["hp_corr", "hp_psnr", "wald_psnr"]].mean()
    order = per_seed.loc["y"].groupby("method").hp_corr.mean().sort_values().index.tolist()
    fig, axs = plt.subplots(1, 3, figsize=(W2, 62 * MM), sharey=True, gridspec_kw={"wspace": 0.12})
    yy = np.arange(len(order))
    panels = [(axs[0], "y", "hp_corr", "(a) HP-corr"),
              (axs[2], "y", "wald_psnr", "(c) Wald fidelity (dB)")]
    for ax, ref, k, title in panels:
        g = per_seed.loc[ref].groupby("method")[k]
        mu, lo, hi = g.mean()[order], g.min()[order], g.max()[order]
        ax.hlines(yy, lo, hi, color=BLUE, lw=1.2)
        ax.plot(mu, yy, "o", color=BLUE, ms=3.5)
        ax.set_title(title, loc="left")
        ax.grid(axis="x")
    for ref, c, lab, dy in (("y", BLUE, "vs corrected y", -0.15), ("r_reg", ORANGE, "vs as-delivered r", 0.15)):
        mu = per_seed.loc[ref].groupby("method").hp_psnr.mean()[order]
        axs[1].plot(mu, yy + dy, "o", color=c, ms=3.5, label=lab)
    axs[1].set_title("(b) HP-PSNR (dB)", loc="left")
    axs[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, handletextpad=0.2,
                  columnspacing=0.8)
    axs[1].grid(axis="x")
    axs[0].set_yticks(yy, [E2_LABELS[m] for m in order])
    fig.text(0.0, -0.1, "Dots: mean over 3 seeds; bars: range over seeds. * exploratory, post hoc; weight chosen on "
             "validation sites.", fontsize=6, color=MUTED)
    save(fig, "fig5_e2")


def fig_qualitative():
    import json
    q = np.load(BENCH / "runs" / "p6_e2" / "analysis" / "qualitative.npz")
    info = json.loads((BENCH / "runs" / "p6_e2" / "analysis" / "qualitative.json").read_text())
    cols = [("x_bicubic", "S2 20 m\n(bicubic)"), ("r_raw", "reference\nas released"), ("y", "reference\ncorrected"),
            ("wald", "Wald"), ("naive_l1", "Naive-L1"), ("reg_l1", "Reg-L1"), ("reg_hp_w30", "Reg-HP\n(λW = 30)")]
    bands = [(3, "B8A"), (4, "B11"), (5, "B12")]
    fig, axs = plt.subplots(len(bands), len(cols), figsize=(W2, W2 * len(bands) / len(cols) + 4 * MM),
                            gridspec_kw={"wspace": 0.04, "hspace": 0.04})
    for i, (b, bname) in enumerate(bands):
        lo, hi = np.percentile(q["y"][b], [2, 98])  # one stretch per row, from the corrected reference
        for j, (k, title) in enumerate(cols):
            ax = axs[i, j]
            ax.imshow(q[k][b], cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if i == 0:
                ax.set_title(title, fontsize=6.3)
            if j == 0:
                ax.set_ylabel(bname, fontsize=7)
    fig.text(0.0, 0.02, f"Tile {info['tile']}, 640 m × 640 m crop; one grey stretch per band (2–98 % of the corrected "
             "reference).", fontsize=6, color=MUTED)
    save(fig, "fig6_qualitative")


def fig_overview():
    """One row per error: SEN2NEON step that causes it -> what we measured -> correction -> where it is tested."""
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(W2, 78 * MM))
    ax.set_xlim(-1.5, 131)
    ax.set_ylim(0, 78)
    ax.axis("off")
    xs, w = [0.3, 34.3, 68.3, 102.3], 27
    heads = ["SEN2NEON step", "Error found (E0)", "Correction", "Tested in"]
    hcol = [BLUE, ORANGE, AQUA, MUTED]
    for x, t, c in zip(xs, heads, hcol):
        ax.text(x + w / 2, 74, t, ha="center", va="center", fontsize=7.2, fontweight="bold", color=c)
    rows = [
        ("co-registration of\nthe airborne image\nto Sentinel-2", "tile shift,\nmedian 0.32 px;\ngrows with |Δt|",
         "shift the reference\nby −(tile shift)", "E1 (simulation)\nE2 (real data)"),
        ("OLS harmonisation\nS2 = m·AVIRIS + n,\nfitted values kept", "contrast reduced\nto ≈ R × S2\n(0.53–0.63 at 20 m)",
         "divide by one\ngain â per tile\nand band", "E1-R (simulation)\nE2 (real data)"),
        ("20 m bands of the\nreleased Sentinel-2\nproduct", "131 tiles\nresampled, not\nnative 20 m",
         "exclude these\ntiles", "all\nexperiments"),
    ]
    h, ys = 17, [52, 30, 8]
    for y, row in zip(ys, rows):
        for x, txt, c in zip(xs, row, hcol):
            ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.6", fc="white",
                                        ec=c, lw=1.2))
            ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=5.7, color=INK, linespacing=1.25)
        for x0, x1 in zip(xs[:-1], xs[1:]):
            ax.annotate("", xy=(x1 - 1.0, y + h / 2), xytext=(x0 + w + 1.0, y + h / 2),
                        arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.3", color=INK, lw=1.4))
    ax.text(65, 2.5, "Training study: 10 strategies in simulation, 9 on real data; HP-corr as the gain-invariant metric",
            ha="center", va="center", fontsize=6.2, color=MUTED, style="italic")
    save(fig, "fig1_overview")


if __name__ == "__main__":
    fig_overview()
    fig_qualitative()
    fig_geometry()
    fig_harmonisation()
    fig_rank_flip()
    fig_e2()
    print("ok")
