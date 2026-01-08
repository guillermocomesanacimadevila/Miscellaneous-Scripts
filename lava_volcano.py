#!/usr/bin/env python3
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch

centromeres = {
    1: 121535434,  2: 92326171,   3: 90504854,   4: 49660117,
    5: 46405641,   6: 58830166,   7: 58054331,   8: 43838887,
    9: 47367679,  10: 39254935,  11: 51644205,  12: 34856694,
    13: 16000000, 14: 16000000,  15: 17000000,  16: 35335801,
    17: 22263006, 18: 15460898,  19: 24681782,  20: 26369569,
    21: 11288129, 22: 13000000
}

results = pd.read_csv(
    "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/lava/ad_scz_age/LAVA_local_rg_bivariate.tsv",
    sep="\t"
)

results["neglog10_p"] = -np.log10(results["p"])

pairs = list(results["pair"].unique())
markers = {pairs[0]: "o", pairs[1]: "s", pairs[2]: "^"}

fdr_sig = results["q_fdr"] < 0.05
bonf_sig = results["sig_paper"].astype(bool)
sig = fdr_sig | bonf_sig

dark_blue = np.array([0.10, 0.15, 0.45, 1.0])
red = np.array([0.92, 0.45, 0.45, 1.0])

cap = 7
above = results["neglog10_p"] > cap

ycap = np.minimum(results["neglog10_p"].to_numpy(), cap)

norm = mcolors.Normalize(vmin=-1, vmax=1)
cmap = matplotlib.colormaps["RdYlBu_r"]
base = cmap(norm(results["rho"].to_numpy()))

t = 1.1 + 0.20 * (np.clip(ycap / cap, 0, 1) ** 0.5)
white = np.array([1, 1, 1, 1.0])
rgba = white * (1 - t)[:, None] + base * t[:, None]
rgba[:, 3] = 1.0

extreme = results.loc[results["neglog10_p"].idxmax()]

def arm_pos_label(row):
    mid = (row["start"] + row["stop"]) / 2
    arm = "p" if mid < centromeres[int(row["chr"])] else "q"
    return f"{int(row['chr'])}{arm}:{row['start']/1e6:.2f}–{row['stop']/1e6:.2f}"

fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)

for pair in pairs:
    idx_base = (results["pair"] == pair) & (~above)
    d = results.loc[idx_base]

    ax.scatter(
        d["rho"],
        np.minimum(d["neglog10_p"], cap),
        c=rgba[idx_base.to_numpy()],
        s=70,
        marker=markers.get(pair, "o"),
        edgecolors="none",
        label=pair.replace("AGE", "LON"),
        zorder=2
    )

    idx_fdr = (results["pair"] == pair) & fdr_sig & (~above)
    d_fdr = results.loc[idx_fdr]
    d_fdr_chr8 = d_fdr[d_fdr["chr"].astype(int) == 8]
    d_fdr_other = d_fdr[d_fdr["chr"].astype(int) != 8]

    if len(d_fdr_other) > 0:
        ax.scatter(
            d_fdr_other["rho"],
            np.minimum(d_fdr_other["neglog10_p"], cap),
            c=[dark_blue],
            s=85,
            marker=markers.get(pair, "o"),
            edgecolors="k",
            linewidths=1.0,
            zorder=6
        )

    if len(d_fdr_chr8) > 0:
        ax.scatter(
            d_fdr_chr8["rho"],
            np.minimum(d_fdr_chr8["neglog10_p"], cap),
            c=[dark_blue],
            s=130,
            marker="*",
            edgecolors="k",
            linewidths=1.0,
            zorder=7
        )

    idx_bonf = (results["pair"] == pair) & bonf_sig & (~above)
    d_bonf = results.loc[idx_bonf]
    d_bonf_chr8 = d_bonf[d_bonf["chr"].astype(int) == 8]
    d_bonf_other = d_bonf[d_bonf["chr"].astype(int) != 8]

    if len(d_bonf_other) > 0:
        ax.scatter(
            d_bonf_other["rho"],
            np.minimum(d_bonf_other["neglog10_p"], cap),
            c=[dark_blue],
            s=130,
            marker=markers.get(pair, "o"),
            edgecolors="k",
            linewidths=1.2,
            zorder=7
        )

    if len(d_bonf_chr8) > 0:
        ax.scatter(
            d_bonf_chr8["rho"],
            np.minimum(d_bonf_chr8["neglog10_p"], cap),
            c=[dark_blue],
            s=190,
            marker="*",
            edgecolors="k",
            linewidths=1.2,
            zorder=8
        )

for pair in pairs:
    idx_cap = (results["pair"] == pair) & above
    d_cap = results.loc[idx_cap]
    if len(d_cap) == 0:
        continue

    ax.scatter(
        d_cap["rho"],
        np.full(len(d_cap), cap),
        c=[dark_blue],
        s=140,
        marker=markers.get(pair, "o"),
        edgecolors="k",
        linewidths=1.2,
        zorder=50,
        clip_on=False
    )

ax.scatter(
    [extreme["rho"]],
    [cap],
    c=[red],
    s=170,
    marker=markers.get(extreme["pair"], "o"),
    edgecolors="k",
    linewidths=1.5,
    zorder=55,
    clip_on=False
)

ax.scatter(
    [extreme["rho"]],
    [cap],
    c="black",
    s=60,
    marker="*",
    edgecolors="none",
    zorder=56,
    clip_on=False
)

xspan = float(results["rho"].max() - results["rho"].min())
dx = 0.02 * xspan if np.isfinite(xspan) and xspan > 0 else 0.03
dy = 0.15
dy20 = 0.35
dy7below = 0.30
dy4up = 0.26

sig_to_label = results.loc[sig & (~above)].copy()
sig_to_label = sig_to_label.sort_values(["neglog10_p", "p"], ascending=[False, True]).head(25)

for _, row in sig_to_label.iterrows():
    ch = int(row["chr"])
    x = float(row["rho"])
    y = float(min(row["neglog10_p"], cap))
    txt = arm_pos_label(row)

    is_chr7p_special = (
        int(row["chr"]) == 7
        and (float(row["start"]) / 1e6) >= 36.49
        and (float(row["start"]) / 1e6) <= 36.53
        and (float(row["stop"]) / 1e6) >= 37.96
        and (float(row["stop"]) / 1e6) <= 38.00
    )

    is_chr7q_special = (
        int(row["chr"]) == 7
        and abs(float(row["start"]) / 1e6 - 93.69) < 0.03
        and abs(float(row["stop"]) / 1e6 - 95.17) < 0.03
    )

    if ch == 6:
        ox, oy = -dx, +dy
        ha, va = "right", "center"

    elif ch == 20:
        ox, oy = 0.0, -dy20
        ha, va = "center", "top"

    elif is_chr7p_special:
        ox, oy = 0.0, -dy7below
        ha, va = "center", "top"

    elif is_chr7q_special:
        ox, oy = +dx * 1.6, +dy * 1.35
        ha, va = "left", "center"

    elif ch == 4:
        ox, oy = 0.0, +dy4up
        ha, va = "center", "bottom"

    elif ch in (8, 13, 17):
        ox, oy = +dx, +dy
        ha, va = "left", "center"

    else:
        ha = "left" if x >= 0 else "right"
        ox = dx if x >= 0 else -dx
        oy = -dy if y > (cap - 0.6) else dy
        va = "center"

    ax.text(
        x + ox,
        y + oy,
        txt,
        fontsize=7,
        fontweight="bold",
        color=(0, 0, 0, 0.85),
        ha=ha,
        va=va,
        zorder=90,
        bbox=dict(
            boxstyle="round,pad=0.18",
            facecolor=(1, 1, 1, 0.75),
            edgecolor=(0, 0, 0, 0.0)
        )
    )

if above.any():
    axins = inset_axes(
        ax,
        width="34%",
        height="32%",
        loc="upper left",
        bbox_to_anchor=(0.40, 0.10, 0.65, 0.85),
        bbox_transform=ax.transAxes,
        borderpad=0
    )

    axins.set_facecolor((1, 1, 1, 0.88))
    axins.patch.set_path_effects([
        pe.SimplePatchShadow(offset=(2, -2), shadow_rgbFace=(0, 0, 0), alpha=0.20),
        pe.Normal()
    ])

    axins.scatter(
        extreme["rho"],
        extreme["neglog10_p"],
        c=[red],
        s=170,
        marker=markers.get(extreme["pair"], "o"),
        edgecolors="k",
        linewidths=1.2,
        zorder=10
    )

    axins.scatter(
        [extreme["rho"]],
        [extreme["neglog10_p"]],
        c="black",
        s=70,
        marker="*",
        edgecolors="none",
        zorder=11
    )

    axins.set_ylabel(r"$-\log_{10}(p)$", fontsize=9)
    axins.tick_params(axis="y", labelsize=8)
    axins.set_xticks([])

    axins.set_ylim(extreme["neglog10_p"] * 0.96, extreme["neglog10_p"] * 1.04)
    axins.set_xlim(extreme["rho"] - 0.05, extreme["rho"] + 0.05)

    axins.text(
        0.5, 0.32,
        arm_pos_label(extreme),
        transform=axins.transAxes,
        ha="center",
        va="center",
        fontsize=7,
        fontweight="bold",
        color=(0, 0, 0, 0.85)
    )

    for sp in axins.spines.values():
        sp.set_linewidth(0.8)

for spine in ax.spines.values():
    spine.set_zorder(0)

ax.set_xlabel(r"Local $r_g$")
ax.set_ylabel(r"$-\log_{10}(p)$")
ax.set_title("")

leg = ax.legend(
    title="Trait pair",
    loc="upper right",
    bbox_to_anchor=(0.99, 0.99),
    frameon=True,
    fancybox=True,
    framealpha=0.92,
    borderpad=0.65,
    labelspacing=0.45,
    handletextpad=0.6,
    handlelength=1.0,
    markerscale=1.2
)

leg.get_title().set_fontweight("bold")
leg.get_title().set_position((6, 0))

frame = leg.get_frame()
frame.set_facecolor((1, 1, 1, 0.92))
frame.set_edgecolor((0, 0, 0, 0.0))
frame.set_linewidth(0.0)
frame.set_path_effects([
    pe.SimplePatchShadow(offset=(2, -2), shadow_rgbFace=(0, 0, 0), alpha=0.18),
    pe.Normal()
])

ax.margins(y=0.01)
ax.set_ylim(0, cap)

plt.show()
