#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

results = pd.read_csv(
    "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/lava/ad_scz_age/LAVA_local_rg_bivariate.tsv",
    sep="\t"
)

results["neglog10_p"] = -np.log10(results["p"])

pairs = list(results["pair"].unique())
markers = {pairs[0]: "o", pairs[1]: "s", pairs[2]: "^"}

fdr_sig = results["q_fdr"] < 0.05
bonf_sig = results["sig_paper"].astype(bool)

dark_blue = np.array([0.10, 0.15, 0.45, 1.0])

cap = 7
above = results["neglog10_p"] > cap

y = np.minimum(results["neglog10_p"].to_numpy(), cap)

norm = mcolors.Normalize(vmin=-1, vmax=1)
cmap = cm.get_cmap("RdYlBu_r")
base = cmap(norm(results["rho"].to_numpy()))

t = 1.1 + 0.20 * (np.clip(y / cap, 0, 1) ** 0.5)
white = np.array([1, 1, 1, 1.0])
rgba = white * (1 - t)[:, None] + base * t[:, None]
rgba[:, 3] = 1.0

plt.figure(figsize=(7, 5))
ax = plt.gca()

for pair in pairs:
    idx_base = (results["pair"] == pair) & (~above)
    d = results.loc[idx_base]

    plt.scatter(
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
    plt.scatter(
        d_fdr["rho"],
        np.minimum(d_fdr["neglog10_p"], cap),
        c=[dark_blue],
        s=85,
        marker=markers.get(pair, "o"),
        edgecolors="k",
        linewidths=1.0,
        zorder=6
    )

    idx_bonf = (results["pair"] == pair) & bonf_sig & (~above)
    d_bonf = results.loc[idx_bonf]
    plt.scatter(
        d_bonf["rho"],
        np.minimum(d_bonf["neglog10_p"], cap),
        c=[dark_blue],
        s=130,
        marker=markers.get(pair, "o"),
        edgecolors="k",
        linewidths=1.2,
        zorder=7
    )

for pair in pairs:
    idx_cap = (results["pair"] == pair) & above
    d_cap = results.loc[idx_cap]
    plt.scatter(
        d_cap["rho"],
        np.full(len(d_cap), cap),
        c=[dark_blue],
        s=130,
        marker=markers.get(pair, "o"),
        edgecolors="k",
        linewidths=1.2,
        zorder=50,
        clip_on=False
    )

for spine in ax.spines.values():
    spine.set_zorder(0)

plt.xlabel(r"Local $r_g$")
plt.ylabel(r"$-\log_{10}(p)$")
plt.title("")

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

plt.margins(y=0.01)
plt.ylim(0, cap)
plt.tight_layout()
plt.show()
