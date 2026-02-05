import glob, re, sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import make_interp_spline
HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[2]
SCR_PATH = PROJECT_ROOT / "Experiments" / "Treatment-Effect" / "scr"
sys.path.insert(0, str(SCR_PATH))
from utils import load_ml_ready, build_treated_events

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "figure.dpi": 320,
    "savefig.dpi": 700,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

def panel_label(ax, label):
    ax.text(-0.06, 1.03, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

def clean(ax):
    ax.grid(False)
    ax.tick_params(axis="both", length=3, width=0.9, color="0.2", labelcolor="0.2")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.spines["left"].set_color("0.2")
    ax.spines["bottom"].set_color("0.2")

def load_horizon_series(base, model="lasso"):
    base = str(base).rstrip("/")
    rows = []
    for f in glob.glob(f"{base}/H*/t_learner_summary.csv"):
        m = re.search(r"/H(\d+)", f)
        if not m:
            continue
        H = int(m.group(1))
        df = pd.read_csv(f)
        if "model" not in df.columns:
            continue
        s = df[df["model"] == model]
        if len(s) == 0:
            continue
        r = s.iloc[0]
        need = ["ATE_DR_weighted", "ATE_DR_ci_low", "ATE_DR_ci_high"]
        if not all(c in df.columns for c in need):
            continue
        rows.append({"H": H, "ATE": float(r["ATE_DR_weighted"]), "lo": float(r["ATE_DR_ci_low"]), "hi": float(r["ATE_DR_ci_high"])})
    out = pd.DataFrame(rows)
    if len(out) == 0:
        found = glob.glob(f"{base}/H*/t_learner_summary.csv")
        raise RuntimeError(f"No usable horizons found in {base}. Files found: {len(found)}")
    return out.sort_values("H").reset_index(drop=True)

def treated_event_counts(df_raw, H_values, L, id_col, day_col, outcome_col, inj_codes=(1, 2)):
    out = []
    for H in H_values:
        t = build_treated_events(
            df_raw,
            L=L,
            H=int(H),
            id_col=id_col,
            day_col=day_col,
            treat_col="is_treated",
            treat_type_col="treatment_type",
            outcome_col=outcome_col,
            inj_codes=inj_codes,
        )
        out.append({"H": int(H), "n_treated_events": int(len(t))})
    return pd.DataFrame(out).sort_values("H").reset_index(drop=True)

def smooth_xy(x, y, k=3, n=1400):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 4:
        return x, y
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    k = min(k, max(1, len(x) - 1))
    xs = np.linspace(x.min(), x.max(), n)
    ys = make_interp_spline(x, y, k=k)(xs)
    return xs, ys

def plot_series(ax, h, ate, lo, hi, label, color, linestyle="-", n=1400, k=3):
    hs, ates = smooth_xy(h, ate, k=k, n=n)
    hs_lo, los = smooth_xy(h, lo, k=k, n=n)
    hs_hi, his = smooth_xy(h, hi, k=k, n=n)
    ax.fill_between(hs, los, his, color=color, alpha=0.12, linewidth=0, zorder=2)
    ax.plot(hs, ates, color=color, lw=2.9, linestyle=linestyle, label=label, zorder=4)
    ax.scatter(h, ate, color=color, s=16, zorder=5, linewidths=0)

L = 3
ID_COL = "id"
DAY_COL = "day"
MODEL = "lasso"
ADD_H0 = True
N_SMOOTH = 1600
K_SPLINE = 3

DATA = "../../Experiments/Treatment-Effect/Data/Symptomtrackingdata_csv-cleaned_with_vars_ml_ready.csv"
base_W = "../../Experiments/Treatment-Effect/outputs/Results_TLearner_W_t"
base_S = "../../Experiments/Treatment-Effect/outputs/Results_TLearner_S_t"

resW = load_horizon_series(base_W, model=MODEL)
resS = load_horizon_series(base_S, model=MODEL)

h_core = np.intersect1d(resW["H"].to_numpy(dtype=int), resS["H"].to_numpy(dtype=int))
h_core = np.sort(h_core)

df_raw = load_ml_ready(DATA)
cntT = treated_event_counts(
    df_raw=df_raw,
    H_values=h_core,
    L=L,
    id_col=ID_COL,
    day_col=DAY_COL,
    outcome_col="W_t",
    inj_codes=(1, 2),
)
nBars_core = cntT.set_index("H").reindex(h_core)["n_treated_events"].to_numpy(dtype=float)

resWc = resW.set_index("H").reindex(h_core).reset_index()
resSc = resS.set_index("H").reindex(h_core).reset_index()

h = h_core.astype(float)
nBars = nBars_core.astype(float)

if ADD_H0:
    h = np.insert(h, 0, 0.0)
    nBars = np.insert(nBars, 0, 0.0)
    resWc = pd.concat([pd.DataFrame([{"H":0,"ATE":0.0,"lo":0.0,"hi":0.0}]), resWc], ignore_index=True)
    resSc = pd.concat([pd.DataFrame([{"H":0,"ATE":0.0,"lo":0.0,"hi":0.0}]), resSc], ignore_index=True)

fig = plt.figure(figsize=(7.35, 4.10), facecolor="white")
ax = fig.add_subplot(111, facecolor="white")
panel_label(ax, "A")
ax.set_title("Causal effect dynamics (DR-ATE)", loc="left", pad=7)

ax2 = ax.twinx()
ax2.bar(h, nBars, width=0.82, alpha=0.075, color="0.1", zorder=0)
ax2.set_ylabel("Treated events (N)", color="0.25")
ax2.tick_params(axis="y", labelsize=7, length=3, width=0.8, colors="0.25")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.spines["bottom"].set_visible(False)
ax2.set_facecolor("none")

ax.axhline(0, lw=1.05, alpha=0.55, color="0.15", zorder=1)

plot_series(ax, h, resWc["ATE"].to_numpy(), resWc["lo"].to_numpy(), resWc["hi"].to_numpy(), r"$W_t$", "tab:blue", linestyle="-", n=N_SMOOTH, k=K_SPLINE)
plot_series(ax, h, resSc["ATE"].to_numpy(), resSc["lo"].to_numpy(), resSc["hi"].to_numpy(), r"$S_t$", "tab:orange", linestyle=(0, (5, 3)), n=N_SMOOTH, k=K_SPLINE)

if ADD_H0:
    ax.axvline(0, color="0.18", lw=1.05, linestyle=(0, (3, 2)), alpha=0.85, zorder=3)
    ax.scatter([0], [0], s=82, color="white", edgecolor="black", zorder=8)
    ax.scatter([0], [0], s=230, facecolors="none", edgecolors="black", linewidth=1.2, zorder=7)
    ax.annotate(
        "treatment (t0)",
        xy=(0, 0),
        xytext=(0.62, 0.08),
        textcoords="data",
        fontsize=8,
        color="0.25",
        arrowprops=dict(arrowstyle="-", lw=0.9, color="0.25", shrinkA=0, shrinkB=6),
    )

ax.set_xlabel("Horizon (days)", color="0.2")
ax.set_ylabel("ATE (Δ from baseline)", color="0.2")

xt = [0, 1, 3, 5, 7, 10, 14] if ADD_H0 else [1, 3, 5, 7, 10, 14]
ax.set_xticks(xt)
ax.yaxis.set_major_locator(MaxNLocator(6))

ymin = np.nanmin([resWc["lo"].min(), resSc["lo"].min()])
ymax = np.nanmax([resWc["hi"].max(), resSc["hi"].max()])
pad = 0.10 * (ymax - ymin) if np.isfinite(ymax - ymin) else 0.2
ax.set_ylim(ymin - pad, ymax + pad)

clean(ax)
ax2.grid(False)

leg = ax.legend(
    frameon=False,
    loc="upper left",
    handlelength=3.0,
    handletextpad=0.7,
    borderaxespad=0.5,
    labelspacing=0.35,
)
for t in leg.get_texts():
    t.set_color("0.05")

ax.text(
    0.01, 0.02,
    "Lines: DR-ATE. Ribbons: 95% CI. Bars: treated episode count.",
    transform=ax.transAxes, fontsize=8, color="0.35"
)

plt.tight_layout()
plt.savefig("aad.png")
plt.show()
