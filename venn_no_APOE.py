import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

path = "/Users/c24102394/Desktop/Desktop - Cardiff University K427JYM6QV/bivariate/AD_SCZ/AD_SCZ.csv"
df = pd.read_csv(path, sep="\t")

nc1 = float(df.loc[0,"nc1@p9"])
nc2 = float(df.loc[0,"nc2@p9"])
nc12 = float(df.loc[0,"nc12@p9"])

onlyA = nc1
onlyB = nc2
both = nc12

pheno1_id = "AD"
pheno2_id = "SCZ"

colour_map = {
    "AD":  "#67a9cf",
    "SCZ": "#ef8a62",
    "LON": "#91cf60"
}

colA = colour_map.get(pheno1_id, "#cccccc")
colB = colour_map.get(pheno2_id, "#cccccc")
colAB = "#bdbdbd"

scale = 4_450_000.0
scale2 = 100.0
scale3 = 20.0

plt.close("all")
fig, ax = plt.subplots(figsize=(1.8, 1.8), dpi=300)
ax.set_axis_off()

v = venn2(
    subsets=(onlyA * scale, onlyB, both * scale2),
    set_labels=("AD", "SCZ"),
    ax=ax
)

pA = v.get_patch_by_id("10")
pB = v.get_patch_by_id("01")
pAB = v.get_patch_by_id("11")

if pA is not None:
    pA.set_color(colA)
if pB is not None:
    pB.set_color(colB)
if pAB is not None:
    pAB.set_color(colAB)

for lbl in v.subset_labels:
    if lbl is not None:
        lbl.set_text("")

tA = v.get_label_by_id("A")
tB = v.get_label_by_id("B")

if tA is not None:
    x, y = tA.get_position()
    ax.text(x - 0.08, y - 0.15, "(-APOE)", ha="center", va="top", fontsize=7)

if tB is not None:
    tB.set_fontsize(11)

for p in [pA, pB, pAB]:
    if p is not None:
        p.set_edgecolor((1, 1, 1, 0.85))
        p.set_linewidth(6)
        p.set_alpha(0.60)

fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
plt.savefig("/Users/c24102394/Desktop/venns.pdf", bbox_inches="tight", transparent=True, dpi=600)
plt.show()
