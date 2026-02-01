import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

path = "/Users/c24102394/Desktop/bivariate/AD_LON/AD_LON.csv"
df = pd.read_csv(path, sep="\t")

nc1 = float(df.loc[0,"nc1@p9"])
nc2 = float(df.loc[0,"nc2@p9"])
nc12 = float(df.loc[0,"nc12@p9"])

onlyA = nc1
onlyB = nc2
both = nc12
#scale = 4_000.0
scale = 1_450_000.0

plt.close("all")
fig, ax = plt.subplots(figsize=(1.8, 1.8), dpi=300)
ax.set_axis_off()

v = venn2(
     subsets=(onlyA*scale, onlyB, both), # onlyA*scale
     set_labels=("AD", "LON"),
     ax=ax
)

# v = venn2(
#     subsets=(onlyA, onlyB, both),
#     set_labels=("SCZ", "LON"),
#     ax=ax
# )

pA = v.get_patch_by_id("10")
pB = v.get_patch_by_id("01")
pAB = v.get_patch_by_id("11")

if pA is not None:
    pA.set_color("#67a9cf")
if pB is not None:
    pB.set_color("#ef8a62")
if pAB is not None:
    pAB.set_color("#bdbdbd")

for lbl in v.subset_labels:
    if lbl is not None:
        lbl.set_text("")

for p in [pA, pB, pAB]:
    if p is not None:
        p.set_edgecolor((1, 1, 1, 0.85))
        p.set_linewidth(6)
        p.set_alpha(0.60)

fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
plt.savefig("venns.png", bbox_inches="tight", transparent=True)
plt.show()
