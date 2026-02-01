import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

path = "/Users/c24102394/Desktop/bivariate/AD_SCZ/AD_SCZ.csv"
df = pd.read_csv(path, sep="\t")
print(df.columns)

nc1 = float(df.loc[0,"nc1@p9"])
nc2 = float(df.loc[0,"nc2@p9"])
nc12 = float(df.loc[0,"nc12@p9"])
onlyA = nc1
onlyB = nc2
both = nc12
scale = 1450000.0
v = venn2(
    subsets=(onlyA*scale, onlyB, both),
    set_labels=("AD", "SCZ")
)

for lbl in v.subset_labels:
    if lbl is not None:
        lbl.set_text("")

plt.tight_layout()
plt.show()
