#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


# gene set ernrichment

ens_id = pd.read_csv("../ref/ensemble/mart_export.txt", sep="\t", dtype=str)
csv_path = Path("../results/CELLxGENE/AD_SCZ/CELLxGENE_blood_brain_normal_AD_SCZ.csv")

df = pd.read_csv(csv_path)

df["Mean Expression (raw counts)"] = pd.to_numeric(df["Mean Expression (raw counts)"], errors="coerce")
df["Pct Cells Expressing"] = pd.to_numeric(df["Pct Cells Expressing"], errors="coerce")
df = df.dropna(subset=["Mean Expression (raw counts)", "Pct Cells Expressing"])

mapping = ens_id.set_index("Gene stable ID")["Gene name"]

gene_id = df["Gene Symbol"].astype(str).str.replace(r"\.\d+$", "", regex=True)
m = gene_id.str.startswith("ENSG", na=False)

df.loc[m, "Gene Symbol"] = gene_id.loc[m].map(mapping).fillna(df.loc[m, "Gene Symbol"])
df = df.dropna(subset=["Gene Symbol"])

df = df.groupby(["Tissue", "Cell Type", "Gene Symbol"], as_index=False).agg({
    "Cell Count": "max",
    "Mean Expression (raw counts)": "mean",
    "Pct Cells Expressing": "mean"
})
