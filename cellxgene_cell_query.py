#!/usr/bin/env python3
import os 
import cellxgene_census
import pandas as pd

census_version = "2025-11-08"
organism = "Homo sapiens"

with cellxgene_census.open_soma(census_version=census_version) as census:
    obs = cellxgene_census.get_obs(
        census,
        organism=organism,
        value_filter="is_primary_data == True and tissue_general in ['blood','brain'] and disease == 'normal'",
        column_names=["tissue_general", "cell_type"]
    )

cell_types = (
    obs[["tissue_general","cell_type"]]
    .drop_duplicates()
    .sort_values(["tissue_general","cell_type"])
)

cell_types.to_csv("celltypes_blood_brain_normal.csv", index=False)
print("Saved celltypes_blood_brain_normal.csv")
print("n cell type groups:", len(cell_types))
print(cell_types.head())
