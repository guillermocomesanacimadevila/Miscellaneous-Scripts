#!/usr/bin/env python3
import os
import argparse
import cellxgene_census
import pandas as pd
from pathlib import Path
import numpy as np
import scipy.sparse as sp

census_version = "2025-11-08"
organism = "Homo sapiens"

def query_cells(disease_state: str, tissues_file: str, out_dir: str):
    out_dir = "../results/CELLxGENE/ref"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(os.path.expanduser(tissues_file), "r") as f:
        tissues = [ln.strip() for ln in f if ln.strip()]
    tissue_list = "[" + ",".join(f"'{t}'" for t in tissues) + "]"
    disease = disease_state.strip().strip("'").strip('"')
    value_filter = (
        "is_primary_data == True"
        f" and tissue_general in {tissue_list}"
        f" and disease == '{disease}'"
    )
    with cellxgene_census.open_soma(census_version=census_version) as census:
        obs = cellxgene_census.get_obs(
            census,
            organism=organism,
            value_filter=value_filter,
            column_names=["tissue_general", "cell_type"],
        )
    cell_types = (
        obs[["tissue_general", "cell_type"]]
        .drop_duplicates()
        .sort_values(["tissue_general", "cell_type"])
    )
    tag = "_".join(tissues)
    out_csv = out_dir / f"celltypes_{tag}_{disease}.csv"
    cell_types.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")
    print("n cell type groups:", len(cell_types))
    return cell_types

def grab_cell_types(cells_file: str,
                    genes_file: str,
                    tissues_file: str,
                    out_dir: str,
                    pheno1_id: str,
                    pheno2_id: str,
                    disease_state: str
                    ):
    out_dir = f"../results/CELLxGENE/{pheno1_id}_{pheno2_id}/"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    disease = disease_state.strip().strip("'").strip('"')
    with open(os.path.expanduser(tissues_file), "r") as f:
        tissues = [ln.strip() for ln in f if ln.strip()]
    tag = "_".join(tissues)
    cells_ref = f"../results/CELLxGENE/ref/celltypes_{tag}_{disease}.csv"
    pairs = pd.read_csv(cells_ref)
    use_all_cells = (cells_file is None) or (str(cells_file).strip().upper() == "ALL")
    if not use_all_cells:
        with open(os.path.expanduser(cells_file), "r") as f:
            cells = [ln.strip() for ln in f if ln.strip()]
        pairs = pairs[pairs["cell_type"].isin(cells)].reset_index(drop=True)

    tissue_list = "[" + ",".join(f"'{t}'" for t in tissues) + "]"
    with cellxgene_census.open_soma(census_version=census_version) as census:
        var = cellxgene_census.get_var(
            census,
            organism=organism,
            column_names=["feature_name"]
        )
    all_genes = var["feature_name"].astype(str).values
    use_all_genes = (genes_file is None) or (str(genes_file).strip().upper() == "ALL")
    if use_all_genes:
        genes = all_genes
    else:
        with open(os.path.expanduser(genes_file), "r") as f:
            genes = [ln.strip() for ln in f if ln.strip()]
    out_csv = out_dir / f"CELLxGENE_{tag}_{disease}_{pheno1_id}_{pheno2_id}.csv"
    pd.DataFrame(columns=[
        "Tissue",
        "Cell Type",
        "Cell Count",
        "Gene Symbol",
        "Mean Expression (raw counts)",
        "Pct Cells Expressing",
    ]).to_csv(out_csv, index=False)
    with cellxgene_census.open_soma(census_version=census_version) as census:
        for tissue, cell_type in pairs[["tissue_general", "cell_type"]].itertuples(index=False, name=None):
            filt = (
                "is_primary_data == True"
                f" and tissue_general in {tissue_list}"
                f" and disease == '{disease}'"
                f" and tissue_general == '{tissue}'"
                f" and cell_type == '{cell_type}'"
            )
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                obs_value_filter=filt,
                obs_column_names=["tissue_general", "cell_type"],
                var_column_names=["feature_name"],
            )
            X = sp.csr_matrix(adata.X)
            n_cells = adata.n_obs
            if use_all_genes:
                Xg = X
                genes_out = adata.var["feature_name"].astype(str).values
            else:
                varnames = adata.var["feature_name"].astype(str).values
                idx = np.array([np.where(varnames == g)[0][0] for g in genes], dtype=int)
                Xg = X[:, idx]
                genes_out = genes
            mean_expr = np.asarray(Xg.mean(axis=0)).ravel()
            pct_expr = np.asarray((Xg > 0).mean(axis=0)).ravel() * 100
            block = pd.DataFrame({
                "Tissue": tissue,
                "Cell Type": cell_type,
                "Cell Count": n_cells,
                "Gene Symbol": genes_out,
                "Mean Expression (raw counts)": mean_expr,
                "Pct Cells Expressing": pct_expr,
            })
            block.to_csv(out_csv, mode="a", header=False, index=False)
            print(tissue, "|", cell_type, "| cells =", n_cells)
    return str(out_csv)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells_file", default=None)
    parser.add_argument("--genes_file", default=None)
    parser.add_argument("--tissues_file", required=True)
    parser.add_argument("--out_dir")
    parser.add_argument("--pheno1", required=True)
    parser.add_argument("--pheno2", required=True)
    parser.add_argument("--disease_state", required=True)
    args = parser.parse_args()
    query_cells(
        disease_state=args.disease_state,
        tissues_file=args.tissues_file,
        out_dir=args.out_dir
    )
    grab_cell_types(
        cells_file=args.cells_file,
        genes_file=args.genes_file,
        tissues_file=args.tissues_file,
        out_dir=args.out_dir,
        pheno1_id=args.pheno1,
        pheno2_id=args.pheno2,
        disease_state=args.disease_state
    )


if __name__ == "__main__":
    main()