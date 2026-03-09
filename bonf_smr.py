#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from statsmodels.stats.multitest import multipletests

# trait_eSMR.merged.tsv
# FDR correct locus-level SMR
# grab coords from loci => map onto SMR outputs for each trait
# count genes.unique() in each - or number of tests, correct over them
# check HEIDI

base = Path("../outputs/SMR")
defined_loci = Path("../outputs/defined_loci")

def locus_specific_smr_correction(pheno1_prefix: str, pheno2_prefix: str, out_dir: str, qtl_type: str):
    gene_map = pd.read_csv("../ref/ensemble/mart_export.txt", sep="\t", dtype=str)
    gene_map = gene_map[["Gene stable ID", "Gene name"]].drop_duplicates()
    gene_map.columns = ["gene_id_clean", "gene_name"]
    out_dir = Path(out_dir) / f"per_locus/{pheno1_prefix}_{pheno2_prefix}"
    os.makedirs(out_dir, exist_ok=True)
    loci_dir = Path(defined_loci) / f"{pheno1_prefix}_{pheno2_prefix}"
    rows = []
    for dir in os.listdir(loci_dir):
        full_path = loci_dir / dir
        if os.path.isdir(full_path) and dir.startswith("locus_"):
            locus_coords = "_".join(dir.split("_")[-3:])
            rows.append({
                "trait1": pheno1_prefix,
                "trait2": pheno2_prefix,
                "pair": f"{pheno1_prefix}_{pheno2_prefix}",
                "coords": locus_coords,
                "locus": dir
            })

    loci = pd.DataFrame(rows)
    dfs = {}
    for trait in [pheno1_prefix, pheno2_prefix]:
        dfs[trait] = pd.read_csv(
            f"{base}/smr_outputs/{qtl_type}/{trait}/trait_eSMR.merged.tsv",
            sep="\t"
        )

    df_pheno1 = dfs[pheno1_prefix].copy()
    df_pheno2 = dfs[pheno2_prefix].copy()
    df_pheno1["CHR"] = df_pheno1["ProbeChr"].astype(str).str.replace("^chr", "", regex=True)
    df_pheno2["CHR"] = df_pheno2["ProbeChr"].astype(str).str.replace("^chr", "", regex=True)
    df_pheno1["BP"] = pd.to_numeric(df_pheno1["Probe_bp"], errors="coerce")
    df_pheno2["BP"] = pd.to_numeric(df_pheno2["Probe_bp"], errors="coerce")
    locus_hits = []
    for _, row in loci.iterrows():
        chrom, start, end = row["coords"].split("_")
        chrom = chrom.replace("chr", "")
        start = int(start)
        end = int(end)
        df1_locus = df_pheno1[
            (df_pheno1["CHR"] == chrom) &
            (df_pheno1["BP"] >= start) &
            (df_pheno1["BP"] <= end)
        ].copy()

        df2_locus = df_pheno2[
            (df_pheno2["CHR"] == chrom) &
            (df_pheno2["BP"] >= start) &
            (df_pheno2["BP"] <= end)
        ].copy()

        locus_hits.append({
            "locus": row["locus"],
            "coords": row["coords"],
            pheno1_prefix: df1_locus,
            pheno2_prefix: df2_locus
        })

        p_bonf_trait1 = 0.05 / df1_locus.shape[0] if df1_locus.shape[0] > 0 else np.nan
        p_bonf_trait2 = 0.05 / df2_locus.shape[0] if df2_locus.shape[0] > 0 else np.nan
        sig_trait1 = df1_locus[df1_locus["p_SMR"] < p_bonf_trait1].copy() if df1_locus.shape[0] > 0 else pd.DataFrame()
        sig_trait2 = df2_locus[df2_locus["p_SMR"] < p_bonf_trait2].copy() if df2_locus.shape[0] > 0 else pd.DataFrame()

        if not sig_trait1.empty:
            sig_trait1["gene_id_clean"] = sig_trait1["gene_id"].astype(str).str.split(".").str[0]
            sig_trait1 = sig_trait1.merge(gene_map, on="gene_id_clean", how="left")

        if not sig_trait2.empty:
            sig_trait2["gene_id_clean"] = sig_trait2["gene_id"].astype(str).str.split(".").str[0]
            sig_trait2 = sig_trait2.merge(gene_map, on="gene_id_clean", how="left")

        print(f"\n{row['locus']} | {row['coords']}")
        print(f"{pheno1_prefix} n_tests: {df1_locus.shape[0]}")
        print(f"{pheno1_prefix} Bonf threshold: {p_bonf_trait1}")
        if not sig_trait1.empty:
            print(f"{pheno1_prefix} SNPs surpassing Bonf:")
            print(
                sig_trait1[
                    ["topSNP", "gene_name", "gene_id_clean", "probeID", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
                ].to_string(index=False)
            )
        else:
            print(f"{pheno1_prefix} SNPs surpassing Bonf: none")

        print(f"{pheno2_prefix} n_tests: {df2_locus.shape[0]}")
        print(f"{pheno2_prefix} Bonf threshold: {p_bonf_trait2}")
        if not sig_trait2.empty:
            print(f"{pheno2_prefix} SNPs surpassing Bonf:")
            print(
                sig_trait2[
                    ["topSNP", "gene_name", "gene_id_clean", "probeID", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
                ].to_string(index=False)
            )
        else:
            print(f"{pheno2_prefix} SNPs surpassing Bonf: none")
    return locus_hits

if __name__ == "__main__":
    locus_specific_smr_correction(
        pheno1_prefix="AD",
        pheno2_prefix="SCZ",
        out_dir="../results/SMR",
        qtl_type="eQTLs"
    )
