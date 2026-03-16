#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import os
import argparse

def define_loci(pheno1: str,
                pheno2: str,
                pheno1_prefix: str,
                pheno2_prefix: str,
                clump_path: str,
                out_dir: str,
                window: int):

    pheno1_df = pd.read_csv(pheno1, sep="\t")
    pheno2_df = pd.read_csv(pheno2, sep="\t")
    pheno1_df.rename(columns={"POS": "BP"}, inplace=True)
    pheno2_df.rename(columns={"POS": "BP"}, inplace=True)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    clump_path = Path(clump_path)
    rows = []
    for dir in os.listdir(clump_path):
        full_path = clump_path / dir
        if os.path.isdir(full_path) and dir.startswith("locus_"):
            locus_coords = "_".join(dir.split("_")[-3:])
            lead_file = full_path / "lead_snps.tsv"
            lead_df = pd.read_csv(lead_file, sep="\t")
            lead_snp = lead_df["lead_snp"].iloc[0]
            rows.append({
                "trait1": pheno1_prefix,
                "trait2": pheno2_prefix,
                "pair": f"{pheno1_prefix}_{pheno2_prefix}",
                "locus": dir,
                "coords": locus_coords,
                "lead_snp": lead_snp
            })
    lead_snps_df = pd.DataFrame(rows, columns=["trait1", "trait2", "pair", "locus", "coords", "lead_snp"])
    lead_snps_df = lead_snps_df.sort_values("coords").reset_index(drop=True)
    lead_snps_df.to_csv(out / "lead_snps.tsv", sep="\t", index=False)
    for _, r in lead_snps_df.iterrows():
        snp = r["lead_snp"]
        lead1 = pheno1_df.loc[pheno1_df["SNP"] == snp].iloc[0]
        chr_ = int(lead1["CHR"])
        bp = int(lead1["BP"])
        start = max(0, bp - window)
        end = bp + window
        pheno1_loc = pheno1_df.loc[(pheno1_df["CHR"] == chr_) & (pheno1_df["BP"] >= start) & (pheno1_df["BP"] <= end)].copy()
        pheno2_loc = pheno2_df.loc[(pheno2_df["CHR"] == chr_) & (pheno2_df["BP"] >= start) & (pheno2_df["BP"] <= end)].copy()
        common = pheno1_loc.merge(pheno2_loc[["SNP"]], on="SNP", how="inner")
        pheno1_common = pheno1_loc.loc[pheno1_loc["SNP"].isin(common["SNP"])].copy()
        pheno2_common = pheno2_loc.loc[pheno2_loc["SNP"].isin(common["SNP"])].copy()
        pheno1_common = pheno1_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        pheno2_common = pheno2_common.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)
        locus_dir = out / f"locus_chr{chr_}_{start}_{end}"
        locus_dir.mkdir(parents=True, exist_ok=True)
        pheno1_common.to_csv(locus_dir / f"gwas_{pheno1_prefix}.ldgwas.tsv", sep="\t", index=False)
        pheno2_common.to_csv(locus_dir / f"gwas_{pheno2_prefix}.ldgwas.tsv", sep="\t", index=False)
        pd.DataFrame({"SNP": common["SNP"]}).sort_values("SNP").to_csv(locus_dir / "common_snps.tsv", sep="\t", index=False)
    return lead_snps_df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno1", required=True)
    ap.add_argument("--pheno2", required=True)
    ap.add_argument("--pheno1_prefix", required=True)
    ap.add_argument("--pheno2_prefix", required=True)
    ap.add_argument("--clump_path", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--window", type=int, default=500000)
    args = ap.parse_args()
    define_loci(args.pheno1,
                args.pheno2,
                args.pheno1_prefix,
                args.pheno2_prefix,
                args.clump_path,
                args.out_dir,
                args.window)

if __name__ == "__main__":
    main()