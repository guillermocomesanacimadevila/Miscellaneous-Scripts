#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def overlap(pheno1_id: str, pheno2_id: str, pheno1_path: str, pheno2_path: str, out_path: str):
    pheno1_path = Path(pheno1_path)
    pheno2_path = Path(pheno2_path)
    
    if pheno1_path.stat().st_size == 0 or pheno1_path.read_text().strip() == "":
        pheno1_df = pd.DataFrame(columns=["SNP"])
    else:
        pheno1_df = pd.read_csv(pheno1_path, sep="\t")

    if pheno2_path.stat().st_size == 0 or pheno2_path.read_text().strip() == "":
        pheno2_df = pd.DataFrame(columns=["SNP"])
    else:
        pheno2_df = pd.read_csv(pheno2_path, sep="\t")

    if pheno1_df.empty or pheno2_df.empty:
        merged = pd.DataFrame(columns=["SNP"])
    else:
        merged = pheno1_df.merge(
            pheno2_df,
            how="inner",
            on="SNP",
            suffixes=(f"_{pheno1_id}", f"_{pheno2_id}")
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno1_id", required=True)
    ap.add_argument("--pheno2_id", required=True)
    ap.add_argument("--pheno1_path", required=True)
    ap.add_argument("--pheno2_path", required=True)
    ap.add_argument("--out_path", required=True)
    args = ap.parse_args()
    overlap(args.pheno1_id, args.pheno2_id, args.pheno1_path, args.pheno2_path, args.out_path)

if __name__ == "__main__":
    main()