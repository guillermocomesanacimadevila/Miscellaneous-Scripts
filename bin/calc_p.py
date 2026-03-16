#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import norm

# p = 2(1 - phi(|Z|))
# z = rg / se
# p = 1 - erf(abs(z) / sqrt(2))

def add_pvals(df: pd.DataFrame, rg_col: str, se_col: str) -> pd.DataFrame:
    df = df.copy()
    df[rg_col] = pd.to_numeric(df[rg_col], errors="coerce")
    df[se_col] = pd.to_numeric(df[se_col], errors="coerce")
    df["Z"] = df[rg_col] / df[se_col]
    df["P"] = 2.0 * norm.sf(np.abs(df["Z"]))
    return df

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno1_prefix", required=True)
    p.add_argument("--pheno2_prefix", required=True)
    p.add_argument("--rg_column_name", required=True)
    p.add_argument("--se_column_name", required=True)
    p.add_argument("--results", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(
        args.results,
        sep=r"\s+",
        engine="python"
    )
    df = add_pvals(df, rg_col=args.rg_column_name, se_col=args.se_column_name)
    out_path = os.path.join(args.out_dir, f"{args.pheno1_prefix}-{args.pheno2_prefix}.cors.pvals.tsv")
    df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
