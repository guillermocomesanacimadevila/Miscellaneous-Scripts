#!/usr/bin/env python3
import os
import gzip
import pandas as pd
import argparse
from pathlib import Path

# compute z-score
# change POS to BP
# gzip df

def reformat(sumstats, beta_col:str, se_col:str, pos_col:str, pheno_id:str, out_dir) -> None:
    df = pd.read_csv(sumstats, sep="\t", engine="python")
    df = df.copy()
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    df["Z"] = df[beta_col] / df[se_col]
    df["CHR"] = df["CHR"].astype(int)
    df[pos_col] = df[pos_col].astype(int)
    df.rename(
        columns={
            pos_col: "BP"
        },
        inplace=True
    )
    df.to_csv(
        out_dir / f"{pheno_id}_mixer_ready.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip"
    )

def main():
    args = argparse.ArgumentParser()
    args.add_argument("--sumstats", required=True)
    args.add_argument("--beta_col", required=True)
    args.add_argument("--se_col", required=True)
    args.add_argument("--pos_col", required=True)
    args.add_argument("--pheno_id", required=True)
    args.add_argument("--out_dir", required=True)
    args = args.parse_args()
    beta_col = args.beta_col
    se_col = args.se_col
    pos_col = args.pos_col
    pheno_id = args.pheno_id
    out_dir = args.out_dir
    reformat(args.sumstats, beta_col, se_col, pos_col, pheno_id, out_dir)

if __name__ == "__main__":
    main()
