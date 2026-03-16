#!/usr/bin/env python3
import argparse
import pandas as pd
import polars as pl
import os
from pathlib import Path

# Add N col (--N) arg
def preprocess_xqtl(qtl_dir: str,
                    qtl_type: str,
                    n: int,
                    qtl_dataset: str):

    qtl_dir = Path(qtl_dir)
    ls = []
    for file in os.listdir(qtl_dir):
        if file.startswith(f"{qtl_type}_") and qtl_dataset in file:
            ls.append(file)

    # now look at each chr
    for file in ls:
        print(f"Starting preprocessing {file}")
        df = pd.read_csv(qtl_dir / file, sep="\t")
        df["N"] = n
        if qtl_dataset == "BrainMeta":
            df.rename(columns={
                "Chr": "CHR",
                "Freq": "FRQ",
                "b": "BETA",
                "p": "P"
            }, inplace=True)
        df = df[[
            "Gene",
            "SNP",
            "CHR",
            "BP",
            "A1",
            "A2",
            "N",
            "FRQ",
            "BETA",
            "SE",
            "P"
        ]]
        df.to_csv(qtl_dir / file, sep="\t", index=False)
        print(f"Finished preprocessing {file}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("qtl_dir", help="directory containing xQTL files")
    parser.add_argument("qtl_type", help="xQTL type (eQTL / sQTL)")
    parser.add_argument("n", type=int, help="sample size")
    parser.add_argument("qtl_dataset", help="dataset name (e.g. BrainMeta)")
    args = parser.parse_args()
    preprocess_xqtl(
        qtl_dir=args.qtl_dir,
        qtl_type=args.qtl_type,
        n=args.n,
        qtl_dataset=args.qtl_dataset
    )

if __name__ == "__main__":
    main()