#!/usr/bin/env python3

import pandas as pd
import argparse

VALID = {"A","C","G","T"}

def prep_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    needed = ["SNP","A1","A2","BETA","SE"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"]   = pd.to_numeric(df["SE"], errors="coerce")
    df = df[df["SNP"].notna()]
    df = df[df["A1"].isin(VALID) & df["A2"].isin(VALID)]
    df = df[df["SE"].notna() & (df["SE"] > 0) & df["BETA"].notna()]
    if "FRQ" in df.columns and "FREQ" not in df.columns:
        df["FREQ"] = pd.to_numeric(df["FRQ"], errors="coerce")
    elif "FREQ" in df.columns:
        df["FREQ"] = pd.to_numeric(df["FREQ"], errors="coerce")
    elif "MAF" in df.columns:
        df["FREQ"] = pd.to_numeric(df["MAF"], errors="coerce")

    df["Z"] = df["BETA"] / df["SE"]
    df = df.sort_values("SNP")
    df = df.drop_duplicates(subset=["SNP"], keep="first")
    return df

def main():
    parser = argparse.ArgumentParser(description="Prepare LAVA sumstats TSV: clean, ensure FREQ and Z.")
    parser.add_argument("-i","--input", required=True, help="Input GWAS summary stats TSV file")
    parser.add_argument("-o","--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    df = pd.read_csv(args.input, sep="\t")
    df = prep_data(df)
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
