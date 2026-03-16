#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def read_sumstats(path):
    df = pd.read_csv(path, sep="\t", dtype={"CHR": "Int64", "POS": "Int64"})
    df = df[df["CHR"].between(1, 22)]
    df["build"] = "GRCh37"
    return df

def remove_mhc(df):
    m = ~((df["CHR"] == 6) & (df["POS"] >= 25119106) & (df["POS"] <= 33854733))
    return df[m].copy()

def prepare_trait(df, prefix):
    cols = ["SNP", "A1", "A2", "P", "CHR", "POS"]
    extra = [c for c in ["BETA", "SE", "N"] if c in df.columns]
    df = df[cols + extra].copy()
    df = remove_mhc(df)
    df = df.drop_duplicates(subset=["SNP"])
    df = df[df["P"].between(0, 1)]
    df = df.rename(columns={c: f"{c}_{prefix}" for c in df.columns if c not in ["SNP", "CHR", "POS"]})
    return df

def harmonise(t1, t2, p1, p2):
    m = t1.merge(t2, on=["SNP", "CHR", "POS"], how="inner")
    same = (m[f"A1_{p1}"] == m[f"A1_{p2}"]) & (m[f"A2_{p1}"] == m[f"A2_{p2}"])
    flip = (m[f"A1_{p1}"] == m[f"A2_{p2}"]) & (m[f"A2_{p1}"] == m[f"A1_{p2}"])
    keep = same | flip
    m = m[keep].copy()
    b2 = f"BETA_{p2}"
    if b2 in m.columns:
        m.loc[flip, b2] = -m.loc[flip, b2]
    m.loc[flip, f"A1_{p2}"] = m.loc[flip, f"A1_{p1}"]
    m.loc[flip, f"A2_{p2}"] = m.loc[flip, f"A2_{p1}"]
    m["build"] = "GRCh37"
    return m

def build_cfdr_input(m, p1, p2):
    out = m[[f"SNP", f"P_{p2}", f"P_{p1}"]].copy()
    out = out.rename(columns={f"P_{p2}": "p1", f"P_{p1}": "p2"})
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trait1", required=True)
    p.add_argument("--trait2", required=True)
    p.add_argument("--prefix1", required=True)
    p.add_argument("--prefix2", required=True)
    p.add_argument("--out_prefix", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()
    t1 = read_sumstats(Path(args.trait1))
    t2 = read_sumstats(Path(args.trait2))
    t1p = prepare_trait(t1, args.prefix1)
    t2p = prepare_trait(t2, args.prefix2)
    m = harmonise(t1p, t2p, args.prefix1, args.prefix2)
    cfdr = build_cfdr_input(m, args.prefix1, args.prefix2)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    m.to_csv(out_dir / f"{args.out_prefix}.harmonised_{args.prefix1}_{args.prefix2}.tsv", sep="\t", index=False)
    cfdr.to_csv(out_dir / f"{args.out_prefix}.cfdr_input_{args.prefix1}_{args.prefix2}.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
