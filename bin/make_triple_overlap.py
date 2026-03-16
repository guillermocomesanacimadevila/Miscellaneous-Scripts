#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def load_hits(conjfdr_dir, pheno):
    p = Path(conjfdr_dir) / pheno / f"{pheno}_shared_hits.tsv"
    if not p.exists():
        raise FileNotFoundError(str(p))
    return pd.read_csv(p, sep="\t")

def load_clumped_leads(conjfdr_dir, pheno):
    p = Path(conjfdr_dir) / pheno / f"{pheno}_clump.clumped"
    if not p.exists():
        return pd.DataFrame(columns=["SNP"])
    df = pd.read_csv(p, sep=r"\s+")
    if "SNP" not in df.columns:
        return pd.DataFrame(columns=["SNP"])
    return df[["SNP"]].dropna().drop_duplicates()

def filter_by_fdr(df, thresh):
    if "conj_fdr" in df.columns:
        return df[df["conj_fdr"] < thresh].copy()
    cols = [c for c in df.columns if c.startswith("conj_fdr")]
    if len(cols) > 0:
        return df[df[cols[0]] < thresh].copy()
    return df.copy()

def suffix(df, tag):
    m = {}
    for c in df.columns:
        if c != "SNP":
            m[c] = f"{c}_{tag}"
    return df.rename(columns=m)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conjfdr-dir", required=True)
    ap.add_argument("--phenos", nargs=2, required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--thresh", type=float, default=0.05)
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p1, p2 = args.phenos[0], args.phenos[1]
    d1 = filter_by_fdr(load_hits(args.conjfdr_dir, p1), args.thresh)
    d2 = filter_by_fdr(load_hits(args.conjfdr_dir, p2), args.thresh)
    d1 = d1.drop_duplicates(subset=["SNP"]).copy()
    d2 = d2.drop_duplicates(subset=["SNP"]).copy()
    d1s = suffix(d1, p1)
    d2s = suffix(d2, p2)
    merged = d1s.merge(d2s, on="SNP", how="inner")
    snps = merged[["SNP"]].drop_duplicates()
    snps.to_csv(out_dir / "overlap_snps.tsv", sep="\t", index=False)
    merged.to_csv(out_dir / "overlap_merged.tsv", sep="\t", index=False)
    leads1 = load_clumped_leads(args.conjfdr_dir, p1)
    leads2 = load_clumped_leads(args.conjfdr_dir, p2)
    all_leads = pd.concat([leads1, leads2], ignore_index=True).drop_duplicates()
    overlap_leads = all_leads.merge(snps, on="SNP", how="inner")
    overlap_leads.to_csv(out_dir / "overlap_leads.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
