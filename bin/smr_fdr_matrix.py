#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import os

# Automated structure
# * Take SMR results on X eQTL datasets
# * For trait Y -> split into different datasets (polars) - 1 per eQTL dataset
# * For a list of genes (protein coding) pertaining to locus Y_i (CLI arg)
# * Find genes in each eQTL dataset
# * Correct for multiple testing on each dataset (FDR)
# * Save onto .tsv a dataframe (rows = Gene name, columns = eQTL dataset where found, and in each cell -> q-val / FDR corrected p

def bh_qvals(p):
    p = np.asarray(p, dtype=float)
    m = p.size
    o = np.argsort(p)
    q = (p[o] * m) / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[o] = q
    return out

def get_genes(s):
    p = Path(s)
    if p.exists():
        return [x.strip() for x in p.read_text().splitlines() if x.strip()]
    return [x.strip() for x in s.split(",") if x.strip()]

def run_one(df, qtl, genes, pcol, alpha, pheidi_col, pheidi_min):
    d = df[df["qtl_name"] == qtl].copy()
    d[pcol] = pd.to_numeric(d[pcol], errors="coerce")
    d = d[np.isfinite(d[pcol])].copy()
    if pheidi_min is not None:
        if pheidi_col not in d.columns:
            raise ValueError(f"Missing `{pheidi_col}`")
        d[pheidi_col] = pd.to_numeric(d[pheidi_col], errors="coerce")
        d = d[np.isfinite(d[pheidi_col])].copy()
        d = d[d[pheidi_col] >= float(pheidi_min)].copy()
    m = len(d)
    if m == 0:
        return {g: None for g in genes}, m, np.nan
    d["qval"] = bh_qvals(d[pcol].to_numpy())
    d = d.sort_values(pcol, ascending=True).reset_index(drop=True)
    bh_line = (np.arange(1, m + 1) / m) * alpha
    ok = d[pcol].to_numpy() <= bh_line
    p_cut = d.loc[ok].iloc[-1][pcol] if ok.any() else np.nan
    out = {}
    idx = d["index"].astype(str)
    for g in genes:
        gg = d[idx == g]
        out[g] = None if gg.empty else float(gg.sort_values(pcol).iloc[0]["qval"])
    return out, m, float(p_cut) if np.isfinite(p_cut) else np.nan

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--pheno_id", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--pcol", default="p_SMR_multi")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--qtls", default="")
    ap.add_argument("--pheidi_col", default="p_HEIDI")
    ap.add_argument("--pheidi_min", type=float, default=None)
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = out_dir / f"{args.pheno_id}.tsv"
    df = pd.read_csv(args.infile, sep="\t")
    print(f"shape: {df.shape}")
    print(df.columns)
    if "qtl_name" not in df.columns or "index" not in df.columns:
        raise ValueError("need columns: qtl_name, index")
    if args.pcol not in df.columns:
        raise ValueError(f"Missing `{args.pcol}`")
    genes = get_genes(args.genes)
    qtls = [x.strip() for x in args.qtls.split(",") if x.strip()] if args.qtls.strip() else sorted(df["qtl_name"].dropna().astype(str).unique())
    mat = pd.DataFrame(index=genes)
    for qtl in qtls:
        g2q, m, p_cut = run_one(df, qtl, genes, args.pcol, args.alpha, args.pheidi_col, args.pheidi_min)
        mat[qtl] = [g2q[g] for g in genes]
        print(f"{qtl}\tm={m}\tBH_p_cutoff={p_cut}")
    mat.to_csv(out_path, sep="\t", na_rep="NA")
    sig = (mat.apply(pd.to_numeric, errors="coerce") <= args.alpha)
    print("sig_counts_per_gene")
    for g, c in sig.sum(axis=1).sort_values(ascending=False).items():
        print(f"{g}\t{int(c)}")
    print(f"wrote: {out_path}")

if __name__ == "__main__":
    main()