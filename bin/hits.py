#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def read_mat(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df

def sig_mask(df, alpha):
    x = df.apply(pd.to_numeric, errors="coerce")
    return x.le(alpha)

def overlap_table(a, b, alpha, trait_a, trait_b):
    sa = sig_mask(a, alpha)
    sb = sig_mask(b, alpha)
    common_cols = [c for c in a.columns if c in b.columns]
    if not common_cols:
        raise ValueError("no overlapping tissues/columns between matrices")
    sa = sa[common_cols]
    sb = sb[common_cols]
    common_genes = sa.index.intersection(sb.index)
    sa = sa.loc[common_genes]
    sb = sb.loc[common_genes]
    both = sa & sb
    keep = both.any(axis=1)
    both = both.loc[keep]
    rows = []
    for g in both.index:
        tissues = both.columns[both.loc[g]].tolist()
        rows.append({
            "gene": g,
            f"n_tissues_{trait_a}": int(sa.loc[g].sum()),
            f"n_tissues_{trait_b}": int(sb.loc[g].sum()),
            "n_tissues_both": int(both.loc[g].sum()),
            "tissues_both": ",".join(tissues)
        })
    if not rows:
        return pd.DataFrame(columns=[
            "gene",
            f"n_tissues_{trait_a}",
            f"n_tissues_{trait_b}",
            "n_tissues_both",
            "tissues_both"
        ])
    out = pd.DataFrame(rows).sort_values(["n_tissues_both", "gene"], ascending=[False, True])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat_a", required=True)
    ap.add_argument("--mat_b", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--trait_a", default="traitA")
    ap.add_argument("--trait_b", default="traitB")
    args = ap.parse_args()
    a = read_mat(Path(args.mat_a))
    b = read_mat(Path(args.mat_b))
    out = overlap_table(a, b, args.alpha, args.trait_a, args.trait_b)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"wrote: {out_path}")
    print(f"overlap genes: {out.shape[0]}")

if __name__ == "__main__":
    main()
