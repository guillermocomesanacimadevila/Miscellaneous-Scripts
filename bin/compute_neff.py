#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def calc_neff(df_path: Path, cases: float, controls: float) -> pd.DataFrame:
    df = pd.read_csv(df_path, sep="\t", dtype=str, low_memory=False)
    if cases <= 0 or controls <= 0:
        raise ValueError("cases and controls must be > 0")
    n_eff = 4.0 / ((1.0 / cases) + (1.0 / controls))
    df["N"] = int(round(n_eff))
    return df

def main():
    ap = argparse.ArgumentParser(description="Add constant N_eff column to LDSC-ready GWAS")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--pheno_id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cases", required=True, type=float)
    ap.add_argument("--controls", required=True, type=float)
    a = ap.parse_args()
    src = Path(a.inp)
    outdir = Path(a.outdir)
    out = outdir / a.pheno_id / "post-qc" / f"{a.pheno_id}.ldsc_ready_neff.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = calc_neff(src, a.cases, a.controls)
    df.to_csv(out, sep="\t", index=False)
    print(f"written: {out}")

if __name__ == "__main__":
    main()