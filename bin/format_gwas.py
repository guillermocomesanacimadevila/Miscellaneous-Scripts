#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def parse_and_reformat(df: pd.DataFrame,
                       snp_col: str,
                       a1_col: str,
                       a2_col: str,
                       freq_col: str,
                       beta_col: str,
                       se_col: str,
                       p_col: str,
                       n_col: str,
                       out_dir: str,
                       pheno_id: str) -> pd.DataFrame:
    df = df.copy()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.rename(columns={
        snp_col: 'SNP',
        a1_col: 'A1',
        a2_col: 'A2',
        freq_col: 'freq',
        beta_col: 'b',
        se_col: 'se',
        p_col: 'p',
        n_col: 'n',
    }, inplace=True)
    df = df[['SNP','A1','A2','freq','b','se','p','n']].copy()
    df['SNP'] = df['SNP'].astype(str)
    df['A1'] = df['A1'].astype(str).str.upper()
    df['A2'] = df['A2'].astype(str).str.upper()
    for c in ['freq','b','se','p','n']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=['SNP','A1','A2','freq','b','se','p'])
    out_path = out_dir / f"{pheno_id}.smr.ma"
    df.to_csv(out_path, sep='\t', index=False)
    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pheno_id', type=str)
    parser.add_argument('--infile', required=True, type=str)
    parser.add_argument('--out_dir', required=True, type=str)
    parser.add_argument('--snp_col', required=True, type=str)
    parser.add_argument('--a1_col', required=True, type=str)
    parser.add_argument('--a2_col', required=True, type=str)
    parser.add_argument('--freq_col', required=True, type=str)
    parser.add_argument('--beta_col', required=True, type=str)
    parser.add_argument('--se_col', required=True, type=str)
    parser.add_argument('--p_col', required=True, type=str)
    parser.add_argument('--n_col', required=True, type=str)
    parser.add_argument('--sep', default='\t', type=str)
    args = parser.parse_args()
    df = pd.read_csv(Path(args.infile), sep=args.sep)
    parse_and_reformat(
        df=df,
        snp_col=args.snp_col,
        a1_col=args.a1_col,
        a2_col=args.a2_col,
        freq_col=args.freq_col,
        beta_col=args.beta_col,
        se_col=args.se_col,
        p_col=args.p_col,
        n_col=args.n_col,
        out_dir=args.out_dir,
        pheno_id=args.pheno_id
    )

if __name__ == '__main__':
    main()
