#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
import os

# mkdir -p ${workflow.launchDir}/outputs/conjFDR/${trait_pair}/mapped
# from ${pheno1_sumstats} -> map SNP ["CHR"] & ["BP"] onto conjFDR res

# STEPS:
# * from AD df -> extract cols = ["SNP", "CHR", "BP"] - from multi-option array
# * parse conjFDR output for trait pair (X, Y)
# * for conjFDR output -> merge on "SNP"

def assemble_pleio_results(
        pheno1_df: pd.DataFrame,
        conjFDR: pd.DataFrame,
        pheno_id: str,
        out_dir: str,
        chr_col: str,
        bp_col: str,
        snp_col: str):

    pheno1_df = pheno1_df.copy()
    conjFDR = conjFDR.copy()
    needed_pheno = [snp_col, bp_col, chr_col]
    missing_pheno = [c for c in needed_pheno if c not in pheno1_df.columns]
    if missing_pheno:
        raise ValueError(f"[{pheno_id}] Missing columns in pheno1_df: {missing_pheno}")

    if snp_col not in conjFDR.columns:
        raise ValueError(f"Missing merge column '{snp_col}' in conjFDR. Available: {list(conjFDR.columns)}")

    pheno1_df = pheno1_df[needed_pheno].dropna(subset=[snp_col]).drop_duplicates(subset=[snp_col])
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    res = conjFDR.merge(pheno1_df, on=snp_col, how="left", validate="m:1")
    res.to_csv(out_dir / "mapped_results.tsv", sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1", required=True, help="Pheno1 sumstats file (tsv/csv). Must contain SNP/CHR/BP (or POS).")
    parser.add_argument("--conjfdr", required=True, help="conjFDR results file (tsv/csv). Must contain SNP column.")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--pheno-id", required=True, help="Phenotype ID label (used for error messages)")
    parser.add_argument("--snp-col", default="SNP", help="SNP column name (default: SNP)")
    parser.add_argument("--chr-col", default="CHR", help="Chromosome column name (default: CHR)")
    parser.add_argument("--bp-col", default="BP", help="Base-pair/position column name (default: BP)")
    args = parser.parse_args()
    pheno1_path = Path(args.pheno1)
    conjfdr_path = Path(args.conjfdr)
    sep_pheno = "\t" if pheno1_path.suffix.lower() in [".tsv", ".txt", ".gz"] else ","
    sep_conj = "\t" if conjfdr_path.suffix.lower() in [".tsv", ".txt", ".gz"] else ","
    pheno1_df = pd.read_csv(pheno1_path, sep=sep_pheno, dtype=str, low_memory=False)
    conj_df = pd.read_csv(conjfdr_path, sep=sep_conj, dtype=str, low_memory=False)
    assemble_pleio_results(
        pheno1_df=pheno1_df,
        conjFDR=conj_df,
        pheno_id=args.pheno_id,
        out_dir=args.outdir,
        chr_col=args.chr_col,
        bp_col=args.bp_col,
        snp_col=args.snp_col,
    )

if __name__ == "__main__":
    main()
