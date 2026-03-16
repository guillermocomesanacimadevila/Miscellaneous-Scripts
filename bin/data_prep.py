#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
import os

# steps:
# * parse argument = lead SNP ID -> MAYBE ON THE POST CLUMPING FILE WE SAVED IT
# * MAP +/- 250 KB AROUND IT FOR EACH TRAIT (locus_0 = AD, SCZ, locus_1 = AD/SCZ, locus_2 = AD/SCZ/LON)
# * Compute Z-score
# * Input = post-ldsc sumstats
# * we need to ensure (nf): ../../outputs/SuSiE/locus_0/AD/ld

def assemble_for_susie(
        pheno_sumstats: pd.DataFrame,
        lead_snp: pd.DataFrame,
        out_dir: str,
        pheno_id: str,
        locus_id: str,
        window_bp: int):
    pheno_sumstats = pheno_sumstats.copy()
    lead_snp = lead_snp.copy()
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if "BP" not in pheno_sumstats.columns:
        pheno_sumstats.rename(
            columns={
                "POS": "BP"
            }, inplace=True)
    pheno_sumstats["Z"] = pheno_sumstats["BETA"] / pheno_sumstats["SE"]
    lead = list(lead_snp["SNP"])
    lead = lead[0]
    # map lead snp in pheno_sumstats -> gather +/- 250Kb -> with respective SNP/CHR/BP/BETA/P/SE -> save as ${locus_id}_${pheno_id}_susie_ready.tsv to out_dir
    hit = pheno_sumstats.loc[
        pheno_sumstats["SNP"] == lead,
        ["CHR", "BP"]
    ]
    lead_chr = hit["CHR"].iloc[0]
    lead_bp = int(hit["BP"].iloc[0])
    start = lead_bp - window_bp
    end = lead_bp + window_bp
    region = pheno_sumstats.loc[
        (pheno_sumstats["CHR"] == lead_chr) &
        (pheno_sumstats["BP"] >= start) &
        (pheno_sumstats["BP"] <= end)
        ]
    out_path = out_dir / f"{locus_id}/{pheno_id}/{locus_id}_{pheno_id}_susie_ready.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    region[["SNP", "CHR", "A1", "A2", "BP", "N", "BETA", "SE", "P", "Z"]].to_csv(
        out_path, sep="\t", index=False
    )
    return out_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--lead_snp", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--pheno_id", required=True)
    parser.add_argument("--locus_id", required=True)
    parser.add_argument("--window_bp", type=int, default=500_000)
    args = parser.parse_args()
    pheno_sumstats = pd.read_csv(args.sumstats, sep="\t")
    lead_snp = pd.read_csv(args.lead_snp, sep="\t")
    out_path = assemble_for_susie(
        pheno_sumstats=pheno_sumstats,
        lead_snp=lead_snp,
        out_dir=args.out_dir,
        pheno_id=args.pheno_id,
        locus_id=args.locus_id,
        window_bp=args.window_bp
    )

if __name__ == "__main__":
    main()