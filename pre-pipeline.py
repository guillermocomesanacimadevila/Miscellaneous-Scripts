#!/usr/bin/env python3
import argparse
from bed_reader import open_bed
import polars as pl
import numpy as np
from pathlib import Path
from tqdm import tqdm

# plink 1.9/2.0 format for bed / bim / fam files
# https://plink.readthedocs.io/en/latest/plink_fmt/
# .fam: FID IID PID MID SEX PHENOTYPE
# .bim: CHR SNP_ID BP A1 A2
# .bed: Categorical GT (0/1/2)
# -> MAF file: SNP MAF

def read_fam(fam_file):
    return pl.read_csv(
        fam_file,
        separator=" ",
        has_header=False,
        new_columns=["FID", "IID", "PAD", "MAD", "SEX", "PHENO"]
    )

def read_bim(bim_file):
    return pl.read_csv(
        bim_file,
        separator="\t",
        has_header=False,
        new_columns=["CHR", "SNP", "CM", "BP", "A1", "A2"]
    ).with_row_index("idx")

def read_sumstats(sumstats_file):
    df = pl.read_csv(
        sumstats_file,
        separator="\t",
        has_header=True,
    ).select(
        pl.col("SNP"),
        pl.col("Beta").cast(pl.Float64),
        pl.col("SE").cast(pl.Float64)
    )
    return df

def read_freq(freq_file):
    return pl.read_csv(
        freq_file,
        separator="\t",
        has_header=True,
    ).select(
        pl.col("SNP"),
        pl.col("MAF").cast(pl.Float64),
    )

def build_snp_params(bim_file, sumstats_file, frq_file):
    df = (
        bim_file.select(["SNP", "idx"])
        .join(sumstats_file.select(["SNP", "Beta", "SE"]), on="SNP", how="left")
        .join(frq_file.select(["SNP", "MAF"]), on="SNP", how="left")
    )
    df = df.with_columns([
        pl.col("Beta").fill_null(0.0).cast(pl.Float64).alias("Beta"),
        pl.col("SE").fill_null(0.0).cast(pl.Float64).alias("SE"),
        pl.col("MAF").fill_null(0.0).cast(pl.Float64).alias("MAF"),
    ])
    # precompute s^2 and MAF * 2 -> for PRS and PRS SE calcs
    df = df.with_columns([
        pl.col("SE").pow(2).alias("se2"),
        (pl.col("MAF") * 2).alias("mu"),
    ])
    beta = df["Beta"].to_numpy()
    se2 = df["se2"].to_numpy()
    mu = df["mu"].to_numpy()
    # checks
    maf0 = df["MAF"][0]
    se0 = df["SE"][0]
    print("MAF (e.g. entry - index [0]):", maf0)
    print("2*MAF (e.g. entry - index [0]):", 2 * maf0)
    print("\nSE (e.g. entry - index [0]):", se0)
    print("SE**2 (e.g. entry - index [0]):", se0 ** 2)
    return beta, se2, mu

def iter_genotype_blocks(bed_prefix, block_snps): # block_snps + arg
    bed = open_bed(bed_prefix + ".bed")
    M = bed.sid_count # 27,452
    n_blocks = (M + block_snps - 1) // block_snps
    blocks = []
    j0 = 0
    with tqdm(total=n_blocks, desc="Reading GTs", unit="block") as pbar:
        while j0 < M: # === takes forever if returned directly from func === #
            j1 = min(j0 + block_snps, M)
            G = bed.read(index=np.s_[:, j0:j1])
            G = G.astype(np.float64)
            G[np.isnan(G)] = 0.0
            blocks.append((j0, j1, G))
            j0 = j1
            pbar.update(1)
            pbar.set_postfix({"SNPs": f"{j1}/{M}"})
    return blocks

def compute_prs(patient_ids, blocks, beta, se2, mu, z=1.96):
    N = len(patient_ids)
    prs = np.zeros(N, dtype=np.float64)
    var = np.zeros(N, dtype=np.float64)
    for j0, j1, G in tqdm(blocks, desc="Computing PRs", unit="block"):
        beta_block = beta[j0:j1]
        se2_block = se2[j0:j1]
        mu_block = mu[j0:j1]
        # (gi,k - 2fk)
        X = G - mu_block
        # PRS,i = sigma(beta,k (gi,k - 2fk))
        prs += np.sum(beta_block * X, axis=1)
        # SE(PRS,i) = sqrt(sigma(((gi,k - 2fk)**2) * SE**2))
        var += np.sum((X * X) * se2_block, axis=1)
    prs_se = np.sqrt(var)
    l95 = prs - z * prs_se
    u95 = prs + z * prs_se
    result = pl.DataFrame({
        "IID": patient_ids,
        "myPRS": prs,
        "myPRS.SE": prs_se,
        "myPRS.L95CI": l95,
        "myPRS.U95CI": u95,
    })
    return result

def main():
    parser = argparse.ArgumentParser(description="Compute PRS + PRS SE + 95% CI")
    parser.add_argument("--bed-prefix", required=True)
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--freq", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--block-snps", type=int, default=5000)
    parser.add_argument("--z", type=float, default=1.96)
    args = parser.parse_args()
    prefix = args.bed_prefix
    bed_path = Path(prefix + ".bed")
    bim_path = Path(prefix + ".bim")
    fam_path = Path(prefix + ".fam")
    sumstats_path = Path(args.sumstats)
    freq_path = Path(args.freq)
    fam = read_fam(str(fam_path))
    bim = read_bim(str(bim_path))
    sumstats = read_sumstats(str(sumstats_path))
    freq = read_freq(str(freq_path))
    beta, se2, mu = build_snp_params(bim, sumstats, freq)
    patient_ids = fam["IID"].to_numpy()
    blocks = iter_genotype_blocks(prefix, args.block_snps)
    result = compute_prs(patient_ids, blocks, beta, se2, mu, z=args.z)
    result.write_csv(
        args.out,
        separator="\t",
    )
    print(f"Done! Results written to {args.out}")

if __name__ == "__main__":
    main()

# ad = read_sumstats("/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/AD_SUM_STATISTICS.txt")
# bim = read_bim("/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/BDR_clumped_AD.bim")
# fam = read_fam("/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/BDR_clumped_AD.fam")
# freq = read_freq("/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/FREQ_file.txt")
# bed = "/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/BDR_clumped_AD.bed"
# dat = build_snp_params(bim, ad, freq)
# print(len(dat)) # 3 vectors (1 for betas, one for s**2 and 1 for MAF x 2)
# print(dat)
# print(iter_genotype_blocks(bed, 5000))
# Check to ensure == no leakage
# import pandas as pd
# res = pd.read_csv("/Users/c24102394/Desktop/PRS-var/PRS_VARIABILITY/prs_results.tsv", sep="\t")
# print(len(df["IID"].unique())) # 514 (initial == 514 patients)
