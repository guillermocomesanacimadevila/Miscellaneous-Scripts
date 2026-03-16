#!/usr/bin/env python3
import os
import pandas as pd
import argparse
from pathlib import Path
import subprocess

def check_ld_between_xqtl_smr_and_susie_hits(
        pheno1_id: str,
        pheno2_id: str,
        xqtl: str,
        res_dir: str,
        chrom,
        start,
        end,
        susie_dir: str,
        ref_bfile: str,
        outdir: str):

    res_dir = Path(res_dir)
    susie_dir = Path(susie_dir)
    outdir = Path(outdir)
    locus_id = f"chr{chrom}_{start}_{end}"
    smr_locus_dir = res_dir / f"{pheno1_id}_{pheno2_id}" / xqtl / locus_id
    overlap_file = susie_dir / "overlaps" / f"{pheno1_id}_{pheno2_id}" / f"{pheno1_id}_{pheno2_id}_{locus_id}_mapped.tsv"
    trait1_file = susie_dir / pheno1_id / f"locus_{locus_id}" / f"cs95_{pheno1_id}_L1_mapped.tsv"
    trait2_file = susie_dir / pheno2_id / f"locus_{locus_id}" / f"cs95_{pheno2_id}_L1_mapped.tsv"
    outdir.mkdir(parents=True, exist_ok=True)

    smr_top_snps = []

    if not smr_locus_dir.exists():
        raise FileNotFoundError(f"{smr_locus_dir} not found")

    for file in os.listdir(smr_locus_dir):
        file = smr_locus_dir / file

        if not file.is_file():
            continue
        if file.suffix != ".tsv":
            continue
        if file.stat().st_size == 0:
            continue

        df = pd.read_csv(file, sep="\t")

        if df.empty:
            continue

        if f"topSNP_{pheno1_id}" in df.columns:
            smr_top_snps.extend(df[f"topSNP_{pheno1_id}"].dropna().astype(str).tolist())

        if f"topSNP_{pheno2_id}" in df.columns:
            smr_top_snps.extend(df[f"topSNP_{pheno2_id}"].dropna().astype(str).tolist())

    smr_top_snps = sorted(set(smr_top_snps))

    if len(smr_top_snps) == 0:
        print("No SMR top SNPs found")
        return pd.DataFrame()

    susie_hits = []

    if overlap_file.exists() and overlap_file.stat().st_size > 0:
        overlap_df = pd.read_csv(overlap_file, sep="\t")

        if not overlap_df.empty:
            overlap_df = overlap_df[["SNP"]].dropna().drop_duplicates().copy()
            overlap_df["susie_source"] = "overlap"
            overlap_df["trait"] = "overlap"
            susie_hits = overlap_df.to_dict("records")

    if len(susie_hits) == 0:
        if trait1_file.exists() and trait1_file.stat().st_size > 0:
            df1 = pd.read_csv(trait1_file, sep="\t")
            if not df1.empty:
                df1["PIP"] = pd.to_numeric(df1["PIP"], errors="coerce")
                df1 = df1[df1["PIP"] > 0].copy()
                if not df1.empty:
                    df1 = df1[["SNP"]].dropna().drop_duplicates().copy()
                    df1["susie_source"] = "traitwise_pip_gt_0"
                    df1["trait"] = pheno1_id
                    susie_hits.extend(df1.to_dict("records"))

        if trait2_file.exists() and trait2_file.stat().st_size > 0:
            df2 = pd.read_csv(trait2_file, sep="\t")
            if not df2.empty:
                df2["PIP"] = pd.to_numeric(df2["PIP"], errors="coerce")
                df2 = df2[df2["PIP"] > 0].copy()
                if not df2.empty:
                    df2 = df2[["SNP"]].dropna().drop_duplicates().copy()
                    df2["susie_source"] = "traitwise_pip_gt_0"
                    df2["trait"] = pheno2_id
                    susie_hits.extend(df2.to_dict("records"))

    if len(susie_hits) == 0:
        print("No SuSiE SNPs found")
        return pd.DataFrame()

    print("SMR top SNPs:", smr_top_snps)
    print("SuSiE SNPs:", [x["SNP"] for x in susie_hits])

    ld_results = []

    for smr_snp in smr_top_snps:
        for hit in susie_hits:
            susie_snp = hit["SNP"]
            susie_source = hit["susie_source"]
            trait = hit["trait"]

            if smr_snp == susie_snp:
                ld_df = pd.DataFrame([{
                    "CHR_A": chrom,
                    "SNP_A": smr_snp,
                    "CHR_B": chrom,
                    "SNP_B": susie_snp,
                    "R2": 1.0,
                    "Dprime": 1.0,
                    "smr_snp": smr_snp,
                    "susie_snp": susie_snp,
                    "susie_source": susie_source,
                    "trait": trait
                }])
                ld_results.append(ld_df)
                continue

            plink_prefix = outdir / f"{locus_id}_{smr_snp}_{susie_snp}"

            cmd = [
                "plink",
                "--bfile", str(ref_bfile),
                "--chr", str(chrom),
                "--from-bp", str(start),
                "--to-bp", str(end),
                "--ld", str(smr_snp), str(susie_snp),
                "--out", str(plink_prefix)
            ]

            subprocess.run(cmd, check=False)

            log_file = Path(f"{plink_prefix}.log")
            nosex_file = Path(f"{plink_prefix}.nosex")

            if not log_file.exists():
                continue
            if log_file.stat().st_size == 0:
                if log_file.exists():
                    log_file.unlink()
                if nosex_file.exists():
                    nosex_file.unlink()
                continue

            r2 = None
            dprime = None

            with open(log_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if "R-sq =" in line and "D' =" in line:
                        left, right = line.split("D' =")
                        r2 = float(left.split("R-sq =")[1].strip())
                        dprime = float(right.strip())
                        break

            if r2 is None:
                if log_file.exists():
                    log_file.unlink()
                if nosex_file.exists():
                    nosex_file.unlink()
                continue

            ld_df = pd.DataFrame([{
                "CHR_A": chrom,
                "SNP_A": smr_snp,
                "CHR_B": chrom,
                "SNP_B": susie_snp,
                "R2": r2,
                "Dprime": dprime,
                "smr_snp": smr_snp,
                "susie_snp": susie_snp,
                "susie_source": susie_source,
                "trait": trait
            }])

            ld_results.append(ld_df)

            if log_file.exists():
                log_file.unlink()
            if nosex_file.exists():
                nosex_file.unlink()

    if len(ld_results) == 0:
        print("No LD results found")
        return pd.DataFrame()

    ld_results = pd.concat(ld_results, axis=0, ignore_index=True)
    ld_results = ld_results[[
        "CHR_A", "SNP_A", "CHR_B", "SNP_B",
        "R2", "Dprime",
        "smr_snp", "susie_snp",
        "susie_source", "trait"
    ]]
    ld_results.to_csv(
        outdir / f"{pheno1_id}_{pheno2_id}_{xqtl}_{locus_id}_smr_vs_susie_ld.tsv",
        sep="\t",
        index=False
    )

    for f in outdir.glob(f"{locus_id}_*"):
        if f.suffix in [".log", ".nosex"]:
            f.unlink()

    return ld_results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1_id", required=True)
    parser.add_argument("--pheno2_id", required=True)
    parser.add_argument("--xqtl", required=True)
    parser.add_argument("--res_dir", required=True)
    parser.add_argument("--chr", required=True, type=int)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--susie_dir", required=True)
    parser.add_argument("--ref_bfile", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    check_ld_between_xqtl_smr_and_susie_hits(
        pheno1_id=args.pheno1_id,
        pheno2_id=args.pheno2_id,
        xqtl=args.xqtl,
        res_dir=args.res_dir,
        chrom=args.chr,
        start=args.start,
        end=args.end,
        susie_dir=args.susie_dir,
        ref_bfile=args.ref_bfile,
        outdir=args.outdir
    )

if __name__ == "__main__":
    main()