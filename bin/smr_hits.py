#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def get_shared_smr_hits(pheno1_id, pheno2_id, pheno1_df, pheno2_df, xqtl, chrom, start, end, outdir):
    df1 = pd.read_csv(pheno1_df, sep="\t")
    df2 = pd.read_csv(pheno2_df, sep="\t")
    df1["ProbeChr"] = pd.to_numeric(df1["ProbeChr"], errors="coerce")
    df1["Probe_bp"] = pd.to_numeric(df1["Probe_bp"], errors="coerce")
    df1["p_SMR_multi"] = pd.to_numeric(df1["p_SMR_multi"], errors="coerce")
    df1["p_HEIDI"] = pd.to_numeric(df1["p_HEIDI"], errors="coerce")
    df2["ProbeChr"] = pd.to_numeric(df2["ProbeChr"], errors="coerce")
    df2["Probe_bp"] = pd.to_numeric(df2["Probe_bp"], errors="coerce")
    df2["p_SMR_multi"] = pd.to_numeric(df2["p_SMR_multi"], errors="coerce")
    df2["p_HEIDI"] = pd.to_numeric(df2["p_HEIDI"], errors="coerce")
    locus_id = f"chr{chrom}_{start}_{end}"
    base_outdir = os.path.join(outdir, f"{pheno1_id}_{pheno2_id}", xqtl, locus_id)
    os.makedirs(base_outdir, exist_ok=True)
    shared_qtls = sorted(set(df1["qtl_name"].dropna()) & set(df2["qtl_name"].dropna()))

    all_shared = []
    for qtl in shared_qtls:
        locus1 = df1[
            (df1["qtl_name"] == qtl) &
            (df1["ProbeChr"] == chrom) &
            (df1["Probe_bp"].between(start, end))
        ].copy()

        locus2 = df2[
            (df2["qtl_name"] == qtl) &
            (df2["ProbeChr"] == chrom) &
            (df2["Probe_bp"].between(start, end))
        ].copy()

        if locus1.empty or locus2.empty:
            print(qtl, "empty in one trait")
            continue

        p_bonf1 = 0.05 / locus1.shape[0]
        p_bonf2 = 0.05 / locus2.shape[0]

        sig1 = locus1[
            (locus1["p_SMR_multi"] < p_bonf1) &
            (locus1["p_HEIDI"] > 0.01)
        ].copy()

        sig2 = locus2[
            (locus2["p_SMR_multi"] < p_bonf2) &
            (locus2["p_HEIDI"] > 0.01)
        ].copy()

        if sig1.empty or sig2.empty:
            print(qtl, "no sig hits")
            continue

        shared = sig1.merge(
            sig2,
            on=["qtl_name", "gene_id"],
            suffixes=(f"_{pheno1_id}", f"_{pheno2_id}")
        )

        print(qtl, shared.shape)

        if shared.empty:
            continue

        shared["shared_gene"] = shared[f"index_{pheno1_id}"]
        shared = shared[[
            "qtl_name",
            "shared_gene",
            "gene_id",
            f"topSNP_{pheno1_id}", f"A1_{pheno1_id}", f"A2_{pheno1_id}", f"b_SMR_{pheno1_id}", f"p_SMR_multi_{pheno1_id}", f"p_HEIDI_{pheno1_id}",
            f"topSNP_{pheno2_id}", f"A1_{pheno2_id}", f"A2_{pheno2_id}", f"b_SMR_{pheno2_id}", f"p_SMR_multi_{pheno2_id}", f"p_HEIDI_{pheno2_id}"
        ]].copy()

        shared.to_csv(
            os.path.join(base_outdir, f"{qtl}_{pheno1_id}_{pheno2_id}_shared.tsv"),
            sep="\t",
            index=False
        )

        all_shared.append(shared)

    if all_shared:
        all_shared_df = pd.concat(all_shared, axis=0, ignore_index=True)
        all_shared_df.to_csv(
            os.path.join(base_outdir, f"{locus_id}_{pheno1_id}_{pheno2_id}_ALL_shared.tsv"),
            sep="\t",
            index=False
        )
        return all_shared_df

    print("No shared hits found.")
    return pd.DataFrame()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1_id", required=True)
    parser.add_argument("--pheno2_id", required=True)
    parser.add_argument("--pheno1_df", required=True)
    parser.add_argument("--pheno2_df", required=True)
    parser.add_argument("--xqtl", required=True)
    parser.add_argument("--chr", required=True, type=int)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    get_shared_smr_hits(
        pheno1_id=args.pheno1_id,
        pheno2_id=args.pheno2_id,
        pheno1_df=args.pheno1_df,
        pheno2_df=args.pheno2_df,
        xqtl=args.xqtl,
        chrom=args.chr,
        start=args.start,
        end=args.end,
        outdir=args.outdir
    )