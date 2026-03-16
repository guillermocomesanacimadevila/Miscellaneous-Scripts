#!/usr/bin/env python3
import os
import pandas as pd
import argparse
from pathlib import Path
from statsmodels.stats.multitest import multipletests

def analyse_qtls(
    smr_res_trait1: str,
    smr_res_trait2: str,
    pheno1_prefix: str,
    pheno2_prefix: str,
    gene_list_path: str,
    out_dir: str,
    heidi_thresh: float,
    qtl_type: str
):
    pd.set_option("display.float_format", "{:.10f}".format)
    t1 = pd.read_csv(smr_res_trait1, sep="\t")
    t2 = pd.read_csv(smr_res_trait2, sep="\t")
    required_cols = {
        "qtl_name","index","p_SMR_multi","p_HEIDI","b_SMR","se_SMR","p_SMR","A1","A2",
        "topSNP","topSNP_chr","topSNP_bp","probeID","gene_id","ProbeChr","Probe_bp",
        "Freq","nsnp_HEIDI"
    }
    miss1 = required_cols - set(t1.columns)
    miss2 = required_cols - set(t2.columns)
    if miss1:
        raise ValueError(f"{smr_res_trait1} missing columns: {sorted(miss1)}")
    if miss2:
        raise ValueError(f"{smr_res_trait2} missing columns: {sorted(miss2)}")
    genes = (
        pd.read_csv(gene_list_path, header=None, sep=None, engine="python")
          .iloc[:, 0]
          .astype(str)
          .str.strip()
          .str.upper()
          .tolist()
    )
    genes = [g for g in genes if g]
    genes_set = set(genes)
    qtls = sorted(set(t1["qtl_name"].unique()) & set(t2["qtl_name"].unique()))
    base_out = Path(out_dir) / "SMR" / "res" / str(qtl_type) / f"{pheno1_prefix}_{pheno2_prefix}"
    base_out.mkdir(parents=True, exist_ok=True)
    keep_cols = [
        "GENE","index","gene_id","qtl_name","probeID","ProbeChr","Probe_bp","topSNP",
        "topSNP_chr","topSNP_bp","A1","A2","Freq","b_SMR","se_SMR","p_SMR","p_SMR_multi",
        "q_SMR_multi","p_HEIDI","nsnp_HEIDI","trait"
    ]

    def write_raw_with_q(df, qtl, trait_prefix, out_path):
        sub = df[df["qtl_name"] == qtl].copy()
        if sub.shape[0] == 0:
            pd.DataFrame(columns=list(df.columns) + ["GENE", "q_SMR_multi", "trait"]).to_csv(out_path, sep="\t", index=False)
            return
        sub["GENE"] = sub["index"].astype(str).str.strip().str.upper()
        p = pd.to_numeric(sub["p_SMR_multi"], errors="coerce")
        q = pd.Series([float("nan")] * sub.shape[0], index=sub.index, dtype="float64")
        mask = p.notna()
        if mask.any():
            q.loc[mask] = multipletests(p.loc[mask].to_numpy(), method="fdr_bh")[1]
        sub["q_SMR_multi"] = q
        sub["trait"] = trait_prefix
        sub.to_csv(out_path, sep="\t", index=False)

    def prep_for_overlap(df, qtl):
        df = df[df["qtl_name"] == qtl].copy()
        if df.shape[0] == 0:
            return df
        df["GENE"] = df["index"].astype(str).str.strip().str.upper()
        df = df[df["GENE"].isin(genes_set)].copy()
        df["p_SMR_multi_num"] = pd.to_numeric(df["p_SMR_multi"], errors="coerce")
        df = df[df["p_SMR_multi_num"].notna()].copy()
        if df.shape[0] == 0:
            return df
        df["q_SMR_multi"] = multipletests(df["p_SMR_multi_num"].to_numpy(), method="fdr_bh")[1]
        df = df.sort_values("p_SMR_multi_num").drop_duplicates("GENE", keep="first")
        return df

    for qtl in qtls:
        qtl_dir = base_out / qtl
        qtl_dir.mkdir(parents=True, exist_ok=True)
        write_raw_with_q(t1, qtl, pheno1_prefix, qtl_dir / f"{pheno1_prefix}_raw_with_q.tsv")
        write_raw_with_q(t2, qtl, pheno2_prefix, qtl_dir / f"{pheno2_prefix}_raw_with_q.tsv")
        a = prep_for_overlap(t1, qtl)
        s = prep_for_overlap(t2, qtl)
        if a.shape[0] == 0 or s.shape[0] == 0:
            pd.DataFrame(columns=keep_cols).to_csv(qtl_dir / "genes_passed_fdr.tsv", sep="\t", index=False)
            pd.DataFrame(columns=keep_cols).to_csv(qtl_dir / "genes_final.tsv", sep="\t", index=False)
            Path(qtl_dir / f"{qtl_type}_genes.txt").write_text("")
            print(f"[{qtl_type} | {qtl}] passed_fdr=0 final=0")
            continue
        a["trait"] = pheno1_prefix
        s["trait"] = pheno2_prefix
        m = a.merge(s, on="GENE", how="inner", suffixes=("_A", "_S"))
        keep_fdr = (
            (pd.to_numeric(m["q_SMR_multi_A"], errors="coerce") < 0.05) &
            (pd.to_numeric(m["q_SMR_multi_S"], errors="coerce") < 0.05) &
            ((pd.to_numeric(m["p_HEIDI_A"], errors="coerce") > 0) | (pd.to_numeric(m["p_HEIDI_S"], errors="coerce") > 0))
        )
        passed = m[keep_fdr].copy()
        keep_final = (
            keep_fdr &
            (pd.to_numeric(m["p_HEIDI_A"], errors="coerce") > float(heidi_thresh)) &
            (pd.to_numeric(m["p_HEIDI_S"], errors="coerce") > float(heidi_thresh))
        )
        finalm = m[keep_final].copy()

        def to_long(mm):
            if mm.shape[0] == 0:
                return pd.DataFrame(columns=keep_cols)
            a_cols = {c+"_A": c for c in ["index","gene_id","qtl_name","probeID","ProbeChr","Probe_bp","topSNP","topSNP_chr","topSNP_bp","A1","A2","Freq","b_SMR","se_SMR","p_SMR","p_SMR_multi","q_SMR_multi","p_HEIDI","nsnp_HEIDI"]}
            s_cols = {c+"_S": c for c in ["index","gene_id","qtl_name","probeID","ProbeChr","Probe_bp","topSNP","topSNP_chr","topSNP_bp","A1","A2","Freq","b_SMR","se_SMR","p_SMR","p_SMR_multi","q_SMR_multi","p_HEIDI","nsnp_HEIDI"]}
            out_a = mm[["GENE"] + list(a_cols.keys())].rename(columns=a_cols).copy()
            out_s = mm[["GENE"] + list(s_cols.keys())].rename(columns=s_cols).copy()
            out_a["trait"] = pheno1_prefix
            out_s["trait"] = pheno2_prefix
            out = pd.concat([out_a, out_s], ignore_index=True)
            return out[keep_cols].sort_values(["GENE","trait"])

        passed_fdr = to_long(passed)
        final = to_long(finalm)
        passed_fdr.to_csv(qtl_dir / "genes_passed_fdr.tsv", sep="\t", index=False)
        final.to_csv(qtl_dir / "genes_final.tsv", sep="\t", index=False)
        final["GENE"].drop_duplicates().to_csv(qtl_dir / f"{qtl_type}_genes.txt", index=False, header=False)
        print(f"[{qtl_type} | {qtl}] passed_fdr={passed_fdr['GENE'].nunique()} final={final['GENE'].nunique()}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smr_res_trait1", required=True)
    p.add_argument("--smr_res_trait2", required=True)
    p.add_argument("--pheno1_prefix", required=True)
    p.add_argument("--pheno2_prefix", required=True)
    p.add_argument("--gene_list", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--heidi_thresh", type=float, default=0.01)
    p.add_argument("--qtl_type", required=True)
    a = p.parse_args()
    analyse_qtls(
        a.smr_res_trait1,
        a.smr_res_trait2,
        a.pheno1_prefix,
        a.pheno2_prefix,
        a.gene_list,
        a.out_dir,
        a.heidi_thresh,
        a.qtl_type
    )

if __name__ == "__main__":
    main()