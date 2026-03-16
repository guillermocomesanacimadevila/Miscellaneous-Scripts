#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
import os

# loop through SMR/res/....
# for file there
# compile all genes.final
# compile all genes_fdr pass
# save their qtl name of interest
# save onto outputs/QTL_hits/s or e_qtl
# 3 files within each -> final files -> all genes_qtl_fdr.tsv / genes_qtl_causal.tsv / gene_list_all.txt

def compile_qtl_hits(res_root:str,out_root:str,pheno_pair:str):
    res_root=Path(res_root); out_root=Path(out_root)
    for qtl_type in ["eQTL","sQTL"]:
        base=res_root/qtl_type/pheno_pair
        if not base.exists():
            continue
        out_dir=out_root/qtl_type/pheno_pair
        out_dir.mkdir(parents=True,exist_ok=True)
        fdr_frames=[]; final_frames=[]
        for ds_dir in sorted([p for p in base.iterdir() if p.is_dir()]):
            dataset=ds_dir.name
            fdr_path=ds_dir/"genes_passed_fdr.tsv"
            final_path=ds_dir/"genes_final.tsv"
            if fdr_path.exists():
                df=pd.read_csv(fdr_path,sep="\t")
                if df.shape[0]:
                    df["qtl_dataset"]=dataset
                    cols=["qtl_dataset"]+[c for c in df.columns if c!="qtl_dataset"]
                    fdr_frames.append(df[cols])
            if final_path.exists():
                df=pd.read_csv(final_path,sep="\t")
                if df.shape[0]:
                    df["qtl_dataset"]=dataset
                    cols=["qtl_dataset"]+[c for c in df.columns if c!="qtl_dataset"]
                    final_frames.append(df[cols])
        fdr_all=pd.concat(fdr_frames,ignore_index=True) if fdr_frames else pd.DataFrame()
        final_all=pd.concat(final_frames,ignore_index=True) if final_frames else pd.DataFrame()
        genes_union=set()
        if fdr_all.shape[0] and "GENE" in fdr_all.columns:
            genes_union |= set(fdr_all["GENE"].astype(str).str.strip().str.upper())
        if final_all.shape[0] and "GENE" in final_all.columns:
            genes_union |= set(final_all["GENE"].astype(str).str.strip().str.upper())
        fdr_all.to_csv(out_dir/"all_genes_qtl_fdr.tsv",sep="\t",index=False)
        final_all.to_csv(out_dir/"genes_qtl_causal.tsv",sep="\t",index=False)
        pd.Series(sorted(g for g in genes_union if g),dtype=str).to_csv(out_dir/"gene_list_all.txt",index=False,header=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--res_root",default="../outputs/SMR/res")
    p.add_argument("--out_root",default="../outputs/QTL_hits")
    p.add_argument("--pheno_pair",default="AD_SCZ")
    a=p.parse_args()
    compile_qtl_hits(a.res_root,a.out_root,a.pheno_pair)

if __name__=="__main__":
    main()

# python bin/compile_qtl_hits.py \
#   --res_root outputs/SMR/res \
#   --out_root outputs/QTL_hits \
#   --pheno_pair AD_SCZ