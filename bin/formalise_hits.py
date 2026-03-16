#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def formalise_hits(pheno1_id:str,pheno2_id:str,qtl_type:str,res_dir:str,p_heidi:float,fdr:float):
    base=Path(res_dir)/qtl_type/f"{pheno1_id}_{pheno2_id}"
    if not base.exists(): raise FileNotFoundError(str(base))
    keep_base=["qtl_name","probeID","ProbeChr","Probe_bp","topSNP","topSNP_chr","topSNP_bp","A1","A2","Freq","b_SMR","se_SMR","p_SMR","p_SMR_multi","q_SMR_multi","p_HEIDI","nsnp_HEIDI","gene_id"]
    keep_a=["GENE"]+[f"{c}_A" for c in keep_base]
    keep_s=[f"{c}_S" for c in keep_base]
    tidy_cols=["GENE"]+[f"{c}_A" for c in keep_base]+[f"{c}_S" for c in keep_base]
    for qtl_dir in sorted([p for p in base.iterdir() if p.is_dir()]):
        f1=qtl_dir/f"{pheno1_id}_raw_with_q.tsv"
        f2=qtl_dir/f"{pheno2_id}_raw_with_q.tsv"
        if (not f1.exists()) or (not f2.exists()):
            continue
        t1=pd.read_csv(f1,sep="\t")
        t2=pd.read_csv(f2,sep="\t")
        t1["GENE"]=t1["index"].astype(str).str.strip().str.upper()
        t2["GENE"]=t2["index"].astype(str).str.strip().str.upper()
        t1=t1.dropna(subset=["GENE"]).copy()
        t2=t2.dropna(subset=["GENE"]).copy()
        overlap=sorted(set(t1["GENE"])&set(t2["GENE"]))
        if len(overlap)==0:
            (qtl_dir/f"FINAL_{pheno1_id}_{pheno2_id}_{qtl_dir.name}_hits_past_heidi.tsv").write_text("")
            (qtl_dir/f"FINAL_{pheno1_id}_{pheno2_id}_{qtl_dir.name}_hits_past_fdr_only.tsv").write_text("")
            continue
        a=t1[t1["GENE"].isin(overlap)].copy()
        b=t2[t2["GENE"].isin(overlap)].copy()
        a["q_SMR_multi"]=pd.to_numeric(a.get("q_SMR_multi"),errors="coerce")
        b["q_SMR_multi"]=pd.to_numeric(b.get("q_SMR_multi"),errors="coerce")
        a["p_HEIDI"]=pd.to_numeric(a.get("p_HEIDI"),errors="coerce")
        b["p_HEIDI"]=pd.to_numeric(b.get("p_HEIDI"),errors="coerce")
        a=a.sort_values(["GENE","q_SMR_multi"],na_position="last").drop_duplicates("GENE",keep="first")
        b=b.sort_values(["GENE","q_SMR_multi"],na_position="last").drop_duplicates("GENE",keep="first")
        m=a.merge(b,on="GENE",how="inner",suffixes=("_A","_S"))
        keep_fdr=(pd.to_numeric(m["q_SMR_multi_A"],errors="coerce")<float(fdr))&(pd.to_numeric(m["q_SMR_multi_S"],errors="coerce")<float(fdr))
        keep_heidi_both=keep_fdr&(pd.to_numeric(m["p_HEIDI_A"],errors="coerce")>float(p_heidi))&(pd.to_numeric(m["p_HEIDI_S"],errors="coerce")>float(p_heidi))
        keep_heidi_any=keep_fdr&((pd.to_numeric(m["p_HEIDI_A"],errors="coerce")>float(p_heidi))|(pd.to_numeric(m["p_HEIDI_S"],errors="coerce")>float(p_heidi)))
        fdr_only=m[keep_heidi_any].copy()
        heidi_both=m[keep_heidi_both].copy()
        for df in (fdr_only,heidi_both):
            cols=[c for c in tidy_cols if c in df.columns]
            df=df[cols].copy()
        fdr_only=fdr_only[[c for c in tidy_cols if c in fdr_only.columns]].copy()
        heidi_both=heidi_both[[c for c in tidy_cols if c in heidi_both.columns]].copy()
        fdr_only.to_csv(qtl_dir/f"FINAL_{pheno1_id}_{pheno2_id}_{qtl_dir.name}_hits_past_fdr_only.tsv",sep="\t",index=False)
        heidi_both.to_csv(qtl_dir/f"FINAL_{pheno1_id}_{pheno2_id}_{qtl_dir.name}_hits_past_heidi.tsv",sep="\t",index=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pheno1_id",required=True)
    p.add_argument("--pheno2_id",required=True)
    p.add_argument("--qtl_type",required=True)
    p.add_argument("--res_dir",required=True)
    p.add_argument("--p_heidi",type=float,default=0.01)
    p.add_argument("--fdr",type=float,default=0.05)
    a=p.parse_args()
    formalise_hits(a.pheno1_id,a.pheno2_id,a.qtl_type,a.res_dir,a.p_heidi,a.fdr)

if __name__=="__main__":
    main()