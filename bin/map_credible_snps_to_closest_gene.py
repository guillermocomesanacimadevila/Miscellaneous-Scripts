#!/usr/bin/env python3
import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
import os

# parse magma output given pheno_prefix X
# check pos of overlapping snps in credible sets
# or individual non-overlapping ->
# create a col ("nearest gene" -> gene ID)
# for each locus file
# if locus file starts like so / inside L{L}/
# find for each of those -> find cs95_{pheno_id}_L{L}.tsv
# pd.read_csv("cs95_{pheno_id}_L{L}.tsv")
# for easch snp look at teh locus coords and chr on the dir name -> bp and check on the gencode gtf
# create a column "nearest gene" ->  and map the gene name from teh gft
# save onto -> base_res/susie/res/{pheno_id}/locus_coords/cs95_{pheno_id}_L{L}_mapped.tsv

def load_gencode_pc_gtf(gtf_path: str):
    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=["seqname","source","feature","start","end","score","strand","frame","attributes"],
        low_memory=False
    )
    gtf = gtf[gtf["feature"].eq("gene")].copy()
    attrs = gtf["attributes"].astype(str).str.strip().str.split(";")
    parsed = []
    for row in attrs:
        d = {}
        for part in row:
            part = part.strip()
            if part and " " in part:
                k, v = part.split(" ", 1)
                d[k] = v.strip().strip('"')
        parsed.append(d)
    gtf["gene_type"] = [x.get("gene_type", x.get("gene_biotype")) for x in parsed]
    gtf["gene_name"] = [x.get("gene_name") for x in parsed]
    gtf = gtf[gtf["gene_type"].eq("protein_coding")].copy()
    gtf["gene_name"] = gtf["gene_name"].astype(str).str.strip()
    gtf = gtf[gtf["gene_name"].ne("") & gtf["gene_name"].str.lower().ne("nan")].copy()
    gtf["gene_name"] = gtf["gene_name"].replace({"FAM63B": "MINDY1", "FAM63A": "MINDY2"})
    gtf["seqname"] = gtf["seqname"].astype(str).str.replace("^chr", "", regex=True)
    gtf["start"] = pd.to_numeric(gtf["start"], errors="coerce")
    gtf["end"] = pd.to_numeric(gtf["end"], errors="coerce")
    gtf = gtf[["seqname","start","end","gene_name"]].dropna().sort_values(["seqname","start","end"]).reset_index(drop=True)
    return gtf

def annotate_bp_nearest_gene(gtf_pc: pd.DataFrame, chrom: str, bps: np.ndarray):
    g_chr = gtf_pc[gtf_pc["seqname"].eq(str(chrom))].reset_index(drop=True)
    starts = g_chr["start"].to_numpy(np.int64)
    ends = g_chr["end"].to_numpy(np.int64)
    names = g_chr["gene_name"].to_numpy()
    inside = (starts[:, None] <= bps[None, :]) & (bps[None, :] <= ends[:, None])
    has_inside = inside.any(axis=0)
    first_inside_idx = inside.argmax(axis=0)
    dist_left = np.abs(starts[:, None] - bps[None, :])
    dist_right = np.abs(ends[:, None] - bps[None, :])
    nearest_idx = np.minimum(dist_left, dist_right).argmin(axis=0)
    chosen_idx = np.where(has_inside, first_inside_idx, nearest_idx)
    return names[chosen_idx]

def run_univ(base_res: Path, gtf_pc: pd.DataFrame, L: str):
    for pheno_dir in sorted([p for p in (base_res / "SuSiE").iterdir() if p.is_dir() and p.name not in {"overlap", "res"}]):
        pheno_id = pheno_dir.name
        out_root = base_res / "SuSiE" / "res" / pheno_id
        for locus_dir in sorted([p for p in pheno_dir.iterdir() if p.is_dir() and p.name.startswith("locus_")]):
            chrom, lo, hi = re.match(r"locus_chr(\d+|X|Y|MT)_(\d+)_(\d+)$", locus_dir.name).groups()
            cs_path = locus_dir / f"L{L}" / f"cs95_{pheno_id}_L{L}.tsv"
            if cs_path.stat().st_size == 0:
                cs = pd.DataFrame(columns=["SNP", "nearest_gene"])
            else:
                cs = pd.read_csv(cs_path, sep="\t")
                if cs.empty:
                    cs = pd.DataFrame(columns=list(cs.columns) + ["nearest_gene"])
                else:
                    bps = pd.to_numeric(cs["BP"], errors="coerce").to_numpy(np.int64)
                    cs["nearest_gene"] = annotate_bp_nearest_gene(gtf_pc, chrom, bps)
            out_dir = out_root / locus_dir.name
            out_dir.mkdir(parents=True, exist_ok=True)
            cs.to_csv(out_dir / f"cs95_{pheno_id}_L{L}_mapped.tsv", sep="\t", index=False)

def run_biv(base_res: Path, gtf_pc: pd.DataFrame):
    overlap_dir = base_res / "SuSiE" / "overlap"
    out_root = base_res / "SuSiE" / "res" / "overlaps"
    out_root.mkdir(parents=True, exist_ok=True)
    for tsv in sorted(overlap_dir.glob("*_chr*.tsv")):
        m = re.match(r"([A-Za-z0-9]+)_([A-Za-z0-9]+)_chr(\d+|X|Y|MT)_(\d+)_(\d+)\.tsv$", tsv.name)
        pheno1, pheno2, chrom, lo, hi = m.groups()
        if tsv.stat().st_size == 0:
            df = pd.DataFrame(columns=["SNP", "nearest_gene"])
        else:
            df = pd.read_csv(tsv, sep="\t")
            if df.empty:
                df = pd.DataFrame(columns=list(df.columns) + ["nearest_gene"])
            else:
                bps = pd.to_numeric(df[f"BP_{pheno1}"], errors="coerce").to_numpy(np.int64)
                df["nearest_gene"] = annotate_bp_nearest_gene(gtf_pc, chrom, bps)
        pair_dir = out_root / f"{pheno1}_{pheno2}"
        pair_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(pair_dir / tsv.name.replace(".tsv", "_mapped.tsv"), sep="\t", index=False)

def run_one(input_path: str, gtf_pc: pd.DataFrame, out_path: str):
    input_path = Path(input_path)
    out_path = Path(out_path)

    m = re.match(r"([A-Za-z0-9]+)_([A-Za-z0-9]+)_chr(\d+|X|Y|MT)_(\d+)_(\d+)\.tsv$", input_path.name)
    pheno1, pheno2, chrom, lo, hi = m.groups()

    if input_path.stat().st_size == 0:
        df = pd.DataFrame(columns=["SNP", "nearest_gene"])
    else:
        df = pd.read_csv(input_path, sep="\t")
        if df.empty:
            df = pd.DataFrame(columns=list(df.columns) + ["nearest_gene"])
        else:
            bps = pd.to_numeric(df[f"BP_{pheno1}"], errors="coerce").to_numpy(np.int64)
            df["nearest_gene"] = annotate_bp_nearest_gene(gtf_pc, chrom, bps)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_res", default=None)
    parser.add_argument("--gtf", default="../ref/GENCODE/gencode.v37lift37.annotation.gtf")
    parser.add_argument("--L", default=None)
    parser.add_argument("--input", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    gtf_pc = load_gencode_pc_gtf(args.gtf)

    if args.input is not None and args.out is not None:
        run_one(args.input, gtf_pc, args.out)
    else:
        base_res = Path(args.base_res)
        run_univ(base_res, gtf_pc, args.L)
        run_biv(base_res, gtf_pc)

if __name__ == "__main__":
    main()