#!/usr/bin/env python3
import argparse
import pandas as pd
import subprocess
from pathlib import Path

# grab output of conjFDR
# grab snp hits,
# two dataframes (one AD, SCZ, - triple overlap = + LON)
# so 2 dfs for each snp based on pairwise traits
# for snp in sumstats_trait1:
#   clump r2 >= 0.6 per 1000Kb
# define lead SNP of locus
# around lead snp keep for each trait
# r2 >= 0.1 with lead over +/- 250Kb
# for each conjFDR SNP hit
# If +/- 500kb of other SNPs that were hits
# Calculate r^2 for those, if r^2 > 0.6 then consider them the same locus and don´t clump them separately into distinct loci
# Saving files -> for each locus save it into pheno_1-pheno_2/ directory and then create another dir -> locus_id (i.e. pheno1_id_chr_bpstart_bp_end) and pheno2_id_chr_bpstart_bp_end

def read_tsv(p):
    return pd.read_csv(p, sep="\t")

def write_tsv(df, p):
    df.to_csv(p, sep="\t", index=False)

def ensure_cols(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError("Missing required columns: " + ",".join(missing))

def make_score_pairwise(df, score_col):
    ensure_cols(df, ["SNP", score_col])
    out = df[["SNP", score_col]].copy()
    out.columns = ["SNP", "P"]
    return out

def make_score_triple(df, c1, c2, how):
    ensure_cols(df, ["SNP", c1, c2])
    x = df[[c1, c2]].copy()
    if how == "max":
        p = x.max(axis=1)
    else:
        p = x.min(axis=1)
    out = pd.DataFrame({"SNP": df["SNP"], "P": p})
    return out

def run(cmd, cwd=None):
    subprocess.run(cmd, check=True, cwd=cwd)

def plink_clump(plink, ref_bfile, clump_in, out_prefix, r2, kb):
    cmd = [
        plink,
        "--bfile", ref_bfile,
        "--clump", str(clump_in),
        "--clump-field", "P",
        "--clump-p1", "1",
        "--clump-p2", "1",
        "--clump-r2", str(r2),
        "--clump-kb", str(kb),
        "--out", str(out_prefix),
    ]
    run(cmd)

def read_leads(clumped_path):
    df = pd.read_csv(clumped_path, sep=r"\s+")
    ensure_cols(df, ["SNP", "CHR", "BP"])
    df = df[["SNP", "CHR", "BP"]].dropna().drop_duplicates()
    df = df.reset_index(drop=True)
    df["locus_id"] = range(len(df))
    return df

def plink_ld(plink, ref_bfile, lead_snp, out_prefix, kb):
    cmd = [
        plink,
        "--bfile", ref_bfile,
        "--r2",
        "--ld-snp", lead_snp,
        "--ld-window-kb", str(kb),
        "--ld-window", "999999",
        "--ld-window-r2", "0",
        "--out", str(out_prefix),
    ]
    run(cmd)

def plink_ld_list(plink, ref_bfile, snp_list_path, out_prefix, kb, r2_min):
    cmd = [
        plink,
        "--bfile", ref_bfile,
        "--r2",
        "--ld-snp-list", str(snp_list_path),
        "--ld-window-kb", str(kb),
        "--ld-window", "999999",
        "--ld-window-r2", str(r2_min),
        "--out", str(out_prefix),
    ]
    run(cmd)

def parse_ld_file(ld_path, r2_min):
    df = pd.read_csv(ld_path, sep=r"\s+")
    ensure_cols(df, ["SNP_A", "SNP_B", "CHR_A", "BP_A", "CHR_B", "BP_B", "R2"])
    df = df[df["R2"] >= r2_min].copy()
    if df.empty:
        return None, None
    snps = pd.concat([df["SNP_A"], df["SNP_B"]]).dropna().drop_duplicates()
    chr_vals = pd.concat([df["CHR_A"], df["CHR_B"]]).dropna().unique()
    chr_val = chr_vals[0]
    bp_min = int(min(df["BP_A"].min(), df["BP_B"].min()))
    bp_max = int(max(df["BP_A"].max(), df["BP_B"].max()))
    return snps, (chr_val, bp_min, bp_max)

def extract_sumstats(sumstats_path, snps):
    df = read_tsv(sumstats_path)
    ensure_cols(df, ["SNP"])
    keep = set(snps.tolist())
    out = df[df["SNP"].isin(keep)].copy()
    return out

def uf_find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def uf_union(parent, rank, a, b):
    ra = uf_find(parent, a)
    rb = uf_find(parent, b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    if rank[ra] == rank[rb]:
        rank[ra] += 1

def merge_leads_by_ld(leads, hits, score_col, plink, ref_bfile, out_dir, kb, r2_min):
    leads = leads.copy()
    score_map = hits.set_index("SNP")[score_col].to_dict()
    leads["score"] = leads["SNP"].map(score_map)
    snp_list = out_dir / "lead_snps.list"
    snp_list.write_text("\n".join(leads["SNP"].astype(str).tolist()) + "\n")
    pref = out_dir / "lead_snps_ld"
    plink_ld_list(plink, ref_bfile, snp_list, pref, kb=kb, r2_min=r2_min)
    ld_path = Path(str(pref) + ".ld")
    if not ld_path.exists():
        leads = leads[["SNP", "CHR", "BP"]].drop_duplicates().reset_index(drop=True)
        leads["locus_id"] = range(len(leads))
        return leads
    ld = pd.read_csv(ld_path, sep=r"\s+")
    ensure_cols(ld, ["SNP_A", "SNP_B", "R2", "CHR_A", "CHR_B", "BP_A", "BP_B"])
    ld = ld[(ld["CHR_A"] == ld["CHR_B"]) & ((ld["BP_A"] - ld["BP_B"]).abs() <= kb * 1000)]
    if ld.empty:
        leads = leads[["SNP", "CHR", "BP"]].drop_duplicates().reset_index(drop=True)
        leads["locus_id"] = range(len(leads))
        return leads
    snps = leads["SNP"].astype(str).tolist()
    parent = {s: s for s in snps}
    rank = {s: 0 for s in snps}
    for _, r in ld.iterrows():
        if float(r["R2"]) >= r2_min:
            a = str(r["SNP_A"])
            b = str(r["SNP_B"])
            if a in parent and b in parent:
                uf_union(parent, rank, a, b)
    groups = {}
    for s in snps:
        root = uf_find(parent, s)
        if root not in groups:
            groups[root] = []
        groups[root].append(s)
    reps = []
    for root, gs in groups.items():
        sub = leads[leads["SNP"].isin(gs)].copy()
        sub["score"] = pd.to_numeric(sub["score"], errors="coerce")
        sub = sub.sort_values(["score", "BP"], ascending=[True, True])
        reps.append(sub.iloc[0][["SNP", "CHR", "BP"]])
    reps = pd.DataFrame(reps).drop_duplicates().reset_index(drop=True)
    reps["locus_id"] = range(len(reps))
    return reps

def locus_dirname(pheno1, pheno2, chr_val, start, end, locus_id):
    return f"locus_{locus_id}_{pheno1}_{pheno2}_chr{chr_val}_{start}_{end}"

def subset_hits(hits, snps):
    ensure_cols(hits, ["SNP"])
    keep = set(snps.tolist())
    return hits[hits["SNP"].isin(keep)].copy()

def pick_lead(df, col, label):
    if col not in df.columns:
        return None
    x = pd.to_numeric(df[col], errors="coerce")
    if x.isna().all():
        return None
    i = x.idxmin()
    return {"label": label, "score_col": col, "lead_snp": str(df.loc[i, "SNP"]), "lead_score": float(x.loc[i])}

def write_leads(df_hits, out_path, mode, args):
    rows = []
    if mode == "pairwise":
        r = pick_lead(df_hits, args.score_col, "lead_pairwise_conj")
        if r is not None:
            rows.append(r)
        r = pick_lead(df_hits, "cfdr12", "lead_pairwise_cfdr12")
        if r is not None:
            rows.append(r)
        r = pick_lead(df_hits, "cfdr21", "lead_pairwise_cfdr21")
        if r is not None:
            rows.append(r)
    else:
        r = pick_lead(df_hits, "P", "lead_triple_score")
        if r is not None:
            rows.append(r)
        r = pick_lead(df_hits, args.triple_col1, "lead_AD_LON")
        if r is not None:
            rows.append(r)
        r = pick_lead(df_hits, args.triple_col2, "lead_SCZ_LON")
        if r is not None:
            rows.append(r)
    if rows:
        write_tsv(pd.DataFrame(rows), out_path)
    else:
        write_tsv(pd.DataFrame(columns=["label","score_col","lead_snp","lead_score"]), out_path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pairwise", "triple"], required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pheno1", required=True)
    ap.add_argument("--pheno2", required=True)
    ap.add_argument("--plink", default="plink")
    ap.add_argument("--ref-bfile", required=True)
    ap.add_argument("--clump-r2", type=float, default=0.6)
    ap.add_argument("--clump-kb", type=int, default=1000)
    ap.add_argument("--ld-r2", type=float, default=0.1)
    ap.add_argument("--ld-kb", type=int, default=250)
    ap.add_argument("--merge-kb", type=int, default=500)
    ap.add_argument("--score-col", default="conj_fdr")
    ap.add_argument("--triple-col1", default="conj_fdr_AD_LON")
    ap.add_argument("--triple-col2", default="conj_fdr_SCZ_LON")
    ap.add_argument("--triple-how", choices=["max", "min"], default="max")
    ap.add_argument("--sumstats", nargs="*", default=[])
    ap.add_argument("--sumstats-names", nargs="*", default=[])
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_dir = out_dir / f"{args.pheno1}-{args.pheno2}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    hits = read_tsv(args.hits).drop_duplicates(subset=["SNP"]).copy()
    if args.mode == "pairwise":
        clump_df = make_score_pairwise(hits, args.score_col)
        score_col = args.score_col
    else:
        clump_df = make_score_triple(hits, args.triple_col1, args.triple_col2, args.triple_how)
        hits = hits.merge(clump_df, on="SNP", how="left")
        score_col = "P"
    clump_in = pair_dir / "clump_input.tsv"
    write_tsv(clump_df, clump_in)
    clump_out = pair_dir / "clump"
    plink_clump(args.plink, args.ref_bfile, clump_in, clump_out, args.clump_r2, args.clump_kb)
    clumped_path = Path(str(clump_out) + ".clumped")
    leads = read_leads(clumped_path)
    leads = merge_leads_by_ld(
        leads=leads,
        hits=hits,
        score_col=score_col,
        plink=args.plink,
        ref_bfile=args.ref_bfile,
        out_dir=pair_dir,
        kb=args.merge_kb,
        r2_min=args.clump_r2
    )
    write_tsv(leads, pair_dir / "lead_snps.tsv")
    loci_rows = []
    for _, r in leads.iterrows():
        locus_id = int(r["locus_id"])
        lead = str(r["SNP"])
        pref = pair_dir / f"tmp_locus_{locus_id}"
        plink_ld(args.plink, args.ref_bfile, lead, pref, args.ld_kb)
        ld_path = Path(str(pref) + ".ld")
        snps, coords = parse_ld_file(ld_path, args.ld_r2)
        if snps is None:
            continue
        chr_val, start, end = coords
        loc_dir = pair_dir / locus_dirname(args.pheno1, args.pheno2, chr_val, start, end, locus_id)
        loc_dir.mkdir(parents=True, exist_ok=True)
        loci_rows.append({
            "locus_id": locus_id,
            "lead_snp": lead,
            "chr": chr_val,
            "start": start,
            "end": end,
            "n_snps": int(len(snps))
        })
        write_tsv(pd.DataFrame({"SNP": snps}), loc_dir / "snps.tsv")
        loc_hits = subset_hits(hits, snps)
        write_tsv(loc_hits, loc_dir / "hits_annot.tsv")
        write_leads(loc_hits, loc_dir / "lead_snps.tsv", args.mode, args)
        if len(args.sumstats) == len(args.sumstats_names) and len(args.sumstats) > 0:
            for p, name in zip(args.sumstats, args.sumstats_names):
                loc = extract_sumstats(p, snps)
                write_tsv(loc, loc_dir / f"gwas_{name}.tsv")
    if len(loci_rows) > 0:
        write_tsv(pd.DataFrame(loci_rows), pair_dir / "locus_coords.tsv")
    else:
        write_tsv(pd.DataFrame(columns=["locus_id","lead_snp","chr","start","end","n_snps"]), pair_dir / "locus_coords.tsv")

if __name__ == "__main__":
    main()
