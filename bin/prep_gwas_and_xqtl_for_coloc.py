#!/usr/bin/env python3
import argparse
import pandas as pd
import polars as pl
import os
from pathlib import Path

# => map lead snp & gene from correspondant SMR hit
# => map it on xQTL - and GWAS => 250kb up and down from it on GWAS / xQTL
# => on xQTL make sure Gene == ADAM10 => Then overlap xQTL SNPs with the GWAS ones
# => save both

# grab snp from the ld files
# results/ld_results
# AD_SCZ_eQTL_chr11_46732038_47732038_smr_vs_susie_ld.tsv
# if file != exists
# print => Yowza! No xQTL hit for that combo! Please specify a genuine SMR hit

def check_smr_hit(pheno1_id: str,
                  pheno2_id: str,
                  ld_dir: str,
                  qtl_type: str,
                  chr,
                  start_locus,
                  end_locus):
    ld_dir = Path(ld_dir)
    target_file = f"{pheno1_id}_{pheno2_id}_{qtl_type}_chr{chr}_{start_locus}_{end_locus}_smr_vs_susie_ld.tsv"
    for file in os.listdir(ld_dir):
        if file == target_file:
            df = pd.read_csv(ld_dir / file, sep="\t")
            print("Found your SMR hit!")
            return df

    print("Yowza! No xQTL hit for that combo. Please specify a genuine SMR hit")
    return None

# need a way to check the qtl dataset
# dataset specific!
# save preproessed to =>
# results/trait1_trait2/qtl_type/data_prep/qtl_dataset/locus_id/...QTL - trait1...sumstats . trait2...sumstats

def assemble_coloc_inputs(pheno1_id: str,
                          pheno2_id: str,
                          ld_dir: str,
                          qtl_type: str,
                          qtl_dataset: str,
                          gwas_loci_dir: str,
                          qtl_loci_dir: str,
                          chr,
                          smr_dir: str,
                          start_locus,
                          end_locus,
                          window: int,
                          out_dir: str):
    ld_dir = Path(ld_dir)
    qtl_loci_dir = Path(qtl_loci_dir)
    gwas_loci_dir = Path(gwas_loci_dir)
    smr_dir = Path(smr_dir)
    out_dir = Path(out_dir) / f"{pheno1_id}_{pheno2_id}/{qtl_type}/data_prep/{qtl_dataset}/chr{chr}_{start_locus}_{end_locus}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # grab snp dataset
    df = check_smr_hit(
        pheno1_id=pheno1_id,
        pheno2_id=pheno2_id,
        ld_dir=ld_dir,
        qtl_type=qtl_type,
        chr=chr,
        start_locus=start_locus,
        end_locus=end_locus
    )

    if df is None:
        return None

    print("SMR hit checked successfully")

    # map gene onto SMR hit file
    smr_res = pd.read_csv(
        smr_dir / f"{pheno1_id}_{pheno2_id}/{qtl_type}/chr{chr}_{start_locus}_{end_locus}/{qtl_type}_{qtl_dataset}_{pheno1_id}_{pheno2_id}_shared.tsv",
        sep="\t"
    )

    # shared_gene
    gene_map_1 = dict(zip(smr_res[f"topSNP_{pheno1_id}"].astype(str), smr_res["shared_gene"]))
    gene_map_2 = dict(zip(smr_res[f"topSNP_{pheno2_id}"].astype(str), smr_res["shared_gene"]))
    df["Gene"] = df["SNP_A"].astype(str).map(gene_map_1)
    df["Gene"] = df["Gene"].fillna(df["SNP_A"].astype(str).map(gene_map_2))

    # grab lead SNP (i.e. SMR hit)
    df_hits = df[["SNP_A", "Gene"]].dropna().drop_duplicates().copy()  # all SMR top xQTL for that locus

    # gwas dfs for that corresponding locus
    pheno1_df = pd.read_csv(
        gwas_loci_dir / f"{pheno1_id}_{pheno2_id}/locus_chr{chr}_{start_locus}_{end_locus}/gwas_{pheno1_id}.ldgwas.tsv",
        sep="\t"
    )
    pheno2_df = pd.read_csv(
        gwas_loci_dir / f"{pheno1_id}_{pheno2_id}/locus_chr{chr}_{start_locus}_{end_locus}/gwas_{pheno2_id}.ldgwas.tsv",
        sep="\t"
    )
    qtl_df = pd.read_csv(
        qtl_loci_dir / f"{qtl_type}_{qtl_dataset}_chr{chr}.txt",
        sep="\t"
    )  # sQTL_BrainMeta_chr15.txt

    # now for SMR hit lead => map xkb up and down from it => and make sure that the gene is the same as mapped one
    # also, overlap SNPs between xQTL for that CHR and both GWAS for traits 1 and 2

    qtl_df["BP"] = pd.to_numeric(qtl_df["BP"], errors="coerce")
    pheno1_df["BP"] = pd.to_numeric(pheno1_df["BP"], errors="coerce")
    pheno2_df["BP"] = pd.to_numeric(pheno2_df["BP"], errors="coerce")
    saved = []
    for _, row in df_hits.iterrows():
        smr_snp = row["SNP_A"]
        gene = row["Gene"]
        print(f"Working on SMR hit {smr_snp} => Gene: {gene}")
        qtl_snp_row = qtl_df[qtl_df["SNP"].astype(str) == str(smr_snp)].copy()

        if qtl_snp_row.empty:
            print(f"Could not find {smr_snp} in xQTL dataset")
            continue

        smr_bp = int(qtl_snp_row["BP"].iloc[0])
        region_start = smr_bp - window
        region_end = smr_bp + window
        print(f"Window for {smr_snp}: chr{chr}:{region_start}-{region_end}")

        qtl_sub = qtl_df[
            (qtl_df["BP"] >= region_start) &
            (qtl_df["BP"] <= region_end) &
            (qtl_df["Gene"].astype(str) == str(gene))
        ].copy()

        if qtl_sub.empty:
            print(f"No xQTL rows found for {gene} around {smr_snp}")
            continue

        # overlaps
        pheno1_sub = pheno1_df[
            (pheno1_df["BP"] >= region_start) &
            (pheno1_df["BP"] <= region_end)
        ].copy()

        pheno2_sub = pheno2_df[
            (pheno2_df["BP"] >= region_start) &
            (pheno2_df["BP"] <= region_end)
        ].copy()

        shared_snps_1 = set(qtl_sub["SNP"].astype(str)).intersection(set(pheno1_sub["SNP"].astype(str)))
        qtl_pheno1 = qtl_sub[qtl_sub["SNP"].astype(str).isin(shared_snps_1)].copy()
        pheno1_overlap = pheno1_sub[pheno1_sub["SNP"].astype(str).isin(shared_snps_1)].copy()
        shared_snps_2 = set(qtl_sub["SNP"].astype(str)).intersection(set(pheno2_sub["SNP"].astype(str)))
        qtl_pheno2 = qtl_sub[qtl_sub["SNP"].astype(str).isin(shared_snps_2)].copy()
        pheno2_overlap = pheno2_sub[pheno2_sub["SNP"].astype(str).isin(shared_snps_2)].copy()

        qtl_pheno1["smr_snp"] = smr_snp
        qtl_pheno1["Gene"] = gene
        pheno1_overlap["smr_snp"] = smr_snp
        pheno1_overlap["Gene"] = gene

        qtl_pheno2["smr_snp"] = smr_snp
        qtl_pheno2["Gene"] = gene
        pheno2_overlap["smr_snp"] = smr_snp
        pheno2_overlap["Gene"] = gene

        # save
        qtl_pheno1.to_csv(out_dir / f"{qtl_type}_{qtl_dataset}_{pheno1_id}_{gene}_{smr_snp}_sumstats.tsv", sep="\t", index=False)
        pheno1_overlap.to_csv(out_dir / f"gwas_{pheno1_id}_{gene}_{smr_snp}_sumstats.tsv", sep="\t", index=False)
        qtl_pheno2.to_csv(out_dir / f"{qtl_type}_{qtl_dataset}_{pheno2_id}_{gene}_{smr_snp}_sumstats.tsv", sep="\t", index=False)
        pheno2_overlap.to_csv(out_dir / f"gwas_{pheno2_id}_{gene}_{smr_snp}_sumstats.tsv", sep="\t", index=False)
        print(f"Saved coloc inputs for {smr_snp} | {gene}")
        saved.append({
            "smr_snp": smr_snp,
            "Gene": gene,
            "region_start": region_start,
            "region_end": region_end,
            f"n_qtl_{pheno1_id}": qtl_pheno1.shape[0],
            f"n_gwas_{pheno1_id}": pheno1_overlap.shape[0],
            f"n_qtl_{pheno2_id}": qtl_pheno2.shape[0],
            f"n_gwas_{pheno2_id}": pheno2_overlap.shape[0]
        })

    summary_df = pd.DataFrame(saved)
    summary_df.to_csv(
        out_dir / f"{pheno1_id}_{pheno2_id}_{qtl_type}_{qtl_dataset}_coloc_summary.tsv",
        sep="\t",
        index=False
    )
    return summary_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1_id", required=True)
    parser.add_argument("--pheno2_id", required=True)
    parser.add_argument("--ld_dir", required=True)
    parser.add_argument("--qtl_type", required=True)
    parser.add_argument("--qtl_dataset", required=True)
    parser.add_argument("--gwas_loci_dir", required=True)
    parser.add_argument("--qtl_loci_dir", required=True)
    parser.add_argument("--chr", required=True, type=int)
    parser.add_argument("--smr_dir", required=True)
    parser.add_argument("--start_locus", required=True, type=int)
    parser.add_argument("--end_locus", required=True, type=int)
    parser.add_argument("--window", required=True, type=int)
    parser.add_argument("--out_dir", required=True)
    args = parser.parse_args()
    assemble_coloc_inputs(
        pheno1_id=args.pheno1_id,
        pheno2_id=args.pheno2_id,
        ld_dir=args.ld_dir,
        qtl_type=args.qtl_type,
        qtl_dataset=args.qtl_dataset,
        gwas_loci_dir=args.gwas_loci_dir,
        qtl_loci_dir=args.qtl_loci_dir,
        chr=args.chr,
        smr_dir=args.smr_dir,
        start_locus=args.start_locus,
        end_locus=args.end_locus,
        window=args.window,
        out_dir=args.out_dir
    )

if __name__ == "__main__":
    main()