#!/usr/bin/env python3
import polars as pl
import pandas as pd
import argparse
import os
from pathlib import Path

# steps
# results/defined_loci/locus...
# grab coordinates
# look in ref/sc_eQTL/dataset/tissue/  - grab locus at those coords
# save onto /results/sc_eQTL_coloc/trait1_trait2/prep/locus_coords_tissue_cell_type
# then run coloc.R -> /results/sc_eQTL_coloc/trait1_trait2/res/locus...celltype..._coloc_results.tsv

# UPDATE:
# * => FINE-MAPPED SNPS MAP ONTO sc-eQTL => save as df with FDR correction
# * => Adapt it to Fujita´s sc-eQTLs (i.e. if dataset == ...)

ref_base = Path("../ref/sc_eQTLs")
base = Path("../results")
susie_dir = Path("../outputs/susie/res")

def prep_data_for_coloc(pheno1_prefix: str,
                        pheno2_prefix: str,
                        dataset: str,
                        out_dir: str,
                        cell_type: list[str]
                        ):
    loci_dir = Path(f"../outputs/defined_loci/{pheno1_prefix}_{pheno2_prefix}")
    if out_dir is None:
        out_dir = f"../results/sc_eQTL_coloc/{pheno1_prefix}_{pheno2_prefix}/{dataset}"
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for dir in os.listdir(loci_dir):
        full_path = loci_dir / dir
        if os.path.isdir(full_path) and dir.startswith("locus_"):
            locus_coords = "_".join(dir.split("_")[-3:])
            rows.append({
                "trait1": pheno1_prefix,
                "trait2": pheno2_prefix,
                "pair": f"{pheno1_prefix}_{pheno2_prefix}",
                "coords": locus_coords,
                "locus": dir
            })
    loci = pd.DataFrame(rows)

    selected_cell_types = [c.strip() for c in cell_type]
    if "ALL" in [c.upper() for c in selected_cell_types]:
        selected_cell_types = None

    # map those coords to the sc_eQTL dataset
    # one1k => just on 1 cell type fora given eQTL dataset
    if dataset.lower() == "onek1":
        sc_eqtl = pl.read_csv(next((ref_base / dataset.lower()).glob("*_eqtl_table.tsv")), separator="\t")

        if selected_cell_types is not None and "CELL_ID" in sc_eqtl.columns:
            sc_eqtl = sc_eqtl.filter(pl.col("CELL_ID").is_in(selected_cell_types))

        # grab coords for each locus => then reformat coloc standard
        for _, row in loci.iterrows():
            chrom, start, end = row["coords"].split("_")
            chrom = chrom.replace("chr", "")
            locus_eqtl = sc_eqtl.filter(
                (pl.col("CHR").cast(pl.Utf8) == chrom) &
                (pl.col("POS") >= int(start)) &
                (pl.col("POS") <= int(end))
            )

            if locus_eqtl.height == 0:
                continue

            if "CELL_ID" in locus_eqtl.columns:
                cell_ids = locus_eqtl["CELL_ID"].unique().to_list()
            else:
                cell_ids = ["onek1"]

            for cid in cell_ids:
                sub = locus_eqtl.filter(pl.col("CELL_ID") == cid) if "CELL_ID" in locus_eqtl.columns else locus_eqtl

                sub = sub.rename({
                    "POS": "BP",
                    "RSID": "SNP",
                    "P_VALUE": "P"
                })

                # rename cols to LDSC/coloc-like format
                # SNP A1 A2 FRQ N BETA SE P CHR BP
                sub = sub.with_columns(
                    pl.lit(982).alias("N"), # onek1
                    pl.col("A2_FREQ_ONEK1K").alias("FRQ")
                )

                base_cols = ["SNP", "A1", "A2", "FRQ", "N", "P", "CHR", "BP"]
                effect_cols = [c for c in ["BETA", "SE"] if c in sub.columns]
                sub = sub.select(base_cols + effect_cols)

                # ===> IMPORTANT - for duplicate SNPs -> pick the one with the strongest signal
                sub = (
                    sub
                    .sort("P")
                    .unique(subset="SNP")
                )
                locus_out = out_dir / f"{row['locus']}_{cid}_sc_eqtl.tsv"
                sub.write_csv(locus_out, separator="\t")

    elif dataset.lower() == "bray":
        sc_files = sorted((ref_base / dataset.lower()).glob("*_cis_eQTL_nominal.tsv")) + \
                   sorted((ref_base / dataset.lower()).glob("*_cis_eQTL_nominal.tsv.gz"))

        for sc_file in sc_files:
            cell_id = sc_file.name.replace("_cis_eQTL_nominal.tsv.gz", "").replace("_cis_eQTL_nominal.tsv", "")

            if selected_cell_types is not None and cell_id not in selected_cell_types:
                continue

            sc_eqtl = pl.read_csv(sc_file, separator="\t")

            # grab coords for each locus => then reformat coloc standard
            for _, row in loci.iterrows():
                chrom, start, end = row["coords"].split("_")
                chrom = chrom.replace("chr", "")
                locus_eqtl = sc_eqtl.filter(
                    (pl.col("CHROM").cast(pl.Utf8) == chrom) &
                    (pl.col("POS") >= int(start)) &
                    (pl.col("POS") <= int(end))
                )

                if locus_eqtl.height == 0:
                    continue

                # rename cols to LDSC/coloc-like format
                # SNP A1 A2 FRQ N BETA SE P CHR BP
                locus_eqtl = locus_eqtl.with_columns(
                    pl.lit(134).alias("N"), # bray et al.
                    pl.col("slope").alias("BETA"),
                    pl.col("AF").alias("FRQ"),
                    pl.col("slope_se").alias("SE"),
                    pl.col("pval_nominal").alias("P"),
                    pl.col("POS").alias("BP"),
                    pl.col("REF").alias("A1"),
                    pl.col("ALT").alias("A2"),
                    pl.col("CHROM").alias("CHR")
                )

                base_cols = ["SNP", "A1", "A2", "FRQ", "N", "P", "CHR", "BP"]
                effect_cols = [c for c in ["BETA", "SE"] if c in locus_eqtl.columns]
                locus_eqtl = locus_eqtl.select(base_cols + effect_cols)

                # ===> IMPORTANT - for duplicate SNPs -> pick the one with the strongest signal
                locus_eqtl = (
                    locus_eqtl
                    .sort("P")
                    .unique(subset="SNP")
                )
                locus_out = out_dir / f"{row['locus']}_{cell_id}_sc_eqtl.tsv"
                locus_eqtl.write_csv(locus_out, separator="\t")

    # Map fine-mapped SNPs to respective sc-eQTL hits
    trait_pairs = [p for p in os.listdir(susie_dir / "overlaps") if (susie_dir / "overlaps" / p).is_dir()]
    mapped_base_dir = Path("../results/sc_eQTL_coloc")
    for pair in trait_pairs:
        trait1, trait2 = pair.split("_")
        loci_dir = Path(f"../outputs/defined_loci/{pair}")
        overlap_dir = susie_dir / "overlaps" / pair
        pair_out_dir = mapped_base_dir / pair / "mapped"
        shared_out_dir = pair_out_dir / "shared"
        trait_specific_out_dir = pair_out_dir / "trait_specific"
        os.makedirs(shared_out_dir, exist_ok=True)
        os.makedirs(trait_specific_out_dir, exist_ok=True)
        loci = [d for d in os.listdir(loci_dir) if (loci_dir / d).is_dir() and d.startswith("locus_")]
        for locus in loci:
            locus_coords = "_".join(locus.split("_")[-3:])
            overlap_file = overlap_dir / f"{pair}_{locus_coords}_mapped.tsv"
            shared_snps = []
            if overlap_file.exists():
                overlap_df = pl.read_csv(overlap_file, separator="\t")
                if overlap_df.height > 0:
                    if "SNP" in overlap_df.columns:
                        shared_snps = overlap_df["SNP"].drop_nulls().to_list()
                    elif "RSID" in overlap_df.columns:
                        shared_snps = overlap_df["RSID"].drop_nulls().to_list()
                    elif "rsid" in overlap_df.columns:
                        shared_snps = overlap_df["rsid"].drop_nulls().to_list()

            if dataset.lower() == "onek1":
                sc_eqtl = pl.read_csv(next((ref_base / dataset.lower()).glob("*_eqtl_table.tsv")), separator="\t")

                if selected_cell_types is not None and "CELL_ID" in sc_eqtl.columns:
                    sc_eqtl = sc_eqtl.filter(pl.col("CELL_ID").is_in(selected_cell_types))

                if len(shared_snps) > 0:
                    shared_snps = list(set(shared_snps))
                    mapped_shared = sc_eqtl.filter(pl.col("RSID").is_in(shared_snps))

                    if mapped_shared.height > 0:
                        if "CELL_ID" in mapped_shared.columns:
                            cell_ids = mapped_shared["CELL_ID"].unique().to_list()
                        else:
                            cell_ids = ["onek1"]

                        for cid in cell_ids:
                            sub = mapped_shared.filter(pl.col("CELL_ID") == cid) if "CELL_ID" in mapped_shared.columns else mapped_shared
                            out_file = shared_out_dir / f"{locus}_{pair}_{cid}_shared_sc_eqtl.tsv"
                            sub.write_csv(out_file, separator="\t")

                else:
                    for trait in [trait1, trait2]:
                        cs_file = susie_dir / trait / locus / f"cs95_{trait}_L1_mapped.tsv"
                        if not cs_file.exists():
                            continue

                        cs_df = pl.read_csv(cs_file, separator="\t")
                        if cs_df.height == 0:
                            continue

                        trait_snps = []
                        if "SNP" in cs_df.columns:
                            trait_snps = cs_df["SNP"].drop_nulls().to_list()
                        elif "RSID" in cs_df.columns:
                            trait_snps = cs_df["RSID"].drop_nulls().to_list()
                        elif "rsid" in cs_df.columns:
                            trait_snps = cs_df["rsid"].drop_nulls().to_list()

                        trait_snps = list(set(trait_snps))
                        if len(trait_snps) == 0:
                            continue

                        mapped_trait = sc_eqtl.filter(pl.col("RSID").is_in(trait_snps))
                        if mapped_trait.height == 0:
                            continue

                        if "CELL_ID" in mapped_trait.columns:
                            cell_ids = mapped_trait["CELL_ID"].unique().to_list()
                        else:
                            cell_ids = ["onek1"]

                        for cid in cell_ids:
                            sub = mapped_trait.filter(pl.col("CELL_ID") == cid) if "CELL_ID" in mapped_trait.columns else mapped_trait
                            out_file = trait_specific_out_dir / f"{locus}_{trait}_{cid}_trait_specific_sc_eqtl.tsv"
                            sub.write_csv(out_file, separator="\t")

            elif dataset.lower() == "bray":
                sc_files = sorted((ref_base / dataset.lower()).glob("*_cis_eQTL_nominal.tsv")) + \
                           sorted((ref_base / dataset.lower()).glob("*_cis_eQTL_nominal.tsv.gz"))

                for sc_file in sc_files:
                    cell_id = sc_file.name.replace("_cis_eQTL_nominal.tsv.gz", "").replace("_cis_eQTL_nominal.tsv", "")

                    if selected_cell_types is not None and cell_id not in selected_cell_types:
                        continue

                    sc_eqtl = pl.read_csv(sc_file, separator="\t")

                    if len(shared_snps) > 0:
                        shared_snps = list(set(shared_snps))
                        mapped_shared = sc_eqtl.filter(pl.col("SNP").is_in(shared_snps))

                        if mapped_shared.height > 0:
                            out_file = shared_out_dir / f"{locus}_{pair}_{cell_id}_shared_sc_eqtl.tsv"
                            mapped_shared.write_csv(out_file, separator="\t")

                    else:
                        for trait in [trait1, trait2]:
                            cs_file = susie_dir / trait / locus / f"cs95_{trait}_L1_mapped.tsv"
                            if not cs_file.exists():
                                continue

                            cs_df = pl.read_csv(cs_file, separator="\t")
                            if cs_df.height == 0:
                                continue

                            trait_snps = []
                            if "SNP" in cs_df.columns:
                                trait_snps = cs_df["SNP"].drop_nulls().to_list()
                            elif "RSID" in cs_df.columns:
                                trait_snps = cs_df["RSID"].drop_nulls().to_list()
                            elif "rsid" in cs_df.columns:
                                trait_snps = cs_df["rsid"].drop_nulls().to_list()

                            trait_snps = list(set(trait_snps))
                            if len(trait_snps) == 0:
                                continue

                            mapped_trait = sc_eqtl.filter(pl.col("SNP").is_in(trait_snps))
                            if mapped_trait.height == 0:
                                continue

                            out_file = trait_specific_out_dir / f"{locus}_{trait}_{cell_id}_trait_specific_sc_eqtl.tsv"
                            mapped_trait.write_csv(out_file, separator="\t")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1_prefix", required=True, type=str)
    parser.add_argument("--pheno2_prefix", required=True, type=str)
    parser.add_argument("--dataset", required=True, type=str)
    parser.add_argument("--out_dir", required=False, type=str, default=None)
    parser.add_argument("--cell_type", nargs="+", required=False, default=["ALL"], help="Cell type(s) to process, e.g. MG OPC Endo-Peri, or ALL")
    args = parser.parse_args()
    prep_data_for_coloc(
        pheno1_prefix=args.pheno1_prefix,
        pheno2_prefix=args.pheno2_prefix,
        dataset=args.dataset,
        out_dir=args.out_dir,
        cell_type=args.cell_type
    )

if __name__ == "__main__":
    main()

