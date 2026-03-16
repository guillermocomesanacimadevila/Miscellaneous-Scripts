#!/usr/bin/env python3
import argparse
import os
import gzip
import shutil
from pathlib import Path
import pandas as pd

base_dir = Path("../results")
loci_dir = Path("../results/defined_loci")

def prep_for_finemap(pheno1_id: str, pheno2_id: str):
    dir = Path(loci_dir / f"{pheno1_id}_{pheno2_id}")
    out_dir = Path(base_dir / "finemap" / f"{pheno1_id}_{pheno2_id}")
    out_dir.mkdir(parents=True, exist_ok=True)
    for d in os.listdir(dir):
        full_path = dir / d
        if os.path.isdir(full_path) and d.startswith("locus_"):
            locus_coords = "_".join(d.split("_")[-3:])
            p1 = None
            p2 = None
            for file in os.listdir(full_path):
                fp = full_path / file
                if file.startswith(f"gwas_{pheno1_id}") and file.endswith(".tsv"):
                    p1 = pd.read_csv(fp, sep="\t")
                if file.startswith(f"gwas_{pheno2_id}") and file.endswith(".tsv"):
                    p2 = pd.read_csv(fp, sep="\t")
            locus_out = out_dir / locus_coords
            locus_out.mkdir(parents=True, exist_ok=True)
            ld_gz = full_path / "ld_gwas" / "locus.ld.gz"
            if ld_gz.exists():
                with gzip.open(ld_gz, "rt") as fin, open(locus_out / "locus.ld", "wt") as fout:
                    shutil.copyfileobj(fin, fout)
            if p1 is not None:
                p1 = p1.rename(columns={
                    "SNP": "rsid",
                    "CHR": "chromosome",
                    "A1": "allele1",
                    "BP": "position",
                    "A2": "allele2",
                    "FRQ": "maf",
                    "BETA": "beta",
                    "SE": "se"
                })
                p1["maf"] = p1["maf"].astype(float)
                p1["maf"] = p1["maf"].where(p1["maf"] <= 0.5, 1.0 - p1["maf"])
                p1 = p1[["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "N"]]
                n1 = int(p1["N"].median())
                p1 = p1.drop(columns=["N"])
                p1.to_csv(locus_out / f"{pheno1_id}.z", sep=" ", index=False)
                with open(locus_out / f"{pheno1_id}.master", "w") as f:
                    f.write("z;ld;n_samples\n")
                    f.write(f"{pheno1_id}.z;locus.ld;{n1}\n")
            if p2 is not None:
                p2 = p2.rename(columns={
                    "SNP": "rsid",
                    "CHR": "chromosome",
                    "A1": "allele1",
                    "BP": "position",
                    "A2": "allele2",
                    "FRQ": "maf",
                    "BETA": "beta",
                    "SE": "se"
                })
                p2["maf"] = p2["maf"].astype(float)
                p2["maf"] = p2["maf"].where(p2["maf"] <= 0.5, 1.0 - p2["maf"])
                p2 = p2[["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "N"]]
                n2 = int(p2["N"].median())
                p2 = p2.drop(columns=["N"])
                p2.to_csv(locus_out / f"{pheno2_id}.z", sep=" ", index=False)
                with open(locus_out / f"{pheno2_id}.master", "w") as f:
                    f.write("z;ld;n_samples\n")
                    f.write(f"{pheno2_id}.z;locus.ld;{n2}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno1_id", required=True)
    parser.add_argument("--pheno2_id", required=True)
    args = parser.parse_args()
    prep_for_finemap(args.pheno1_id, args.pheno2_id)

if __name__ == "__main__":
    main()