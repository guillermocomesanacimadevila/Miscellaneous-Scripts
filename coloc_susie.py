import polars as pl
import pandas as pd
from pathlib import Path
import subprocess
from drugmr.utils import filter_mr_targets
from drugmr import paths
import argparse


def get_ld_matrix(
        file: str,
        id: str,
        ref_prefix: str,
        maf: float,
        out_dir: str):

    script = "bin/ld_matrix.sh"
    cmd = ["bash", script, str(file), id, str(ref_prefix), str(maf), str(out_dir)]
    print(f"[TRACKING] Running LD matrix for {id}")
    subprocess.run(cmd, check=True)


def coloc_susie(
        pqtl_dataset: str,
        protein_id: str,
        pheno_id: str,
        gwas_file: str,
        pqtl_file: str,
        n_cases: int,
        n_controls: int,
        local_results_dir: str = "results"):

    script = "bin/coloc_susie.R"
    cmd = ["Rscript", script, pqtl_dataset, protein_id, pheno_id, str(gwas_file), str(pqtl_file), str(n_cases), str(n_controls), str(local_results_dir)]
    print(f"[TRACKING] Running COLOC SuSiE for {protein_id}")
    subprocess.run(cmd, check=True)


def grab_ld_matrix_for_cis_mr_hits_and_coloc_susie(
        pqtl_dataset: str,
        pheno_id: str,
        pqtl_dir: str,
        maf: float,
        ref_bfile: str,
        n_cases: int,
        n_controls: int,
        local_results_dir: str = "results"):

    pqtl_res = pl.read_csv(paths.mr_out(pqtl_dataset, pheno_id, local_results_dir), separator="\t")
    targets = filter_mr_targets(pqtl_res)
    print(f"[TRACKING] Proteins passing cis-MR filters: {len(targets)}")
    pqtl_dir = Path(pqtl_dir)
    out_dir = paths.coloc_susie_out(pqtl_dataset, pheno_id, local_results_dir).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []

    for protein in targets:
        protein_dir = pqtl_dir / protein
        gwas_file = protein_dir / "gwas.parquet"
        pqtl_file = protein_dir / "pqtl.parquet"
        protein_file = out_dir / f"{pheno_id}_{protein}_coloc_susie.tsv"

        try:
            # LD matrices are the expensive step (plink over the ref panel + a big
            # fread in R) - only run ld_matrix.sh again if not already there
            gwas_ld = protein_dir / f"{pheno_id}_ld_matrix.ld"
            pqtl_ld = protein_dir / f"{protein}_ld_matrix.ld"

            if gwas_ld.exists():
                print(f"[TRACKING] GWAS LD matrix already exists for {protein}, skipping")
            else:
                get_ld_matrix(
                    file=str(gwas_file),
                    id=pheno_id,
                    ref_prefix=ref_bfile,
                    maf=maf,
                    out_dir=str(protein_dir),
                )

            if pqtl_ld.exists():
                print(f"[TRACKING] pQTL LD matrix already exists for {protein}, skipping")
            else:
                get_ld_matrix(
                    file=str(pqtl_file),
                    id=protein,
                    ref_prefix=ref_bfile,
                    maf=maf,
                    out_dir=str(protein_dir),
                )

            coloc_susie(
                pqtl_dataset=pqtl_dataset,
                protein_id=protein,
                pheno_id=pheno_id,
                gwas_file=str(gwas_file),
                pqtl_file=str(pqtl_file),
                n_cases=n_cases,
                n_controls=n_controls,
                local_results_dir=local_results_dir,
            )
            results.append(pd.read_csv(protein_file, sep="\t"))
            protein_file.unlink()
        except subprocess.CalledProcessError as error:
            # this is a sensitivity check, not a filter (see coloc_targets.py's
            # pairwise_coloc() for the actual gating COLOC step) - susieR's LD/sumstats
            # consistency check is far stricter than coloc.abf's and genuinely does fail
            # on individual loci (e.g. reference-panel LD not matching the sumstats'
            # sample size), so 1 protein failing must not take down the whole run
            print(f"[CONCERN] COLOC SuSiE failed for {protein}, skipping - {error}")

    if not results:
        print("[CONCERN] No COLOC SuSiE results were generated.")
        return

    master = pd.concat(results, ignore_index=True)
    out_file = paths.coloc_susie_out(pqtl_dataset, pheno_id, local_results_dir)
    master.to_csv(out_file, sep="\t", index=False)
    print(f"[DONE] Saved master COLOC SuSiE table: {out_file}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pheno_id", required=True)
    p.add_argument("--pqtl_dataset", required=True)
    p.add_argument("--pqtl_dir", required=True)
    p.add_argument("--maf", type=float, required=True)
    p.add_argument("--ref_bfile", required=True)
    p.add_argument("--n_cases", type=int, required=True)
    p.add_argument("--n_controls", type=int, required=True)
    p.add_argument("--local_results_dir", default="results")
    args = p.parse_args()

    grab_ld_matrix_for_cis_mr_hits_and_coloc_susie(
        pqtl_dataset=args.pqtl_dataset,
        pheno_id=args.pheno_id,
        pqtl_dir=args.pqtl_dir,
        maf=args.maf,
        ref_bfile=args.ref_bfile,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        local_results_dir=args.local_results_dir,
    )


if __name__ == "__main__":
    main()