#!/usr/bin/env python3
import argparse
import os
import time
from pathlib import Path
import subprocess

# ref/xQTLs/eQTLs/eQTL_eQTLGen
# sQTL_BrainMeta
# eQTL_eQTLGen

def check_format(xqtl: str, ref_dir: str, qtl_dataset: str):
    target_dir = Path(ref_dir) / f"{xqtl}s/{qtl_dataset}"
    files = [f for f in os.listdir(target_dir) if not f.startswith(".")]
    smr_ext = (".besd", ".epi", ".esi", ".summary")
    if all(file.endswith(smr_ext) for file in files):
        print("SMR format! Needs adjustment!")
        time.sleep(2)
        print("No worries bro, I got this!")
        return True
    print("Not pure SMR format, moving on.")
    return False

def format_smr_qtls_for_coloc(xqtl: str,
        chroms: list, # 11, 15, 8
        ref_dir: str,
        qtl_dataset: str,
        out_dir: str):
    # step 1 => define chroms of interest
    ref_dir = Path(ref_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    qtl_dir = Path(ref_dir) / f"{xqtl}s/{qtl_dataset}"

    # define outcome = need for transformation?
    outcome = check_format(xqtl, ref_dir, qtl_dataset)

    if outcome == True:
        for chrom in chroms:
            besd_files = list(qtl_dir.glob(f"*chr{chrom}.besd"))
            besd_prefix = besd_files[0].with_suffix("")
            chr_out = out_dir / f"{qtl_dataset}_chr{chrom}"
            cmd = [
                "smr",
                "--beqtl-summary", str(besd_prefix),
                "--query", "1",
                "--chr", str(chrom),
                "--out", str(chr_out)
            ]
            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)
        print("Done formatting SMR QTL dataset for requested chromosomes only.")
        return
    print("Dataset already looks like text / sumstats format, nothing to do.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xqtl", help="xQTL type (eQTL / sQTL / pQTL)")
    parser.add_argument("chroms", nargs="+", type=int, help="chromosomes")
    parser.add_argument("ref_dir", help="reference directory")
    parser.add_argument("qtl_dataset", help="qtl dataset")
    parser.add_argument("out_dir", help="output directory")
    args = parser.parse_args()
    format_smr_qtls_for_coloc(
        xqtl=args.xqtl,
        chroms=args.chroms,
        ref_dir=args.ref_dir,
        qtl_dataset=args.qtl_dataset,
        out_dir=args.out_dir
    )

if __name__ == "__main__":
    main()