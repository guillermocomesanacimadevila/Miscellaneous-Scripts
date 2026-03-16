#!/usr/bin/env bash
set -euo pipefail

pheno1_id="${1:-}"
pheno2_id="${2:-}"
loci_dir="${3:-}"
ref_prefix="${4:-}"

[ -n "${pheno1_id}" ] && [ -n "${pheno2_id}" ] && [ -n "${loci_dir}" ] && [ -n "${ref_prefix}" ] || {
  echo "usage: bin/gen_ld_matrix.sh PHENO1 PHENO2 LOCI_DIR REF_PREFIX"
  exit 1
}

task_root="$(pwd)"
ref_prefix="${task_root}/${ref_prefix}"
find "${loci_dir}" -type d -name 'locus_chr*_*_*' | while read -r locus_dir; do
  locus_name="$(basename "${locus_dir}")"
  chr_var="$(echo "${locus_name}" | cut -d'_' -f2 | sed 's/^chr//')"
  cd "${locus_dir}"
  awk -F'\t' 'NR>1{print $1}' "gwas_${pheno1_id}.ldgwas.tsv" | sort -u > "${pheno1_id}.snps"
  awk -F'\t' 'NR>1{print $1}' "gwas_${pheno2_id}.ldgwas.tsv" | sort -u > "${pheno2_id}.snps"
  comm -12 "${pheno1_id}.snps" "${pheno2_id}.snps" > gwas.snps.intersection
  mkdir -p ld_gwas

  plink \
    --bfile "${ref_prefix}.${chr_var}" \
    --extract gwas.snps.intersection \
    --maf 0.01 \
    --geno 0.02 \
    --make-bed \
    --out ld_gwas/locus_tmp

  cut -f2 ld_gwas/locus_tmp.bim > ld_gwas/snps_in_ld_order.txt

  plink \
    --bfile ld_gwas/locus_tmp \
    --r square gz \
    --out ld_gwas/locus

  python - <<PY
import pandas as pd
from pathlib import Path

order = pd.read_csv("ld_gwas/snps_in_ld_order.txt", header=None)[0].astype(str).tolist()

def reorder(f):
    f = Path(f)
    df = pd.read_csv(f, sep="\\t", dtype={"SNP": str})
    df = df.set_index("SNP")
    keep = [x for x in order if x in df.index]
    out = f.with_name(f.name.replace(".ldgwas.tsv", ".ldorder.tsv"))
    df.loc[keep].reset_index().to_csv(out, sep="\\t", index=False)

reorder("gwas_${pheno1_id}.ldgwas.tsv")
reorder("gwas_${pheno2_id}.ldgwas.tsv")
PY

  cd "${task_root}"
done