#!/usr/bin/env bash
set -euo pipefail

# file (sumstats parquet, must have a SNP column)
# id
# ref_prefix (plink bfile prefix, e.g. dat/ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL)
# maf
# outdir

file="$1"
id="$2"
ref_prefix="$3"
maf="$4"
outdir="$5"

mkdir -p "$outdir"

snplist_file="${outdir}/${id}.snplist.txt"
ld_out="${outdir}/${id}_ld_matrix"

echo "[TRACKING] Extracting SNP list for ${id}"
python3 -c "
import pyarrow.parquet as pq
snps = pq.read_table('${file}', columns=['SNP']).column('SNP').to_pylist()
open('${snplist_file}', 'w').write('\n'.join(snps) + '\n')
"

if [ ! -s "$snplist_file" ]; then
    echo "[CONCERN] No SNPs found in ${file}"
    exit 1
fi

echo "[TRACKING] Running LD matrix for ${id}"
plink \
    --bfile "$ref_prefix" \
    --extract "$snplist_file" \
    --maf "$maf" \
    --r square \
    --write-snplist \
    --out "$ld_out"

echo "[DONE] Saved LD matrix: ${ld_out}.ld (SNP order: ${ld_out}.snplist)"