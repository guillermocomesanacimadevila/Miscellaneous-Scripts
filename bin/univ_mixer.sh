#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $0 <USER> <SUMSTATS.gz> <PHENO> [--mixer-root PATH] [--mixer-ref PATH] [--outdir PATH] [--rep INT]

Defaults:
  --mixer-root /scratch/<USER>/mixer_dropbox
  --mixer-ref <mixer-root>/reference/ldsc/1000G_EUR_Phase3_plink
  --outdir /scratch/<USER>/mixer_runs/results/univariate/<PHENO>
  --rep 1

Example:
  $0 c...... /scratch/c....../mixer_ready/SCZ_mixer_ready.tsv.gz SCZ
EOF
}

if [[ $# -lt 3 ]]; then usage; exit 1; fi

USER="$1"
SUMSTATS="$2"
PHENO="$3"
shift 3

REP=1
MIXER_ROOT="/scratch/${USER}/mixer_dropbox"
MIXER_REF=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mixer-root) MIXER_ROOT="$2"; shift 2 ;;
    --mixer-ref) MIXER_REF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --rep) REP="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${MIXER_REF}" ]]; then
  MIXER_REF="${MIXER_ROOT}/reference/ldsc/1000G_EUR_Phase3_plink"
fi
if [[ -z "${OUTDIR}" ]]; then
  OUTDIR="/scratch/${USER}/mixer_runs/results/univariate/${PHENO}"
fi

if command -v module >/dev/null 2>&1; then
  module purge
  module load singularity
fi

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER="singularity"
else
  echo "ERROR: apptainer/singularity not found on this machine" >&2
  exit 127
fi

MIXER_SIF="${MIXER_ROOT}/singularity/v2.2.1/mixer.sif"

if [[ "${OUTDIR}" == /scratch/* ]]; then
  BIND_SCRATCH="--bind /scratch:/scratch"
else
  BIND_SCRATCH=""
fi

MIXER_PY="${CONTAINER} exec ${BIND_SCRATCH} --bind ${PWD}:${PWD} ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

for f in "${SUMSTATS}" "${MIXER_SIF}" "${MIXER_REF}/1000G.EUR.QC.@.bim" "${MIXER_REF}/1000G.EUR.QC.@.run4.ld" "${MIXER_REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty file: ${f}" >&2; exit 2; }
done

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "==> [${PHENO}] split_sumstats"
${MIXER_PY} split_sumstats \
  --trait1-file "${SUMSTATS}" \
  --out "${PHENO}.chr@.sumstats.gz"

echo "==> [${PHENO}] fit1 (rep${REP})"
${MIXER_PY} fit1 \
  --trait1-file "${SUMSTATS}" \
  --bim-file "${MIXER_REF}/1000G.EUR.QC.@.bim" \
  --ld-file "${MIXER_REF}/1000G.EUR.QC.@.run4.ld" \
  --extract "${MIXER_REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps" \
  --seed "${REP}" \
  --out "${PHENO}.rep${REP}"

echo "==> [${PHENO}] test1 (rep${REP})"
${MIXER_PY} test1 \
  --trait1-file "${SUMSTATS}" \
  --bim-file "${MIXER_REF}/1000G.EUR.QC.@.bim" \
  --ld-file "${MIXER_REF}/1000G.EUR.QC.@.run4.ld" \
  --extract "${MIXER_REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps" \
  --load-params-file "${PHENO}.rep${REP}.fit.json" \
  --out "${PHENO}.rep${REP}"

echo "==> [${PHENO}] outputs present?"
ls -lh "${PHENO}.rep${REP}.fit.json" "${PHENO}.rep${REP}.json" "${PHENO}.rep${REP}.log" 2>/dev/null || true

echo "==> [${PHENO}] sanity: analysis fields (no jq)"
python3 - <<PY
import json
for fn in ["${PHENO}.rep${REP}.fit.json","${PHENO}.rep${REP}.json"]:
    try:
        with open(fn) as f:
            d=json.load(f)
        print(fn, "analysis =", d.get("analysis"))
    except FileNotFoundError:
        print(fn, "MISSING")
PY
