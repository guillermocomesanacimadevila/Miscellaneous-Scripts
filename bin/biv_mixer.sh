#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $0 <USER> <SUM1.gz> <SUM2.gz> <PHENO1> <PHENO2> [--rep INT] [--univdir PATH] [--outdir PATH] [--mixer-root PATH] [--mixer-ref PATH]

Example:
  $0 c.c24102394 /scratch/c.c24102394/mixer_ready/AD_mixer_ready.tsv.gz /scratch/c.c24102394/mixer_ready/SCZ_mixer_ready.tsv.gz AD SCZ --rep 1
EOF
}

if [[ $# -lt 5 ]]; then usage; exit 1; fi

USER="$1"; SUM1="$2"; SUM2="$3"; PHENO1="$4"; PHENO2="$5"
shift 5

REP=1
MIXER_ROOT="/scratch/${USER}/mixer_dropbox"
MIXER_REF=""
UNIVDIR="/scratch/${USER}/mixer_runs/results/univariate"
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rep) REP="$2"; shift 2 ;;
    --univdir) UNIVDIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --mixer-root) MIXER_ROOT="$2"; shift 2 ;;
    --mixer-ref) MIXER_REF="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${MIXER_REF}" ]]; then
  MIXER_REF="${MIXER_ROOT}/reference/ldsc/1000G_EUR_Phase3_plink"
fi

PAIR="${PHENO1}_${PHENO2}"

if [[ -z "${OUTDIR}" ]]; then
  OUTDIR="/scratch/${USER}/mixer_runs/results/bivariate/${PAIR}"
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

if [[ "${OUTDIR}" == /scratch/* ]] || [[ "${UNIVDIR}" == /scratch/* ]] || [[ "${SUM1}" == /scratch/* ]] || [[ "${SUM2}" == /scratch/* ]]; then
  BIND_SCRATCH="--bind /scratch:/scratch"
else
  BIND_SCRATCH=""
fi

MIXER_PY="${CONTAINER} exec ${BIND_SCRATCH} --bind ${PWD}:${PWD} ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
MIXER_FIG="${CONTAINER} exec ${BIND_SCRATCH} --bind ${PWD}:${PWD} ${MIXER_SIF} python /tools/mixer/precimed/mixer_figures.py"

U1="${UNIVDIR}/${PHENO1}/${PHENO1}.rep${REP}.json"
U2="${UNIVDIR}/${PHENO2}/${PHENO2}.rep${REP}.json"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

for f in \
  "${SUM1}" "${SUM2}" "${U1}" "${U2}" \
  "${MIXER_SIF}" \
  "${MIXER_REF}/1000G.EUR.QC.@.bim" \
  "${MIXER_REF}/1000G.EUR.QC.@.run4.ld" \
  "${MIXER_REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps"
do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty: ${f}" >&2; exit 2; }
done

echo "==> sanity: univariate jsons must be analysis==univariate"
python3 - <<PY
import json, sys
for fn in ["${U1}","${U2}"]:
    with open(fn) as f:
        d=json.load(f)
    if d.get("analysis") != "univariate":
        print("ERROR:", fn, "analysis =", d.get("analysis"), file=sys.stderr)
        sys.exit(3)
    print("OK:", fn, "analysis = univariate")
PY

echo "==> [${PAIR}] fit2 (rep${REP})"
${MIXER_PY} fit2 \
  --trait1-file "${SUM1}" \
  --trait2-file "${SUM2}" \
  --trait1-params-file "${U1}" \
  --trait2-params-file "${U2}" \
  --out "${PAIR}.fit.rep${REP}" \
  --extract "${MIXER_REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps" \
  --bim-file "${MIXER_REF}/1000G.EUR.QC.@.bim" \
  --ld-file "${MIXER_REF}/1000G.EUR.QC.@.run4.ld"

echo "==> [${PAIR}] test2 (rep${REP})"
${MIXER_PY} test2 \
  --trait1-file "${SUM1}" \
  --trait2-file "${SUM2}" \
  --load-params-file "${PAIR}.fit.rep${REP}.json" \
  --out "${PAIR}.apply.rep${REP}" \
  --bim-file "${MIXER_REF}/1000G.EUR.QC.@.bim" \
  --ld-file "${MIXER_REF}/1000G.EUR.QC.@.run4.ld"

echo "==> [${PAIR}] figures (point_estimate only)"
${MIXER_FIG} two \
  --json-fit "${PAIR}.fit.rep${REP}.json" \
  --json-test "${PAIR}.apply.rep${REP}.json" \
  --out "${PAIR}" \
  --statistic point_estimate || true

echo "==> [${PAIR}] outputs present?"
ls -lh "${PAIR}.fit.rep${REP}.json" "${PAIR}.apply.rep${REP}.json" \
      "${PAIR}.fit.rep${REP}.log" "${PAIR}.apply.rep${REP}.log" \
      ${PAIR}*.pdf 2>/dev/null || true
