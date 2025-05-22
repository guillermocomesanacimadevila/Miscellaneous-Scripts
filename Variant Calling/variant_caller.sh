#!/bin/bash

#============================================#
# Variant Calling Pipeline with Coverage
# Author: Guillermo Comesaña Cimadevila
# Date: $(date +%F)
#============================================#

set -euo pipefail

echo "============================================"
echo "      Variant Calling Pipeline Launcher      "
echo "============================================"

# --- Conda Setup --- #
ENV_NAME="varcall_env"

if ! command -v conda &> /dev/null; then
    echo "[ERROR] Conda not found. Please install Miniconda or Anaconda."
    exit 1
fi

if ! conda env list | grep -q "$ENV_NAME"; then
    echo "[INFO] Creating Conda environment '$ENV_NAME'..."
    conda create -y -n "$ENV_NAME" bwa samtools bcftools
else
    echo "[INFO] Using existing Conda environment '$ENV_NAME'"
fi

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# --- User Input --- #
read -rp "Enter path to reference genome (FASTA): " REF
read -rp "Enter path to read 1 (FASTQ): " READ1
read -rp "Enter path to read 2 (FASTQ, leave empty for single-end): " READ2
read -rp "Enter output prefix (e.g. sample1): " OUTPREFIX

# --- Input Validation --- #
for file in "$REF" "$READ1"; do
    if [ ! -f "$file" ]; then
        echo "[ERROR] File not found: $file"
        exit 1
    fi
done

if [ -n "$READ2" ] && [ ! -f "$READ2" ]; then
    echo "[ERROR] File not found: $READ2"
    exit 1
fi

# --- Logging --- #
LOGFILE="varcall_${OUTPREFIX}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOGFILE") 2>&1

# --- Index Reference (if needed) --- #
if [ ! -f "${REF}.bwt" ]; then
    echo "[INFO] Indexing reference genome..."
    bwa index "$REF"
else
    echo "[INFO] Reference already indexed."
fi

# --- Alignment --- #
echo "[INFO] Aligning reads with BWA..."
if [ -z "$READ2" ]; then
    bwa mem "$REF" "$READ1" > "${OUTPREFIX}.sam"
else
    bwa mem "$REF" "$READ1" "$READ2" > "${OUTPREFIX}.sam"
fi

# --- SAM to Sorted BAM --- #
echo "[INFO] Converting SAM to sorted BAM..."
samtools view -bS "${OUTPREFIX}.sam" | samtools sort -o "${OUTPREFIX}.sorted.bam"
samtools index "${OUTPREFIX}.sorted.bam"

# --- Coverage and Depth --- #
echo "[INFO] Calculating per-base depth..."
samtools depth -a "${OUTPREFIX}.sorted.bam" > "${OUTPREFIX}.depth.txt"

echo "[INFO] Calculating average depth and breadth of coverage..."
TOTAL_BASES=$(awk '{sum += 1} END {print sum}' "${OUTPREFIX}.depth.txt")
TOTAL_COVERED=$(awk '$3 > 0 {sum += 1} END {print sum}' "${OUTPREFIX}.depth.txt")
TOTAL_DEPTH=$(awk '{sum += $3} END {print sum}' "${OUTPREFIX}.depth.txt")

if [ "$TOTAL_BASES" -gt 0 ]; then
    AVG_DEPTH=$(echo "scale=2; $TOTAL_DEPTH / $TOTAL_BASES" | bc)
    COVERAGE_PERCENT=$(echo "scale=2; ($TOTAL_COVERED / $TOTAL_BASES) * 100" | bc)
else
    AVG_DEPTH=0
    COVERAGE_PERCENT=0
fi

echo "[RESULT] Average depth: ${AVG_DEPTH}x"
echo "[RESULT] Breadth of coverage: ${COVERAGE_PERCENT}%"

# Save summary
echo -e "Average_Depth\tBreadth_of_Coverage(%)" > "${OUTPREFIX}.coverage_summary.txt"
echo -e "${AVG_DEPTH}\t${COVERAGE_PERCENT}" >> "${OUTPREFIX}.coverage_summary.txt"

# --- Variant Calling --- #
echo "[INFO] Calling variants with bcftools..."
bcftools mpileup -Ou -f "$REF" "${OUTPREFIX}.sorted.bam" | \
    bcftools call -mv -Ob -o "${OUTPREFIX}.bcf"
bcftools view "${OUTPREFIX}.bcf" > "${OUTPREFIX}.vcf"

# --- Done --- #
echo "============================================"
echo "      Pipeline completed successfully!       "
echo "============================================"
echo "  Output files:"
echo "    - Sorted BAM: ${OUTPREFIX}.sorted.bam"
echo "    - Depth file: ${OUTPREFIX}.depth.txt"
echo "    - Coverage summary: ${OUTPREFIX}.coverage_summary.txt"
echo "    - VCF: ${OUTPREFIX}.vcf"
echo "    - Log: ${LOGFILE}"
