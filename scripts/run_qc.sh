#!/bin/bash
set -euo pipefail

# =========================
# Paths
# =========================

# Project on external drive
PROJECT_D="/mnt/d/GIS/yeast_metabolomics_2017"

# Local temporary work folder inside WSL
BASE="$HOME/yeast_qc_tmp"

# Input/output inside local temp folder
RAW_LOCAL="$BASE/raw_fastq"
QC_LOCAL="$BASE/qc"
FASTQC_OUT="$QC_LOCAL/fastqc_raw"
LOGDIR="$BASE/logs"
LOG="$LOGDIR/qc.log"

# Final destination on /mnt/d
QC_FINAL="$PROJECT_D/qc"

# =========================
# Parallel settings
# 10 total threads = 5 jobs x 2 threads
# =========================
PARALLEL_JOBS=5
FASTQC_THREADS=2

mkdir -p "$RAW_LOCAL" "$FASTQC_OUT" "$LOGDIR" "$QC_FINAL"

# =========================
# Timer
# =========================
TOTAL_START=$(date +%s)

echo "=== QC START ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "PROJECT_D=$PROJECT_D" | tee -a "$LOG"
echo "RAW_LOCAL=$RAW_LOCAL" | tee -a "$LOG"
echo "FASTQC_OUT=$FASTQC_OUT" | tee -a "$LOG"
echo "QC_FINAL=$QC_FINAL" | tee -a "$LOG"
echo "PARALLEL_JOBS=$PARALLEL_JOBS" | tee -a "$LOG"
echo "FASTQC_THREADS=$FASTQC_THREADS" | tee -a "$LOG"
echo | tee -a "$LOG"

# =========================
# Check source FASTQ files
# =========================
if ! compgen -G "$PROJECT_D/raw_fastq/*.fastq.gz" > /dev/null; then
    echo "ERROR: No FASTQ files found in $PROJECT_D/raw_fastq" | tee -a "$LOG"
    exit 1
fi

# =========================
# Copy FASTQ files locally
# =========================
echo "Copying FASTQ files to local workspace..." | tee -a "$LOG"
cp "$PROJECT_D"/raw_fastq/*.fastq.gz "$RAW_LOCAL"/

# =========================
# Run FastQC in parallel
# =========================
echo "Running FastQC in parallel..." | tee -a "$LOG"

export FASTQC_OUT FASTQC_THREADS LOG

find "$RAW_LOCAL" -maxdepth 1 -type f -name "*.fastq.gz" | sort | \
xargs -n 1 -P "$PARALLEL_JOBS" -I {} bash -c '
    fq="$1"
    name=$(basename "$fq")
    echo "[START] $name" >> "$LOG"
    fastqc "$fq" -o "$FASTQC_OUT" --threads "$FASTQC_THREADS" >> "$LOG" 2>&1
    echo "[DONE]  $name" >> "$LOG"
' _ {}

FASTQC_END=$(date +%s)
FASTQC_MIN=$(( (FASTQC_END - TOTAL_START) / 60 ))

echo "FastQC finished in ${FASTQC_MIN} minute(s)" | tee -a "$LOG"

# =========================
# Run MultiQC
# =========================
MULTIQC_START=$(date +%s)

echo "Running MultiQC..." | tee -a "$LOG"
multiqc "$FASTQC_OUT" -o "$QC_LOCAL" >> "$LOG" 2>&1

MULTIQC_END=$(date +%s)
MULTIQC_MIN=$(( (MULTIQC_END - MULTIQC_START) / 60 ))

echo "MultiQC finished in ${MULTIQC_MIN} minute(s)" | tee -a "$LOG"

# =========================
# Move results back to /mnt/d
# =========================
echo "Moving results back to $QC_FINAL ..." | tee -a "$LOG"

rm -rf "$QC_FINAL/fastqc_raw"
mkdir -p "$QC_FINAL"

mv "$FASTQC_OUT" "$QC_FINAL/"
mv "$QC_LOCAL"/multiqc_report.html "$QC_FINAL/"
mv "$QC_LOCAL"/multiqc_data "$QC_FINAL/"
cp "$LOG" "$QC_FINAL/qc.log"

# =========================
# Clean local temp files
# =========================
echo "Cleaning local temporary files..." | tee -a "$LOG"
rm -rf "$BASE"

TOTAL_END=$(date +%s)
TOTAL_MIN=$(( (TOTAL_END - TOTAL_START) / 60 ))

echo "=== QC DONE ==="
echo "Total QC time: ${TOTAL_MIN} minute(s)"
echo "Final results saved in:"
echo "$QC_FINAL"
echo "$QC_FINAL/multiqc_report.html"
