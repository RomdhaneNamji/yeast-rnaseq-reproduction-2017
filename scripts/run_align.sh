#!/bin/bash
set -euo pipefail

# =========================
# Paths
# =========================

# Project on external drive
PROJECT_D="/mnt/d/GIS/yeast_metabolomics_2017"

# Local temporary workspace in WSL
BASE="$HOME/yeast_align_tmp"
RAW_LOCAL="$BASE/raw_fastq"
ALIGN_LOCAL="$BASE/alignments"
LOGDIR="$BASE/logs"

# Final destination on /mnt/d
ALIGN_FINAL="$PROJECT_D/alignments"

# Reference on /mnt/d
INDEX_PREFIX="$PROJECT_D/ref/index/sacCer3"

# Metadata file
METADATA="$PROJECT_D/metadata/sample_metadata.csv"

mkdir -p "$RAW_LOCAL" "$ALIGN_LOCAL" "$LOGDIR" "$ALIGN_FINAL"

# =========================
# Parallel settings
# 10 total threads = 2 jobs x 5 threads
# =========================
PARALLEL_JOBS=2
HISAT2_THREADS=5
SORT_THREADS=2

# =========================
# Timer
# =========================
TOTAL_START=$(date +%s)

echo "=== ALIGNMENT START ==="
date
echo "PROJECT_D=$PROJECT_D"
echo "RAW_LOCAL=$RAW_LOCAL"
echo "ALIGN_LOCAL=$ALIGN_LOCAL"
echo "ALIGN_FINAL=$ALIGN_FINAL"
echo "INDEX_PREFIX=$INDEX_PREFIX"
echo "METADATA=$METADATA"
echo

# =========================
# Checks
# =========================
if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    echo "ERROR: HISAT2 index not found: ${INDEX_PREFIX}.1.ht2"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Metadata file not found: $METADATA"
    exit 1
fi

# =========================
# Copy FASTQ locally
# =========================
echo "Copying FASTQ files to local workspace..."
cp "$PROJECT_D"/raw_fastq/*.fastq.gz "$RAW_LOCAL"/

# =========================
# Run alignment in parallel
# sample_metadata.csv columns:
# sample,condition,read1,read2
# =========================
echo "Running HISAT2 alignments..."

tail -n +2 "$METADATA" | \
xargs -I {} -P "$PARALLEL_JOBS" bash -c '
    line="$1"

    sample=$(echo "$line" | cut -d, -f1)
    read1=$(echo "$line" | cut -d, -f3 | xargs basename)
    read2=$(echo "$line" | cut -d, -f4 | xargs basename)

    r1="'"$RAW_LOCAL"'/$read1"
    r2="'"$RAW_LOCAL"'/$read2"

    echo "[START] $sample"

    hisat2 \
      -x "'"$INDEX_PREFIX"'" \
      -1 "$r1" \
      -2 "$r2" \
      -p "'"$HISAT2_THREADS"'" \
      --summary-file "'"$ALIGN_LOCAL"'/${sample}.hisat2.summary.txt" \
      2> "'"$ALIGN_LOCAL"'/${sample}.hisat2.log" \
    | samtools view -@ "'"$SORT_THREADS"'" -bS - \
    | samtools sort -@ "'"$SORT_THREADS"'" -o "'"$ALIGN_LOCAL"'/${sample}.sorted.bam" -

    samtools index -@ "'"$SORT_THREADS"'" "'"$ALIGN_LOCAL"'/${sample}.sorted.bam"
    samtools flagstat -@ "'"$SORT_THREADS"'" "'"$ALIGN_LOCAL"'/${sample}.sorted.bam" > "'"$ALIGN_LOCAL"'/${sample}.flagstat.txt"

    echo "[DONE]  $sample"
' _ {}

# =========================
# Move results back to /mnt/d
# =========================
echo
echo "Moving BAM files and logs back to $ALIGN_FINAL ..."
mkdir -p "$ALIGN_FINAL"

mv "$ALIGN_LOCAL"/*.sorted.bam "$ALIGN_FINAL"/
mv "$ALIGN_LOCAL"/*.sorted.bam.bai "$ALIGN_FINAL"/
mv "$ALIGN_LOCAL"/*.hisat2.log "$ALIGN_FINAL"/
mv "$ALIGN_LOCAL"/*.hisat2.summary.txt "$ALIGN_FINAL"/
mv "$ALIGN_LOCAL"/*.flagstat.txt "$ALIGN_FINAL"/

# =========================
# Clean local temp files
# =========================
echo "Cleaning local temporary files..."
rm -rf "$BASE"

TOTAL_END=$(date +%s)
TOTAL_MIN=$(( (TOTAL_END - TOTAL_START) / 60 ))

echo
echo "=== ALIGNMENT DONE ==="
date
echo "Total alignment time: ${TOTAL_MIN} minute(s)"
echo "Final results saved in:"
echo "$ALIGN_FINAL"
