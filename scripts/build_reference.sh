#!/bin/bash
set -euo pipefail


#===================================
# Paths
#====================================


PROJECT="/mnt/d/GIS/yeast_metabolomics_2017"
REFDIR="$PROJECT/ref"
DOWNLOADS="$REFDIR/downloads"
INDEXDIR="$REFDIR/index"
LOGDIR="$REFDIR/logs"
LOG="$LOGDIR/build_reference.log"

FASTA="$DOWNLOADS/sacCer3.fa"
ANNOT="$DOWNLOADS/sacCer3.gff3"

mkdir -p "$INDEXDIR" "$LOGDIR"

START=$(date +%s)

echo "=== BUILD REFERENCE START ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "FASTA=$FASTA" | tee -a "$LOG"
echo "ANNOT=$ANNOT" | tee -a "$LOG"
echo | tee -a "$LOG"

# Check files exist
if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found: $FASTA" | tee -a "$LOG"
    exit 1
fi

if [[ ! -f "$ANNOT" ]]; then
    echo "ERROR: Annotation file not found: $ANNOT" | tee -a "$LOG"
    exit 1
fi

# =================================
# Build HISAT2 index
# =================================


echo "Building HISAT2 index ..." | tee -a "$LOG"
hisat2-build "$FASTA" "$INDEXDIR/sacCer3" >> "$LOG" 2>&1


END=$(date +%s)
MINUTES=$(( (END - START) / 60 ))

echo | tee -a "$LOG"
echo "=== BUILD REFERENCE DONE ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "Time: ${MINUTES} minute(s)" | tee -a "$LOG"
echo "Index prefix: $INDEXDIR/sacCer3" | tee -a "$LOG"
