#!/bin/bash

# File with one SRR accession per line
INPUT="./metadata/srr_runs.txt"

# Output folder for FASTQ files
OUTDIR="./raw_fastq"

# Log folder
LOGDIR="./logs"

mkdir -p "$OUTDIR" "$LOGDIR"

TOTAL_START=$(date +%s)

echo "Starting ENA download process..."
echo "Input file: $INPUT"
echo "Output folder: $OUTDIR"
echo

while read -r run; do
    # Skip empty lines
    if [[ -z "$run" ]]; then
        continue
    fi

    echo "Processing $run"

    # Skip if both files already exist
    if [[ -f "$OUTDIR/${run}_1.fastq.gz" && -f "$OUTDIR/${run}_2.fastq.gz" ]]; then
        echo "$run already downloaded, skipping"
        echo
        continue
    fi

    SAMPLE_START=$(date +%s)

    # Get ENA file paths
    line=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp&format=tsv")
    entry=$(echo "$line" | tail -n 1)
    fastqs=$(echo "$entry" | cut -f2)

    url1=$(echo "$fastqs" | cut -d';' -f1)
    url2=$(echo "$fastqs" | cut -d';' -f2)

    # Check that both URLs were found
    if [[ -z "$url1" || -z "$url2" ]]; then
        echo "ERROR: could not get FASTQ links for $run" | tee -a "$LOGDIR/ena_errors.log"
        echo
        continue
    fi

    echo "  Downloading read 1 ..."
    curl -L "https://${url1}" -o "$OUTDIR/${run}_1.fastq.gz" >> "$LOGDIR/download.log" 2>&1

    echo "  Downloading read 2 ..."
    curl -L "https://${url2}" -o "$OUTDIR/${run}_2.fastq.gz" >> "$LOGDIR/download.log" 2>&1

    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))

    echo "$run DONE in ${SAMPLE_TIME} seconds"
    echo

done < "$INPUT"

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

echo "All downloads finished."
echo "Total time: ${TOTAL_TIME} seconds"
