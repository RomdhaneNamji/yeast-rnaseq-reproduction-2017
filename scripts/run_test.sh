#!/bin/bash

run=SRR5482575
mkdir -p raw_fastq

# Get the ENA table row
line=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" | tail -n 1)

echo "ENA line:"
echo "$line"
echo

# Extract only the fastq_ftp column (2nd column)
fastqs=$(echo "$line" | cut -f2)

echo "FASTQ paths:"
echo "$fastqs"
echo

# Split into read1 and read2
url1=$(echo "$fastqs" | cut -d';' -f1)
url2=$(echo "$fastqs" | cut -d';' -f2)

echo "URL1 = $url1"
echo "URL2 = $url2"
echo

# Download
curl -L "https://${url1}" -o "raw_fastq/${run}_1.fastq.gz"
curl -L "https://${url2}" -o "raw_fastq/${run}_2.fastq.gz"


