#!/bin/bash
set -euo pipefail

PROJECT="/mnt/d/GIS/yeast_metabolomics_2017"
ALIGNDIR="$PROJECT/alignments"
COUNTDIR="$PROJECT/counts/strandedness_test"
ANNOT="$PROJECT/ref/downloads/sacCer3.gff3"

mkdir -p "$COUNTDIR"

for s in 0 1 2; do
    echo "Testing -s $s ..."
    featureCounts \
      -T 10 \
      -F GFF \
      -t gene \
      -g gene_id \
      -p \
      -s "$s" \
      -a "$ANNOT" \
      -o "$COUNTDIR/gene_counts_s${s}.txt" \
      "$ALIGNDIR"/*.sorted.bam
done
