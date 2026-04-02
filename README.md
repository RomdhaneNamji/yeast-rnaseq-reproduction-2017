# Yeast RNA-seq reproduction

This repository contains a reproduction workflow for a yeast colony biofilm RNA-seq study.

## Project goals
- reproduce the main transcriptional differences between aerial and root cells
- perform quality control, alignment, counting, differential expression, and GO analysis
- compare reproduced results with the published study

## Workflow
1. Download RNA-seq data
2. Run FastQC and MultiQC
3. Align reads with HISAT2
4. Count reads with featureCounts
5. Perform differential expression with DESeq2
6. Run GO enrichment analysis

## Main tools
- HISAT2
- samtools
- featureCounts
- DESeq2
- R

## Main findings
- aerial cells showed stress, glucose starvation, and sporulation-related signatures
- root cells showed translation- and transport-related signatures
- the reproduced biological pattern was consistent with the published study

## Notes
Raw FASTQ and BAM files are not included in this repository because of file size.
