#!/bin/bash
# =======================================================================================
# featureCounts Script
# =======================================================================================
# This script runs featureCounts on all sorted BAM files in the alignment directory to 
# generate a gene-level count matrix for downstream differential expression analysis.
# It automatically detects whether the data is paired-end or single-end, selects the 
# appropriate mode, and uses a GTF or GFF3 annotation file for read summarization.
#
# Usage:
# ./counts.sh [NUM_THREADS]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# =======================================================================================
set -euo pipefail

# Define input and output directories
COUNTS_INPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/alignment/bam"
COUNTS_OUTPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/counts"
LOGS_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/logs"
REF_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/annotation"

# Make sure output/log directory exists
mkdir -p "$COUNTS_OUTPUT"
mkdir -p "$LOGS_DIR"

# Check annotation file exists and is not empty (GTF or GFF3)
if [ -s "${REF_DIR}/annotation.gtf" ]; then
    ANNOTATION="${REF_DIR}/annotation.gtf"
elif [ -s "${REF_DIR}/annotation.gff3" ]; then
    ANNOTATION="${REF_DIR}/annotation.gff3"
else
    echo "ERROR: annotation.gtf or annotation.gff3 is missing or empty"
    exit 1
fi

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

# Pick the first BAM to test pairedness
TEST_BAM=$(ls "$COUNTS_INPUT"/*.sorted.bam | head -n 1)

# Count alignments with paired flag
PAIRED_COUNT=$(samtools view -c -f 1 "$TEST_BAM")

# Run featureCounts
if [ "$PAIRED_COUNT" -gt 0 ]; then
    featureCounts -p -T "$THREADS" -a "$ANNOTATION" -t exon -g gene_id -M -O -o "$COUNTS_OUTPUT/alignment_counts.txt" "$COUNTS_INPUT"/*.sorted.bam \
        > "$LOGS_DIR/featureCounts.log" 2>&1
else
    featureCounts -T "$THREADS" -a "$ANNOTATION" -t exon -g gene_id  -M -O -o "$COUNTS_OUTPUT/alignment_counts.txt" "$COUNTS_INPUT"/*.sorted.bam \
        > "$LOGS_DIR/featureCounts.log" 2>&1
fi
