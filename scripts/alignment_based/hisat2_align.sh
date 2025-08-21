#!/bin/bash
# =======================================================================================
# HISAT2 Alignment Script
# =======================================================================================
# This script m# Maps NGS reads against a single reference genome using HISAT2.
# Index will be built from reference genome (reference/reference_genome) if pre-built index 
# is not in reference/genome_index folder. Each sample uses multiple threads if specified.
#
# Usage:
# ./hisat2_align.sh [NUM_THREADS]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# =======================================================================================
set -euo pipefail

# Directories
HISAT_INPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/trimmed"
HISAT_OUTPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/alignment/bam/"
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"
LOGS_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/logs"
HISAT_INDEX="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/genome_index"
HISAT_GENOME="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/reference_genome"

# Make sure output directory exists
mkdir -p "$HISAT_OUTPUT"
mkdir -p "$LOGS_DIR"
mkdir -p "$HISAT_INDEX"

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

# Check samples.txt exists and is not empty
if [ ! -s "${CONFIG_DIR}/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# Check if index files exist. If not, build the index (first check if reference genome exists)
if [ -z "$(ls -A "$HISAT_INDEX")" ]; then
    echo "HISAT2 index missing. Building index..."
    if [ -z "$(ls -A "$HISAT_GENOME")" ]; then
        echo "ERROR: Reference genome folder is empty."
        exit 1
    fi
    # Take first file in the genome folder
    GENOME_FILE=$(ls -1 "$HISAT_GENOME" | head -n 1)
    GENOME_FILE="$HISAT_GENOME/$GENOME_FILE"
    echo "Building HISAT2 index from $GENOME_FILE..."
    mkdir -p "$HISAT_INDEX"
    hisat2-build "$GENOME_FILE" "$HISAT_INDEX/index_file"
    echo "Index built successfully."
fi

# Loop through samples
while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]; do
    R1="${HISAT_INPUT}/${SAMPLE}_1_paired.fastq.gz"
    R2="${HISAT_INPUT}/${SAMPLE}_2_paired.fastq.gz"
    SAMPLE_LOG="$LOGS_DIR/${SAMPLE}_hisat2.log"

    # Check if input files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "ERROR: Missing FASTQ for $SAMPLE"
        exit 1
    fi

    # Check if input files are corrupted
    if ! gzip -t "$R1" || ! gzip -t "$R2"; then
        echo "ERROR: Corrupted FASTQ for $SAMPLE"
        exit 1
    fi

    # Run HISAT2, convert to BAM, log HISAT2 messages
    echo "Running HISAT2 for $SAMPLE..." 
    hisat2 -p "$THREADS" \
        -x "$HISAT_INDEX/index_file" \
        -1 "$R1" -2 "$R2" \
        -S /dev/stdout 2> "$SAMPLE_LOG" | \
    samtools sort -@ "$THREADS" -o "$HISAT_OUTPUT/${SAMPLE}.sorted.bam"

done < "$CONFIG_DIR/samples.txt"