#!/bin/bash
# ===================================================================================
# FastQC + MultiQC Pipeline Script
# ===================================================================================
# This script performs quality control on RNA-seq FASTQ files. It supports sequential
# or parallel execution based on the number of threads. MultiQC is run at the end to
# aggregate all FastQC reports.
#
# Usage:
# ./quality_control.sh [NUM_THREADS]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# ===================================================================================

set -euo pipefail

# Define input and output directories
INPUT_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/raw"
OUTPUT_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/qc"
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Check samples.txt exists and is not empty
if [ ! -s "${CONFIG_DIR}/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

if [ "$THREADS" -gt 1 ]; then
    # Run in parallel
    echo "Running FastQC in parallel using $THREADS threads..."

    # Decide how many samples to run simultaneously
    TOTAL_CORES=$THREADS
    SAMPLES_AT_ONCE=1
    if [ "$THREADS" -gt 3 ]; then
        SAMPLES_AT_ONCE=2
    fi

    # Export variables so parallel jobs can access them
    export INPUT_DIR OUTPUT_DIR

    # Run FastQC in parallel
    parallel -j "$SAMPLES_AT_ONCE" -a "$CONFIG_DIR/samples.txt" '
        SAMPLE={}
        # Find each sample in input directory 
        R1="${INPUT_DIR}/${SAMPLE}_1.fastq.gz"
        R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"

        # Check if files exist
        if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
            echo "ERROR: Missing FASTQ for $SAMPLE"
            exit 1
        fi

        # Check if files are corrupted
        if ! gzip -t "$R1" || ! gzip -t "$R2"; then
            echo "ERROR: Corrupted FASTQ for $SAMPLE"
            exit 1
        fi

        # Guarantees at least 1 thread per job
        FASTQC_THREADS=$((TOTAL_CORES / SAMPLES_AT_ONCE))
        [ $FASTQC_THREADS -lt 1 ] && FASTQC_THREADS=1

        # Run FastQC on both files - parallel processing
        fastqc -o "$OUTPUT_DIR" -t "$FASTQC_THREADS" "$R1" "$R2"
    '

else
    # Run sequentially
    echo "Running FastQC sequentially..."
    # Loop through samples
    while read SAMPLE; do
        echo "Processing $SAMPLE..."

        # Find each sample in input directory 
        R1="${INPUT_DIR}/${SAMPLE}_1.fastq.gz"
        R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"

        # Check if files exist
        if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
            echo "ERROR: Missing FASTQ for $SAMPLE"
            exit 1
        fi

        # Check if files are corrupted
        if ! gzip -t "$R1" || ! gzip -t "$R2"; then
            echo "ERROR: Corrupted FASTQ for $SAMPLE"
            exit 1
        fi

        # Run FastQC on both files - single core
        fastqc -o "$OUTPUT_DIR" -t 1 "$R1" "$R2"

    done < "$CONFIG_DIR/samples.txt"
fi

# Run MultiQC to aggregate QC results
echo "Running MultiQC to summarize FastQC reports..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR"