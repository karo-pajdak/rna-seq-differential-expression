#!/bin/bash
# ===================================================================================
# FastQC + MultiQC Pipeline Script
# ===================================================================================
# This script performs quality control on RNA-seq FASTQ files. It supports sequential
# or parallel execution based on the number of threads. MultiQC is run at the end to
# aggregate all FastQC reports.
#
# Usage:
# ./quality_control.sh [NUM_THREADS] [INPUT_DIR] [OUTPUT_DIR]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# ===================================================================================

set -euo pipefail

# Define input and output directories
INPUT_DIR=${2:-"$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/raw"}
OUTPUT_DIR=${3:-"$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/qc/raw"}
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"
LOGS_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/logs"

# Make sure output/log directory exists
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGS_DIR"

# Check samples.txt exists and is not empty
if [ ! -s "${CONFIG_DIR}/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

# Timestamp for logs
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

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
    export INPUT_DIR OUTPUT_DIR LOGS_DIR TOTAL_CORES SAMPLES_AT_ONCE TIMESTAMP

    # Run FastQC in parallel
    parallel -j "$SAMPLES_AT_ONCE" -a "$CONFIG_DIR/samples.txt" '
        SAMPLE={}

        # Inline suffix detection to differentiate between trimmed and raw
        if [ -f "${INPUT_DIR}/${SAMPLE}_1_paired.fastq.gz" ]; then
            SUFFIX1="_1_paired.fastq.gz"
            SUFFIX2="_2_paired.fastq.gz"
        else
            SUFFIX1="_1.fastq.gz"
            SUFFIX2="_2.fastq.gz"
        fi

        R1="${INPUT_DIR}/${SAMPLE}${SUFFIX1}"
        R2="${INPUT_DIR}/${SAMPLE}${SUFFIX2}"
        LOG_FILE="${LOGS_DIR}/${SAMPLE}_fastqc_${TIMESTAMP}.log"

        if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
            echo "ERROR: Missing FASTQ for $SAMPLE"
            echo "  Looking for: $R1"
            echo "  Looking for: $R2"
            exit 1
        fi

        if ! gzip -t "$R1" || ! gzip -t "$R2"; then
            echo "ERROR: Corrupted FASTQ for $SAMPLE"
            exit 1
        fi

        # Guarantees at least 1 thread per job
        FASTQC_THREADS=$((TOTAL_CORES / SAMPLES_AT_ONCE))
        [ $FASTQC_THREADS -lt 1 ] && FASTQC_THREADS=1

        echo "Processing $SAMPLE with $FASTQC_THREADS threads..."
        fastqc -o "$OUTPUT_DIR" -t "$FASTQC_THREADS" "$R1" "$R2" > "$LOG_FILE" 2>&1
        echo "Finished $SAMPLE. Log saved to $LOG_FILE"
    '

else
    # Run sequentially
    echo "Running FastQC sequentially..."
    # Loop through samples
    while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]; do
        echo "Processing $SAMPLE..."

        # Inline suffix detection to differentiate between trimmed and raw
        if [ -f "${INPUT_DIR}/${SAMPLE}_1_paired.fastq.gz" ]; then
            SUFFIX1="_1_paired.fastq.gz"
            SUFFIX2="_2_paired.fastq.gz"
        else
            SUFFIX1="_1.fastq.gz"
            SUFFIX2="_2.fastq.gz"
        fi

        # Find each sample in input directory 
        R1="${INPUT_DIR}/${SAMPLE}${SUFFIX1}"
        R2="${INPUT_DIR}/${SAMPLE}${SUFFIX2}"
        LOG_FILE="${LOGS_DIR}/${SAMPLE}_fastqc_${TIMESTAMP}.log"

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
        fastqc -o "$OUTPUT_DIR" -t 1 "$R1" "$R2" > "$LOG_FILE" 2>&1
        echo "Finished $SAMPLE. Log saved to $LOG_FILE"

    done < "$CONFIG_DIR/samples.txt"
fi

# Run MultiQC to aggregate QC results
echo "Running MultiQC to summarize FastQC reports..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" > "$LOGS_DIR/multiqc_$(basename "$OUTPUT_DIR")_${TIMESTAMP}.log" 2>&1