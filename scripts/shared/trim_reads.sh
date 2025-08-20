#!/bin/bash
# =======================================================================================
# Trimmomatic Pipeline Script
# =======================================================================================
# This script cleans raw RNA-seq FASTQ files by removing adapters and low-quality bases
# using Trimmomatic. It checks that input files exist and are not corrupted, then saves 
#trimmed paired and unpaired reads to an output folder. Logs are created for each sample, 
#and the script can run sequentially or in parallel on multiple cores.
#
# Usage:
# ./trim_reads.sh [NUM_THREADS]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# =======================================================================================

set -euo pipefail

# Directories
TRIM_INPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/raw"
TRIM_OUTPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/trimmed"
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"

# Adapter file
ADAPTERS="${CONFIG_DIR}/TruSeq3-PE.fa"

# Make sure output directory exists
mkdir -p "$TRIM_OUTPUT"

# Check samples.txt exists and is not empty
if [ ! -s "${CONFIG_DIR}/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

if [ "$THREADS" -gt 1 ]; then
    # Run in parallel
    echo "Running Trimmomatic in parallel using $THREADS threads..."

    # Decide how many samples to run simultaneously
    TOTAL_CORES=$THREADS
    SAMPLES_AT_ONCE=1
    if [ "$THREADS" -gt 3 ]; then
        SAMPLES_AT_ONCE=2
    fi

    # Export variables so parallel jobs can access them
    export TRIM_INPUT TRIM_OUTPUT CONFIG_DIR ADAPTERS

    # Run Trimmomatic in parallel
    parallel -j "$SAMPLES_AT_ONCE" -a "$CONFIG_DIR/samples.txt" '
        SAMPLE={}
        # Find each sample in input directory 
        RAW_R1="${TRIM_INPUT}/${SAMPLE}_1.fastq.gz"
        RAW_R2="${TRIM_INPUT}/${SAMPLE}_2.fastq.gz"

        # Output filenames
        PAIRED_R1="${TRIM_OUTPUT}/${SAMPLE}_1_paired.fastq.gz"
        PAIRED_R2="${TRIM_OUTPUT}/${SAMPLE}_2_paired.fastq.gz"
        UNPAIRED_R1="${TRIM_OUTPUT}/${SAMPLE}_1_unpaired.fastq.gz"
        UNPAIRED_R2="${TRIM_OUTPUT}/${SAMPLE}_2_unpaired.fastq.gz"
        LOG_FILE="${TRIM_OUTPUT}/${SAMPLE}_trimmomatic.log"

        # Check if input files exist
        if [ ! -f "$RAW_R1" ] || [ ! -f "$RAW_R2" ]; then
            echo "ERROR: Missing FASTQ for $SAMPLE"
            exit 1
        fi

        # Check if input files are corrupted
        if ! gzip -t "$RAW_R1" || ! gzip -t "$RAW_R2"; then
            echo "ERROR: Corrupted FASTQ for $SAMPLE"
            exit 1
        fi

        # Run Trimmomatic
        trimmomatic PE "$RAW_R1" "$RAW_R2" \
            "$PAIRED_R1" "$UNPAIRED_R1" \
            "$PAIRED_R2" "$UNPAIRED_R2" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:36 \
            > "$LOG_FILE" 2>&1

        echo "Finished $SAMPLE. Log saved to $LOG_FILE"
    '

else
    # Run Trimmmomatic sequentially
    echo "Running Trimmomatic sequentially..."
# Loop through samples
    while read SAMPLE; do
        RAW_R1="${TRIM_INPUT}/${SAMPLE}_1.fastq.gz"
        RAW_R2="${TRIM_INPUT}/${SAMPLE}_2.fastq.gz"

        # Output filenames
        PAIRED_R1="${TRIM_OUTPUT}/${SAMPLE}_1_paired.fastq.gz"
        PAIRED_R2="${TRIM_OUTPUT}/${SAMPLE}_2_paired.fastq.gz"
        UNPAIRED_R1="${TRIM_OUTPUT}/${SAMPLE}_1_unpaired.fastq.gz"
        UNPAIRED_R2="${TRIM_OUTPUT}/${SAMPLE}_2_unpaired.fastq.gz"
        LOG_FILE="${TRIM_OUTPUT}/${SAMPLE}_trimmomatic.log"

        # Check if input files exist
        if [ ! -f "$RAW_R1" ] || [ ! -f "$RAW_R2" ]; then
            echo "ERROR: Missing FASTQ for $SAMPLE"
            exit 1
        fi

        # Check if input files are corrupted
        if ! gzip -t "$RAW_R1" || ! gzip -t "$RAW_R2"; then
            echo "ERROR: Corrupted FASTQ for $SAMPLE"
            exit 1
        fi

        # Run Trimmomatic (paired-end)
        trimmomatic PE "$RAW_R1" "$RAW_R2" \
            "$PAIRED_R1" "$UNPAIRED_R1" \
            "$PAIRED_R2" "$UNPAIRED_R2" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:36 \
            > "$LOG_FILE" 2>&1

        echo "Finished $SAMPLE. Log saved to $LOG_FILE"
    done < "$CONFIG_DIR/samples.txt"
fi

echo "All samples processed successfully."