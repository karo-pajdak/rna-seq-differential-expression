#!/bin/bash
# =======================================================================================
# salmon Script
# =======================================================================================
# This script 
#
# Usage:
# ./salmon.sh [NUM_THREADS]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# =======================================================================================

set -euo pipefail

# Define input and output directories
SALMON_INPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/trimmed"
SALMON_OUTPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/salmon"
SALMON_INDEX="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/transcriptome_index"
SALMON_GENOME="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/reference_transcripts"
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"
LOGS_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/logs"

# Make sure output directory exists
mkdir -p "$SALMON_OUTPUT"
mkdir -p "$LOGS_DIR"
mkdir -p "$SALMON_INDEX"

# User specifies how many cores used; default = 1 (Sequential)
THREADS=${1:-1}

# Check samples.txt exists and is not empty
if [ ! -s "${CONFIG_DIR}/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# Check if index files exist. If not, build the index (first check if reference transcriptome exists)
if [ -z "$(ls -A "$SALMON_INDEX")" ]; then
    echo "Salmon index missing. Building index..."
    if [ -z "$(ls -A "$SALMON_GENOME")" ]; then
        echo "ERROR: Reference transcriptome folder is empty."
        exit 1
    fi
    # Take first file in the genome folder
    TRANSCRIPTOME_FILE=$(ls -1 "$SALMON_GENOME" | head -n 1)
    TRANSCRIPTOME_FILE="$SALMON_GENOME/$TRANSCRIPTOME_FILE"
    echo "Building Salmon index from $TRANSCRIPTOME_FILE..."
    mkdir -p "$SALMON_INDEX"
    salmon index -t "$TRANSCRIPTOME_FILE" -i "$SALMON_INDEX/index_file"
    echo "Index built successfully."
fi