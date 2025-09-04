#!/bin/bash
# =======================================================================================
# Salmon Quantification Script
# =======================================================================================
# This script quantifies RNA-seq reads at the transcript level using Salmon.
# Builds a transcriptome-only or decoy-aware index if missing, then quantifies all samples
# listed in config/samples.txt. Each sample can use multiple threads if specified.
#
# Usage:
# ./salmon.sh [NUM_THREADS] [USE_DECOY]
# - NUM_THREADS: optional, default = 1 (sequential). Set >1 for parallel mode.
# - USE_DECOY: optional, default = false. Set true to build decoy-aware index.
# =======================================================================================
set -euo pipefail

# Directories
SALMON_INPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/trimmed"
SALMON_OUTPUT="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/results/salmon"
SALMON_INDEX="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/transcriptome_index"
SALMON_TRANSCRIPTOME_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/reference_transcripts"
SALMON_GENOME_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/reference_genome"
SALMON_DECOY="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/decoy"
CONFIG_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/config"
LOGS_DIR="$HOME/BioinformaticsPortfolio/rna-seq-differential-expression/logs"

mkdir -p "$SALMON_OUTPUT" "$LOGS_DIR" "$SALMON_INDEX" "$SALMON_DECOY"

# User inputs
THREADS=${1:-1}
USE_DECOY=${2:-false}   # Set true to build decoy-aware index

# Check samples.txt
if [ ! -s "$CONFIG_DIR/samples.txt" ]; then
    echo "ERROR: samples.txt is missing or empty"
    exit 1
fi

# Pick transcriptome FASTA
TRANSCRIPTOME_FILE=$(ls -1 "$SALMON_TRANSCRIPTOME_DIR"/*.fa* | head -n 1)
if [ -z "$TRANSCRIPTOME_FILE" ]; then
    echo "ERROR: No transcriptome FASTA found in $SALMON_TRANSCRIPTOME_DIR"
    exit 1
fi

# Pick genome FASTA if decoy mode
GENOME_FILE=$(ls -1 "$SALMON_GENOME_DIR"/*.fa* | head -n 1)

# Build index if missing
if [ -z "$(ls -A "$SALMON_INDEX" 2>/dev/null)" ]; then
    if [ "$USE_DECOY" == "true" ]; then
        if [ -z "$GENOME_FILE" ]; then
            echo "ERROR: Decoy mode selected but genome FASTA not found in $SALMON_GENOME_DIR"
            exit 1
        fi
        echo "Building decoy-aware Salmon index..."
        grep "^>" "$GENOME_FILE" | cut -d " " -f1 | sed 's/>//' > "$SALMON_DECOY/decoys.txt"
        GENTROME_FILE="$SALMON_GENOME_DIR/gentrome.fa"
        cat "$TRANSCRIPTOME_FILE" "$GENOME_FILE" > "$GENTROME_FILE"
        salmon index -t "$GENTROME_FILE" -d "$SALMON_DECOY/decoys.txt" -i "$SALMON_INDEX" --gencode -p "$THREADS"
        echo "Decoy-aware index built successfully."
    else
        echo "Building transcriptome-only Salmon index..."
        salmon index -t "$TRANSCRIPTOME_FILE" -i "$SALMON_INDEX" --gencode -p "$THREADS"
        echo "Transcriptome-only index built successfully."
        echo "Note: A decoy-aware index could be built for higher accuracy (need 100+ GB memory)."
    fi
else
    echo "Salmon index already exists at $SALMON_INDEX"
fi

# Quantification loop
while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]; do
    R1="${SALMON_INPUT}/${SAMPLE}_1_paired.fastq.gz"
    R2="${SALMON_INPUT}/${SAMPLE}_2_paired.fastq.gz"
    SAMPLE_LOG="$LOGS_DIR/${SAMPLE}_salmon.log"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "ERROR: Missing FASTQ for $SAMPLE"
        exit 1
    fi

    if ! gzip -t "$R1" || ! gzip -t "$R2"; then
        echo "ERROR: Corrupted FASTQ for $SAMPLE"
        exit 1
    fi

    echo "Running Salmon for $SAMPLE..."
    salmon quant -i "$SALMON_INDEX" -l A \
        -1 "$R1" -2 "$R2" \
        --validateMappings --seqBias --gcBias --posBias \
        -p "$THREADS" -o "$SALMON_OUTPUT/$SAMPLE" \
        > "$SAMPLE_LOG" 2>&1

done < "$CONFIG_DIR/samples.txt"