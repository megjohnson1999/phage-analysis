#!/bin/bash

# Test case execution script for phage-analysis pipeline
set -e

# Create output directory if it doesn't exist
mkdir -p test_data/results

# Run the snakemake workflow with the test configuration
echo "Running phage-analysis pipeline with test data..."

# Check if we want to run with or without Reneo
if [ "$1" == "--skip-reneo" ]; then
    echo "Skipping Reneo step, using direct assembly filtering..."
    CONFIG="test_data/test_config.yaml"
    # Create a temporary config file with Reneo disabled
    TMP_CONFIG=$(mktemp)
    cat $CONFIG | sed 's/use_reneo: true/use_reneo: false/' > $TMP_CONFIG
    
    snakemake --use-conda \
        --configfile $TMP_CONFIG \
        --cores 4 \
        --printshellcmds \
        --rerun-incomplete \
        --keep-going \
        -p test_data/results/01_filtered_assembly/filtered_assembly_1KB.fasta
        
    rm $TMP_CONFIG
else
    echo "Using Reneo for assembly binning..."
    snakemake --use-conda \
        --configfile test_data/test_config.yaml \
        --cores 4 \
        --printshellcmds \
        --rerun-incomplete \
        --keep-going \
        -p test_data/results/01_reneo_output/genomes_and_unresolved_edges.fasta
fi

echo "Test complete!"