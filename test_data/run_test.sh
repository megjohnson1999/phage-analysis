#!/bin/bash

# Test case execution script for phage-analysis pipeline
set -e

# Create output directory if it doesn't exist
mkdir -p test_data/results

# Run the snakemake workflow with the test configuration
echo "Running phage-analysis pipeline with test data..."
snakemake --use-conda \
    --configfile test_data/test_config.yaml \
    --cores 4 \
    --printshellcmds \
    --rerun-incomplete \
    --keep-going \
    -n  # Dry run first to check execution plan

echo "To execute the actual workflow, remove the -n flag"
echo "Test complete!"