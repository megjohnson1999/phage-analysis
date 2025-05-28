#!/bin/bash

# Local test script for phage analysis pipeline
# This script runs the workflow with local conda environments

cd "$(dirname "$0")/.."

echo "Starting phage analysis pipeline test with local conda environments..."
echo "Working directory: $(pwd)"

# Check if snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "Error: snakemake not found. Please install snakemake first:"
    echo "conda install -c conda-forge snakemake"
    exit 1
fi

# Create output directory
mkdir -p test_data/results_local

# Run the workflow with local config
echo "Running workflow with local conda environments..."
snakemake --configfile test_data/test_config_local.yaml \
    --use-conda \
    --conda-create-envs-only \
    --cores 2 \
    --directory workflow/ \
    --snakefile workflow/Snakefile

echo "Conda environments created. Now running the actual workflow..."

# Run just the protein splitting part to test the issue
snakemake --configfile test_data/test_config_local.yaml \
    --use-conda \
    --cores 2 \
    --directory workflow/ \
    --snakefile workflow/Snakefile \
    test_data/results_local/03_split_proteins_by_phage/split_protein_list.txt

echo "Test completed!"