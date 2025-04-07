#!/bin/bash

# Test case execution script for phage-analysis pipeline
set -e

# Create mock database files (empty files just for testing)
echo "Setting up mock databases..."
for DB in mmseqs2 checkv genomad iphop phabox2 phacts vcontact3; do
    touch test_data/db/$DB/db_mock.txt
done

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