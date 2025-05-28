#!/bin/bash
# Script to test phage-specific PHACTS analysis within the Snakemake workflow
# Usage: bash scripts/test_phage_specific_phacts.sh <proteins_faa_file>

set -e

# Check if an input file was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <protein_file.faa>"
    echo "e.g., $0 no_reneo_full/03_split_proteins/proteins.part_001.faa"
    exit 1
fi

INPUT_FILE="$1"
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# Determine the working directory of the workflow 
WORKFLOW_DIR="$(cd "$(dirname "$0")/../workflow" && pwd)"
echo "=== Phage-specific PHACTS analysis test ==="
echo "Using workflow directory: ${WORKFLOW_DIR}"
echo "Using input file: ${INPUT_FILE}"

# Create a temporary config file using settings from the test config
WORKFLOW_TEST_CONFIG="$WORKFLOW_DIR/../test_data/test_config.yaml"
TMP_CONFIG=$(mktemp)

# Copy the test config and just update the output directory
cat "$WORKFLOW_TEST_CONFIG" > "$TMP_CONFIG"
echo "output_dir: \"$(dirname "$INPUT_FILE")/phage_specific_test\"" >> "$TMP_CONFIG"

# Create the output directory
mkdir -p "$(dirname "$INPUT_FILE")/phage_specific_test/03_split_proteins_by_phage"
mkdir -p "$(dirname "$INPUT_FILE")/phage_specific_test/03_orf_predictions"
mkdir -p "$(dirname "$INPUT_FILE")/phage_specific_test/logs"

# Copy the input file to the required location
cp "$INPUT_FILE" "$(dirname "$INPUT_FILE")/phage_specific_test/03_orf_predictions/proteins.faa"

echo ""
echo "=== Running phage-specific protein splitting ==="
# Run the specific rules for phage-specific PHACTS analysis
cd "$WORKFLOW_DIR" && snakemake \
    --snakefile "$WORKFLOW_DIR/Snakefile" \
    --configfile $TMP_CONFIG \
    --use-conda \
    --conda-frontend mamba \
    --cores 4 \
    --printshellcmds \
    --keep-going \
    "$(dirname "$INPUT_FILE")/phage_specific_test/03_split_proteins_by_phage/split_protein_list.txt"

# Check the results
echo ""
echo "=== Checking split results ==="
if [ -f "$(dirname "$INPUT_FILE")/phage_specific_test/03_split_proteins_by_phage/split_protein_list.txt" ]; then
    echo "Split protein list created successfully."
    
    SPLIT_COUNT=$(ls -1 "$(dirname "$INPUT_FILE")/phage_specific_test/03_split_proteins_by_phage/"*.faa 2>/dev/null | wc -l)
    echo "Number of split phage protein files: $SPLIT_COUNT"
    
    # List a few example files
    echo "Example split files:"
    ls -1 "$(dirname "$INPUT_FILE")/phage_specific_test/03_split_proteins_by_phage/"*.faa 2>/dev/null | head -n 5
    
    # Run the rest of the PHACTS analysis on the split files
    echo ""
    echo "=== Running phage-specific PHACTS prediction ==="
    cd "$WORKFLOW_DIR" && snakemake \
        --snakefile "$WORKFLOW_DIR/Snakefile" \
        --configfile $TMP_CONFIG \
        --use-conda \
        --conda-frontend mamba \
        --cores 4 \
        --printshellcmds \
        --keep-going \
        "$(dirname "$INPUT_FILE")/phage_specific_test/03_phacts_results_by_phage/phacts_predictions_compiled.tsv"
    
    # Show the final results
    echo ""
    echo "=== Final PHACTS results ==="
    if [ -f "$(dirname "$INPUT_FILE")/phage_specific_test/03_phacts_results_by_phage/phacts_predictions_compiled.tsv" ]; then
        echo "PHACTS predictions file created successfully."
        echo "Contents of predictions file:"
        cat "$(dirname "$INPUT_FILE")/phage_specific_test/03_phacts_results_by_phage/phacts_predictions_compiled.tsv"
    else
        echo "Error: PHACTS predictions file was not created."
    fi
else
    echo "Error: Split protein list was not created."
fi

# Clean up temporary files but keep the results directory
rm $TMP_CONFIG
echo ""
echo "=== Test complete ==="
echo "You can find the results in: $(dirname "$INPUT_FILE")/phage_specific_test/"
echo ""
echo "If you want to run the test with a different input file, clean up first with:"
echo "rm -rf $(dirname "$INPUT_FILE")/phage_specific_test"
echo ""
echo "View PHACTS predictions with:"
echo "cat $(dirname "$INPUT_FILE")/phage_specific_test/03_phacts_results_by_phage/phacts_predictions_compiled.tsv"
echo ""