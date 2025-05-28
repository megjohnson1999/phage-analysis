#!/bin/bash
# Script to test splitting proteins by phage ID

# Set up paths
CURRENT_DIR="$(pwd)"
WORKFLOW_DIR="$(cd "$(dirname "$0")/../workflow" && pwd)"
SCRIPTS_DIR="${WORKFLOW_DIR}/scripts"
TEST_DIR="${CURRENT_DIR}/test_phage_splitting"

# Create test directory
mkdir -p "${TEST_DIR}"
cd "${TEST_DIR}"

# Check if an input file was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_protein_file.faa>"
    echo "e.g., $0 /path/to/proteins.part_095.faa"
    exit 1
fi

INPUT_FILE="$1"
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# Activate conda environment (activate the environment that has biopython)
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate phacts || conda activate base

# Create a symbolic link to the input file
ln -sf "${INPUT_FILE}" "${TEST_DIR}/proteins.faa"

# Run the splitting script
echo "Running splitting script..."
python "${SCRIPTS_DIR}/split_proteins_by_phage.py" \
    --input "${TEST_DIR}/proteins.faa" \
    --output-dir "${TEST_DIR}/split_proteins" \
    --split-list "${TEST_DIR}/split_list.txt" \
    --show-examples

# Show results
echo ""
echo "Files created in ${TEST_DIR}/split_proteins/:"
ls -la "${TEST_DIR}/split_proteins/"

echo ""
echo "Split list contents:"
cat "${TEST_DIR}/split_list.txt"

echo ""
echo "Test complete. You can examine the split files in ${TEST_DIR}/split_proteins/"