#!/bin/bash

# Wrapper script for Reneo that handles expected failures
# Reneo is known to fail at certain steps but still produces the required output files

set -o pipefail

# Parse command line arguments
INPUT_GRAPH=""
READS_DIR=""
OUTPUT_DIR=""
THREADS=""
MINLENGTH=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_GRAPH="$2"
            shift 2
            ;;
        --reads)
            READS_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --minlength)
            MINLENGTH="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$INPUT_GRAPH" ] || [ -z "$READS_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Missing required parameters"
    echo "Usage: $0 --input <graph> --reads <dir> --output <dir> [--threads <n>] [--minlength <n>]"
    exit 1
fi

# Set defaults if not provided
THREADS=${THREADS:-24}
MINLENGTH=${MINLENGTH:-1000}

echo "Running Reneo with the following parameters:"
echo "  Input graph: $INPUT_GRAPH"
echo "  Reads directory: $READS_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo "  Minimum length: $MINLENGTH"

# Run Reneo, capturing the exit code
echo "Starting Reneo run..."
reneo run --input "$INPUT_GRAPH" \
    --reads "$READS_DIR" \
    --minlength "$MINLENGTH" \
    --output "$OUTPUT_DIR" \
    --threads "$THREADS" 2>&1 | tee reneo_output.log

RENEO_EXIT_CODE=${PIPESTATUS[0]}

echo "Reneo exited with code: $RENEO_EXIT_CODE"

# Check if the expected output file exists
EXPECTED_OUTPUT="$OUTPUT_DIR/genomes_and_unresolved_edges.fasta"

if [ -f "$EXPECTED_OUTPUT" ]; then
    echo "SUCCESS: Expected output file found at $EXPECTED_OUTPUT"
    
    # Check if file is not empty
    if [ -s "$EXPECTED_OUTPUT" ]; then
        echo "Output file contains data ($(wc -l < "$EXPECTED_OUTPUT") lines)"
        echo "Reneo completed successfully (output files generated despite exit code)"
        exit 0
    else
        echo "ERROR: Output file exists but is empty"
        exit 1
    fi
else
    echo "ERROR: Expected output file not found at $EXPECTED_OUTPUT"
    echo "Reneo failed to produce the required output"
    
    # List what files were created
    echo "Files in output directory:"
    ls -la "$OUTPUT_DIR/" || echo "Could not list output directory"
    
    exit 1
fi