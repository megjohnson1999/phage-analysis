#!/bin/bash

# This script tests PHACTS installation and import issues
# Run this script directly from your project root with:
# bash scripts/test_phacts_install.sh OUTPUT_DIR
# Where OUTPUT_DIR is the same as your workflow output directory

set -e  # Exit on error

if [ $# -lt 1 ]; then
    echo "Usage: $0 OUTPUT_DIR"
    echo "Example: $0 /path/to/output"
    exit 1
fi

OUTPUT_DIR="$1"
PHACTS_DIR="${OUTPUT_DIR}/db/phacts/PHACTS"
PHACTS_PATH="${PHACTS_DIR}/phacts.py"
LOG_FILE="${OUTPUT_DIR}/phacts_test_log.txt"

echo "==== PHACTS Installation Test ====" | tee "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "PHACTS directory: $PHACTS_DIR" | tee -a "$LOG_FILE"
echo "PHACTS script: $PHACTS_PATH" | tee -a "$LOG_FILE"

# Check if directories exist
echo -n "Output directory exists: " | tee -a "$LOG_FILE"
if [ -d "$OUTPUT_DIR" ]; then echo "Yes" | tee -a "$LOG_FILE"; else echo "No" | tee -a "$LOG_FILE"; fi

echo -n "PHACTS directory exists: " | tee -a "$LOG_FILE"
if [ -d "$PHACTS_DIR" ]; then echo "Yes" | tee -a "$LOG_FILE"; else echo "No" | tee -a "$LOG_FILE"; fi

echo -n "PHACTS script exists: " | tee -a "$LOG_FILE"
if [ -f "$PHACTS_PATH" ]; then echo "Yes" | tee -a "$LOG_FILE"; else echo "No" | tee -a "$LOG_FILE"; fi

# If PHACTS directory doesn't exist, create it
if [ ! -d "$PHACTS_DIR" ]; then
    echo "PHACTS not found. Installing PHACTS..." | tee -a "$LOG_FILE"
    mkdir -p "${OUTPUT_DIR}/db/phacts"
    cd "${OUTPUT_DIR}/db/phacts"
    git clone https://github.com/deprekate/PHACTS.git | tee -a "$LOG_FILE"
    cd PHACTS
    # Create __init__.py file to make it a proper Python package
    touch __init__.py
    cd ..
    touch .installed
fi

# Create __init__.py if it doesn't exist
if [ ! -f "$PHACTS_DIR/__init__.py" ]; then
    echo "Creating missing __init__.py file" | tee -a "$LOG_FILE"
    touch "$PHACTS_DIR/__init__.py"
fi

# Create a test input file
TEST_DIR="${OUTPUT_DIR}/phacts_test"
mkdir -p "$TEST_DIR"
TEST_FASTA="${TEST_DIR}/test.fasta"

# Create a simple test fasta file
cat > "$TEST_FASTA" << EOF
>test_sequence
MTHSVDSVVGYDISGIYLSGHKAFLDSLASVISDGGTTKSEKVLQRALTQYTGGLGYMVVDL
KVEPTKLAAMIDISLLALENCKAACKKAGINYERPGMLLVDIGTGSGSLSAGEIRKAGYDHV
TVIDIEKDRSWALLDELKQWSGYQVVRVSDEEDAAKLFHEAGIDFLGGIDLTARTGEMLSLA
AEFADVKVYFVQTEQEAQALFKDGATFHSVNLPARLSRFTSLTLPSVPVGTLIHRLSQDSKL
EOF

echo -e "\n==== Python Environment Test ====" | tee -a "$LOG_FILE"
echo "Python version:" | tee -a "$LOG_FILE"
python --version 2>&1 | tee -a "$LOG_FILE"

echo -e "\n==== PHACTS Import Test ====" | tee -a "$LOG_FILE"
# Run Python script to test import
python - << EOF | tee -a "$LOG_FILE"
import sys
import os

print(f"Current working directory: {os.getcwd()}")
print(f"Python executable: {sys.executable}")

# Add PHACTS to Python path
phacts_dir = "${PHACTS_DIR}"
phacts_parent = os.path.dirname("${PHACTS_DIR}")
print(f"Adding to Python path: {phacts_dir}")
print(f"Adding to Python path: {phacts_parent}")
sys.path.insert(0, phacts_dir)
sys.path.insert(0, phacts_parent)

print("Python path:")
for path in sys.path:
    print(f"  - {path}")

print("\nTrying to import phacts...")
try:
    import phacts
    print("IMPORT SUCCESSFUL!")
    print(f"phacts location: {phacts.__file__}")
except ImportError as e:
    print(f"IMPORT FAILED: {e}")
    
print("\nTrying to import from directory...")
try:
    import PHACTS.phacts
    print("IMPORT SUCCESSFUL with PHACTS.phacts!")
except ImportError as e:
    print(f"IMPORT FAILED with PHACTS.phacts: {e}")
EOF

echo -e "\n==== PHACTS Execution Test ====" | tee -a "$LOG_FILE"

# Set PYTHONPATH explicitly
export PYTHONPATH="${PHACTS_DIR}:$(dirname ${PHACTS_DIR}):$PYTHONPATH"
echo "PYTHONPATH: $PYTHONPATH" | tee -a "$LOG_FILE"

# Try to run PHACTS directly
echo -e "\nAttempting to run PHACTS directly..." | tee -a "$LOG_FILE"
COMMAND="python \"$PHACTS_PATH\" \"$TEST_FASTA\" -o \"$TEST_DIR\""
echo "Command: $COMMAND" | tee -a "$LOG_FILE"
eval $COMMAND >> "$LOG_FILE" 2>&1 || echo "Direct execution failed" | tee -a "$LOG_FILE"

# Try to run as a module
echo -e "\nAttempting to run PHACTS as a module..." | tee -a "$LOG_FILE"
cd "$(dirname ${PHACTS_DIR})"
COMMAND="python -m PHACTS.phacts \"$TEST_FASTA\" -o \"$TEST_DIR\""
echo "Command: $COMMAND" | tee -a "$LOG_FILE"
eval $COMMAND >> "$LOG_FILE" 2>&1 || echo "Module execution failed" | tee -a "$LOG_FILE"

# Check if any output was created
echo -e "\nChecking for output files..." | tee -a "$LOG_FILE"
ls -la "$TEST_DIR" | tee -a "$LOG_FILE"

echo -e "\n==== Test Complete ====" | tee -a "$LOG_FILE"
echo "See full logs at: $LOG_FILE" | tee -a "$LOG_FILE"