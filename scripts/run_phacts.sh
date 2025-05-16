#!/bin/bash
#
# PHACTS wrapper script for the phage-analysis pipeline
# This script handles the proper environment setup and path resolution for PHACTS
#

# Usage information
function show_usage {
    echo "Usage: $0 INPUT_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "PHACTS Wrapper Script for phage lifestyle prediction"
    echo ""
    echo "Required arguments:"
    echo "  INPUT_FILE               Input protein FASTA file"
    echo "  -o, --output OUTPUT_DIR  Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -h, --help               Show this help message"
    echo "  -d, --debug              Enable debug mode (extra logging)"
    echo "  -p, --phacts PATH        Explicitly set PHACTS path"
    echo ""
}

# Defaults
DEBUG=false
INPUT_FILE=""
OUTPUT_DIR=""
CUSTOM_PHACTS_PATH=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--phacts)
            CUSTOM_PHACTS_PATH="$2"
            shift 2
            ;;
        -d|--debug)
            DEBUG=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            # First non-flag argument is the input file
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
            else
                echo "Error: Unexpected argument: $1"
                show_usage
                exit 1
            fi
            shift
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: Input file is required"
    show_usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory (-o) is required"
    show_usage
    exit 1
fi

# Create log file in output directory
mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/phacts_run.log"
touch "$LOG_FILE"

# Helper function for logging
function log {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Debug logging
if [[ "$DEBUG" == true ]]; then
    log "DEBUG: Input file: $INPUT_FILE"
    log "DEBUG: Output directory: $OUTPUT_DIR"
    log "DEBUG: Running on $(hostname)"
    log "DEBUG: Python version: $(python --version 2>&1)"
fi

# Find PHACTS installation
function find_phacts {
    # 1. Try custom path if provided
    if [[ -n "$CUSTOM_PHACTS_PATH" ]]; then
        if [[ -f "$CUSTOM_PHACTS_PATH" ]]; then
            log "Using custom PHACTS path: $CUSTOM_PHACTS_PATH"
            echo "$CUSTOM_PHACTS_PATH"
            return 0
        else
            log "Warning: Custom PHACTS path not found: $CUSTOM_PHACTS_PATH"
        fi
    fi
    
    # 2. Try path relative to this script (from workflow installation)
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
    WORKFLOW_PHACTS="$SCRIPT_DIR/../output/db/phacts/PHACTS/phacts.py"
    if [[ -f "$WORKFLOW_PHACTS" ]]; then
        log "Using workflow-installed PHACTS: $WORKFLOW_PHACTS"
        echo "$WORKFLOW_PHACTS"
        return 0
    fi
    
    # 3. Check common installation locations
    COMMON_LOCATIONS=(
        "$HOME/Software/PHACTS/PHACTS/phacts.py"
        "$HOME/software/PHACTS/phacts.py"
        "$HOME/Softwares/PHACTS/phacts.py"
        "$HOME/softwares/PHACTS/phacts.py"
        "$HOME/phacts/phacts.py"
        "/usr/local/PHACTS/phacts.py"
    )
    
    for loc in "${COMMON_LOCATIONS[@]}"; do
        if [[ -f "$loc" ]]; then
            log "Found PHACTS in common location: $loc"
            echo "$loc"
            return 0
        fi
    done
    
    # 4. Try to find in PATH
    PHACTS_IN_PATH=$(which phacts.py 2>/dev/null || true)
    if [[ -n "$PHACTS_IN_PATH" ]]; then
        log "Found PHACTS in PATH: $PHACTS_IN_PATH"
        echo "$PHACTS_IN_PATH"
        return 0
    fi
    
    # No PHACTS found
    log "ERROR: Could not find PHACTS installation"
    log "Please install PHACTS using scripts/install_phacts.sh"
    return 1
}

# Find PHACTS and get its directory
PHACTS_PATH=$(find_phacts)
if [[ $? -ne 0 ]]; then
    exit 1
fi

PHACTS_DIR="$(dirname "$PHACTS_PATH")"
PHACTS_PARENT="$(dirname "$PHACTS_DIR")"

# Ensure PHACTS directory has __init__.py (required for module imports)
if [[ ! -f "$PHACTS_DIR/__init__.py" ]]; then
    log "Creating missing __init__.py file in PHACTS directory"
    touch "$PHACTS_DIR/__init__.py"
fi

# Set up proper PYTHONPATH for PHACTS
export PYTHONPATH="$PHACTS_DIR:$PHACTS_PARENT:$PYTHONPATH"
if [[ "$DEBUG" == true ]]; then
    log "DEBUG: PYTHONPATH=$PYTHONPATH"
fi

# Run PHACTS with proper environment
log "Running PHACTS on input file: $INPUT_FILE"
python "$PHACTS_PATH" "$INPUT_FILE" -o "$OUTPUT_DIR" >> "$LOG_FILE" 2>&1

# Check exit status
if [[ $? -eq 0 ]]; then
    # Check if output was created
    if [[ -f "$OUTPUT_DIR/prediction.txt" ]]; then
        log "PHACTS prediction completed successfully"
        exit 0
    else
        log "Warning: PHACTS ran but did not produce a prediction.txt file"
        exit 2
    fi
else
    log "Error: PHACTS execution failed. See $LOG_FILE for details"
    
    # Try alternative execution method as fallback
    log "Attempting alternative execution method..."
    cd "$PHACTS_PARENT" && python -m PHACTS.phacts "$INPUT_FILE" -o "$OUTPUT_DIR" >> "$LOG_FILE" 2>&1
    
    if [[ $? -eq 0 ]] && [[ -f "$OUTPUT_DIR/prediction.txt" ]]; then
        log "PHACTS prediction completed successfully with alternative method"
        exit 0
    else
        log "All execution attempts failed"
        exit 1
    fi
fi