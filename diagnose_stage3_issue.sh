#!/bin/bash

################################################################################
# Stage 3 Diagnostic Script
#
# This script investigates issues found during Stage 3 verification
################################################################################

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_dir>"
    exit 1
fi

OUTDIR="$1"

echo "============================================================"
echo "Stage 3 Diagnostic Analysis"
echo "============================================================"
echo ""

################################################################################
# 1. Check what rules ran vs didn't run
################################################################################
echo "[1] Checking which Snakemake rules completed..."
echo ""

# Look for checkpoint/completion markers
echo "Checkpoint files:"
find "$OUTDIR" -name ".*.done" -o -name "*_ready" -o -name "*_done" 2>/dev/null | sort

echo ""
echo "Recent files created (last 10):"
find "$OUTDIR" -type f -printf '%T+ %p\n' 2>/dev/null | sort -r | head -10

echo ""

################################################################################
# 2. Check if this output dir has data from previous runs
################################################################################
echo "[2] Checking for mixed data from different runs..."
echo ""

# Check iPhop sequence IDs
if [ -f "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" ]; then
    echo "First 3 sequence IDs from iPhop:"
    tail -n +2 "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" | cut -f1 | head -3
    echo ""
    echo "Last 3 sequence IDs from iPhop:"
    tail -n +2 "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" | cut -f1 | tail -3
    echo ""
fi

# Check MMseqs2 sequence IDs
if [ -f "$OUTDIR/03_genomic_info/mmseqs_taxonomy.tsv" ]; then
    echo "First 3 sequence IDs from MMseqs2:"
    tail -n +2 "$OUTDIR/03_genomic_info/mmseqs_taxonomy.tsv" | cut -f1 | head -3
    echo ""
    echo "Last 3 sequence IDs from MMseqs2:"
    tail -n +2 "$OUTDIR/03_genomic_info/mmseqs_taxonomy.tsv" | cut -f1 | tail -3
    echo ""
fi

# Check BACPHLIP (should be correct)
if [ -f "$OUTDIR/03_genomic_info/bacphlip_lifestyle.tsv" ]; then
    echo "First 3 sequence IDs from BACPHLIP:"
    tail -n +2 "$OUTDIR/03_genomic_info/bacphlip_lifestyle.tsv" | cut -f1 | head -3
    echo ""
fi

################################################################################
# 3. Check for the actual input file that was used
################################################################################
echo "[3] Looking for input sequence file that was used..."
echo ""

# Check for split sequences (these are created from input)
if [ -f "$OUTDIR/03_split_seqs/split_file_list.txt" ]; then
    echo "Split sequence files were created:"
    head -5 "$OUTDIR/03_split_seqs/split_file_list.txt"
    echo ""

    # Check first split file
    first_split=$(head -1 "$OUTDIR/03_split_seqs/split_file_list.txt")
    if [ -f "$first_split" ]; then
        echo "Sequences in first split file:"
        grep -c "^>" "$first_split"
        echo ""
        echo "First sequence ID in split file:"
        grep "^>" "$first_split" | head -1
        echo ""
    fi
else
    echo "No split_file_list.txt found - splits may not have been created"
fi

################################################################################
# 4. Check for Prodigal and CheckV output directories
################################################################################
echo "[4] Checking for missing tool directories..."
echo ""

dirs_to_check=(
    "03_genomic_info"
    "03_checkv_final"
    "03_prodigal_output"
)

for dir in "${dirs_to_check[@]}"; do
    if [ -d "$OUTDIR/$dir" ]; then
        file_count=$(find "$OUTDIR/$dir" -type f 2>/dev/null | wc -l)
        echo "✓ $dir/ exists ($file_count files)"
        ls -lh "$OUTDIR/$dir/" | head -10
    else
        echo "✗ $dir/ does NOT exist"
    fi
    echo ""
done

################################################################################
# 5. Check log files for specific failures
################################################################################
echo "[5] Checking specific tool log files..."
echo ""

# CheckV log
if [ -f "$OUTDIR/logs/checkv_final.log" ] || [ -f "$OUTDIR/03_checkv_final/checkv.log" ]; then
    echo "=== CheckV log ==="
    checkv_log=$(find "$OUTDIR" -name "*checkv*.log" 2>/dev/null | head -1)
    if [ -n "$checkv_log" ]; then
        echo "Log file: $checkv_log"
        tail -20 "$checkv_log"
    fi
    echo ""
fi

# Prodigal log
if [ -f "$OUTDIR/logs/prodigal.log" ]; then
    echo "=== Prodigal log ==="
    tail -20 "$OUTDIR/logs/prodigal.log"
    echo ""
fi

# Phabox2 log
phabox_log=$(find "$OUTDIR" -name "*phabox*.log" 2>/dev/null | head -1)
if [ -n "$phabox_log" ]; then
    echo "=== Phabox2 log (last 30 lines) ==="
    tail -30 "$phabox_log"
    echo ""
fi

################################################################################
# 6. Check the actual input file used by the pipeline
################################################################################
echo "[6] Checking pipeline configuration..."
echo ""

# Look for config file or workflow info
if [ -f "$OUTDIR/.snakemake/metadata/config.yaml" ]; then
    echo "Found Snakemake config:"
    cat "$OUTDIR/.snakemake/metadata/config.yaml" 2>/dev/null | grep -E "start_from|input_|output_dir" || echo "Config not readable"
    echo ""
fi

################################################################################
# Summary
################################################################################
echo "============================================================"
echo "Diagnostic Summary"
echo "============================================================"
echo ""
echo "Based on the findings above, check:"
echo ""
echo "1. Were CheckV and Prodigal rules actually run?"
echo "   - Look for their log files and output directories"
echo ""
echo "2. Are iPhop/MMseqs2 results from a previous run?"
echo "   - Compare sequence IDs with BACPHLIP (which has correct count)"
echo ""
echo "3. What was the actual input file used?"
echo "   - Check split_file_list.txt and compare IDs"
echo ""
echo "4. Why did Phabox2 fail?"
echo "   - Check the Phabox2 log file above"
echo ""
