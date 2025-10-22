#!/bin/bash

################################################################################
# Regenerate Final Contig Summary Table
#
# This script re-runs just the final summary table generation with the fixed
# Python script, without re-running the entire pipeline.
################################################################################

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_dir>"
    echo ""
    echo "Example:"
    echo "  $0 test_bacphlip"
    exit 1
fi

OUTDIR="$1"

# Verify output directory exists
if [ ! -d "$OUTDIR" ]; then
    echo "Error: Output directory does not exist: $OUTDIR"
    exit 1
fi

echo "============================================================"
echo "Regenerating Final Contig Summary Table"
echo "============================================================"
echo "Output directory: $OUTDIR"
echo ""

# Backup old summary if it exists
if [ -f "$OUTDIR/final_contig_summary.tsv" ]; then
    echo "Backing up old summary table..."
    cp "$OUTDIR/final_contig_summary.tsv" "$OUTDIR/final_contig_summary.tsv.backup"
    echo "  Backup saved: $OUTDIR/final_contig_summary.tsv.backup"
    echo ""
fi

# Determine the input sequences file
if [ -f "$OUTDIR/02_clustering/vOTU_repSeqs.fasta" ]; then
    INPUT_SEQS="$OUTDIR/02_clustering/vOTU_repSeqs.fasta"
elif [ -f "$OUTDIR/01_phage_prediction/phageContigs.fasta" ]; then
    INPUT_SEQS="$OUTDIR/01_phage_prediction/phageContigs.fasta"
else
    # Look for any clustered or phage sequences in the output dir
    echo "Searching for input sequences..."
    INPUT_SEQS=$(find "$OUTDIR" -name "vOTU_repSeqs.fasta" -o -name "phageContigs.fasta" | head -1)
    if [ -z "$INPUT_SEQS" ]; then
        echo "Error: Could not find input sequences file in $OUTDIR"
        exit 1
    fi
fi

echo "Input sequences: $INPUT_SEQS"
echo ""

# Build the command
CMD="python workflow/scripts/create_final_contig_table.py"
CMD="$CMD --phage-seqs $INPUT_SEQS"

# Add CheckV results if available
if [ -f "$OUTDIR/03_checkv_final/quality_summary.tsv" ]; then
    CMD="$CMD --checkv-results $OUTDIR/03_checkv_final/quality_summary.tsv"
    echo "Found: CheckV quality results"
fi

# Add consensus taxonomy if available
if [ -f "$OUTDIR/03_genomic_info/consensus_taxonomy.tsv" ]; then
    CMD="$CMD --consensus-taxonomy $OUTDIR/03_genomic_info/consensus_taxonomy.tsv"
    echo "Found: Consensus taxonomy"
fi

# Add lifestyle consensus if available
if [ -f "$OUTDIR/03_genomic_info/lifestyle_consensus.tsv" ]; then
    CMD="$CMD --lifestyle-consensus $OUTDIR/03_genomic_info/lifestyle_consensus.tsv"
    echo "Found: Lifestyle consensus"
fi

# Add iPhop predictions if available
if [ -f "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" ]; then
    CMD="$CMD --iphop-predictions $OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv"
    echo "Found: iPhop host predictions"
fi

# Set output path
CMD="$CMD --output $OUTDIR/final_contig_summary.tsv"

echo ""
echo "Running command:"
echo "$CMD"
echo ""

# Execute the command
eval $CMD

echo ""
echo "============================================================"
echo "Done!"
echo "============================================================"
echo ""
echo "New summary table: $OUTDIR/final_contig_summary.tsv"
if [ -f "$OUTDIR/final_contig_summary.tsv.backup" ]; then
    echo "Old backup:        $OUTDIR/final_contig_summary.tsv.backup"
fi
echo ""
echo "Check the new summary:"
echo "  head $OUTDIR/final_contig_summary.tsv"
echo ""
echo "Compare with backup:"
echo "  diff <(head -1 $OUTDIR/final_contig_summary.tsv.backup) <(head -1 $OUTDIR/final_contig_summary.tsv)"
echo ""
