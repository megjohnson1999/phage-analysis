#!/bin/bash

################################################################################
# Stage 3 Output Verification Script
#
# This script verifies that all expected outputs from Stage 3 (Comprehensive
# Characterization) were generated correctly when starting from clustered
# phage sequences.
#
# Usage:
#   bash verify_stage3_outputs.sh <output_dir> [input_clustered_seqs]
#
# Arguments:
#   output_dir            - Path to pipeline output directory
#   input_clustered_seqs  - (Optional) Path to input clustered sequences FASTA
#
# Example:
#   bash verify_stage3_outputs.sh /path/to/outputs /path/to/vOTU_repSeqs.fasta
################################################################################

set -euo pipefail

# Color codes for output (disable if not supported)
if [[ -t 1 ]]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    BOLD='\033[1m'
    NC='\033[0m' # No Color
else
    GREEN=''
    RED=''
    YELLOW=''
    BLUE=''
    BOLD=''
    NC=''
fi

# Parse arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_dir> [input_clustered_seqs]"
    echo ""
    echo "Arguments:"
    echo "  output_dir            - Path to pipeline output directory"
    echo "  input_clustered_seqs  - (Optional) Path to input clustered sequences FASTA"
    exit 1
fi

OUTDIR="$1"
INPUT_SEQS="${2:-}"

# Verify output directory exists
if [ ! -d "$OUTDIR" ]; then
    echo -e "${RED}Error: Output directory does not exist: $OUTDIR${NC}"
    exit 1
fi

echo -e "${BOLD}============================================================${NC}"
echo -e "${BOLD}Stage 3 Output Verification${NC}"
echo -e "${BOLD}============================================================${NC}"
echo -e "Output directory: ${BLUE}$OUTDIR${NC}"
if [ -n "$INPUT_SEQS" ]; then
    echo -e "Input sequences:  ${BLUE}$INPUT_SEQS${NC}"
fi
echo ""

# Track overall status
MISSING_FILES=0
EMPTY_FILES=0
WARNINGS=0
TOTAL_CHECKS=0

################################################################################
# 1. Check all expected output files exist
################################################################################
echo -e "${BOLD}[1/6] Checking for expected output files...${NC}"
echo ""

EXPECTED_FILES=(
    "03_checkv_final/quality_summary.tsv"
    "03_iphop_results/iphop_predictions_compiled.tsv"
    "03_orf_predictions/proteins.faa"
    "03_orf_predictions/genes.fna"
    "03_genomic_info/mmseqs_taxonomy.tsv"
    "03_genomic_info/phabox_output/taxonomy.tsv"
    "03_genomic_info/phabox_output/lifestyle.tsv"
    "03_genomic_info/bacphlip_lifestyle.tsv"
    "03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"
    "03_genomic_info/consensus_taxonomy.tsv"
    "03_genomic_info/lifestyle_consensus.tsv"
    "final_contig_summary.tsv"
)

for file in "${EXPECTED_FILES[@]}"; do
    TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
    filepath="$OUTDIR/$file"

    if [ -f "$filepath" ]; then
        # Get file size (works on both Linux and macOS)
        if command -v stat >/dev/null 2>&1; then
            size=$(stat -f%z "$filepath" 2>/dev/null || stat -c%s "$filepath" 2>/dev/null)
        else
            size=$(wc -c < "$filepath")
        fi
        echo -e "${GREEN}✓${NC} $file (${size} bytes)"
    else
        echo -e "${RED}✗ MISSING${NC}: $file"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

# Check vContact3 directory
TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
if [ -d "$OUTDIR/03_genomic_info/vc3_output" ]; then
    file_count=$(find "$OUTDIR/03_genomic_info/vc3_output" -type f | wc -l)
    echo -e "${GREEN}✓${NC} 03_genomic_info/vc3_output/ (directory exists, ${file_count} files)"
else
    echo -e "${RED}✗ MISSING${NC}: 03_genomic_info/vc3_output/"
    MISSING_FILES=$((MISSING_FILES + 1))
fi

echo ""

################################################################################
# 2. Count sequences in each output
################################################################################
echo -e "${BOLD}[2/6] Counting sequences in outputs...${NC}"
echo ""

# Function to count sequences in FASTA
count_fasta() {
    if [ -f "$1" ]; then
        grep -c "^>" "$1" 2>/dev/null || echo "0"
    else
        echo "0"
    fi
}

# Function to count lines in TSV (excluding header)
count_tsv() {
    if [ -f "$1" ]; then
        tail -n +2 "$1" 2>/dev/null | wc -l | tr -d ' '
    else
        echo "0"
    fi
}

# Count input sequences if provided
if [ -n "$INPUT_SEQS" ] && [ -f "$INPUT_SEQS" ]; then
    INPUT_COUNT=$(count_fasta "$INPUT_SEQS")
    echo "Input sequences (clustered reps): $INPUT_COUNT"
else
    INPUT_COUNT=""
    echo "Input sequences: Not provided"
fi

echo ""
echo "Output sequence counts:"
echo "  CheckV quality assessments:       $(count_tsv "$OUTDIR/03_checkv_final/quality_summary.tsv")"
echo "  iPhop host predictions:           $(count_tsv "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv")"
echo "  Predicted proteins (ORFs):        $(count_fasta "$OUTDIR/03_orf_predictions/proteins.faa")"
echo "  Predicted genes:                  $(count_fasta "$OUTDIR/03_orf_predictions/genes.fna")"
echo "  MMseqs2 taxonomy:                 $(count_tsv "$OUTDIR/03_genomic_info/mmseqs_taxonomy.tsv")"
echo "  Phabox2 taxonomy:                 $(count_tsv "$OUTDIR/03_genomic_info/phabox_output/taxonomy.tsv")"
echo "  Phabox2 lifestyle:                $(count_tsv "$OUTDIR/03_genomic_info/phabox_output/lifestyle.tsv")"
echo "  BACPHLIP lifestyle:               $(count_tsv "$OUTDIR/03_genomic_info/bacphlip_lifestyle.tsv")"
echo "  Consensus taxonomy:               $(count_tsv "$OUTDIR/03_genomic_info/consensus_taxonomy.tsv")"
echo "  Lifestyle consensus:              $(count_tsv "$OUTDIR/03_genomic_info/lifestyle_consensus.tsv")"
echo "  Final contig summary:             $(count_tsv "$OUTDIR/final_contig_summary.tsv")"

echo ""

################################################################################
# 3. Check for empty files or failures
################################################################################
echo -e "${BOLD}[3/6] Checking for empty or problematic files...${NC}"
echo ""

CRITICAL_FILES=(
    "03_checkv_final/quality_summary.tsv"
    "03_iphop_results/iphop_predictions_compiled.tsv"
    "03_orf_predictions/proteins.faa"
    "final_contig_summary.tsv"
)

for file in "${CRITICAL_FILES[@]}"; do
    TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
    filepath="$OUTDIR/$file"

    if [ -f "$filepath" ]; then
        lines=$(wc -l < "$filepath" 2>/dev/null | tr -d ' ')
        if [ "$lines" -le 1 ]; then
            echo -e "${YELLOW}⚠ WARNING${NC}: $file appears empty or has only header (${lines} lines)"
            EMPTY_FILES=$((EMPTY_FILES + 1))
        else
            echo -e "${GREEN}✓${NC} $file has data (${lines} lines)"
        fi
    else
        echo -e "${RED}✗${NC} $file does not exist"
    fi
done

echo ""

################################################################################
# 4. Sample the data quality
################################################################################
echo -e "${BOLD}[4/6] Sampling data quality...${NC}"
echo ""

# CheckV completeness distribution
if [ -f "$OUTDIR/03_checkv_final/quality_summary.tsv" ]; then
    echo "CheckV quality assessment distribution:"
    tail -n +2 "$OUTDIR/03_checkv_final/quality_summary.tsv" | \
        cut -f7 | sort | uniq -c | \
        awk '{printf "  %-20s %s\n", $2, $1}'
    echo ""
fi

# iPhop prediction success rate
if [ -f "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" ]; then
    echo "iPhop host prediction success:"
    tail -n +2 "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" | \
        awk -F'\t' '{
            if ($2 != "" && $2 != "NA" && $2 != "No host prediction")
                print "  Has prediction"
            else
                print "  No prediction"
        }' | sort | uniq -c | awk '{printf "  %-20s %s sequences\n", $2" "$3, $1}'
    echo ""
fi

# Lifestyle consensus
if [ -f "$OUTDIR/03_genomic_info/lifestyle_consensus.tsv" ]; then
    echo "Lifestyle prediction distribution:"
    tail -n +2 "$OUTDIR/03_genomic_info/lifestyle_consensus.tsv" | \
        cut -f2 | sort | uniq -c | \
        awk '{printf "  %-20s %s\n", $2, $1}'
    echo ""
fi

# Taxonomy family level (top 10)
if [ -f "$OUTDIR/03_genomic_info/consensus_taxonomy.tsv" ]; then
    echo "Top 10 phage families (consensus taxonomy):"
    tail -n +2 "$OUTDIR/03_genomic_info/consensus_taxonomy.tsv" | \
        cut -f5 | grep -v "^$" | sort | uniq -c | sort -rn | head -10 | \
        awk '{printf "  %-40s %s\n", $2, $1}'
    echo ""
fi

################################################################################
# 5. Check log files for errors
################################################################################
echo -e "${BOLD}[5/6] Checking for errors in log files...${NC}"
echo ""

# Find log files in common locations
LOG_DIRS=(
    "$OUTDIR"
    "$OUTDIR/.snakemake/log"
    "$OUTDIR/logs"
)

found_logs=0
error_logs=0

for logdir in "${LOG_DIRS[@]}"; do
    if [ -d "$logdir" ]; then
        while IFS= read -r -d '' logfile; do
            found_logs=$((found_logs + 1))
            if grep -qi "error\|failed\|exception" "$logfile" 2>/dev/null; then
                echo -e "${YELLOW}⚠ Found potential errors in${NC}: $(basename "$logfile")"
                grep -i "error\|failed\|exception" "$logfile" 2>/dev/null | head -3 | sed 's/^/    /'
                echo ""
                error_logs=$((error_logs + 1))
                WARNINGS=$((WARNINGS + 1))
            fi
        done < <(find "$logdir" -maxdepth 2 -type f \( -name "*.log" -o -name "*.err" \) -print0 2>/dev/null)
    fi
done

if [ $found_logs -eq 0 ]; then
    echo "No log files found in standard locations"
elif [ $error_logs -eq 0 ]; then
    echo -e "${GREEN}✓${NC} No errors detected in $found_logs log files"
fi

echo ""

################################################################################
# 6. Verify sequence ID consistency
################################################################################
echo -e "${BOLD}[6/6] Verifying sequence ID consistency...${NC}"
echo ""

# Create temp directory for ID files
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Extract sequence IDs from different files
if [ -f "$OUTDIR/03_checkv_final/quality_summary.tsv" ]; then
    tail -n +2 "$OUTDIR/03_checkv_final/quality_summary.tsv" | cut -f1 | sort > "$TMPDIR/checkv_ids.txt"
fi

if [ -f "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" ]; then
    tail -n +2 "$OUTDIR/03_iphop_results/iphop_predictions_compiled.tsv" | cut -f1 | sort > "$TMPDIR/iphop_ids.txt"
fi

if [ -f "$OUTDIR/final_contig_summary.tsv" ]; then
    tail -n +2 "$OUTDIR/final_contig_summary.tsv" | cut -f1 | sort > "$TMPDIR/summary_ids.txt"
fi

# Compare counts
checkv_count=$(wc -l < "$TMPDIR/checkv_ids.txt" 2>/dev/null | tr -d ' ')
iphop_count=$(wc -l < "$TMPDIR/iphop_ids.txt" 2>/dev/null | tr -d ' ')
summary_count=$(wc -l < "$TMPDIR/summary_ids.txt" 2>/dev/null | tr -d ' ')

echo "Sequence counts by file:"
echo "  CheckV quality:     $checkv_count"
echo "  iPhop predictions:  $iphop_count"
echo "  Final summary:      $summary_count"
echo ""

# Check if counts match
TOTAL_CHECKS=$((TOTAL_CHECKS + 1))
if [ "$checkv_count" == "$summary_count" ] && [ "$checkv_count" == "$iphop_count" ]; then
    echo -e "${GREEN}✓${NC} All sequence counts are consistent"

    # Check if IDs actually match (not just counts)
    if diff "$TMPDIR/checkv_ids.txt" "$TMPDIR/summary_ids.txt" >/dev/null 2>&1; then
        echo -e "${GREEN}✓${NC} All sequence IDs match across outputs"
    else
        echo -e "${YELLOW}⚠ WARNING${NC}: Sequence IDs differ between files (but counts match)"
        WARNINGS=$((WARNINGS + 1))
    fi
else
    echo -e "${RED}✗${NC} Sequence count mismatch detected!"
    WARNINGS=$((WARNINGS + 1))

    # Show which sequences might be missing
    if [ -s "$TMPDIR/checkv_ids.txt" ] && [ -s "$TMPDIR/summary_ids.txt" ]; then
        missing=$(comm -23 "$TMPDIR/checkv_ids.txt" "$TMPDIR/summary_ids.txt" | wc -l | tr -d ' ')
        if [ "$missing" -gt 0 ]; then
            echo "  Sequences in CheckV but missing from summary: $missing"
            echo "  First 5 missing IDs:"
            comm -23 "$TMPDIR/checkv_ids.txt" "$TMPDIR/summary_ids.txt" | head -5 | sed 's/^/    /'
        fi
    fi
fi

echo ""

################################################################################
# Final Summary
################################################################################
echo -e "${BOLD}============================================================${NC}"
echo -e "${BOLD}Verification Summary${NC}"
echo -e "${BOLD}============================================================${NC}"
echo ""

if [ $MISSING_FILES -eq 0 ] && [ $EMPTY_FILES -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}${BOLD}✓ ALL CHECKS PASSED!${NC}"
    echo ""
    echo "Your Stage 3 outputs appear complete and consistent."
    EXIT_CODE=0
else
    echo -e "${YELLOW}⚠ Some issues were detected:${NC}"
    echo "  Missing files:      $MISSING_FILES"
    echo "  Empty files:        $EMPTY_FILES"
    echo "  Warnings:           $WARNINGS"
    echo ""

    if [ $MISSING_FILES -gt 0 ]; then
        echo -e "${RED}Action needed:${NC} Some expected output files are missing."
        echo "  - Check if the pipeline completed successfully"
        echo "  - Review Snakemake logs for failed rules"
        EXIT_CODE=1
    elif [ $EMPTY_FILES -gt 0 ]; then
        echo -e "${YELLOW}Action needed:${NC} Some output files are empty."
        echo "  - This might indicate tool failures or no data passing filters"
        echo "  - Check individual tool log files for details"
        EXIT_CODE=1
    else
        echo -e "${YELLOW}Note:${NC} Minor warnings detected but critical files are present."
        echo "  - Review warnings above for details"
        EXIT_CODE=0
    fi
fi

echo ""
echo "For detailed results, examine:"
echo "  • Final summary table: $OUTDIR/final_contig_summary.tsv"
if [ -f "$OUTDIR/Pipeline_Summary_Report.html" ]; then
    echo "  • HTML report: $OUTDIR/Pipeline_Summary_Report.html"
fi
echo ""

exit $EXIT_CODE
