# Phage-Specific PHACTS Analysis

This branch adds a phage-specific approach to PHACTS lifestyle prediction, addressing an issue where the original workflow could assign the same prediction to multiple phages.

## Problem Being Solved

In the original workflow:
1. Prodigal predicts proteins from all phages, creating a single `proteins.faa` file
2. This file is split into arbitrary batches (by sequence count, not phage ID)
3. PHACTS makes a single lifestyle prediction for each batch
4. The output aggregation assigns this single prediction to all phage IDs found in the batch

This is problematic because a single PHACTS prediction is applied to proteins from multiple phages, which may have different lifestyles.

## Solution

This branch adds a new set of rules that:

1. Splits proteins by phage ID, ensuring each file contains proteins from only a single phage
2. Runs PHACTS on each phage-specific file, producing accurate per-phage predictions
3. Aggregates results with a simple 1:1 mapping (filename to phage ID)

## New Files

- `workflow/scripts/split_proteins_by_phage.py`: Python script that extracts phage IDs from protein headers and splits files accordingly
- `workflow/rules/split_proteins_by_phage.smk`: Snakemake rules for phage-specific PHACTS analysis
- `scripts/test_phage_splitting.sh`: Utility script to test the phage ID extraction

## Header Format and Phage ID Extraction

The script assumes protein headers follow a format like:
```
>read_13355_7 # 4873 # 5514 # 1 # ID=8902_7;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.352
```

The phage ID is extracted as the part before the last underscore and number (e.g., `read_13355` from `read_13355_7`).

## Usage

The workflow will now run both the original PHACTS analysis (for backward compatibility) and the new phage-specific analysis. The new results will be available in:

```
{output_dir}/03_phacts_results_by_phage/phacts_predictions_compiled.tsv
```

You can test the phage ID extraction on a sample file using:
```bash
./scripts/test_phage_splitting.sh /path/to/your/proteins.part_xxx.faa
```

## Benefits

- Each PHACTS prediction is made using only proteins from a single phage
- Predictions accurately reflect the lifestyle of specific phages
- No ambiguity in associating predictions with phages
- Results are more biologically meaningful