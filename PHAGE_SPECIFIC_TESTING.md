# Phage-Specific PHACTS Analysis Testing Guide

## Overview

This document provides a comprehensive guide for testing the phage-specific PHACTS analysis functionality that splits proteins by phage ID before running PHACTS predictions.

## What This Feature Does

The phage-specific PHACTS analysis:
1. Takes protein FASTA files as input
2. Splits them by phage ID using the header format
3. Runs PHACTS prediction on each phage's proteins separately
4. Compiles the results into a single TSV file

This approach improves the accuracy of PHACTS lifestyle predictions by ensuring that each prediction is based only on proteins from a single phage.

## Testing via Snakemake (Recommended)

The easiest way to test this feature is to use the integrated test rules in Snakemake:

```bash
# Navigate to the workflow directory
cd workflow/

# Run the test with a protein file
snakemake --use-conda \
    --cores 4 \
    test_phage_specific/path/to/your/protein/file.faa/test_report.txt
```

For example:
```bash
snakemake --use-conda \
    --cores 4 \
    test_phage_specific/no_reneo_full/03_split_proteins/proteins.part_001.faa/test_report.txt
```

### Test Output Location

The test creates output in the `workflow/test_phage_specific/` directory, with a subdirectory specific to your input file path. For example:

```
workflow/test_phage_specific/no_reneo_full/03_split_proteins/proteins.part_001.faa/
```

This directory will contain:
- Split protein files
- PHACTS prediction results
- Log files
- A summary test report

### Running Specific Test Steps

If you want to run only specific parts of the test:

```bash
# Just the splitting step
snakemake --use-conda \
    test_phage_specific/path/to/your/protein/file.faa/split_list.txt

# Splitting and PHACTS prediction
snakemake --use-conda \
    test_phage_specific/path/to/your/protein/file.faa/phacts_results/phacts_predictions_compiled.tsv
```

## Running in the Full Workflow

The phage-specific PHACTS analysis is integrated into the main workflow. To run it as part of the full pipeline:

```bash
# From the workflow directory
snakemake --use-conda \
    --cores 4 \
    {output_dir}/03_phacts_results_by_phage/phacts_predictions_compiled.tsv
```

Where `{output_dir}` is your configured output directory.

## Debugging and Logs

Log files are stored in:
- Test mode: `workflow/test_phage_specific/{input_file}/logs/`
- Regular workflow: `{output_dir}/logs/`

## Expected Output Format

The final predictions file will have the format:

```
phage_id    lifestyle    probability
contig_123  temperate    0.92
contig_456  virulent     0.88
...
```

## Clean Up

To clean up the test files:

```bash
rm -rf workflow/test_phage_specific/
```