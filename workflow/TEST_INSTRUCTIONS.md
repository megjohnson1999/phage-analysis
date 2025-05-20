# Testing Phage-Specific PHACTS Analysis

This document provides simplified instructions for testing the phage-specific PHACTS analysis.

## How to Run the Test

Use the following command to test with any protein FASTA file:

```bash
# Navigate to the workflow directory
cd workflow/

# Run the test with an absolute path to your protein file
snakemake --configfile ../test_data/test_config.yaml --use-conda --cores 4 \
  --config input_path=/full/path/to/your/protein_file.faa \
  test_phage_specific/results/test_report.txt
```

For example:

```bash
# Example with full path
snakemake --configfile ../test_data/test_config.yaml --use-conda --cores 4 \
  --config input_path=/home/megan.j/no_reneo_full/03_split_proteins/proteins.part_001.faa \
  test_phage_specific/results/test_report.txt
```

## What the Test Does

1. Takes your protein FASTA file as input
2. Splits it by phage ID
3. Runs PHACTS on each phage-specific file
4. Compiles the results

## Viewing Results

After the test completes, you can view:

- The test report: `test_phage_specific/results/test_report.txt`
- The predictions: `test_phage_specific/results/phacts_results/phacts_predictions_compiled.tsv`
- Log files: `test_phage_specific/results/logs/`

## Cleaning Up

To clean up the test files:

```bash
rm -rf test_phage_specific/
```