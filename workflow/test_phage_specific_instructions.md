# Testing Phage-Specific PHACTS Analysis

This document provides instructions for testing the phage-specific PHACTS analysis using Snakemake rules.

## Testing with a Specific Protein File

You can test the phage-specific splitting and PHACTS analysis with any protein FASTA file using the following command:

```bash
# Navigate to the workflow directory
cd workflow/

# Run the test rule with your protein file
snakemake --use-conda \
    --cores 4 \
    test_phage_specific/path/to/your/protein/file.faa/test_report.txt
```

Replace `path/to/your/protein/file.faa` with the path to your protein FASTA file.

### Example

```bash
# Test with a protein file from a previous run
snakemake --use-conda \
    --cores 4 \
    test_phage_specific/no_reneo_full/03_split_proteins/proteins.part_001.faa/test_report.txt
```

## What the Test Does

The test performs the following steps:

1. Splits the input protein file by phage ID
2. Runs PHACTS prediction on each phage-specific file
3. Aggregates the results into a single TSV file
4. Generates a test report with summary information

## Test Output

The test creates the following directory structure:

```
test_phage_specific/
└── path/to/your/protein/file.faa/
    ├── logs/                           # Log files
    ├── split_list.txt                  # List of split files
    ├── split_proteins/                 # Directory with phage-specific protein files
    │   ├── phage1.faa
    │   ├── phage2.faa
    │   └── ...
    ├── phacts_results/
    │   ├── phacts_predictions_compiled.tsv  # Final predictions
    │   └── tmp/                            # Temporary files
    │       ├── phage1/
    │       ├── phage2/
    │       └── ...
    └── test_report.txt                 # Summary report
```

## Viewing Results

The main results file is:
```
test_phage_specific/path/to/your/protein/file.faa/phacts_results/phacts_predictions_compiled.tsv
```

The summary report is:
```
test_phage_specific/path/to/your/protein/file.faa/test_report.txt
```

## Running Individual Test Steps

If you want to run specific test steps:

```bash
# Just run the splitting
snakemake --use-conda \
    test_phage_specific/path/to/your/protein/file.faa/split_list.txt

# Run the splitting and PHACTS prediction
snakemake --use-conda \
    test_phage_specific/path/to/your/protein/file.faa/phacts_results/phacts_predictions_compiled.tsv
```

## Cleanup

To clean up test files:

```bash
rm -rf workflow/test_phage_specific/
```