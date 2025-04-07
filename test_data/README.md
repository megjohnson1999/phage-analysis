# Phage Analysis Pipeline Test Data

This directory contains a minimal test dataset for validating the phage analysis pipeline functionality.

## Test Dataset Content

- **Assemblies**: Single test assembly with 5 contigs:
  - 1 complete phage contig (30kb)
  - 1 incomplete phage contig (5kb)
  - 1 phage-related element (3kb)
  - 1 bacterial contig (8kb)
  - 1 short contig below size threshold (800bp)

- **Reads**: Simulated paired-end reads for the test sample
  - Limited read pairs for testing only

- **Mock Databases**: Empty placeholder files for testing the workflow structure

## Running the Test

The test is designed to validate the workflow structure without requiring full databases:

```bash
# Perform a dry run to check execution plan
./test_data/run_test.sh

# For actual execution (after removing -n flag):
# ./test_data/run_test.sh
```

## Expected Results

- Pipeline should identify 3 contigs as potential phages
- Clustering should generate at least one viral OTU
- Short contig should be filtered out
- Bacterial contig should be classified as non-viral

## Notes

- This is a minimal test for workflow validation
- For actual analysis, real reference databases are required
- Modify test_data/test_config.yaml to adjust parameters