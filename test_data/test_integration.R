#!/usr/bin/env Rscript
# Test script to verify the phage prediction integration works

# Set working directory to test_data
setwd(dirname(sys.frame(1)$ofile))

# Source the main script (modify path as needed)
source("../workflow/scripts/01_phagePrediction.R")

# Test with mock data
process_phage_predictions(
  phold_file = "mock_phold_results.tsv",
  jager_file = "mock_jaeger_results.tsv", 
  genomad_file = "mock_genomad_results.tsv",
  checkv_file = "mock_checkv_results.tsv",
  output_table = "test_phage_predictions.tsv",
  output_ids = "test_contig_ids.txt"
)

cat("Test completed successfully!\n")
cat("Output files created:\n")
cat("- test_phage_predictions.tsv\n") 
cat("- test_contig_ids.txt\n")