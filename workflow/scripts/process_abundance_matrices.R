#!/usr/bin/env Rscript
#
# Process CoverM abundance data into wide-format matrices
#
# Creates TPM and count matrices filtered by coverage threshold
# Based on logic from psoConstructoin_genomicInfoFormating.Rmd

suppressPackageStartupMessages({
  library(tidyverse)
  library(vroom)
  library(reshape2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: process_abundance_matrices.R <coverm_input> <tpm_output> <count_output> <coverage_threshold>")
}

coverm_file <- args[1]
tpm_output <- args[2]
count_output <- args[3]
coverage_threshold <- as.numeric(args[4])

cat("========================================\n")
cat("Processing Abundance Matrices\n")
cat("========================================\n\n")

cat("Parameters:\n")
cat(sprintf("  CoverM input:        %s\n", coverm_file))
cat(sprintf("  TPM output:          %s\n", tpm_output))
cat(sprintf("  Count output:        %s\n", count_output))
cat(sprintf("  Coverage threshold:  %.2f\n", coverage_threshold))
cat("\n")

# Load CoverM output (long format)
cat("Loading CoverM results...\n")
abd_cov <- vroom(coverm_file,
                 col_names = c("sampleID", "contigID", "rpkm", "tpm",
                               "count", "variance", "mean",
                               "covered_fraction", "covered_bases"),
                 skip = 1,  # Skip header row from CoverM output
                 show_col_types = FALSE)

cat(sprintf("  Loaded %d rows\n", nrow(abd_cov)))
cat(sprintf("  Samples: %d\n", length(unique(abd_cov$sampleID))))
cat(sprintf("  Contigs: %d\n", length(unique(abd_cov$contigID))))
cat("\n")

# Clean sample names (remove common suffixes)
abd_cov <- abd_cov %>%
  mutate(sampleID = gsub("_stats$", "", sampleID))

# Process TPM matrix
cat("Creating TPM abundance matrix...\n")
cat(sprintf("  Filtering by coverage >= %.2f\n", coverage_threshold))

abd_cov_filt_tpm <- abd_cov %>%
  mutate(newTPM = ifelse(covered_fraction >= coverage_threshold, tpm, 0))

# Reshape to wide format (contigs × samples)
tpm_matrix <- dcast(abd_cov_filt_tpm, contigID ~ sampleID, value.var = "newTPM")

# Remove contigs with 0 total abundance
tpm_matrix <- tpm_matrix %>%
  mutate(suma = rowSums(across(where(is.numeric)))) %>%
  filter(suma > 0) %>%
  select(-suma)

cat(sprintf("  TPM matrix: %d contigs × %d samples\n", nrow(tpm_matrix), ncol(tpm_matrix) - 1))
cat(sprintf("  Contigs removed (0 abundance): %d\n",
            length(unique(abd_cov$contigID)) - nrow(tpm_matrix)))

# Save TPM matrix
write_tsv(tpm_matrix, tpm_output)
cat(sprintf("  Saved: %s\n", tpm_output))
cat("\n")

# Process Count matrix
cat("Creating count abundance matrix...\n")
cat(sprintf("  Filtering by coverage >= %.2f\n", coverage_threshold))

abd_cov_filt_count <- abd_cov %>%
  mutate(newCount = ifelse(covered_fraction >= coverage_threshold, count, 0))

# Reshape to wide format (contigs × samples)
count_matrix <- dcast(abd_cov_filt_count, contigID ~ sampleID, value.var = "newCount")

# Remove contigs with 0 total abundance
count_matrix <- count_matrix %>%
  mutate(suma = rowSums(across(where(is.numeric)))) %>%
  filter(suma > 0) %>%
  select(-suma)

cat(sprintf("  Count matrix: %d contigs × %d samples\n", nrow(count_matrix), ncol(count_matrix) - 1))
cat(sprintf("  Contigs removed (0 abundance): %d\n",
            length(unique(abd_cov$contigID)) - nrow(count_matrix)))

# Save count matrix
write_tsv(count_matrix, count_output)
cat(sprintf("  Saved: %s\n", count_output))
cat("\n")

cat("========================================\n")
cat("Abundance matrices created successfully!\n")
cat("========================================\n")
