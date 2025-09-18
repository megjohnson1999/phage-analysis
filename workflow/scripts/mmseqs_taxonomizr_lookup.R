#!/usr/bin/env Rscript
# MMSeqs2 Taxonomic Lookup using taxonomizr
# 
# This script replicates the original R script logic for processing mmseqs2 results
# with taxonomizr, maintaining the exact same workflow and logic.

library(tidyverse)
library(taxonomizr)
library(vroom)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript mmseqs_taxonomizr_lookup.R <input_mmseqs> <taxonomizr_db> <output_file>")
}

input_file <- args[1]
taxonomizr_db <- args[2]
output_file <- args[3]

cat("Processing mmseqs2 results with taxonomizr...\n")
cat("Input file:", input_file, "\n")
cat("Taxonomizr database:", taxonomizr_db, "\n")
cat("Output file:", output_file, "\n")

# Load mmseqs2 results
# Check if file has headers or not
first_line <- readLines(input_file, n = 1)
if (grepl("^(query|contig|edge|virus_comp|scaffold|NODE)", first_line)) {
    # File has headers
    mmseqs_topHit <- vroom(file = input_file, col_names = TRUE)
} else {
    # File doesn't have headers - use the same column names as original script
    mmseqs_topHit <- vroom(file = input_file, col_names = FALSE)
    colnames(mmseqs_topHit) <- c("query", "target", "evalue", "pident", "fident", "nident", "mismatch", "qcov", "tcov", "qstart", "qend", "qlen", "tstart", "tend", "tlen", "alnlen", "bits", "qheader", "theader", "taxid", "taxname", "taxlineage")
}

cat("Loaded", nrow(mmseqs_topHit), "mmseqs2 results\n")

# Replicate the exact filtering logic from the original script
mmseqs_best_hits_by_bitscore <- mmseqs_topHit %>%
  filter(grepl("d_Viruses", taxlineage, ignore.case = TRUE)) %>%
  filter(pident >= 70) %>%
  group_by(query) %>%
  slice(which.max(bits)) %>%
  ungroup() %>%
  select(query, taxid)

cat("Found", nrow(mmseqs_best_hits_by_bitscore), "viral hits with pident >= 70\n")

if (nrow(mmseqs_best_hits_by_bitscore) == 0) {
    cat("Warning: No viral hits found. Creating empty output file.\n")
    # Create empty dataframe with proper structure
    empty_result <- data.frame(
        contigID = character(0),
        blast_superkingdom = character(0),
        blast_phylum = character(0),
        blast_class = character(0),
        blast_order = character(0),
        blast_family = character(0),
        blast_genus = character(0),
        blast_species = character(0)
    )
    vroom_write(empty_result, output_file)
    quit(status = 0)
}

# Use taxonomizr exactly as in the original script
cat("Looking up taxonomic lineages using taxonomizr...\n")
linage <- as.data.frame(getTaxonomy(mmseqs_best_hits_by_bitscore$taxid, taxonomizr_db))

# Check for NA values (taxonomic IDs not found in database)
na_count <- sum(is.na(linage$superkingdom))
if (na_count > 0) {
    cat("Warning:", na_count, "taxonomic IDs could not be resolved\n")
}

# Add taxonomy information exactly as in the original script
mmseqs_best_hits_by_bitscore$blast_superkingdom <- linage$superkingdom
mmseqs_best_hits_by_bitscore$blast_phylum  <- linage$phylum
mmseqs_best_hits_by_bitscore$blast_class  <- linage$class
mmseqs_best_hits_by_bitscore$blast_order  <- linage$order
mmseqs_best_hits_by_bitscore$blast_family  <- linage$family
mmseqs_best_hits_by_bitscore$blast_genus  <- linage$genus
mmseqs_best_hits_by_bitscore$blast_species  <- linage$species

# Rename query column to contigID to match expected format
mmseqs_best_hits_by_bitscore <- mmseqs_best_hits_by_bitscore %>%
    select(-taxid) %>%
    rename(contigID = query)

# Adding superkingdom to the unclassified (matching original script logic)
mmseqs_best_hits_by_bitscore <- mmseqs_best_hits_by_bitscore %>%
    mutate(blast_superkingdom = ifelse(is.na(blast_superkingdom), "Viruses", blast_superkingdom))

# Log summary statistics
cat("Taxonomy lookup completed:\n")
cat("  - Total contigs:", nrow(mmseqs_best_hits_by_bitscore), "\n")
for (level in c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) {
    col_name <- paste0("blast_", level)
    non_null <- sum(!is.na(mmseqs_best_hits_by_bitscore[[col_name]]))
    percentage <- round(100 * non_null / nrow(mmseqs_best_hits_by_bitscore), 1)
    cat("  -", level, ":", non_null, "assignments (", percentage, "%)\n")
}

# Save results
vroom_write(mmseqs_best_hits_by_bitscore, output_file)
cat("Results saved to:", output_file, "\n")
cat("MMSeqs2 taxonomizr lookup completed successfully!\n")