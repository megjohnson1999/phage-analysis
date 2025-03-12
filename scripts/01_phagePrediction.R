# Load necessary libraries
packages <- c("dplyr", "vroom", "reshape2")
options(repos = c(CRAN = "https://cran.r-project.org"))
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(vroom)
library(reshape2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript script.R <phold_file> <jager_file> <genomad_file> <checkv_file> <output_dir>")
}

phold_file <- args[1]
jager_file <- args[2]
genomad_file <- args[3]
checkv_file <- args[4]
output_dir <- args[5]

# Read and process pholdVirPassed data
pholdVirPassed <- vroom(phold_file, delim = "\t", col_names = TRUE)
colnames(pholdVirPassed)[7] <- "category"

# Filter pholdVirPassed
pholdVirPassed_filt <- pholdVirPassed %>%
  filter(product != "hypothetical protein") %>%
  filter(category != "other") %>%
  filter(category != "moron, auxiliary metabolic gene and host takeover") %>%
  filter(category != "transcription regulation") %>%
  filter(category != "unknown function") %>%
  filter(category != "DNA, RNA and nucleotide metabolism")

# Group and count phage-related proteins
phrogsPerContig <- pholdVirPassed_filt %>%
  group_by(contig_id, category) %>%
  dplyr::count(contig_id, category) %>%
  dcast(contig_id ~ category, value.var = "n")

phrogsPerContig[is.na(phrogsPerContig)] <- 0
colnames(phrogsPerContig) <- c("contig_id", "Connector", "Head_Packaging", "Integration_Excision", "Lysis", "Tail")

# Add functional diversity
phrogsPerContig <- phrogsPerContig %>%
  rowwise() %>%
  mutate(
    Functionaldiversity = sum(c_across(-contig_id) > 0)
  ) %>%
  ungroup()

# Read and process jager data
jager <- vroom(jager_file, delim = "\t", col_names = TRUE)
phages_jager <- jager %>%
  filter(prediction == "Phage") %>%
  filter(realiability_score >= 0.6)

# Read and process genomad data
genomad <- vroom(genomad_file, delim = "\t", col_names = TRUE)
colnames(genomad) <- c("contig_id", "genomad_Length", "genomad_topology", 
                       "genomad_coordinates", "genomad_n_genes", "genomad_genetic_code", 
                       "genomad_virus_score", "genomad_fdr", "genomad_n_hallmarks","genomad_marker_enrichment", "genomad_taxonomy")

# Read and process checkv data
checkv <- vroom(checkv_file, delim = "\t", col_names = TRUE)

# Create phagePredictedContigs table
phagePredictedContigs <- checkv[, c("contig_id", "contig_length", "completeness")] %>%
  left_join(phrogsPerContig[, c("contig_id", "Functionaldiversity")], by = "contig_id") %>%
  left_join(phages_jager[, c("contig_id", "realiability_score")], by = "contig_id") %>%
  left_join(genomad[, c("contig_id", "genomad_topology", "genomad_virus_score")], by = "contig_id") %>%
  mutate(
    Functionaldiversity = ifelse(is.na(Functionaldiversity), 0, Functionaldiversity),
    is_phage = case_when(
      Functionaldiversity >= 3 ~ TRUE,
      Functionaldiversity < 3 & genomad_topology %in% c("DTR", "ITR", "Provirus") ~ TRUE,
      Functionaldiversity >= 1 & Functionaldiversity < 3 & (!is.na(realiability_score) | !is.na(genomad_virus_score)) ~ TRUE,
      Functionaldiversity < 1 & !is.na(realiability_score) & !is.na(genomad_virus_score) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(is_phage == TRUE) %>%
  select(-is_phage)

# Extract contig_id vector
contig_ids <- phagePredictedContigs$contig_id

# Save outputs
output_table_path <- file.path(output_dir, "phagePredictedContigs.tsv")
output_vector_path <- file.path(output_dir, "contig_ids.txt")

vroom::vroom_write(phagePredictedContigs, output_table_path, delim = "\t")
writeLines(contig_ids, output_vector_path)

cat("Outputs saved to:", output_dir, "\n")