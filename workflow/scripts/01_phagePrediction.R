#!/usr/bin/env Rscript
# Integrate phage prediction results from multiple tools

# Load necessary libraries
packages <- c("dplyr", "vroom", "reshape2", "tidyr", "logger")
options(repos = c(CRAN = "https://cran.r-project.org"))
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(reshape2)
  library(tidyr)
  library(logger)
})

# Set up logging
log_info("Starting phage prediction integration")

# Function to process inputs
process_phage_predictions <- function(phold_file, jager_file, genomad_file, checkv_file, output_table, output_ids) {
  log_info("Reading input files:")
  log_info(paste("phold file:", phold_file))
  log_info(paste("jager file:", jager_file))
  log_info(paste("genomad file:", genomad_file))
  log_info(paste("checkv file:", checkv_file))
  
  # Read and process pholdVirPassed data
  log_info("Processing phold data")
  pholdVirPassed <- vroom(phold_file, delim = "\t", col_names = TRUE)
  log_info(paste("PHOLD raw data rows:", nrow(pholdVirPassed)))
  colnames(pholdVirPassed)[7] <- "category"
  
  # Filter pholdVirPassed
  pholdVirPassed_filt <- pholdVirPassed %>%
    filter(product != "hypothetical protein") %>%
    filter(category != "other") %>%
    filter(category != "moron, auxiliary metabolic gene and host takeover") %>%
    filter(category != "transcription regulation") %>%
    filter(category != "unknown function") %>%
    filter(category != "DNA, RNA and nucleotide metabolism")
  
  log_info(paste("PHOLD filtered data rows:", nrow(pholdVirPassed_filt)))
  
  # Group and count phage-related proteins
  # Handle empty PHOLD results gracefully
  if (nrow(pholdVirPassed_filt) == 0) {
    log_info("No PHOLD predictions found - creating empty structure")
    phrogsPerContig <- data.frame(
      contig_id = character(0),
      Connector = numeric(0),
      Head_Packaging = numeric(0),
      Integration_Excision = numeric(0),
      Lysis = numeric(0),
      Tail = numeric(0)
    )
  } else {
    phrogsPerContig <- pholdVirPassed_filt %>%
      group_by(contig_id, category) %>%
      dplyr::count(contig_id, category) %>%
      dcast(contig_id ~ category, value.var = "n")
    
    phrogsPerContig[is.na(phrogsPerContig)] <- 0
    colnames(phrogsPerContig) <- c("contig_id", "Connector", "Head_Packaging", "Integration_Excision", "Lysis", "Tail")
  }
  
  # Add functional diversity
  if (nrow(phrogsPerContig) == 0) {
    phrogsPerContig$Functionaldiversity <- numeric(0)
  } else {
    phrogsPerContig <- phrogsPerContig %>%
      rowwise() %>%
      mutate(
        Functionaldiversity = sum(c_across(-contig_id) > 0)
      ) %>%
      ungroup()
  }
  
  # Read and process jager data
  log_info("Processing jaeger data")
  jager <- vroom(jager_file, delim = "\t", col_names = TRUE)
  phages_jager <- jager %>%
    filter(prediction == "Phage") %>%
    filter(realiability_score >= 0.6)
  
  # Read and process genomad data
  log_info("Processing genomad data")
  genomad <- vroom(genomad_file, delim = "\t", col_names = TRUE)
  colnames(genomad) <- c("contig_id", "genomad_Length", "genomad_topology", 
                         "genomad_coordinates", "genomad_n_genes", "genomad_genetic_code", 
                         "genomad_virus_score", "genomad_fdr", "genomad_n_hallmarks","genomad_marker_enrichment", "genomad_taxonomy")
  
  # Read and process checkv data
  log_info("Processing checkv data")
  checkv <- vroom(checkv_file, delim = "\t", col_names = TRUE)
  
  # Create phagePredictedContigs table
  log_info("Integrating results to identify phage contigs")
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
  log_info(paste("Found", length(contig_ids), "phage contigs"))
  
  # Save outputs
  log_info(paste("Writing outputs to:", output_table, "and", output_ids))
  vroom::vroom_write(phagePredictedContigs, output_table, delim = "\t")
  writeLines(contig_ids, output_ids)
  
  log_info("Phage prediction integration completed")
}

# Handle snakemake or command line execution
if (exists("snakemake")) {
  log_appender(appender_file(snakemake@log[[1]]))
  
  process_phage_predictions(
    snakemake@input$phold,
    snakemake@input$jaeger,
    snakemake@input$genomad,
    snakemake@input$checkv,
    snakemake@output$predictions,
    snakemake@output$contig_ids
  )
} else {
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
  
  # Create output paths
  output_table_path <- file.path(output_dir, "phagePredictedContigs.tsv")
  output_vector_path <- file.path(output_dir, "contig_ids.txt")
  
  # Process inputs
  process_phage_predictions(
    phold_file,
    jager_file,
    genomad_file, 
    checkv_file,
    output_table_path,
    output_vector_path
  )
}