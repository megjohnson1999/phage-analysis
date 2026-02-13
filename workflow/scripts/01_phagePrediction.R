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
process_phage_predictions <- function(phold_file, jager_file, genomad_file, checkv_file,
                                       output_table, output_ids, output_viral_summary,
                                       mmseqs_lca_file = NULL) {
  log_info("Reading input files:")
  log_info(paste("phold file:", phold_file))
  log_info(paste("jager file:", jager_file))
  log_info(paste("genomad file:", genomad_file))
  log_info(paste("checkv file:", checkv_file))
  if (!is.null(mmseqs_lca_file)) {
    log_info(paste("mmseqs_lca file:", mmseqs_lca_file))
  }

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

  # Read and process jager data (keep ALL predictions, not just phages)
  log_info("Processing jaeger data")
  jager <- vroom(jager_file, delim = "\t", col_names = TRUE)

  # Create filtered version for phage calling (score >= 0.6 and prediction == "Phage")
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

  # Read mmseqs LCA data if provided
  mmseqs_lca <- NULL
  if (!is.null(mmseqs_lca_file) && file.exists(mmseqs_lca_file)) {
    log_info("Processing mmseqs LCA data")
    # The filtered_lca.tsv has no header - 10 columns
    mmseqs_lca <- vroom(mmseqs_lca_file, delim = "\t", col_names = FALSE)
    colnames(mmseqs_lca) <- c("contig_id", "taxid", "rank", "taxname",
                              "col5", "col6", "col7", "col8", "col9", "col10")
    mmseqs_lca <- mmseqs_lca %>% select(contig_id, taxid, taxname)
    colnames(mmseqs_lca) <- c("contig_id", "mmseqs_taxid", "mmseqs_lca")
    log_info(paste("MMseqs LCA data rows:", nrow(mmseqs_lca)))
  }

  # ============================================================================
  # Create COMPREHENSIVE viral contigs summary (ALL viral contigs)
  # ============================================================================
  log_info("Creating comprehensive viral contigs summary")

  # Start with CheckV as the base (has all viral contigs)
  viral_summary <- checkv %>%
    select(contig_id, contig_length, checkv_quality, completeness, contamination, provirus)

  # Add mmseqs LCA info if available
  if (!is.null(mmseqs_lca)) {
    viral_summary <- viral_summary %>%
      left_join(mmseqs_lca, by = "contig_id")
  }

  # Add PHOLD functional diversity
  viral_summary <- viral_summary %>%
    left_join(phrogsPerContig %>% select(contig_id, Functionaldiversity), by = "contig_id") %>%
    mutate(Functionaldiversity = ifelse(is.na(Functionaldiversity), 0, Functionaldiversity))

  # Add ALL Jaeger results (not just filtered phages)
  jaeger_summary <- jager %>%
    select(contig_id, prediction, realiability_score) %>%
    rename(jaeger_prediction = prediction, jaeger_score = realiability_score)

  viral_summary <- viral_summary %>%
    left_join(jaeger_summary, by = "contig_id")

  # Add genomad results
  genomad_summary <- genomad %>%
    select(contig_id, genomad_topology, genomad_virus_score, genomad_taxonomy)

  viral_summary <- viral_summary %>%
    left_join(genomad_summary, by = "contig_id")

  # Add is_phage flag and prediction_rule
  viral_summary <- viral_summary %>%
    mutate(
      # Check if Jaeger called it a phage with sufficient score
      jaeger_phage = !is.na(jaeger_score) & jaeger_prediction == "Phage" & jaeger_score >= 0.6,
      # Determine is_phage status
      is_phage = case_when(
        Functionaldiversity >= 3 ~ TRUE,
        Functionaldiversity < 3 & genomad_topology %in% c("DTR", "ITR", "Provirus") ~ TRUE,
        Functionaldiversity >= 1 & Functionaldiversity < 3 & (jaeger_phage | !is.na(genomad_virus_score)) ~ TRUE,
        Functionaldiversity < 1 & jaeger_phage & !is.na(genomad_virus_score) ~ TRUE,
        TRUE ~ FALSE
      ),
      # Create prediction rule explaining why it passed or failed
      prediction_rule = case_when(
        Functionaldiversity >= 3 ~ "high_functional_diversity",
        Functionaldiversity < 3 & genomad_topology %in% c("DTR", "ITR", "Provirus") ~ "special_topology",
        Functionaldiversity >= 1 & Functionaldiversity < 3 & (jaeger_phage | !is.na(genomad_virus_score)) ~ "moderate_diversity_with_tool_support",
        Functionaldiversity < 1 & jaeger_phage & !is.na(genomad_virus_score) ~ "dual_tool_agreement",
        # Exclusion reasons
        Functionaldiversity < 1 & !jaeger_phage & is.na(genomad_virus_score) ~ "low_diversity_no_tool_support",
        Functionaldiversity >= 1 & Functionaldiversity < 3 & !jaeger_phage & is.na(genomad_virus_score) ~ "moderate_diversity_no_tool_support",
        Functionaldiversity < 1 & (jaeger_phage | !is.na(genomad_virus_score)) ~ "low_diversity_single_tool_only",
        TRUE ~ "no_criteria_met"
      )
    ) %>%
    select(-jaeger_phage)  # Remove helper column

  # Reorder columns for clarity
  col_order <- c("contig_id", "contig_length", "is_phage", "prediction_rule",
                 "checkv_quality", "completeness", "contamination", "provirus",
                 "Functionaldiversity", "jaeger_prediction", "jaeger_score",
                 "genomad_topology", "genomad_virus_score", "genomad_taxonomy")

  # Add mmseqs columns if present
  if (!is.null(mmseqs_lca)) {
    col_order <- c(col_order[1:2], "mmseqs_taxid", "mmseqs_lca", col_order[3:length(col_order)])
  }

  # Select only columns that exist
  col_order <- col_order[col_order %in% colnames(viral_summary)]
  viral_summary <- viral_summary %>% select(all_of(col_order))

  log_info(paste("Total viral contigs:", nrow(viral_summary)))
  log_info(paste("Phage contigs:", sum(viral_summary$is_phage)))
  log_info(paste("Non-phage viral contigs:", sum(!viral_summary$is_phage)))

  # Write viral contigs summary
  log_info(paste("Writing viral contigs summary to:", output_viral_summary))
  vroom::vroom_write(viral_summary, output_viral_summary, delim = "\t")

  # ============================================================================
  # Create filtered phage predictions table (existing behavior)
  # ============================================================================
  log_info("Creating filtered phage predictions table")

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

  # Get optional mmseqs_lca input (may not be present in all entry points)
  mmseqs_lca_file <- NULL
  if ("mmseqs_lca" %in% names(snakemake@input)) {
    mmseqs_lca_file <- snakemake@input$mmseqs_lca
  }

  process_phage_predictions(
    snakemake@input$phold,
    snakemake@input$jaeger,
    snakemake@input$genomad,
    snakemake@input$checkv,
    snakemake@output$predictions,
    snakemake@output$contig_ids,
    snakemake@output$viral_summary,
    mmseqs_lca_file
  )
} else {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: Rscript script.R <phold_file> <jager_file> <genomad_file> <checkv_file> <output_dir> [mmseqs_lca_file]")
  }

  phold_file <- args[1]
  jager_file <- args[2]
  genomad_file <- args[3]
  checkv_file <- args[4]
  output_dir <- args[5]
  mmseqs_lca_file <- if (length(args) >= 6) args[6] else NULL

  # Create output paths
  output_table_path <- file.path(output_dir, "phagePredictedContigs.tsv")
  output_vector_path <- file.path(output_dir, "contig_ids.txt")
  output_viral_summary_path <- file.path(output_dir, "viral_contigs_summary.tsv")

  # Process inputs
  process_phage_predictions(
    phold_file,
    jager_file,
    genomad_file,
    checkv_file,
    output_table_path,
    output_vector_path,
    output_viral_summary_path,
    mmseqs_lca_file
  )
}
