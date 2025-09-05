#!/usr/bin/env Rscript
# Taxonomic Consensus Script
# 
# This script replicates the taxonomic consensus logic from the original R script,
# using taxonomizr for mmseqs2 results and combining with other tools exactly as
# in the original workflow.

library(tidyverse)
library(taxonomizr)
library(vroom)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Usage: Rscript taxonomic_consensus.R <mmseqs_file> <phabox_taxonomy> <phabox_lifestyle> <vc3_dir> <taxonomizr_db> <output_taxonomy> <output_summary>")
}

mmseqs_file <- args[1]
phabox_taxonomy_file <- args[2]
phabox_lifestyle_file <- args[3]
vc3_dir <- args[4]
taxonomizr_db <- args[5]
output_taxonomy <- args[6]
output_summary <- args[7]

cat("Starting taxonomic consensus analysis...\n")
cat("MMseqs file:", mmseqs_file, "\n")
cat("Phabox taxonomy:", phabox_taxonomy_file, "\n")
cat("Phabox lifestyle:", phabox_lifestyle_file, "\n")
cat("vContact3 dir:", vc3_dir, "\n")
cat("Taxonomizr database:", taxonomizr_db, "\n")

# Helper function from original script
`%notin%` <- function(x,y) !(x %in% y)

# ===== MMSeqs + taxonomizr processing (replicating original R script) =====
cat("Processing MMSeqs results with taxonomizr...\n")

# Load mmseqs2 results
mmseqs_topHit <- vroom(file = mmseqs_file, col_names = TRUE, show_col_types = FALSE)

if (nrow(mmseqs_topHit) == 0) {
    cat("Warning: No mmseqs2 results found\n")
    final_mmseqs <- data.frame(
        contigID = character(0),
        blast_superkingdom = character(0),
        blast_phylum = character(0),
        blast_class = character(0),
        blast_order = character(0),
        blast_family = character(0),
        blast_genus = character(0),
        blast_species = character(0)
    )
} else {
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
        final_mmseqs <- data.frame(
            contigID = character(0),
            blast_superkingdom = character(0),
            blast_phylum = character(0),
            blast_class = character(0),
            blast_order = character(0),
            blast_family = character(0),
            blast_genus = character(0),
            blast_species = character(0)
        )
    } else {
        # Use taxonomizr exactly as in the original script
        cat("Looking up taxonomic lineages using taxonomizr...\n")
        linage <- as.data.frame(getTaxonomy(mmseqs_best_hits_by_bitscore$taxid, taxonomizr_db))
        
        # Add taxonomy information exactly as in the original script
        mmseqs_best_hits_by_bitscore$blast_superkingdom <- linage$superkingdom
        mmseqs_best_hits_by_bitscore$blast_phylum  <- linage$phylum
        mmseqs_best_hits_by_bitscore$blast_class  <- linage$class
        mmseqs_best_hits_by_bitscore$blast_order  <- linage$order
        mmseqs_best_hits_by_bitscore$blast_family  <- linage$family
        mmseqs_best_hits_by_bitscore$blast_genus  <- linage$genus
        mmseqs_best_hits_by_bitscore$blast_species  <- linage$species
        
        final_mmseqs <- mmseqs_best_hits_by_bitscore %>%
            select(-taxid) %>%
            rename(contigID = query)
        
        # Adding superkingdom to the unclassified (matching original script logic)
        final_mmseqs <- final_mmseqs %>%
            mutate(blast_superkingdom = ifelse(is.na(blast_superkingdom), "Viruses", blast_superkingdom))
    }
}

cat("MMSeqs taxonomizr lookup complete:", nrow(final_mmseqs), "contigs\n")

# ===== Load Phabox results =====
cat("Loading Phabox2 results...\n")

# Check for new Phabox2 format first
phabox_summary_file <- file.path(dirname(phabox_taxonomy_file), "final_prediction", "final_prediction_summary.tsv")

if (file.exists(phabox_summary_file)) {
    cat("Found Phabox2 summary file format\n")
    phabox_df <- vroom(phabox_summary_file, show_col_types = FALSE)
    
    # Filter for viral predictions only
    phabox_viral <- phabox_df %>% filter(Pred == "virus")
    
    if (nrow(phabox_viral) == 0) {
        final_phagcn <- data.frame(
            contigID = character(0),
            phagcn_superkingdom = character(0),
            phagcn_phylum = character(0),
            phagcn_class = character(0),
            phagcn_order = character(0),
            phagcn_family = character(0),
            phagcn_genus = character(0)
        )
    } else {
        # Parse Phabox2 lineage (similar to original PhaGCN parsing)
        extract_taxonomy <- function(tax_string) {
          parts <- str_split(tax_string, ";", simplify = TRUE)
          levels <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
          taxonomy <- setNames(vector("list", length(levels)), levels)
          taxonomy[] <- lapply(taxonomy, function(x) NA)
          
          for (part in parts) {
            if (str_detect(part, ":")) {
              split_part <- str_split(part, ":", simplify = TRUE)
              if (length(split_part) >= 2) {
                key <- tolower(split_part[1])
                value <- split_part[2]
                if (key %in% names(taxonomy)) {
                  taxonomy[[key]] <- value
                }
              }
            }
          }
          
          return(taxonomy)
        }
        
        final_phagcn <- phabox_viral %>%
          mutate(
            TaxonomyData = lapply(Lineage, extract_taxonomy),
            phagcn_superkingdom = sapply(TaxonomyData, `[[`, "superkingdom"),
            phagcn_phylum = sapply(TaxonomyData, `[[`, "phylum"),
            phagcn_class = sapply(TaxonomyData, `[[`, "class"),
            phagcn_order = sapply(TaxonomyData, `[[`, "order"),
            phagcn_family = sapply(TaxonomyData, `[[`, "family"),
            phagcn_genus = sapply(TaxonomyData, `[[`, "genus"),
            phagcn_species = sapply(TaxonomyData, `[[`, "species")
          ) %>%
          select(Accession, phagcn_superkingdom, phagcn_phylum, phagcn_class, phagcn_order, phagcn_family, phagcn_genus) %>%
          rename(contigID = Accession)
    }
} else {
    cat("Using legacy Phabox2 format (if files exist)\n")
    # Create empty dataframe for consistency
    final_phagcn <- data.frame(
        contigID = character(0),
        phagcn_superkingdom = character(0),
        phagcn_phylum = character(0),
        phagcn_class = character(0),
        phagcn_order = character(0),
        phagcn_family = character(0),
        phagcn_genus = character(0)
    )
}

cat("Phabox2 results loaded:", nrow(final_phagcn), "contigs\n")

# ===== Load vContact3 results (replicating original VC3 processing) =====
cat("Loading vContact3 results...\n")

vc3_file <- file.path(vc3_dir, "final_assignments.csv")
if (!file.exists(vc3_file)) {
    vc3_file <- file.path(vc3_dir, "exports", "final_assignments.csv")
}

if (!file.exists(vc3_file)) {
    cat("Warning: vContact3 results not found\n")
    final_vc3 <- data.frame(
        contigID = character(0),
        vc3_superkingdom = character(0),
        vc3_phylum = character(0),
        vc3_class = character(0),
        vc3_order = character(0),
        vc3_family = character(0),
        vc3_subfamily = character(0),
        vc3_genus = character(0)
    )
} else {
    vc3 <- vroom(file = vc3_file, show_col_types = FALSE)
    
    # Replicate original VC3 processing exactly
    vc3 <- vc3 %>%
      filter(grepl("contig", GenomeName) == TRUE |
               grepl("edge", GenomeName) == TRUE | 
               grepl("virus_comp", GenomeName) == TRUE) %>%
      select(GenomeName, `realm (prediction)`,
             `phylum (prediction)`, `class (prediction)`, `order (prediction)`,
             `family (prediction)`, `subfamily (prediction)`, `genus (prediction)`) 
    
    colnames(vc3) <- c("contigID", "vc3_realm", "vc3_phylum", "vc3_class", "vc3_order", "vc3_family", "vc3_subfamily", "vc3_genus")
    
    final_vc3 <- vc3 %>%
      filter(vc3_realm != "No Realm") %>%
      filter(vc3_realm != "No prediction") %>%
      mutate(across(everything(), ~ifelse(grepl("novel", .) == TRUE, NA, .))) %>%
      mutate(vc3_superkingdom = "Viruses",
             vc3_phylum = ifelse(grepl("Phixviricota", vc3_phylum) == TRUE,
                                 "Phixviricota", vc3_phylum),
             vc3_class = ifelse(grepl("Malgrandaviricetes", vc3_class) == TRUE,
                                 "Malgrandaviricetes", vc3_class)) %>%
      select(contigID, vc3_superkingdom, vc3_phylum, vc3_class, vc3_order, vc3_family, vc3_subfamily, vc3_genus)
}

cat("vContact3 results loaded:", nrow(final_vc3), "contigs\n")

# ===== Taxonomic merging and consensus (replicating original R script) =====
cat("Creating taxonomic consensus...\n")

# Merge all datasets by contigID (exactly as in original script)
merged_taxonomies <- final_mmseqs %>%
  full_join(final_phagcn, by = "contigID") %>%
  full_join(final_vc3, by = "contigID")

if (nrow(merged_taxonomies) == 0) {
    cat("Warning: No taxonomic data to process\n")
    consensus_taxonomy <- data.frame(
        contigID = character(0),
        superkingdom = character(0),
        phylum = character(0),
        class = character(0),
        order = character(0),
        family = character(0),
        genus = character(0),
        species = character(0)
    )
} else {
    # Apply the exact consensus logic from the original R script
    consensus_taxonomy <- merged_taxonomies %>%
      # Step 1: Determine the consensus superkingdom
      mutate(superkingdom = coalesce(blast_superkingdom, 
                                     phagcn_superkingdom, 
                                     vc3_superkingdom)) %>%
      mutate(superkingdom = ifelse(is.na(superkingdom), "Viruses", superkingdom)) %>%
      # Step 2: Define NA-aware hierarchical checks exactly as in original
      mutate(
        # PHYLUM
        phylum = coalesce(
          ifelse(
            (blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom)),
            blast_phylum, NA
          ),
          ifelse(
            (phagcn_superkingdom == superkingdom) | (is.na(phagcn_superkingdom) & is.na(superkingdom)),
            phagcn_phylum, NA
          ),
          ifelse(
            (vc3_superkingdom == superkingdom) | (is.na(vc3_superkingdom) & is.na(superkingdom)),
            vc3_phylum, NA
          )
        ),
        # CLASS
        class = coalesce(
          ifelse(
            ((blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom))) &
            ((blast_phylum == phylum) | (is.na(blast_phylum) & is.na(phylum))),
            blast_class, NA
          ),
          ifelse(
            ((phagcn_superkingdom == superkingdom) | (is.na(phagcn_superkingdom) & is.na(superkingdom))) &
            ((phagcn_phylum == phylum) | (is.na(phagcn_phylum) & is.na(phylum))),
            phagcn_class, NA
          ),
          ifelse(
            ((vc3_superkingdom == superkingdom) | (is.na(vc3_superkingdom) & is.na(superkingdom))) &
            ((vc3_phylum == phylum) | (is.na(vc3_phylum) & is.na(phylum))),
            vc3_class, NA
          )
        ),
        # ORDER
        order = coalesce(
          ifelse(
            ((blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom))) &
            ((blast_phylum == phylum) | (is.na(blast_phylum) & is.na(phylum))) &
            ((blast_class == class) | (is.na(blast_class) & is.na(class))),
            blast_order, NA
          ),
          ifelse(
            ((phagcn_superkingdom == superkingdom) | (is.na(phagcn_superkingdom) & is.na(superkingdom))) &
            ((phagcn_phylum == phylum) | (is.na(phagcn_phylum) & is.na(phylum))) &
            ((phagcn_class == class) | (is.na(phagcn_class) & is.na(class))),
            phagcn_order, NA
          ),
          ifelse(
            ((vc3_superkingdom == superkingdom) | (is.na(vc3_superkingdom) & is.na(superkingdom))) &
            ((vc3_phylum == phylum) | (is.na(vc3_phylum) & is.na(phylum))) &
            ((vc3_class == class) | (is.na(vc3_class) & is.na(class))),
            vc3_order, NA
          )
        ),
        # FAMILY
        family = coalesce(
          ifelse(
            ((blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom))) &
            ((blast_phylum == phylum) | (is.na(blast_phylum) & is.na(phylum))) &
            ((blast_class == class) | (is.na(blast_class) & is.na(class))) &
            ((blast_order == order) | (is.na(blast_order) & is.na(order))),
            blast_family, NA
          ),
          ifelse(
            ((phagcn_superkingdom == superkingdom) | (is.na(phagcn_superkingdom) & is.na(superkingdom))) &
            ((phagcn_phylum == phylum) | (is.na(phagcn_phylum) & is.na(phylum))) &
            ((phagcn_class == class) | (is.na(phagcn_class) & is.na(class))) &
            ((phagcn_order == order) | (is.na(phagcn_order) & is.na(order))),
            phagcn_family, NA
          ),
          ifelse(
            ((vc3_superkingdom == superkingdom) | (is.na(vc3_superkingdom) & is.na(superkingdom))) &
            ((vc3_phylum == phylum) | (is.na(vc3_phylum) & is.na(phylum))) &
            ((vc3_class == class) | (is.na(vc3_class) & is.na(class))) &
            ((vc3_order == order) | (is.na(vc3_order) & is.na(order))),
            vc3_family, NA
          )
        ),
        # GENUS
        genus = coalesce(
          ifelse(
            ((blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom))) &
            ((blast_phylum == phylum) | (is.na(blast_phylum) & is.na(phylum))) &
            ((blast_class == class) | (is.na(blast_class) & is.na(class))) &
            ((blast_order == order) | (is.na(blast_order) & is.na(order))) &
            ((blast_family == family) | (is.na(blast_family) & is.na(family))),
            blast_genus, NA
          ),
          ifelse(
            ((phagcn_superkingdom == superkingdom) | (is.na(phagcn_superkingdom) & is.na(superkingdom))) &
            ((phagcn_phylum == phylum) | (is.na(phagcn_phylum) & is.na(phylum))) &
            ((phagcn_class == class) | (is.na(phagcn_class) & is.na(class))) &
            ((phagcn_order == order) | (is.na(phagcn_order) & is.na(order))) &
            ((phagcn_family == family) | (is.na(phagcn_family) & is.na(family))),
            phagcn_genus, NA
          ),
          ifelse(
            ((vc3_superkingdom == superkingdom) | (is.na(vc3_superkingdom) & is.na(superkingdom))) &
            ((vc3_phylum == phylum) | (is.na(vc3_phylum) & is.na(phylum))) &
            ((vc3_class == class) | (is.na(vc3_class) & is.na(class))) &
            ((vc3_order == order) | (is.na(vc3_order) & is.na(order))) &
            ((vc3_family == family) | (is.na(vc3_family) & is.na(family))),
            vc3_genus, NA
          )
        ),
        # SPECIES (only blast in this simplified version)
        species = ifelse(
          ((blast_superkingdom == superkingdom) | (is.na(blast_superkingdom) & is.na(superkingdom))) &
          ((blast_phylum == phylum) | (is.na(blast_phylum) & is.na(phylum))) &
          ((blast_class == class) | (is.na(blast_class) & is.na(class))) &
          ((blast_order == order) | (is.na(blast_order) & is.na(order))) &
          ((blast_family == family) | (is.na(blast_family) & is.na(family))) &
          ((blast_genus == genus) | (is.na(blast_genus) & is.na(genus))),
          blast_species, NA
        )
      ) %>%
      # Final selection of columns
      select(contigID, superkingdom, phylum, class, order, family, genus, species)
}

cat("Consensus taxonomy created for", nrow(consensus_taxonomy), "contigs\n")

# Save consensus taxonomy
vroom_write(consensus_taxonomy, output_taxonomy)
cat("Saved consensus taxonomy to:", output_taxonomy, "\n")

# Create summary statistics (simple version)
summary_stats <- list(
    total_contigs = nrow(consensus_taxonomy),
    tool_contributions = list(
        mmseqs2 = nrow(final_mmseqs),
        phabox2 = nrow(final_phagcn),
        vcontact3 = nrow(final_vc3)
    ),
    taxonomy_coverage = list()
)

# Calculate coverage for each taxonomic level
for (level in c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) {
    if (nrow(consensus_taxonomy) > 0) {
        non_null <- sum(!is.na(consensus_taxonomy[[level]]))
        summary_stats$taxonomy_coverage[[level]] <- list(
            count = non_null,
            percentage = round(100 * non_null / nrow(consensus_taxonomy), 2)
        )
    } else {
        summary_stats$taxonomy_coverage[[level]] <- list(count = 0, percentage = 0)
    }
}

# Save summary as JSON
library(jsonlite)
write_json(summary_stats, output_summary, pretty = TRUE)
cat("Saved consensus summary to:", output_summary, "\n")

cat("Taxonomic consensus completed successfully!\n")