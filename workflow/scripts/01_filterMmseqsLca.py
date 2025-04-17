#!/usr/bin/env python3
# Filter mmseqs2 taxonomy results to select viral contigs

import pandas as pd
import argparse
from Bio import SeqIO
import os
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("filter_mmseqs")

def safe_first_taxid(lineage):
    """Safely extract the first taxid from the lineage string."""
    if pd.isnull(lineage):
        return None
    parts = lineage.split(";")
    try:
        return int(parts[0])
    except (IndexError, ValueError):
        logger.warning(f"Invalid lineage entry: {lineage}")
        return None

def filter_mmseqs_output(input_file, output_file, ids_output_file, fasta_file, missing_contigs_output_file):
    """
    Filter mmseqs2 taxonomy results to select viral contigs.
    
    Args:
        input_file: Path to the mmseqs2 LCA table
        output_file: Path for the filtered output table
        ids_output_file: Path for the filtered contig IDs
        fasta_file: Path to the contigs FASTA file
        missing_contigs_output_file: Path for missing contig IDs
    """
    logger.info(f"Loading mmseqs2 table from {input_file}")
    
    # Load the table
    df = pd.read_csv(input_file, sep="\t", header=None, names=[
        "contig_id", "tax_id", "rank", "label", "fragments", "assigned_fragments",
        "label_fragments", "support", "lineage"
    ])
    
    logger.info(f"Loaded {len(df)} entries from mmseqs2 table")

    # Condition 1: First taxID in lineage is 10239 (virus)
    df["first_taxid_in_lineage"] = df["lineage"].apply(safe_first_taxid)
    condition_virus = df["first_taxid_in_lineage"] == 10239

    # Condition 2: taxID is 1
    condition_taxid_one = df["tax_id"] == 1

    # Condition 3: taxID is 0
    condition_taxid_zero = df["tax_id"] == 0

    # Condition 4: Entries not in the table (from FASTA file)
    logger.info(f"Reading FASTA file: {fasta_file}")
    contig_ids_in_table = set(df["contig_id"])
    contig_ids_in_fasta = {record.id for record in SeqIO.parse(fasta_file, "fasta")}
    missing_contigs = contig_ids_in_fasta - contig_ids_in_table
    
    logger.info(f"Found {len(missing_contigs)} contigs in FASTA but not in mmseqs2 table")

    # Save missing contigs to a file
    os.makedirs(os.path.dirname(missing_contigs_output_file), exist_ok=True)
    with open(missing_contigs_output_file, "w") as mc_file:
        for contig_id in sorted(missing_contigs):
            mc_file.write(f"{contig_id}\n")

    # Apply conditions and filter the dataframe
    filtered_df = df[condition_virus | condition_taxid_one | condition_taxid_zero]
    logger.info(f"Filtered to {len(filtered_df)} entries after applying conditions")
    
    passing_contig_ids = set(filtered_df["contig_id"]).union(missing_contigs)
    logger.info(f"Total {len(passing_contig_ids)} contigs passing all filters")

    # Save filtered table
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    filtered_df.to_csv(output_file, sep="\t", index=False, header=False)
    
    # Save passing contig IDs
    os.makedirs(os.path.dirname(ids_output_file), exist_ok=True)
    with open(ids_output_file, "w") as ids_file:
        for contig_id in sorted(passing_contig_ids):
            ids_file.write(f"{contig_id}\n")
    
    logger.info("Filtering completed successfully")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter mmseqs2 easy-taxonomy output based on specific criteria.")
    parser.add_argument("--mmseqs_LCA_table", required=True, help="Path to the input mmseqs2 output table.")
    parser.add_argument("--o_filtered_LCA_table", required=True, help="Path to the output filtered table.")
    parser.add_argument("--o_passing_contig_ids", required=True, help="Path to the output file for filtered contig IDs.")
    parser.add_argument("--contigs", required=True, help="Path to the FASTA file of contigs.")
    parser.add_argument("--o_missing_contig_ids", required=True, help="Path to the output file for missing contig IDs.")

    args = parser.parse_args()
    
    # Get inputs from snakemake
    if 'snakemake' in globals():
        snakemake_logger = logging.getLogger("snakemake.filter_mmseqs")
        logger.handlers = snakemake_logger.handlers
        logger.setLevel(snakemake_logger.level)
        
        filter_mmseqs_output(
            snakemake.input.lca_table,
            snakemake.output.filtered_lca,
            snakemake.output.passing_ids,
            snakemake.input.contigs,
            snakemake.output.missing_ids
        )
    else:
        # Command line execution
        filter_mmseqs_output(
            args.mmseqs_LCA_table,
            args.o_filtered_LCA_table,
            args.o_passing_contig_ids,
            args.contigs,
            args.o_missing_contig_ids
        )
