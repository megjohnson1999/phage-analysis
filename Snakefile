"""
Phage Analysis Pipeline

A Snakemake workflow for phage prediction, clustering, and characterization.
"""

configfile: "config/config.yaml"

# Import rules
include: "workflow/rules/01_prediction.smk"
include: "workflow/rules/02_clustering.smk"
include: "workflow/rules/03_analysis.smk"

# Default target rule
rule all:
    input:
        # Final output files depending on whether clustering is enabled
        lambda wildcards: 
            expand(
                "{outdir}/03_iphop_results/iphop_predictions_compiled.tsv",
                outdir=config["output_dir"]
            ) +
            expand(
                "{outdir}/03_phacts_results/phacts_predictions_compiled.tsv",
                outdir=config["output_dir"]
            ) +
            expand(
                "{outdir}/03_genomic_info/mmseqs_taxonomy.tsv",
                outdir=config["output_dir"]
            ) +
            (["{outdir}/02_clustering/vOTU_repSeqs.fasta".format(
                outdir=config["output_dir"]
             )] if config.get("do_clustering", True) else [])