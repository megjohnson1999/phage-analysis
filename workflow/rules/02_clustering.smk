"""
Rules for clustering phage contigs into vOTUs.
This step is optional and can be skipped based on config settings.
"""

# Determine if clustering should be skipped
def skip_clustering(wildcards):
    return not config.get("do_clustering", True)

# Define the phage input file - either from clustering or from prediction
def get_phage_input(wildcards):
    if config.get("do_clustering", True):
        return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    else:
        return f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"

# 1. Cluster phage contigs with vclust
rule cluster_phages:
    input:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
    output:
        clustering_dir = directory(f"{config['output_dir']}/02_clustering/clustering_results"),
        filtered_fasta = f"{config['output_dir']}/02_clustering/clustering_results/filtered.fasta",
        clusters = f"{config['output_dir']}/02_clustering/clusters.tsv"
    log:
        f"{config['output_dir']}/logs/cluster_phages.log"
    conda:
        config["conda_envs"]["vclust"]
    resources:
        mem_mb = config["resources"]["clustering"]["mem_mb"],
        threads = config["resources"]["clustering"]["threads"],
        time = config["resources"]["clustering"]["time"]
    shell:
        """
        # Filter contigs by length
        seqkit seq -m {config[resources][vclust][min_length]} {input.phage_contigs} > {output.filtered_fasta}
        
        # Run vclust for clustering
        vclust --in {output.filtered_fasta} \
            --out $(dirname {output.clusters}) \
            --id {config[resources][vclust][identity]} \
            --cov {config[resources][vclust][coverage]} \
            --threads {resources.threads} > {log} 2>&1
            
        # Move clusters file to expected output
        mv $(dirname {output.clusters})/clusters.tsv {output.clusters}
        """

# 2. Extract vOTU representative sequences
rule extract_votu_representatives:
    input:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta",
        clusters = f"{config['output_dir']}/02_clustering/clusters.tsv"
    output:
        rep_seqs = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_votu_representatives.log"
    conda:
        "bioconda::seqkit"
    shell:
        """
        # Extract representative sequence IDs from clusters file
        awk '{{print $1}}' {input.clusters} > temp_rep_seqs.txt
        
        # Extract representative sequences
        seqkit grep -f temp_rep_seqs.txt {input.phage_contigs} > {output.rep_seqs} 2> {log}
        
        # Clean up
        rm temp_rep_seqs.txt
        """