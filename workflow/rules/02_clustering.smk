"""
Rules for clustering phage contigs into vOTUs.
This step is optional and can be skipped based on config settings.
"""

# Determine if clustering should be skipped
def skip_clustering(wildcards):
    return not config.get("do_clustering", True)

# Define the phage input file - handles different entry points and clustering options
def get_phage_input(wildcards):
    # If starting from clustered sequences, always use that file (via symlink)
    if config.get("start_from") == "clustered_sequences":
        return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    # Otherwise, use clustering output if enabled, or prediction output if not
    elif config.get("do_clustering", True):
        return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    else:
        return f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"

# 1. Cluster phage contigs with vclust
rule cluster_phages:
    input:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
    output:
        clustering_dir = directory(f"{config['output_dir']}/02_clustering/clustering_results"),
        clusters = f"{config['output_dir']}/02_clustering/clusters.tsv",
        rep_seqs_list = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.tsv"
    log:
        f"{config['output_dir']}/logs/cluster_phages.log"
    conda:
        config["conda_envs"]["vclust"]
    threads: 24
    shell:
        """
        # Create output directory
        mkdir -p {output.clustering_dir}
        
        # Check input sequences
        NUM_SEQS=$(grep -c "^>" {input.phage_contigs} || true)
        echo "Number of input sequences: $NUM_SEQS" > {log} 2>&1
        
        if [ "$NUM_SEQS" -eq 0 ]; then
            echo "WARNING: No input sequences. Creating empty output files." >> {log} 2>&1
            echo -e "sequence_id\tcluster_id" > {output.clusters}
            touch {output.rep_seqs_list}
            exit 0
        fi
        
        # Filter genome pairs with at least 20 common 25-mers and 95% identity
        echo "Running vclust prefilter..." >> {log} 2>&1
        vclust prefilter -i {input.phage_contigs} \
            -o {output.clustering_dir}/vclust_fltr.txt \
            --min-kmers 20 \
            --min-ident 0.95 >> {log} 2>&1
            
        # Check if prefilter produced output
        if [ ! -s {output.clustering_dir}/vclust_fltr.txt ]; then
            echo "WARNING: vclust prefilter produced no output. Creating single-sequence clusters." >> {log} 2>&1
            # Create a clusters file where each sequence is its own cluster
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1"\t"$1}}' > {output.clusters}
            # For single-sequence clusters, all sequences are representatives
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1}}' > {output.rep_seqs_list}
            exit 0
        fi
            
        # Align genome pairs filtered by the prefiltering command
        echo "Running vclust align..." >> {log} 2>&1
        vclust align -i {input.phage_contigs} \
            -o {output.clustering_dir}/vclust_ani.tsv \
            --filter {output.clustering_dir}/vclust_fltr.txt >> {log} 2>&1
            
        # Check if alignment produced output
        if [ ! -s {output.clustering_dir}/vclust_ani.tsv ]; then
            echo "WARNING: vclust align produced no output. Creating single-sequence clusters." >> {log} 2>&1
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1"\t"$1}}' > {output.clusters}
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1}}' > {output.rep_seqs_list}
            exit 0
        fi
            
        # Cluster genomes based on ANI similarity measure
        echo "Running vclust cluster..." >> {log} 2>&1
        vclust cluster -i {output.clustering_dir}/vclust_ani.tsv \
            -o {output.clustering_dir}/vclust_clusters.tsv \
            --algorithm leiden \
            --metric ani \
            --ids {output.clustering_dir}/vclust_ani.ids.tsv \
            --ani {config[params][vclust][identity]} \
            --qcov {config[params][vclust][coverage]} \
            --rcov {config[params][vclust][coverage]} >> {log} 2>&1 || {{
                echo "WARNING: vclust cluster failed or produced no clusters." >> {log} 2>&1
                touch {output.clustering_dir}/vclust_clusters.tsv
            }}
            
        # Check if clustering produced any clusters
        if [ ! -s {output.clustering_dir}/vclust_clusters.tsv ]; then
            echo "WARNING: No clusters formed. Creating single-sequence clusters for all input sequences." >> {log} 2>&1
            # Create a clusters file where each sequence is its own cluster
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1"\t"$1}}' > {output.clusters}
            # All sequences are representatives
            grep "^>" {input.phage_contigs} | sed 's/>//g' | awk '{{print $1}}' > {output.rep_seqs_list}
        else
            # Extract lengths of the sequences
            echo "Extracting sequence lengths..." >> {log} 2>&1
            seqkit fx2tab --length --name {input.phage_contigs} > {output.clustering_dir}/seq_lengths.tsv 2>> {log}
            
            # Sort the clusters and lengths data files before joining
            sort -k1,1 {output.clustering_dir}/vclust_clusters.tsv > {output.clustering_dir}/vclust_sorted_clusters.tsv
            sort -k1,1 {output.clustering_dir}/seq_lengths.tsv > {output.clustering_dir}/sorted_seq_lengths.tsv
            
            # Join the sorted files
            join -1 1 -2 1 -t $'\t' {output.clustering_dir}/vclust_sorted_clusters.tsv {output.clustering_dir}/sorted_seq_lengths.tsv > {output.clustering_dir}/joined_Clusters_length.tsv
            
            # Use awk to find the longest sequence in each cluster
            echo "Selecting longest sequence per cluster as representative..." >> {log} 2>&1
            awk -F'\t' '{{
                if (!($2 in max_len) || max_len[$2] < $3) {{
                    max_len[$2] = $3;
                    max_seq[$2] = $1;
                }}
            }} END {{
                for (c in max_seq)
                    print max_seq[c];
            }}' {output.clustering_dir}/joined_Clusters_length.tsv > {output.rep_seqs_list}
            
            # Copy clusters file to expected output location
            cp {output.clustering_dir}/vclust_clusters.tsv {output.clusters}
        fi
        
        echo "Clustering complete. Found $(wc -l < {output.rep_seqs_list}) vOTU representatives." >> {log} 2>&1
        """

# 2. Extract vOTU representative sequences
rule extract_votu_representatives:
    input:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta",
        rep_seqs_list = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.tsv"
    output:
        rep_seqs = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_votu_representatives.log"
    conda:
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Check if representative sequences list has data
        if [ ! -s {input.rep_seqs_list} ]; then
            echo "WARNING: No representative sequences found. Using all phage contigs as representatives." > {log} 2>&1
            cp {input.phage_contigs} {output.rep_seqs}
            exit 0
        fi
        
        # Extract representative sequences using the pre-calculated list
        echo "Extracting $(wc -l < {input.rep_seqs_list}) vOTU representative sequences..." > {log} 2>&1
        seqkit grep -f {input.rep_seqs_list} {input.phage_contigs} > {output.rep_seqs} 2>> {log}
        
        # Check if extraction was successful
        if [ ! -s {output.rep_seqs} ]; then
            echo "WARNING: No sequences extracted. Using all phage contigs as representatives." >> {log} 2>&1
            cp {input.phage_contigs} {output.rep_seqs}
        fi
        
        echo "Extracted $(grep -c '^>' {output.rep_seqs} || echo 0) sequences." >> {log} 2>&1
        """
