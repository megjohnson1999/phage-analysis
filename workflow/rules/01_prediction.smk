"""
Rules for phage prediction from metagenomic assemblies.
"""

# 1. Run Reneo for binning
rule reneo_binning:
    input:
        assembly = config["assembly_file"],
        reads_dir = config["reads_dir"]
    output:
        f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges.fasta"
    log:
        f"{config['output_dir']}/logs/reneo_binning.log"
    conda:
        config["conda_envs"]["reneo"]
    threads: 24
    shell:
        """
        mkdir -p {config[output_dir]}/01_reneo_output

        # Run Reneo for binning
        reneo run --input {input.assembly} \
            --reads {input.reads_dir} \
            --minlength 1000 \
            --output {config[output_dir]}/01_reneo_output \
            --threads {threads} > {log} 2>&1
        """

# 1b. Filter contigs by length (1KB)
rule contig_length_filter:
    input:
        f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges.fasta"
    output:
        f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta"
    log:
        f"{config['output_dir']}/logs/contig_length_filter.log"
    conda:
        config["conda_envs"]["seqkit"]
    threads: 8
    shell:
        """   
        seqkit seq --min-len 1000 -g \
            "{input}" > "{output}"
        """

# 1c. Filter assembly directly (if not using Reneo)
rule direct_contig_filter:
    input:
        assembly = config["assembly_file"]
    output:
        filtered_assembly = f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta"
    log:
        f"{config['output_dir']}/logs/direct_contig_filter.log"
    conda:
        config["conda_envs"]["seqkit"]
    threads: 8
    shell:
        """
        mkdir -p {config[output_dir]}/01_filtered_assembly
        seqkit seq --min-len 1000 -g \
            "{input.assembly}" > "{output.filtered_assembly}"
        """

# 2. Run mmseqs2 for taxonomy assignment
rule mmseqs_taxonomy:
    input:
        # Use the appropriate filtered contigs based on whether Reneo is enabled
        filtered_contigs = lambda wildcards: 
            f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if config.get("use_reneo", True) else 
            f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta"
    output:
        lca_table = f"{config['output_dir']}/01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv"
    log:
        f"{config['output_dir']}/logs/mmseqs_taxonomy.log"
    conda:
        config["conda_envs"]["mmseqs2"]
    threads: 24
    shell:
        """
        # Create temporary directory
        TMP_DIR=$(mktemp -d)
        
        # Run mmseqs2 for taxonomy assignment
        mmseqs easy-taxonomy {input.filtered_contigs} \
            {config[databases][mmseqs2][db]} \
            $(dirname {output.lca_table}) \
            $TMP_DIR \
            --min-length 30 \
            -e 1e-15 \
            --search-type 2 \
            -s 4.0 \
            --shuffle 0 \
            --lca-mode 2 \
            -a \
            --tax-lineage 2 \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
            --threads {threads} \
            --split-mode 0 \
            --orf-filter 1 \
            > {log} 2>&1
            
        # Clean up
        rm -rf $TMP_DIR
        """

# 3. Filter mmseqs2 results for viral contigs
rule filter_mmseqs_lca:
    input:
        lca_table = f"{config['output_dir']}/01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv",
        contigs = lambda wildcards: 
            f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if config.get("use_reneo", True) else 
            f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta"
    output:
        filtered_lca = f"{config['output_dir']}/01_filtered_mmseqs/filtered_lca.tsv",
        passing_ids = f"{config['output_dir']}/01_filtered_mmseqs/passing_contig_ids.txt",
        missing_ids = f"{config['output_dir']}/01_filtered_mmseqs/missing_contig_ids.txt"
    log:
        f"{config['output_dir']}/logs/filter_mmseqs_lca.log"
    conda:
        "python:3.9"
    script:
        "../scripts/01_filterMmseqsLca.py"

# 3b. Extract passing viral contigs
rule extract_viral_contigs:
    input:
        contigs = lambda wildcards: 
            f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if config.get("use_reneo", True) else 
            f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta",
        passing_ids = f"{config['output_dir']}/01_filtered_mmseqs/passing_contig_ids.txt"
    output:
        viral_contigs = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_viral_contigs.log"
    conda:
        "bioconda::seqkit"
    shell:
        """
        # Extract viral contigs from the filtered assembly
        seqkit grep -f {input.passing_ids} {input.contigs} > {output.viral_contigs} 2> {log}
        """

# 4. Run Jaeger for phage prediction
rule jaeger_prediction:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        results = directory(f"{config['output_dir']}/01_jaeger_output"),
        predictions = f"{config['output_dir']}/01_jaeger_output/final_predictions_scored.tsv"
    log:
        f"{config['output_dir']}/logs/jaeger_prediction.log"
    conda:
        config["conda_envs"]["jaeger"]
    shell:
        """
        Jaeger -i {input.assembly} -o {output.results} \
            -s 2.5 \
            --fsize 1000 \
            --stride 1000 > {log} 2>&1
        """

# 5. Run GeNomad for viral prediction
rule genomad_prediction:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        results = directory(f"{config['output_dir']}/01_genomad_output"),
        virus_summary = f"{config['output_dir']}/01_genomad_output/summary/virus_summary.tsv"
    log:
        f"{config['output_dir']}/logs/genomad_prediction.log"
    conda:
        config["conda_envs"]["genomad"]
    threads: 24
    shell:
        """
        genomad end-to-end --min-score 0.6 \
            --cleanup \
            --threads {threads} \
            {input.assembly} \
            {output.results} \
            {config[databases][genomad][db]} > {log} 2>&1
        """

# 6. Run Phold for protein annotation
rule phold_prediction:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        results = directory(f"{config['output_dir']}/01_phold_output"),
        predictions = f"{config['output_dir']}/01_phold_output/phold_per_cds_predictions.tsv"
    log:
        f"{config['output_dir']}/logs/phold_prediction.log"
    conda:
        config["conda_envs"]["phold"]
    threads: 24
    shell:
        """
        # Run phold
        phold run -i {input.assembly} \
            -o {output.results} \
            -d {config[databases][phold][db]} \
            -t {threads} --cpu --force > {log} 2>&1
        """

# 7. Run CheckV for quality assessment
rule checkv_assessment:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        results = directory(f"{config['output_dir']}/01_checkv_output"),
        quality_summary = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    log:
        f"{config['output_dir']}/logs/checkv_assessment.log"
    conda:
        config["conda_envs"]["checkv"]
    threads: 24
    shell:
        """
        # Run CheckV for viral quality assessment
        checkv end_to_end {input.assembly} \
            {output.results} \
            -d {config[databases][checkv][db]} \
            -t {threads} > {log} 2>&1
        """

# 8. Integrate phage prediction results
rule integrate_phage_predictions:
    input:
        phold = f"{config['output_dir']}/01_phold_output/phold_per_cds_predictions.tsv",
        jaeger = f"{config['output_dir']}/01_jaeger_output/final_predictions_scored.tsv",
        genomad = f"{config['output_dir']}/01_genomad_output/summary/virus_summary.tsv",
        checkv = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    output:
        predictions = f"{config['output_dir']}/01_phage_predictions/phagePredictedContigs.tsv",
        contig_ids = f"{config['output_dir']}/01_phage_predictions/contig_ids.txt"
    log:
        f"{config['output_dir']}/logs/integrate_phage_predictions.log"
    conda:
        config["conda_envs"]["r"]
    script:
        "../scripts/01_phagePrediction.R"

# 9. Extract phage contigs
rule extract_phage_contigs:
    input:
        contigs = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta",
        contig_ids = f"{config['output_dir']}/01_phage_predictions/contig_ids.txt"
    output:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_phage_contigs.log"
    conda:
        "bioconda::seqkit"
    shell:
        """
        # Extract phage contigs from viral contigs
        seqkit grep -f {input.contig_ids} {input.contigs} > {output.phage_contigs} 2> {log}
        """
