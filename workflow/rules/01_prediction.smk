"""
Rules for phage prediction from metagenomic assemblies.
"""

# 1. Run Reneo for binning
rule reneo_binning:
    input:
        assembly = config["assembly_file"],
        reads_dir = config["reads_dir"]
    output:
        directory(f"{config['output_dir']}/01_reneo_output")
    log:
        f"{config['output_dir']}/logs/reneo_binning.log"
    conda:
        config["conda_envs"]["reneo"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["prediction"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        # Run Reneo for binning
        reneo run -a {input.assembly} -r {input.reads_dir} -o {output} \
            -t {resources.threads} > {log} 2>&1
        """

# 2. Run mmseqs2 for taxonomy assignment
rule mmseqs_taxonomy:
    input:
        assembly = config["assembly_file"]
    output:
        lca_table = f"{config['output_dir']}/01_mmseqs_output/lca.tsv"
    log:
        f"{config['output_dir']}/logs/mmseqs_taxonomy.log"
    conda:
        config["conda_envs"]["mmseqs2"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["prediction"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        # Create temporary directory
        TMP_DIR=$(mktemp -d)
        
        # Run mmseqs2 for taxonomy assignment
        mmseqs easy-taxonomy {input.assembly} {config[resources][mmseqs2][db]} \
            $(dirname {output.lca_table}) $TMP_DIR \
            --threads {resources.threads} \
            --lca-ranks species,genus,family,order,class,phylum,superkingdom \
            > {log} 2>&1
            
        # Clean up
        rm -rf $TMP_DIR
        """

# 3. Filter mmseqs2 results for viral contigs
rule filter_mmseqs_lca:
    input:
        lca_table = f"{config['output_dir']}/01_mmseqs_output/lca.tsv",
        contigs = config["assembly_file"]
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

# 4. Run Jaeger for phage prediction
rule jaeger_prediction:
    input:
        assembly = config["assembly_file"]
    output:
        results = directory(f"{config['output_dir']}/01_jaeger_output"),
        predictions = f"{config['output_dir']}/01_jaeger_output/final_predictions_scored.tsv"
    log:
        f"{config['output_dir']}/logs/jaeger_prediction.log"
    conda:
        config["conda_envs"]["jaeger"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["jaeger"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        jaeger.py -i {input.assembly} -o {output.results} \
            -t {resources.threads} > {log} 2>&1
        """

# 5. Run GeNomad for viral prediction
rule genomad_prediction:
    input:
        assembly = config["assembly_file"]
    output:
        results = directory(f"{config['output_dir']}/01_genomad_output"),
        virus_summary = f"{config['output_dir']}/01_genomad_output/summary/virus_summary.tsv"
    log:
        f"{config['output_dir']}/logs/genomad_prediction.log"
    conda:
        config["conda_envs"]["genomad"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["genomad"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        genomad end-to-end {input.assembly} {output.results} \
            $(dirname {input.assembly})/genomad_db \
            --threads {resources.threads} > {log} 2>&1
        """

# 6. Run Phold for protein annotation
rule phold_prediction:
    input:
        assembly = config["assembly_file"]
    output:
        results = directory(f"{config['output_dir']}/01_phold_output"),
        predictions = f"{config['output_dir']}/01_phold_output/phold_per_cds_predictions.tsv"
    log:
        f"{config['output_dir']}/logs/phold_prediction.log"
    conda:
        config["conda_envs"]["phold"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["phold"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        # Run ORF prediction
        prodigal -i {input.assembly} -a {output.results}/proteins.faa -d {output.results}/genes.fna

        # Run phold
        phold predict --orf-file {output.results}/proteins.faa \
            --output-dir {output.results} \
            --threads {resources.threads} > {log} 2>&1
        """

# 7. Run CheckV for quality assessment
rule checkv_assessment:
    input:
        assembly = config["assembly_file"]
    output:
        results = directory(f"{config['output_dir']}/01_checkv_output"),
        quality_summary = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    log:
        f"{config['output_dir']}/logs/checkv_assessment.log"
    conda:
        config["conda_envs"]["checkv"]
    resources:
        mem_mb = config["resources"]["prediction"]["mem_mb"],
        threads = config["resources"]["checkv"]["threads"],
        time = config["resources"]["prediction"]["time"]
    shell:
        """
        # Run CheckV for viral quality assessment
        checkv end-to-end {input.assembly} {output.results} \
            -t {resources.threads} > {log} 2>&1
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
        assembly = config["assembly_file"],
        contig_ids = f"{config['output_dir']}/01_phage_predictions/contig_ids.txt"
    output:
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_phage_contigs.log"
    conda:
        "bioconda::seqkit"
    shell:
        """
        # Extract phage contigs from assembly
        seqkit grep -f {input.contig_ids} {input.assembly} > {output.phage_contigs} 2> {log}
        """