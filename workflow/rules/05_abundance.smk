"""
Rules for calculating per-sample abundance using CoverM

This module provides optional per-sample abundance calculation for phage contigs.
Abundance is calculated by mapping reads from individual samples back to the final phage sequences.
"""

# Rule to calculate per-sample abundance using CoverM
rule calculate_abundance:
    input:
        phage_seqs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta",
        reads_dir = config["reads_dir"]
    output:
        abundance = f"{config['output_dir']}/04_abundance/coverm_results.tsv"
    params:
        output_dir = f"{config['output_dir']}/04_abundance"
    log:
        f"{config['output_dir']}/logs/abundance/calculate_abundance.log"
    conda:
        "../envs/coverm.yaml"
    threads: 24
    resources:
        mem_mb=100000,
        runtime=1440  # 24 hours
    shell:
        """
        echo "========================================" > {log}
        echo "Calculating Per-Sample Abundance" >> {log}
        echo "========================================" >> {log}
        echo "" >> {log}

        # Create output directory
        mkdir -p {params.output_dir}

        # Find all FASTQ files in the reads directory
        READS=({input.reads_dir}/*.fastq.gz {input.reads_dir}/*.fq.gz {input.reads_dir}/*.fastq {input.reads_dir}/*.fq)

        # Check if we found any read files
        if [ ${{#READS[@]}} -eq 0 ]; then
            echo "ERROR: No FASTQ files found in {input.reads_dir}" >> {log}
            echo "Looking for: *.fastq.gz, *.fq.gz, *.fastq, *.fq" >> {log}
            exit 1
        fi

        echo "Found ${{#READS[@]}} read files" >> {log}
        echo "" >> {log}

        # Run CoverM in genome mode
        # This calculates abundance metrics by mapping reads to phage contigs
        echo "Running CoverM..." >> {log}

        coverm genome \
            --coupled ${{READS[@]}} \
            --genome-fasta-files {input.phage_seqs} \
            --methods rpkm tpm count variance mean covered_fraction covered_bases \
            --min-covered-fraction 0 \
            --threads {threads} \
            --output-format sparse \
            > {output.abundance} 2>> {log}

        if [ $? -eq 0 ]; then
            echo "" >> {log}
            echo "CoverM completed successfully" >> {log}
            echo "Output: {output.abundance}" >> {log}

            # Count samples and contigs
            NUM_SAMPLES=$(tail -n +2 {output.abundance} | cut -f1 | sort -u | wc -l)
            NUM_CONTIGS=$(tail -n +2 {output.abundance} | cut -f2 | sort -u | wc -l)
            echo "  Samples analyzed: $NUM_SAMPLES" >> {log}
            echo "  Contigs in output: $NUM_CONTIGS" >> {log}
        else
            echo "" >> {log}
            echo "ERROR: CoverM failed" >> {log}
            exit 1
        fi
        """


# Rule to process CoverM output into abundance matrices
rule process_abundance_matrices:
    input:
        coverm = f"{config['output_dir']}/04_abundance/coverm_results.tsv"
    output:
        tpm_matrix = f"{config['output_dir']}/04_abundance/tpm_matrix.tsv",
        count_matrix = f"{config['output_dir']}/04_abundance/count_matrix.tsv"
    params:
        coverage_threshold = config.get("abundance_coverage_threshold", 0.75)
    log:
        f"{config['output_dir']}/logs/abundance/process_abundance_matrices.log"
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        mem_mb=50000,
        runtime=120  # 2 hours
    shell:
        """
        Rscript {workflow.basedir}/scripts/process_abundance_matrices.R \
            {input.coverm} \
            {output.tpm_matrix} \
            {output.count_matrix} \
            {params.coverage_threshold} \
            > {log} 2>&1

        echo "" >> {log}
        echo "Abundance matrices created:" >> {log}
        echo "  TPM matrix: {output.tpm_matrix}" >> {log}
        echo "  Count matrix: {output.count_matrix}" >> {log}
        """
