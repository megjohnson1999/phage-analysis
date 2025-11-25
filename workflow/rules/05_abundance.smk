"""
Rules for calculating per-sample abundance using CoverM

This module provides optional per-sample abundance calculation for phage contigs and ORFs.
Abundance is calculated by mapping reads from individual samples back to the final phage sequences.

When calculate_abundance is enabled, this module generates:
- Contig-level abundance (TPM and count matrices)
- ORF-level abundance (TPM and count matrices)
- ORF annotations linking ORFs to contigs and functions

Note: Uses get_phage_input() function defined in 02_clustering.smk
"""

# Rule to calculate per-sample abundance using CoverM
rule calculate_abundance:
    input:
        phage_seqs = get_phage_input,  # Use helper function instead of hardcoded path
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
        set +u  # Disable strict mode for unbound variables

        echo "========================================" > {log}
        echo "Calculating Per-Sample Abundance" >> {log}
        echo "========================================" >> {log}
        echo "" >> {log}

        # Create output directory
        mkdir -p {params.output_dir}

        # Find paired-end read files
        # Look for common paired-end patterns: _R1/_R2, _1/_2, .1/.2
        echo "Searching for paired-end read files..." >> {log}

        # Create arrays for R1 and R2 files
        declare -a R1_FILES
        declare -a R2_FILES

        # Search for R1/R2 pattern (both compressed and uncompressed)
        for r1 in {input.reads_dir}/*_R1*.fastq.gz {input.reads_dir}/*_R1*.fq.gz \
                  {input.reads_dir}/*_1.fastq.gz {input.reads_dir}/*_1.fq.gz \
                  {input.reads_dir}/*.1.fastq.gz {input.reads_dir}/*.1.fq.gz \
                  {input.reads_dir}/*_R1*.fastq {input.reads_dir}/*_R1*.fq \
                  {input.reads_dir}/*_1.fastq {input.reads_dir}/*_1.fq \
                  {input.reads_dir}/*.1.fastq {input.reads_dir}/*.1.fq; do
            if [ -f "$r1" ]; then
                # Try to find corresponding R2 file
                r2="${{r1/_R1/_R2}}"
                r2="${{r2/_1./_2.}}"
                r2="${{r2/.1./.2.}}"
                r2="${{r2/_1_/_2_}}"

                if [ -f "$r2" ]; then
                    R1_FILES+=("$r1")
                    R2_FILES+=("$r2")
                    echo "  Found pair: $(basename $r1) / $(basename $r2)" >> {log}
                fi
            fi
        done

        # Check if we found any paired files
        if [ ${{#R1_FILES[@]}} -eq 0 ]; then
            echo "ERROR: No paired-end FASTQ files found in {input.reads_dir}" >> {log}
            echo "Looked for patterns: _R1/_R2, _1/_2, .1/.2" >> {log}
            exit 1
        fi

        echo "" >> {log}
        echo "Found ${{#R1_FILES[@]}} paired-end samples" >> {log}
        echo "" >> {log}

        # Build interleaved list of R1 R2 pairs for --coupled
        PAIRED_READS=()
        for i in "${{!R1_FILES[@]}}"; do
            PAIRED_READS+=("${{R1_FILES[$i]}}" "${{R2_FILES[$i]}}")
        done

        # Run CoverM in contig mode
        # This calculates abundance metrics for each individual contig
        echo "Running CoverM on ${{#R1_FILES[@]}} sample pairs..." >> {log}

        coverm contig \
            --coupled "${{PAIRED_READS[@]}}" \
            --reference {input.phage_seqs} \
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


# ORF-level abundance calculation
# Rule to calculate per-ORF abundance using CoverM
rule calculate_orf_abundance:
    input:
        orf_seqs = f"{config['output_dir']}/03_orf_predictions/genes.fna",
        reads_dir = config["reads_dir"]
    output:
        abundance = f"{config['output_dir']}/04_abundance/orf_coverm_results.tsv"
    params:
        output_dir = f"{config['output_dir']}/04_abundance"
    log:
        f"{config['output_dir']}/logs/abundance/calculate_orf_abundance.log"
    conda:
        "../envs/coverm.yaml"
    threads: 24
    resources:
        mem_mb=100000,
        runtime=1440  # 24 hours
    shell:
        """
        set +u  # Disable strict mode for unbound variables

        echo "========================================" > {log}
        echo "Calculating Per-Sample ORF Abundance" >> {log}
        echo "========================================" >> {log}
        echo "" >> {log}

        # Create output directory
        mkdir -p {params.output_dir}

        # Check if ORF file exists and has content
        if [ ! -f "{input.orf_seqs}" ] || [ ! -s "{input.orf_seqs}" ]; then
            echo "ERROR: ORF sequences file is empty or missing: {input.orf_seqs}" >> {log}
            exit 1
        fi

        ORF_COUNT=$(grep -c ">" {input.orf_seqs} || echo "0")
        echo "Found $ORF_COUNT ORF sequences" >> {log}
        echo "" >> {log}

        # Find paired-end read files
        # Look for common paired-end patterns: _R1/_R2, _1/_2, .1/.2
        echo "Searching for paired-end read files..." >> {log}

        # Create arrays for R1 and R2 files
        declare -a R1_FILES
        declare -a R2_FILES

        # Search for R1/R2 pattern (both compressed and uncompressed)
        for r1 in {input.reads_dir}/*_R1*.fastq.gz {input.reads_dir}/*_R1*.fq.gz \
                  {input.reads_dir}/*_1.fastq.gz {input.reads_dir}/*_1.fq.gz \
                  {input.reads_dir}/*.1.fastq.gz {input.reads_dir}/*.1.fq.gz \
                  {input.reads_dir}/*_R1*.fastq {input.reads_dir}/*_R1*.fq \
                  {input.reads_dir}/*_1.fastq {input.reads_dir}/*_1.fq \
                  {input.reads_dir}/*.1.fastq {input.reads_dir}/*.1.fq; do
            if [ -f "$r1" ]; then
                # Try to find corresponding R2 file
                r2="${{r1/_R1/_R2}}"
                r2="${{r2/_1./_2.}}"
                r2="${{r2/.1./.2.}}"
                r2="${{r2/_1_/_2_}}"

                if [ -f "$r2" ]; then
                    R1_FILES+=("$r1")
                    R2_FILES+=("$r2")
                    echo "  Found pair: $(basename $r1) / $(basename $r2)" >> {log}
                fi
            fi
        done

        # Check if we found any paired files
        if [ ${{#R1_FILES[@]}} -eq 0 ]; then
            echo "ERROR: No paired-end FASTQ files found in {input.reads_dir}" >> {log}
            echo "Looked for patterns: _R1/_R2, _1/_2, .1/.2" >> {log}
            exit 1
        fi

        echo "" >> {log}
        echo "Found ${{#R1_FILES[@]}} paired-end samples" >> {log}
        echo "" >> {log}

        # Build interleaved list of R1 R2 pairs for --coupled
        PAIRED_READS=()
        for i in "${{!R1_FILES[@]}}"; do
            PAIRED_READS+=("${{R1_FILES[$i]}}" "${{R2_FILES[$i]}}")
        done

        # Run CoverM in contig mode on ORF sequences
        # This calculates abundance metrics for each individual ORF
        echo "Running CoverM on ORF sequences with ${{#R1_FILES[@]}} sample pairs..." >> {log}

        coverm contig \
            --coupled "${{PAIRED_READS[@]}}" \
            --reference {input.orf_seqs} \
            --methods rpkm tpm count variance mean covered_fraction covered_bases \
            --min-covered-fraction 0 \
            --threads {threads} \
            --output-format sparse \
            > {output.abundance} 2>> {log}

        if [ $? -eq 0 ]; then
            echo "" >> {log}
            echo "CoverM ORF abundance completed successfully" >> {log}
            echo "Output: {output.abundance}" >> {log}

            # Count samples and ORFs
            NUM_SAMPLES=$(tail -n +2 {output.abundance} | cut -f1 | sort -u | wc -l)
            NUM_ORFS=$(tail -n +2 {output.abundance} | cut -f2 | sort -u | wc -l)
            echo "  Samples analyzed: $NUM_SAMPLES" >> {log}
            echo "  ORFs in output: $NUM_ORFS" >> {log}
        else
            echo "" >> {log}
            echo "ERROR: CoverM ORF abundance failed" >> {log}
            exit 1
        fi
        """


# Rule to process ORF CoverM output into abundance matrices
rule process_orf_abundance_matrices:
    input:
        coverm = f"{config['output_dir']}/04_abundance/orf_coverm_results.tsv"
    output:
        tpm_matrix = f"{config['output_dir']}/04_abundance/orf_tpm_matrix.tsv",
        count_matrix = f"{config['output_dir']}/04_abundance/orf_count_matrix.tsv"
    params:
        coverage_threshold = config.get("abundance_coverage_threshold", 0.75)
    log:
        f"{config['output_dir']}/logs/abundance/process_orf_abundance_matrices.log"
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
        echo "ORF abundance matrices created:" >> {log}
        echo "  ORF TPM matrix: {output.tpm_matrix}" >> {log}
        echo "  ORF count matrix: {output.count_matrix}" >> {log}
        """


# Rule to create ORF annotations linking ORFs to contigs and functions
rule create_orf_annotations:
    input:
        orf_seqs = f"{config['output_dir']}/03_orf_predictions/genes.fna",
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
    output:
        annotations = f"{config['output_dir']}/04_abundance/orf_annotations.tsv"
    log:
        f"{config['output_dir']}/logs/abundance/create_orf_annotations.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/create_orf_annotations.py \
            --genes {input.orf_seqs} \
            --proteins {input.proteins} \
            --output {output.annotations} \
            > {log} 2>&1
        """
