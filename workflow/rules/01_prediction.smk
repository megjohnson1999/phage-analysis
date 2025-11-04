"""
Rules for phage prediction from metagenomic assemblies.
"""

# Helper function to determine if Reneo should be skipped
def should_skip_reneo(wildcards):
    # Skip Reneo if assembly_graph is empty or not provided
    # or if the file doesn't exist
    if not config.get("assembly_graph") or config.get("assembly_graph") == "":
        return True
    import os
    return not os.path.exists(config.get("assembly_graph", ""))

# Define a flag for skipping Reneo
# Set the global variable in the workflow module
import sys
import os
workflow.globals["use_reneo"] = not should_skip_reneo(None)

# Helper function to determine which assembly file to use
def get_assembly_input(wildcards):
    # If assembly_file is provided and not empty, use it
    if config.get("assembly_file") and config.get("assembly_file") != "":
        return config["assembly_file"]
    # If no assembly_file but assembly_graph exists, Reneo must run first 
    # and we'll use its output as the assembly file
    # This will only be used if workflow.globals["use_reneo"] is True
    if workflow.globals["use_reneo"]:
        return f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges.fasta"
    # Fallback case - this should never happen as the validation checks should catch it
    raise ValueError("Neither assembly_file nor valid assembly_graph were provided")

# Helper function to get input for reneo_binning
def get_reneo_input(wildcards):
    # Re-evaluate whether to use reneo based on current config
    # This is necessary because workflow.globals might not be set during DAG construction
    use_reneo_local = not should_skip_reneo(wildcards)
    
    
    if use_reneo_local:
        result = {
            "assembly_graph": config["assembly_graph"],
            "reads_dir": config["reads_dir"]
        }
        return result
    else:
        # Return dummy inputs that will never be used
        return {
            "assembly_graph": "/dev/null",
            "reads_dir": "/dev/null"
        }

# Helper function to get input for direct_contig_filter
def get_direct_filter_input(wildcards):
    import os
    if not workflow.globals["use_reneo"] and config.get("assembly_file") and config.get("assembly_file") != "":
        # Check if the assembly file actually exists
        if os.path.exists(config["assembly_file"]):
            return {"assembly": config["assembly_file"]}
    # Return dummy input for all other cases
    return {"assembly": "/dev/null"}

# 1. Run Reneo for binning
# Always define this rule, but it will only run when use_reneo is True
rule reneo_binning:
    input:
        unpack(get_reneo_input)
    output:
        f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges.fasta"
    log:
        f"{config['output_dir']}/logs/reneo_binning.log"
    conda:
        config["conda_envs"]["reneo"] if not config.get("conda_base_path") else None
    params:
        conda_env = config["conda_envs"]["reneo"],
        conda_base = config.get("conda_base_path", "")
    threads: 24
    shell:
        """
        if [ "{input.assembly_graph}" = "/dev/null" ]; then
            echo "Skipping Reneo - not configured for this run" > {log}
            touch {output}
            exit 0
        fi
        
        mkdir -p {config[output_dir]}/01_reneo_output

        # Handle conda environment activation
        if [ -n "{params.conda_base}" ]; then
            # Use manually specified conda installation
            echo "Using conda base path: {params.conda_base}" >> {log} 2>&1
            
            # Check if this is a path to an existing environment or just environment name
            if [[ "{params.conda_env}" == /* ]]; then
                # Full path to environment
                source {params.conda_base}/bin/activate
                conda activate {params.conda_env}
            else
                # Just environment name
                source {params.conda_base}/bin/activate
                conda activate {params.conda_env}
            fi
        else
            # Rely on Snakemake's conda handling or current environment
            echo "Using Snakemake's conda environment handling" >> {log} 2>&1
            
            # Check if reneo is available
            if ! command -v reneo &> /dev/null; then
                echo "ERROR: reneo not found in environment. Please ensure:" >> {log}
                echo "1. Reneo is installed in the conda environment specified in config" >> {log}
                echo "2. Or set conda_base_path in config to use an existing Reneo environment" >> {log}
                echo "3. Reneo requires a Gurobi license - see documentation" >> {log}
                exit 1
            fi
        fi

        # Run Reneo using wrapper script that handles expected failures
        # The wrapper script will handle koverage_genomes failures gracefully
        bash {workflow.basedir}/scripts/run_reneo_wrapper.sh \
            --input {input.assembly_graph} \
            --reads {input.reads_dir} \
            --minlength 1000 \
            --output {config[output_dir]}/01_reneo_output \
            --threads {threads} > {log} 2>&1 || {{
                echo "Reneo wrapper reported failure, checking for partial output..." >> {log}
                # Check if we got the output we need despite the failure
                if [ -f "{output}" ] && [ -s "{output}" ]; then
                    echo "Found valid output file despite Reneo failure, continuing..." >> {log}
                    exit 0
                else
                    echo "No valid output found, failing the rule" >> {log}
                    exit 1
                fi
            }}

        # Clean up large temporary directories to save space
        echo "Cleaning up temporary directories..." >> {log} 2>&1
        if [ -d "{config[output_dir]}/01_reneo_output/temp" ]; then
            echo "Removing temp directory (size: $(du -sh {config[output_dir]}/01_reneo_output/temp 2>/dev/null | cut -f1))" >> {log} 2>&1
            rm -rf "{config[output_dir]}/01_reneo_output/temp"
        fi
        if [ -d "{config[output_dir]}/01_reneo_output/work" ]; then
            echo "Removing work directory (size: $(du -sh {config[output_dir]}/01_reneo_output/work 2>/dev/null | cut -f1))" >> {log} 2>&1
            rm -rf "{config[output_dir]}/01_reneo_output/work"
        fi
        if [ -d "{config[output_dir]}/01_reneo_output/.snakemake" ]; then
            echo "Removing .snakemake directory (size: $(du -sh {config[output_dir]}/01_reneo_output/.snakemake 2>/dev/null | cut -f1))" >> {log} 2>&1
            rm -rf "{config[output_dir]}/01_reneo_output/.snakemake"
        fi
        echo "Cleanup complete" >> {log} 2>&1
        """

# 1b. Filter contigs by length (1KB) from Reneo output
# Always define this rule
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
        if [ ! -s "{input}" ]; then
            echo "Input file is empty or missing, creating empty output" > {log}
            touch {output}
            exit 0
        fi
        
        seqkit seq --min-len 1000 -g \
            "{input}" > "{output}"
        """

# 1c. Filter assembly directly (when not using Reneo)
# Always define this rule
rule direct_contig_filter:
    input:
        unpack(get_direct_filter_input)
    output:
        filtered_assembly = f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta"
    log:
        f"{config['output_dir']}/logs/direct_contig_filter.log"
    conda:
        config["conda_envs"]["seqkit"]
    threads: 8
    shell:
        """
        if [ "{input.assembly}" = "/dev/null" ]; then
            echo "Skipping direct filter - not configured for this run" > {log}
            mkdir -p {config[output_dir]}/01_filtered_assembly
            touch {output.filtered_assembly}
            exit 0
        fi
        
        mkdir -p {config[output_dir]}/01_filtered_assembly
        seqkit seq --min-len 1000 -g \
            "{input.assembly}" > "{output.filtered_assembly}"
        """

# 2. Run mmseqs2 for taxonomy assignment
rule mmseqs_taxonomy:
    input:
        # Use the appropriate filtered contigs based on whether Reneo is enabled
        # During DAG construction, we need to specify the exact file that will exist
        filtered_contigs = f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if workflow.globals["use_reneo"] else 
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
        
        # Get the base input filename without extension
        INPUT_BASE=$(basename {input.filtered_contigs} .fasta)

        mkdir -p {config[output_dir]}/01_mmseqs_output
        
        # Run mmseqs2 for taxonomy assignment
        mmseqs easy-taxonomy {input.filtered_contigs} \
            {config[databases][mmseqs2][db]} \
            {config[output_dir]}/01_mmseqs_output/temp_results \
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
        
        # Move the generated LCA file to the expected output location
        mv {config[output_dir]}/01_mmseqs_output/temp_results_lca.tsv {output.lca_table}
            
        # Clean up
        rm -rf $TMP_DIR
        """

# 3. Filter mmseqs2 results for viral contigs
rule filter_mmseqs_lca:
    input:
        lca_table = f"{config['output_dir']}/01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv",
        contigs = f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if workflow.globals["use_reneo"] else 
            f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta"
    output:
        filtered_lca = f"{config['output_dir']}/01_filtered_mmseqs/filtered_lca.tsv",
        passing_ids = f"{config['output_dir']}/01_filtered_mmseqs/passing_contig_ids.txt",
        missing_ids = f"{config['output_dir']}/01_filtered_mmseqs/missing_contig_ids.txt"
    log:
        f"{config['output_dir']}/logs/filter_mmseqs_lca.log"
    conda:
        config["conda_envs"]["python"]
    shell:
        """
        # Use python script to filter mmseqs2 results
        python {workflow.basedir}/scripts/01_filterMmseqsLca.py \
            --mmseqs_LCA_table {input.lca_table} \
            --contigs {input.contigs} \
            --o_filtered_LCA_table {output.filtered_lca} \
            --o_passing_contig_ids {output.passing_ids} \
            --o_missing_contig_ids {output.missing_ids} \
            > {log} 2>&1
        """


# 3b. Extract passing viral contigs
rule extract_viral_contigs:
    input:
        contigs = f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" 
            if workflow.globals["use_reneo"] else 
            f"{config['output_dir']}/01_filtered_assembly/filtered_assembly_1KB.fasta",
        passing_ids = f"{config['output_dir']}/01_filtered_mmseqs/passing_contig_ids.txt"
    output:
        viral_contigs = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    log:
        f"{config['output_dir']}/logs/extract_viral_contigs.log"
    conda:
        config["conda_envs"]["seqkit"]
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
        predictions = f"{config['output_dir']}/01_jaeger_output/passing_Viralcontigs_default_jaeger.tsv"
    log:
        f"{config['output_dir']}/logs/jaeger_prediction.log"
    conda:
        config["conda_envs"]["jaeger"]
    shell:
        """
        mkdir -p {output.results}
        
        Jaeger -i {input.assembly} -o {output.results} \
            -s 2.5 \
            --fsize 1000 \
            --stride 1000 > {log} 2>&1
            
        # Check if the output exists - if not, create an empty file to satisfy Snakemake
        if [ ! -f {output.predictions} ]; then
            echo "Warning: Jaeger did not produce output. Creating empty file." >> {log}
            touch {output.predictions}
        fi
        """

# 5. Run GeNomad for viral prediction
rule genomad_prediction:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        results = directory(f"{config['output_dir']}/01_genomad_output"),
        virus_summary = f"{config['output_dir']}/01_genomad_output/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv"
    log:
        f"{config['output_dir']}/logs/genomad_prediction.log"
    conda:
        config["conda_envs"]["genomad"]
    threads: 24
    shell:
        """
        #genomad download-database {config[databases][genomad][db]}

        genomad end-to-end --min-score 0.6 \
            --cleanup \
            --threads {threads} \
            {input.assembly} \
            {output.results} \
            {config[databases][genomad][db]} > {log} 2>&1
        """

# 6a. Split viral contigs for parallel PHOLD processing
rule split_viral_contigs_for_phold:
    input:
        assembly = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        split_dir = directory(f"{config['output_dir']}/01_phold_split_seqs"),
        split_list = f"{config['output_dir']}/01_phold_split_seqs/split_file_list.txt"
    params:
        # Number of sequences per chunk - adjust based on expected sequence sizes
        chunk_size = 1000
    log:
        f"{config['output_dir']}/logs/split_viral_contigs_for_phold.log"
    conda:
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Create output directory and ensure directory for final output exists
        mkdir -p {output.split_dir}
        mkdir -p {config[output_dir]}/01_phold_output/tmp
        
        # Count total sequences to calculate appropriate chunking
        TOTAL_SEQS=$(seqkit stats -T {input.assembly} | tail -n 1 | cut -f 4)
        echo "Total sequences: $TOTAL_SEQS" > {log} 2>&1
        
        # Check if input file is empty
        if [ "$TOTAL_SEQS" -eq 0 ]; then
            echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
            # Create an empty placeholder file to satisfy workflow
            touch {output.split_dir}/empty.fasta
            echo "{output.split_dir}/empty.fasta" > {output.split_list}
            echo "Created empty placeholder file" >> {log} 2>&1
        else
            # Split FASTA file into chunked files with multiple sequences per file
            seqkit split {input.assembly} -O {output.split_dir} -p {params.chunk_size} >> {log} 2>&1
            
            # Create list of split files - use absolute paths for reliability
            find {output.split_dir} -name "*.fasta" -type f | sort > {output.split_list}
            
            # Report chunking results
            CHUNK_COUNT=$(wc -l < {output.split_list})
            echo "Created $CHUNK_COUNT chunk files from $TOTAL_SEQS sequences" >> {log} 2>&1
            
            # Validate the split files
            for file in $(cat {output.split_list}); do
                if [ ! -s "$file" ]; then
                    echo "Warning: Empty split file detected: $file" >> {log} 2>&1
                fi
            done
        fi
        """

# Helper function to get PHOLD samples
def get_phold_samples():
    # After split_viral_contigs_for_phold is run, this reads the split file list
    split_list = f"{config['output_dir']}/01_phold_split_seqs/split_file_list.txt"
    
    # Make sure the required modules are imported
    import os
    import glob
    
    # Primary method: read from the split file list if it exists
    if os.path.exists(split_list):
        try:
            with open(split_list, "r") as f:
                files = [line.strip() for line in f if line.strip()]
            # Extract filenames without extension and ensure unique values
            samples = [os.path.splitext(os.path.basename(file))[0] for file in files]
            # Filter out any empty strings
            samples = [s for s in samples if s]
            # If we found samples, return them
            if samples:
                return samples
        except Exception as e:
            print(f"Warning: Error reading split file list: {e}")
            # Continue to fallback methods
    
    # Fallback method 1: Use glob to find fasta files directly
    try:
        split_dir = f"{config['output_dir']}/01_phold_split_seqs"
        if os.path.exists(split_dir) and os.path.isdir(split_dir):
            samples = [os.path.splitext(os.path.basename(f))[0] 
                    for f in glob.glob(f"{split_dir}/*.fasta")]
            # Filter out any empty strings
            samples = [s for s in samples if s]
            if samples:
                return samples
    except Exception as e:
        print(f"Warning: Error using glob to find fasta files: {e}")
    
    # Fallback method 2: Check tmp dir for existing predictions
    try:
        tmp_dir = f"{config['output_dir']}/01_phold_output/tmp"
        if os.path.exists(tmp_dir) and os.path.isdir(tmp_dir):
            samples = [d for d in os.listdir(tmp_dir) 
                      if os.path.isdir(os.path.join(tmp_dir, d))]
            if samples:
                return samples
    except Exception as e:
        print(f"Warning: Error checking tmp dir for predictions: {e}")
    
    # If all methods fail, return empty list
    return []

# 6b. Checkpoint to wait for split files to be created
checkpoint wait_for_phold_splits:
    input:
        split_list = f"{config['output_dir']}/01_phold_split_seqs/split_file_list.txt"
    output:
        touch(f"{config['output_dir']}/01_phold_output/.splits_ready")
    shell:
        "mkdir -p $(dirname {output})"

# 6c. Run PHOLD on a chunk of contigs
rule phold_single_prediction:
    input:
        checkpoint = f"{config['output_dir']}/01_phold_output/.splits_ready",
        contig_file = f"{config['output_dir']}/01_phold_split_seqs/{{sample}}.fasta"
    output:
        results_dir = directory(f"{config['output_dir']}/01_phold_output/tmp/{{sample}}"),
        predictions = f"{config['output_dir']}/01_phold_output/tmp/{{sample}}/phold_per_cds_predictions.tsv"
    log:
        f"{config['output_dir']}/logs/phold_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phold"]
    threads: 24
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}
        
        # Run phold on chunk of contigs
        (phold run -i {input.contig_file} \
            -o {output.results_dir} \
            -d {config[databases][phold][db]} \
            -t {threads} --cpu --force > {log} 2>&1) || true
            
        # Check if output file exists, if not create empty file with header
        if [ ! -f "{output.predictions}" ]; then
            echo "WARNING: PHOLD failed to create output for {wildcards.sample}, creating empty output file" >> {log}
            # Create directory if it doesn't exist
            mkdir -p $(dirname {output.predictions})
            # Create empty file with header structure
            echo -e "contig_id\torf_id\tstart\tend\tstrand\taa_length\tcategory\tproduct\thit\tevalue\tidentity" > {output.predictions}
            echo "Created empty PHOLD predictions file with header only" >> {log}
        fi
        """

# 6d. Helper rule to force running all PHOLD predictions
rule run_all_phold_predictions:
    input:
        checkpoint = f"{config['output_dir']}/01_phold_output/.splits_ready",
        # For actual runs, get samples from the split files
        # For dry runs, this will be an empty list, which is fine
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/01_phold_output/tmp/{{sample}}/phold_per_cds_predictions.tsv",
            sample=get_phold_samples()
        )
    output:
        touch(f"{config['output_dir']}/01_phold_output/.all_predictions_done")

# 6e. Aggregate PHOLD results
rule phold_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        # Aggregation only happens after all individual predictions are done
        all_done = f"{config['output_dir']}/01_phold_output/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/01_phold_output/tmp/{{sample}}/phold_per_cds_predictions.tsv",
            sample=get_phold_samples()
        )
    output:
        predictions = f"{config['output_dir']}/01_phold_output/phold_per_cds_predictions.tsv"
    log:
        f"{config['output_dir']}/logs/phold_aggregate_results.log"
    conda:
        config["conda_envs"]["phold"]
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHOLD results" > {log} 2>&1
        
        # Create a temporary directory for processing
        TMP_DIR=$(mktemp -d)
        
        # First create a list of all prediction files
        find "$RESULTS_DIR/tmp" -name "phold_per_cds_predictions.tsv" -type f > "$TMP_DIR/prediction_files.txt"
        
        # Count how many files we found
        FILE_COUNT=$(wc -l < "$TMP_DIR/prediction_files.txt")
        echo "Found $FILE_COUNT prediction files to process" >> {log} 2>&1
        
        if [ "$FILE_COUNT" -gt 0 ]; then
            # Get the first file to extract header
            FIRST_FILE=$(head -n 1 "$TMP_DIR/prediction_files.txt")
            
            # Create the output file with header
            head -n 1 "$FIRST_FILE" > {output.predictions}
            
            # Process files in batches to avoid command line length limits
            # Use a while loop to read the file list instead of xargs
            while read -r pred_file; do
                # Skip header line (first line) from each file
                awk 'NR>1' "$pred_file" >> "$TMP_DIR/aggregated_data.tmp"
            done < "$TMP_DIR/prediction_files.txt"
            
            # Append all data to the output file
            cat "$TMP_DIR/aggregated_data.tmp" >> {output.predictions}
            
            # Count records in final file
            RECORD_COUNT=$(wc -l < {output.predictions})
            RECORD_COUNT=$((RECORD_COUNT - 1))  # Subtract 1 for header
            echo "Successfully compiled PHOLD results with $RECORD_COUNT data records" >> {log} 2>&1
        else
            # Create empty output with header structure
            echo "contig_id\torf_id\tstart\tend\tstrand\taa_length\tcategory\tproduct\thit\tevalue\tidentity" > {output.predictions}
            echo "No PHOLD result files found, created empty file with header" >> {log} 2>&1
        fi
        
        # Clean up temporary directory
        rm -rf "$TMP_DIR"
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
        jaeger = f"{config['output_dir']}/01_jaeger_output/passing_Viralcontigs_default_jaeger.tsv",
        genomad = f"{config['output_dir']}/01_genomad_output/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv",
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
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Extract phage contigs from viral contigs
        seqkit grep -f {input.contig_ids} {input.contigs} > {output.phage_contigs} 2> {log}
        """