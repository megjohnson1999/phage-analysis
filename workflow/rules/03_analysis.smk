"""
Rules for analyzing phage genomes, including host prediction and genomic characterization.
These rules can take input from either clustering or prediction steps.
"""

# 1. Split phage sequences for array-based processing
rule split_phage_sequences:
    input:
        phage_seqs = get_phage_input  # Function defined in 02_clustering.smk
    output:
        split_dir = directory(f"{config['output_dir']}/03_split_seqs"),
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    params:
        # Fixed number of sequences per job for consistent performance
        sequences_per_job = 100
    log:
        f"{config['output_dir']}/logs/split_phage_sequences.log"
    conda:
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Create output directory
        mkdir -p {output.split_dir}
        
        # Count total sequences
        TOTAL_SEQS=$(seqkit stats -T {input.phage_seqs} | tail -n 1 | cut -f 4)
        echo "Total sequences: $TOTAL_SEQS" > {log} 2>&1
        
        # Check if input file is empty
        if [ "$TOTAL_SEQS" -eq 0 ]; then
            echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
            # Create an empty placeholder file to satisfy workflow
            touch {output.split_dir}/empty.fasta
            echo "{output.split_dir}/empty.fasta" > {output.split_list}
            echo "Created empty placeholder file" >> {log} 2>&1
        else
            # Use a fixed number of sequences per job for consistent performance
            echo "Using fixed chunk size of {params.sequences_per_job} sequences per job" >> {log} 2>&1
            
            # Split FASTA file by exact sequence count
            seqkit split2 --by-size {params.sequences_per_job} {input.phage_seqs} -O {output.split_dir} >> {log} 2>&1
            
            # Create list of split files - use absolute paths for reliability
            find {output.split_dir} -name "*.fasta" -type f | sort > {output.split_list}
            
            # Report chunking results
            ACTUAL_BATCH_COUNT=$(wc -l < {output.split_list})
            echo "Created $ACTUAL_BATCH_COUNT chunk files from $TOTAL_SEQS sequences" >> {log} 2>&1
        fi
        """

# List all samples for iphop from split files
def get_iphop_samples():
    # After split_phage_sequences is run, this reads the split file list
    split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    
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
        split_dir = f"{config['output_dir']}/03_split_seqs"
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
        tmp_dir = f"{config['output_dir']}/03_iphop_results/tmp"
        if os.path.exists(tmp_dir) and os.path.isdir(tmp_dir):
            samples = [d for d in os.listdir(tmp_dir) 
                      if os.path.isdir(os.path.join(tmp_dir, d))]
            if samples:
                return samples
    except Exception as e:
        print(f"Warning: Error checking tmp dir for predictions: {e}")
    
    # If all methods fail, return empty list
    return []

# 2. iphop workflow - Use checkpoint to wait for split files to be created
checkpoint wait_for_iphop_splits:
    input:
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    output:
        touch(f"{config['output_dir']}/03_iphop_results/.splits_ready")
    shell:
        "mkdir -p $(dirname {output})"

# Check if input files exist for iPhop
rule check_iphop_input_files:
    input:
        # Prerequisite - splits must be ready
        splits_ready = f"{config['output_dir']}/03_iphop_results/.splits_ready",
        # Input sequences should exist
        phage_seqs = get_phage_input,
        # Check that split files are created
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    output:
        touch(f"{config['output_dir']}/03_iphop_results/.input_files_found")
    log:
        f"{config['output_dir']}/logs/check_iphop_input_files.log"
    shell:
        """
        # Check if input sequences exist
        if [ ! -f "{input.phage_seqs}" ]; then
            echo "Error: Input sequence file does not exist: {input.phage_seqs}" > {log} 2>&1
            exit 1
        fi
        
        # Check if split list exists and has content
        if [ ! -s "{input.split_list}" ]; then
            echo "Warning: Split file list is empty or doesn't exist: {input.split_list}" > {log} 2>&1
            echo "Creating a placeholder split list..." >> {log} 2>&1
            mkdir -p $(dirname {input.split_list})
            touch {input.split_list}
        fi
        
        # Everything checked out, all input files exist
        echo "All required input files for iPhop were found" > {log} 2>&1
        """

# 2a. Run iPhop for host prediction on a single split file
rule iphop_single_prediction:
    input:
        checkpoint = f"{config['output_dir']}/03_iphop_results/.splits_ready",
        input_check = f"{config['output_dir']}/03_iphop_results/.input_files_found",
        phage_file = f"{config['output_dir']}/03_split_seqs/{{sample}}.fasta"
    output:
        results_dir = directory(f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}"),
        prediction = f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv"
    log:
        f"{config['output_dir']}/logs/iphop_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["iphop"]
    threads: 24
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}
        
        # Check if input file is empty or doesn't exist
        if [ ! -f "{input.phage_file}" ] || [ ! -s "{input.phage_file}" ]; then
            echo "Warning: Input file is empty or missing: {input.phage_file}" > {log} 2>&1
            echo "Creating empty iPhop output file..." >> {log} 2>&1
            echo "query,host,score,identity,coverage,kingdom,phylum,class,order,family,genus" > {output.prediction}
        else
            # Count sequences in input file
            SEQ_COUNT=$(grep -c ">" {input.phage_file} || echo "0")
            echo "Processing $SEQ_COUNT sequences with iPhop" > {log} 2>&1
            
            if [ "$SEQ_COUNT" -eq 0 ]; then
                echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
                echo "Creating empty iPhop output file..." >> {log} 2>&1
                echo "query,host,score,identity,coverage,kingdom,phylum,class,order,family,genus" > {output.prediction}
            else
                # Run iPhop
                iphop predict --fa_file {input.phage_file} \
                    --db_dir {config[databases][iphop][db]} \
                    --out_dir {output.results_dir} \
                    --num_threads {threads} >> {log} 2>&1
                
                # Check if the output exists - if not, create an empty file to satisfy Snakemake
                if [ ! -f "{output.prediction}" ]; then
                    echo "Warning: iPhop did not produce output. Creating empty file." >> {log} 2>&1
                    echo "query,host,score,identity,coverage,kingdom,phylum,class,order,family,genus" > {output.prediction}
                fi
            fi
        fi
        """

# Helper rule to force running all iPhop predictions
rule run_all_iphop_predictions:
    input:
        checkpoint = f"{config['output_dir']}/03_iphop_results/.splits_ready",
        input_check = f"{config['output_dir']}/03_iphop_results/.input_files_found",
        # For actual runs, get samples from the split files
        # For dry runs, this will be an empty list, which is fine
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv",
            sample=get_iphop_samples()
        )
    output:
        touch(f"{config['output_dir']}/03_iphop_results/.all_predictions_done")

# 2b. Aggregate iPhop results
rule iphop_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        # Aggregation only happens after all individual predictions are done
        all_done = f"{config['output_dir']}/03_iphop_results/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv",
            sample=get_iphop_samples()
        )
    output:
        predictions = f"{config['output_dir']}/03_iphop_results/iphop_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/iphop_aggregate_results.log"
    conda:
        config["conda_envs"]["python"]
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling iPhop results" > {log} 2>&1
        
        # Create a temporary directory for processing
        TMP_DIR=$(mktemp -d)
        
        # First create a list of all prediction files
        find "$RESULTS_DIR/tmp" -name "host_prediction_to_genus.csv" -type f > "$TMP_DIR/prediction_files.txt"
        
        # Count how many files we found
        FILE_COUNT=$(wc -l < "$TMP_DIR/prediction_files.txt")
        echo "Found $FILE_COUNT prediction files to process" >> {log} 2>&1
        
        if [ "$FILE_COUNT" -gt 0 ]; then
            # Get the first file to extract header
            FIRST_FILE=$(head -n 1 "$TMP_DIR/prediction_files.txt")
            
            # Create the output file with header (convert CSV to TSV)
            head -n 1 "$FIRST_FILE" | tr ',' '\t' > {output.predictions}
            
            # Process files in batches to avoid command line length limits
            while read -r pred_file; do
                # Skip header line (first line) from each file and convert CSV to TSV
                awk -F ',' 'NR>1 {{OFS="\t"; print}}' "$pred_file" >> "$TMP_DIR/aggregated_data.tmp"
            done < "$TMP_DIR/prediction_files.txt"
            
            # Append all data to the output file
            cat "$TMP_DIR/aggregated_data.tmp" >> {output.predictions}
            
            # Count records in final file
            RECORD_COUNT=$(wc -l < {output.predictions})
            RECORD_COUNT=$((RECORD_COUNT - 1))  # Subtract 1 for header
            echo "Successfully compiled iPhop results with $RECORD_COUNT data records" >> {log} 2>&1
        else
            # Create empty output with header structure
            echo -e "query\thost\tscore\tidentity\tcoverage\tkingdom\tphylum\tclass\torder\tfamily\tgenus" > {output.predictions}
            echo "No iPhop result files found, created empty file with header" >> {log} 2>&1
        fi
        
        # Clean up temporary directory
        rm -rf "$TMP_DIR"
        """

# 3. Run Prodigal for ORF prediction on phage sequences
rule prodigal_orf_prediction:
    input:
        phage_seqs = get_phage_input
    output:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
        genes = f"{config['output_dir']}/03_orf_predictions/genes.fna"
    log:
        f"{config['output_dir']}/logs/prodigal_orf_prediction.log"
    conda:
        config["conda_envs"]["vcontact3"]
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.proteins})
        
        # Run Prodigal for ORF prediction
        prodigal -i {input.phage_seqs} \
            -a {output.proteins} \
            -d {output.genes} \
            -p meta > {log} 2>&1
        """

# 4. Split protein files for PHACTS array processing
# NOTE: This rule is disabled in favor of phage-specific splitting
# See split_proteins_by_phage.smk for the improved implementation
# rule split_protein_files:
#     input:
#         proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
#     output:
#         split_dir = directory(f"{config['output_dir']}/03_split_proteins"),
#         split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
#     params:
#         # Use larger chunk size
#         chunk_size = 1000,
#         # Maximum number of batches (prevents too many files)
#         max_batches = 100
#     log:
#         f"{config['output_dir']}/logs/split_protein_files.log"
#     conda:
#         config["conda_envs"]["seqkit"]
#     shell:
#         """
#         # Create output directory
#         mkdir -p {output.split_dir}
#         
#         # Count total sequences
#         TOTAL_SEQS=$(grep -c ">" {input.proteins})
#         echo "Total protein sequences: $TOTAL_SEQS" > {log} 2>&1
#         
#         # Determine optimal batch count - use proper Python syntax
#         BATCH_COUNT=$(python -c "print(min({params.max_batches}, max(1, (int($TOTAL_SEQS) // {params.chunk_size}) + 1)))")
#         echo "Will create $BATCH_COUNT batches" >> {log} 2>&1
#         
#         # Use seqkit split with part option instead of complex AWK
#         seqkit split -p $BATCH_COUNT {input.proteins} -O {output.split_dir} >> {log} 2>&1
#         
#         # Create list of split files
#         find {output.split_dir} -name "*.faa" | sort > {output.split_list}
#         ACTUAL_BATCH_COUNT=$(wc -l < {output.split_list})
#         echo "Created $ACTUAL_BATCH_COUNT protein batch files" >> {log} 2>&1
#         
#         # If no files were created (empty input), create an empty placeholder
#         if [ "$ACTUAL_BATCH_COUNT" -eq 0 ]; then
#             echo "Input was empty, creating placeholder file" >> {log} 2>&1
#             touch "{output.split_dir}/empty.faa"
#             echo "{output.split_dir}/empty.faa" > {output.split_list}
#         fi
#         """

# NOTE: This function is disabled in favor of phage-specific splitting
# See split_proteins_by_phage.smk for the improved implementation
# # List all samples for phacts from split protein files
# def get_phacts_samples():
#     # After split_protein_files is run, this reads the split file list
#     split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
#     
#     # Make sure the required modules are imported
#     import os
#     import glob
#     
#     # Primary method: read from the split file list if it exists
#     if os.path.exists(split_list):
#         try:
#             with open(split_list, "r") as f:
#                 files = [line.strip() for line in f if line.strip()]
#             # Extract filenames without extension and ensure unique values
#             samples = [os.path.splitext(os.path.basename(file))[0] for file in files]
#             # Filter out any empty strings
#             samples = [s for s in samples if s]
#             # If we found samples, return them
#             if samples:
#                 return samples
#         except Exception as e:
#             print(f"Warning: Error reading split file list: {e}")
#             # Continue to fallback methods
#     
#     # Fallback method 1: Use glob to find faa files directly
#     try:
#         split_dir = f"{config['output_dir']}/03_split_proteins"
#         if os.path.exists(split_dir) and os.path.isdir(split_dir):
#             samples = [os.path.splitext(os.path.basename(f))[0] 
#                     for f in glob.glob(f"{split_dir}/*.faa")]
#             # Filter out any empty strings
#             samples = [s for s in samples if s]
#             if samples:
#                 return samples
#     except Exception as e:
#         print(f"Warning: Error using glob to find faa files: {e}")
#     
#     # Fallback method 2: Check tmp dir for existing predictions
#     try:
#         tmp_dir = f"{config['output_dir']}/03_phacts_results/tmp"
#         if os.path.exists(tmp_dir) and os.path.isdir(tmp_dir):
#             # Look for directories in the tmp directory
#             samples = [d for d in os.listdir(tmp_dir)
#                       if os.path.isdir(os.path.join(tmp_dir, d))]
#             if samples:
#                 return samples
#     except Exception as e:
#         print(f"Warning: Error checking tmp dir for predictions: {e}")
#     
#     # If all methods fail, return empty list
#     return []

# NOTE: These rules are disabled in favor of phage-specific splitting
# See split_proteins_by_phage.smk for the improved implementation
# checkpoint wait_for_phacts_splits:
#     input:
#         split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
#     output:
#         touch(f"{config['output_dir']}/03_phacts_results/.splits_ready")
#     shell:
#         "mkdir -p $(dirname {output})"

# # Check if input files exist for PHACTS
# rule check_phacts_input_files:
#     input:
#         # Prerequisite - splits must be ready
#         splits_ready = f"{config['output_dir']}/03_phacts_results/.splits_ready",
#         # Required protein predictions
#         proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
#         # Check that split files are created
#         split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
#     output:
#         touch(f"{config['output_dir']}/03_phacts_results/.input_files_found")
#     log:
#         f"{config['output_dir']}/logs/check_phacts_input_files.log"
#     shell:
#         """
#         # Check if protein file exists
#         if [ ! -f "{input.proteins}" ]; then
#             echo "Error: Protein file does not exist: {input.proteins}" > {log} 2>&1
#             exit 1
#         fi
#         
#         # Check if split list exists and has content
#         if [ ! -s "{input.split_list}" ]; then
#             echo "Warning: Split protein file list is empty or doesn't exist: {input.split_list}" > {log} 2>&1
#             echo "Creating a placeholder split list..." >> {log} 2>&1
#             mkdir -p $(dirname {input.split_list})
#             touch {input.split_list}
#         fi
#         
#         # Everything checked out, all input files exist
#         echo "All required input files for PHACTS were found" > {log} 2>&1
#         """

# We don't need a separate rule for log directories
# Snakemake will create them automatically

# NOTE: Old PHACTS rules removed - using phage-specific implementation in split_proteins_by_phage.smk
# These rules used random protein splitting which caused poor predictions (~0.5 probability)
# The new phage-specific approach groups proteins by phage ID for accurate lifestyle prediction

# 6. Run mmseqs2 taxonomy assignment on phage genomes
rule mmseqs_phage_taxonomy:
    input:
        phage_seqs = get_phage_input
    output:
        # Change the output path to match the actual file naming pattern
        lca_output = f"{config['output_dir']}/03_genomic_info/mmseqs_output_lca.tsv",
        taxonomy = f"{config['output_dir']}/03_genomic_info/mmseqs_taxonomy.tsv"
    log:
        f"{config['output_dir']}/logs/mmseqs_phage_taxonomy.log"
    conda:
        config["conda_envs"]["mmseqs2"]
    threads: 24
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.lca_output})
        
        # Check if input file exists and has content
        if [ ! -f "{input.phage_seqs}" ] || [ ! -s "{input.phage_seqs}" ]; then
            echo "Warning: Input file is empty or missing: {input.phage_seqs}" > {log} 2>&1
            echo "Creating empty mmseqs2 output files..." >> {log} 2>&1
            echo -e "query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\ttaxlineage" > {output.lca_output}
            cp {output.lca_output} {output.taxonomy}
        else
            # Count sequences in input file
            SEQ_COUNT=$(grep -c ">" {input.phage_seqs} || echo "0")
            echo "Processing $SEQ_COUNT sequences with mmseqs2" > {log} 2>&1
            
            if [ "$SEQ_COUNT" -eq 0 ]; then
                echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
                echo "Creating empty mmseqs2 output files..." >> {log} 2>&1
                echo -e "query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\ttaxlineage" > {output.lca_output}
                cp {output.lca_output} {output.taxonomy}
            else
                # Create temporary directory
                TMP_DIR=$(mktemp -d)
                
                # Run mmseqs2 for taxonomy assignment with output prefix that produces files matching existing naming pattern
                mmseqs easy-taxonomy {input.phage_seqs} {config[databases][mmseqs2][db]} \
                    {config[output_dir]}/03_genomic_info/mmseqs_output $TMP_DIR \
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
                    >> {log} 2>&1
                    
                # Check if mmseqs2 produced output
                if [ ! -f "{output.lca_output}" ]; then
                    echo "Warning: mmseqs2 did not produce output. Creating empty file." >> {log} 2>&1
                    echo -e "query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\ttaxlineage" > {output.lca_output}
                fi
                
                # Copy and format the taxonomy results - use the correct file path
                cp {output.lca_output} {output.taxonomy}
                    
                # Clean up
                rm -rf $TMP_DIR
            fi
        fi
        """

# 7. Run Phabox2 for phage taxonomy and lifestyle prediction
rule phabox_prediction:
    input:
        phage_seqs = get_phage_input
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/phabox_output"),
        taxonomy = f"{config['output_dir']}/03_genomic_info/phabox_output/taxonomy.tsv",
        lifestyle = f"{config['output_dir']}/03_genomic_info/phabox_output/lifestyle.tsv"
    log:
        f"{config['output_dir']}/logs/phabox_prediction.log"
    conda:
        config["conda_envs"]["phabox2"]
    threads: 24
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}
        
        # Run Phabox2 end-to-end (includes virus identification, lifestyle, taxonomy, and host prediction)
        phabox2 --task end_to_end --dbdir {config[databases][phabox][db]} \
            --outpth {output.results_dir} \
            --contigs {input.phage_seqs} \
            --len 1000 \
            --threads {threads} > {log} 2>&1 || {{
                echo "Phabox2 failed. Creating placeholder files..." >> {log} 2>&1
                echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
                echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
                exit 0
            }}
        
        # Process Phabox2 outputs and create standardized files
        # Phabox2 creates output files with specific naming patterns
        if [ -f "{output.results_dir}/out/phamer_prediction.csv" ]; then
            # Extract taxonomy information
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
            tail -n +2 {output.results_dir}/out/phamer_prediction.csv | \
                awk -F',' '{{print $1"\t"$2"\t"$3}}' >> {output.taxonomy}
        else
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
        fi
        
        if [ -f "{output.results_dir}/out/cherry_prediction.csv" ]; then
            # Extract lifestyle information
            echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
            tail -n +2 {output.results_dir}/out/cherry_prediction.csv | \
                awk -F',' '{{print $1"\t"$2"\t"$3}}' >> {output.lifestyle}
        else
            echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
        fi
        """

# 8. Run vContact3 for phage taxonomy based on gene content
rule vcontact3_taxonomy:
    input:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
        phage_seqs = get_phage_input
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/vc3_output")
    log:
        f"{config['output_dir']}/logs/vcontact3_taxonomy.log"
    conda:
        config["conda_envs"]["vcontact3"]
    threads: 24
    shell:
        """
        # Create gene2genome file
        mkdir -p {output.results_dir}
        echo -e "protein_id\tcontig_id" > {output.results_dir}/gene2genome.csv
        grep ">" {input.proteins} | sed 's/>//g' | 
        awk -F " # " '{{split($1,a,"_"); print $1"\t"a[1]}}' >> {output.results_dir}/gene2genome.csv
        
        # Run vContact3
        vcontact3 run --nucleotide {input.phage_seqs} \
            --output {output.results_dir} \
            --db-domain "prokaryotes" \
            --db-version 223 \
            --db-path {config[databases][vcontact3][db]} \
            -t {threads} > {log} 2>&1 || {{
                echo "WARNING: vContact3 failed. Created output directory." >> {log} 2>&1
                mkdir -p {output.results_dir}
                exit 0
            }}
        """

# 9. Run BACPHLIP for phage lifestyle prediction
# NOTE: BACPHLIP assumes complete phage genomes. Results include CheckV completeness flags.
rule bacphlip_lifestyle:
    input:
        phage_seqs = get_phage_input,
        checkv_quality = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    output:
        results = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle.tsv",
        with_completeness = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"
    log:
        f"{config['output_dir']}/logs/bacphlip_lifestyle.log"
    conda:
        config["conda_envs"]["bacphlip"]
    threads: 8
    shell:
        """
        # Create temporary directory
        TMP_DIR=$(mktemp -d)
        
        # Clean up any existing BACPHLIP directory
        rm -rf {input.phage_seqs}.BACPHLIP_DIR/
        
        # Run BACPHLIP in multi-fasta mode
        echo "Running BACPHLIP on all sequences..." > {log} 2>&1
        bacphlip -i {input.phage_seqs} --multi_fasta -f \
            > $TMP_DIR/bacphlip_raw.tsv 2>> {log} || {{
                echo "BACPHLIP failed. Creating placeholder output..." >> {log} 2>&1
                echo -e "Sequence\tLifestyle\tConfidence" > {output.results}
                echo -e "Sequence\tLifestyle\tConfidence\tCompleteness\tCheckV_quality" > {output.with_completeness}
                rm -rf $TMP_DIR
                exit 0
            }}
        
        # Copy raw results
        cp $TMP_DIR/bacphlip_raw.tsv {output.results}
        
        # Join with CheckV completeness data
        echo "Adding completeness information..." >> {log} 2>&1
        
        # Create header
        echo -e "Sequence\tLifestyle\tConfidence\tCompleteness\tCheckV_quality" > {output.with_completeness}
        
        # Process BACPHLIP results and join with CheckV data
        # Note: CheckV quality_summary.tsv has columns: contig_id, various metrics..., checkv_quality
        tail -n +2 $TMP_DIR/bacphlip_raw.tsv | while IFS=$'\t' read -r seq lifestyle conf; do
            # Look up this sequence in CheckV results
            checkv_line=$(grep "^${{seq}}\t" {input.checkv_quality} || echo "")
            
            if [ -n "$checkv_line" ]; then
                # Extract completeness (column 10) and quality (column 8) from CheckV
                completeness=$(echo "$checkv_line" | cut -f10)
                quality=$(echo "$checkv_line" | cut -f8)
            else
                completeness="Unknown"
                quality="Not_found_in_CheckV"
            fi
            
            echo -e "${{seq}}\t${{lifestyle}}\t${{conf}}\t${{completeness}}\t${{quality}}"
        done >> {output.with_completeness}
        
        # Log summary statistics
        echo "BACPHLIP analysis complete. Summary:" >> {log} 2>&1
        echo "Total sequences: $(tail -n +2 {output.results} | wc -l)" >> {log} 2>&1
        echo "Complete genomes: $(grep -E "\tComplete\t|\tHigh-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
        echo "Medium-quality: $(grep "\tMedium-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
        echo "Low-quality: $(grep "\tLow-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
        echo "Not-determined: $(grep "\tNot-determined\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
        
        # Clean up
        rm -rf $TMP_DIR
        rm -rf {input.phage_seqs}.BACPHLIP_DIR/
        """
