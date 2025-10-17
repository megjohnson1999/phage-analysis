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

        # Collect both genome-level and genus-level predictions (prioritize genome-level)
        find "$RESULTS_DIR/tmp" -name "Host_prediction_to_genome_m90.csv" -type f > "$TMP_DIR/genome_files.txt"
        find "$RESULTS_DIR/tmp" -name "host_prediction_to_genus.csv" -type f > "$TMP_DIR/genus_files.txt"

        # Count how many files we found
        GENOME_COUNT=$(wc -l < "$TMP_DIR/genome_files.txt")
        GENUS_COUNT=$(wc -l < "$TMP_DIR/genus_files.txt")
        echo "Found $GENOME_COUNT genome-level prediction files and $GENUS_COUNT genus-level prediction files" >> {log} 2>&1

        # Start with genome-level predictions (more specific)
        TOTAL_RECORDS=0
        if [ "$GENOME_COUNT" -gt 0 ]; then
            # Get the first genome file to extract header
            FIRST_FILE=$(head -n 1 "$TMP_DIR/genome_files.txt")

            # Create the output file with header (convert CSV to TSV)
            head -n 1 "$FIRST_FILE" | tr ',' '\t' > {output.predictions}

            # Process genome-level files
            while read -r pred_file; do
                # Skip header line (first line) from each file and convert CSV to TSV
                awk -F ',' 'NR>1 && $1!="" {{OFS="\t"; print}}' "$pred_file" >> "$TMP_DIR/aggregated_data.tmp"
            done < "$TMP_DIR/genome_files.txt"

            # Count genome-level records
            if [ -f "$TMP_DIR/aggregated_data.tmp" ]; then
                GENOME_RECORDS=$(wc -l < "$TMP_DIR/aggregated_data.tmp")
                echo "Added $GENOME_RECORDS genome-level host predictions" >> {log} 2>&1
                TOTAL_RECORDS=$GENOME_RECORDS
            fi

        elif [ "$GENUS_COUNT" -gt 0 ]; then
            # Fallback to genus-level if no genome-level found
            FIRST_FILE=$(head -n 1 "$TMP_DIR/genus_files.txt")
            head -n 1 "$FIRST_FILE" | tr ',' '\t' > {output.predictions}

            # Process genus-level files
            while read -r pred_file; do
                awk -F ',' 'NR>1 && $1!="" {{OFS="\t"; print}}' "$pred_file" >> "$TMP_DIR/aggregated_data.tmp"
            done < "$TMP_DIR/genus_files.txt"

            if [ -f "$TMP_DIR/aggregated_data.tmp" ]; then
                GENUS_RECORDS=$(wc -l < "$TMP_DIR/aggregated_data.tmp")
                echo "Added $GENUS_RECORDS genus-level host predictions" >> {log} 2>&1
                TOTAL_RECORDS=$GENUS_RECORDS
            fi
        fi

        # Append all data to the output file
        if [ -f "$TMP_DIR/aggregated_data.tmp" ] && [ -s "$TMP_DIR/aggregated_data.tmp" ]; then
            cat "$TMP_DIR/aggregated_data.tmp" >> {output.predictions}
            echo "Successfully compiled iPhop results with $TOTAL_RECORDS data records" >> {log} 2>&1
        else
            # Create empty output with header structure if no predictions found
            echo -e "Virus\tHost_genome\tHost_taxonomy\tMain_method\tConfidence_score\tAdditional_methods" > {output.predictions}
            echo "No iPhop predictions found, created empty file with header" >> {log} 2>&1
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
        
        # Check if input file exists and has content
        if [ ! -f "{input.phage_seqs}" ] || [ ! -s "{input.phage_seqs}" ]; then
            echo "Warning: Input file is empty or missing: {input.phage_seqs}" > {log} 2>&1
            echo "Creating empty Prodigal output files..." >> {log} 2>&1
            touch {output.proteins}
            touch {output.genes}
        else
            # Count sequences in input file
            SEQ_COUNT=$(grep -c ">" {input.phage_seqs} || echo "0")
            echo "Processing $SEQ_COUNT sequences with Prodigal" > {log} 2>&1
            
            if [ "$SEQ_COUNT" -eq 0 ]; then
                echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
                echo "Creating empty Prodigal output files..." >> {log} 2>&1
                touch {output.proteins}
                touch {output.genes}
            else
                # Run Prodigal for ORF prediction
                prodigal -i {input.phage_seqs} \
                    -a {output.proteins} \
                    -d {output.genes} \
                    -p meta >> {log} 2>&1
                
                # Check if outputs were created
                if [ ! -f "{output.proteins}" ]; then
                    echo "Warning: Prodigal did not create protein output. Creating empty file." >> {log} 2>&1
                    touch {output.proteins}
                fi
                if [ ! -f "{output.genes}" ]; then
                    echo "Warning: Prodigal did not create gene output. Creating empty file." >> {log} 2>&1
                    touch {output.genes}
                fi
            fi
        fi
        """


# 6. Run mmseqs2 taxonomy assignment on phage genomes
rule mmseqs_phage_taxonomy:
    input:
        phage_seqs = get_phage_input
    output:
        # Change the output path to match the actual file naming pattern
        lca_output = f"{config['output_dir']}/03_genomic_info/mmseqs_output_lca.tsv",
        tophit_output = f"{config['output_dir']}/03_genomic_info/mmseqs_output_tophit_aln",
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
                
                # Copy the tophit_aln file (with taxlineage) for taxonomic consensus
                if [ -f "{output.tophit_output}" ]; then
                    cp {output.tophit_output} {output.taxonomy}
                else
                    echo "Warning: tophit_aln file not found, using LCA output" >> {log} 2>&1
                    cp {output.lca_output} {output.taxonomy}
                fi
                    
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
        
        # Check if input file exists and has content
        if [ ! -f "{input.phage_seqs}" ] || [ ! -s "{input.phage_seqs}" ]; then
            echo "Warning: Input file is empty or missing: {input.phage_seqs}" > {log} 2>&1
            echo "Creating empty Phabox2 output files..." >> {log} 2>&1
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
            echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
        else
            # Count sequences in input file
            SEQ_COUNT=$(grep -c ">" {input.phage_seqs} || echo "0")
            echo "Processing $SEQ_COUNT sequences with Phabox2" > {log} 2>&1
            
            if [ "$SEQ_COUNT" -eq 0 ]; then
                echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
                echo "Creating empty Phabox2 output files..." >> {log} 2>&1
                echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
                echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
            else
                # Run Phabox2 end-to-end (includes virus identification, lifestyle, taxonomy, and host prediction)
                phabox2 --task end_to_end --dbdir {config[databases][phabox][db]} \
                    --outpth {output.results_dir} \
                    --contigs {input.phage_seqs} \
                    --len 1000 \
                    --threads {threads} >> {log} 2>&1 || {{
                        echo "Phabox2 failed. Creating placeholder files..." >> {log} 2>&1
                        echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
                        echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
                        exit 0
                    }}
            fi
        fi
        
        # Process Phabox2 outputs and create standardized files
        # Handle both new format (final_prediction_summary.tsv) and legacy format
        
        # Check for new comprehensive format first
        if [ -f "{output.results_dir}/final_prediction/final_prediction_summary.tsv" ]; then
            echo "Processing Phabox2 new format: final_prediction_summary.tsv" >> {log}
            
            # Extract taxonomy information from new format
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
            tail -n +2 {output.results_dir}/final_prediction/final_prediction_summary.tsv | \
                awk -F'\t' '$3=="virus" {{
                    # Extract taxonomy from Lineage column (column 7)
                    # Convert confidence score from PhaGCNScore (column 8) 
                    lineage = $7; 
                    if (lineage == "-" || lineage == "") lineage = "unclassified";
                    conf = $8;
                    if (conf == "" || conf == "-") conf = "0.0";
                    print $1"\t"lineage"\t"conf
                }}' >> {output.taxonomy}
            
            # Extract lifestyle information from new format  
            echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
            tail -n +2 {output.results_dir}/final_prediction/final_prediction_summary.tsv | \
                awk -F'\t' '$3=="virus" {{
                    # Extract lifestyle from TYPE column (column 11)
                    # Use PhaTYPScore for confidence (column 12)
                    lifestyle = $11;
                    if (lifestyle == "-" || lifestyle == "") lifestyle = "unknown";
                    conf = $12;
                    if (conf == "" || conf == "-") conf = "0.0";
                    print $1"\t"lifestyle"\t"conf
                }}' >> {output.lifestyle}
        
        # Fall back to legacy format processing
        elif [ -f "{output.results_dir}/out/phamer_prediction.csv" ]; then
            echo "Processing Phabox2 legacy format: phamer_prediction.csv" >> {log}
            # Extract taxonomy information
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
            tail -n +2 {output.results_dir}/out/phamer_prediction.csv | \
                awk -F',' '{{print $1"\t"$2"\t"$3}}' >> {output.taxonomy}
        else
            echo "No Phabox2 taxonomy results found, creating empty file" >> {log}
            echo -e "contig_id\ttaxonomy_prediction\tconfidence" > {output.taxonomy}
        fi
        
        # Handle lifestyle predictions - only if we didn't already process new format
        if [ ! -f "{output.results_dir}/final_prediction/final_prediction_summary.tsv" ]; then
            if [ -f "{output.results_dir}/out/cherry_prediction.csv" ]; then
                echo "Processing Phabox2 legacy format: cherry_prediction.csv" >> {log}
                # Extract lifestyle information
                echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
                tail -n +2 {output.results_dir}/out/cherry_prediction.csv | \
                    awk -F',' '{{print $1"\t"$2"\t"$3}}' >> {output.lifestyle}
            else
                echo "No Phabox2 lifestyle results found, creating empty file" >> {log}
                echo -e "contig_id\tlifestyle_prediction\tconfidence" > {output.lifestyle}
            fi
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
        # Create output directory
        mkdir -p {output.results_dir}
        
        # Check if input files exist and have content
        if [ ! -f "{input.phage_seqs}" ] || [ ! -s "{input.phage_seqs}" ] || [ ! -f "{input.proteins}" ] || [ ! -s "{input.proteins}" ]; then
            echo "Warning: Input files are empty or missing" > {log} 2>&1
            echo "phage_seqs: {input.phage_seqs}" >> {log} 2>&1
            echo "proteins: {input.proteins}" >> {log} 2>&1
            echo "Skipping vContact3 analysis..." >> {log} 2>&1
            # Create empty gene2genome file
            echo -e "protein_id\tcontig_id" > {output.results_dir}/gene2genome.csv
        else
            # Count sequences in input files
            SEQ_COUNT=$(grep -c ">" {input.phage_seqs} || echo "0")
            PROT_COUNT=$(grep -c ">" {input.proteins} || echo "0")
            echo "Processing $SEQ_COUNT sequences and $PROT_COUNT proteins with vContact3" > {log} 2>&1
            
            if [ "$SEQ_COUNT" -eq 0 ] || [ "$PROT_COUNT" -eq 0 ]; then
                echo "Warning: No sequences or proteins to process" >> {log} 2>&1
                echo "Skipping vContact3 analysis..." >> {log} 2>&1
                # Create empty gene2genome file
                echo -e "protein_id\tcontig_id" > {output.results_dir}/gene2genome.csv
            else
                # Create gene2genome file
                echo -e "protein_id\tcontig_id" > {output.results_dir}/gene2genome.csv
                grep ">" {input.proteins} | sed 's/>//g' | 
                awk -F " # " '{{split($1,a,"_"); print $1"\t"a[1]}}' >> {output.results_dir}/gene2genome.csv
                
                # Run vContact3
                vcontact3 run --nucleotide {input.phage_seqs} \
                    --output {output.results_dir} \
                    --db-domain "prokaryotes" \
                    --db-version 223 \
                    --db-path {config[databases][vcontact3][db]} \
                    -t {threads} >> {log} 2>&1 || {{
                        echo "WARNING: vContact3 failed. Created output directory." >> {log} 2>&1
                        mkdir -p {output.results_dir}
                        exit 0
                    }}
            fi
        fi
        """

# 9. Run BACPHLIP for phage lifestyle prediction
# NOTE: BACPHLIP assumes complete phage genomes. Results include CheckV completeness flags.
rule bacphlip_lifestyle:
    input:
        phage_seqs = get_phage_input,
        checkv_quality = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/bacphlip_output"),
        results = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle.tsv",
        with_completeness = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"
    log:
        f"{config['output_dir']}/logs/bacphlip_lifestyle.log"
    conda:
        config["conda_envs"]["bacphlip"]
    threads: 8
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}

        # Check if input file exists and has content
        if [ ! -f "{input.phage_seqs}" ] || [ ! -s "{input.phage_seqs}" ]; then
            echo "Warning: Input file is empty or missing: {input.phage_seqs}" > {log} 2>&1
            echo "Creating empty BACPHLIP output files..." >> {log} 2>&1
            echo -e "Sequence\tVirulent\tTemperate" > {output.results}
            echo -e "Sequence\tVirulent\tTemperate\tCompleteness\tCheckV_quality" > {output.with_completeness}
        else
            # Count sequences in input file
            SEQ_COUNT=$(grep -c ">" {input.phage_seqs} || echo "0")
            echo "Processing $SEQ_COUNT sequences with BACPHLIP" > {log} 2>&1

            if [ "$SEQ_COUNT" -eq 0 ]; then
                echo "Warning: Input file contains 0 sequences" >> {log} 2>&1
                echo "Creating empty BACPHLIP output files..." >> {log} 2>&1
                echo -e "Sequence\tVirulent\tTemperate" > {output.results}
                echo -e "Sequence\tVirulent\tTemperate\tCompleteness\tCheckV_quality" > {output.with_completeness}
            else
                # Create temporary directory
                TMP_DIR=$(mktemp -d)

                # Clean up any existing BACPHLIP directories
                rm -rf {input.phage_seqs}.BACPHLIP_DIR/
                rm -rf {input.phage_seqs}.bacphlip/

                # Run BACPHLIP in multi-fasta mode
                echo "Running BACPHLIP on all sequences..." >> {log} 2>&1
                bacphlip -i {input.phage_seqs} --multi_fasta -f \
                    > $TMP_DIR/bacphlip_raw.tsv 2>> {log} || {{
                        echo "BACPHLIP failed. Creating placeholder output..." >> {log} 2>&1
                        echo -e "Sequence\tVirulent\tTemperate" > {output.results}
                        echo -e "Sequence\tVirulent\tTemperate\tCompleteness\tCheckV_quality" > {output.with_completeness}
                        rm -rf $TMP_DIR
                        exit 0
                    }}

                # Copy the auto-generated BACPHLIP directory to our organized output location
                # BACPHLIP creates a directory at {input}.bacphlip with intermediate results
                if [ -d "{input.phage_seqs}.bacphlip" ]; then
                    echo "Copying BACPHLIP output directory to organized location..." >> {log} 2>&1
                    cp -r {input.phage_seqs}.bacphlip/* {output.results_dir}/ 2>> {log} || true
                    # Clean up the auto-generated directory in the clustering folder
                    rm -rf {input.phage_seqs}.bacphlip/
                fi

                # Copy raw results (stdout output from BACPHLIP)
                cp $TMP_DIR/bacphlip_raw.tsv {output.results}

                # Join with CheckV completeness data
                echo "Adding completeness information..." >> {log} 2>&1

                # Create header - keep original format with Virulent and Temperate confidence scores
                echo -e "Sequence\tVirulent\tTemperate\tCompleteness\tCheckV_quality" > {output.with_completeness}

                # Process BACPHLIP results and join with CheckV data
                # BACPHLIP output format: Sequence\tVirulent_confidence\tTemperate_confidence
                tail -n +1 $TMP_DIR/bacphlip_raw.tsv | while IFS=$'\t' read -r seq virulent_conf temperate_conf; do
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

                    echo -e "${{seq}}\t${{virulent_conf}}\t${{temperate_conf}}\t${{completeness}}\t${{quality}}"
                done >> {output.with_completeness}

                # Log summary statistics
                echo "BACPHLIP analysis complete. Summary:" >> {log} 2>&1
                echo "Total sequences: $(wc -l < $TMP_DIR/bacphlip_raw.tsv)" >> {log} 2>&1
                echo "Complete genomes: $(grep -E "\tComplete\t|\tHigh-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
                echo "Medium-quality: $(grep "\tMedium-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
                echo "Low-quality: $(grep "\tLow-quality\t" {output.with_completeness} | wc -l)" >> {log} 2>&1
                echo "Not-determined: $(grep "\tNot-determined\t" {output.with_completeness} | wc -l)" >> {log} 2>&1

                # Clean up
                rm -rf $TMP_DIR
                rm -rf {input.phage_seqs}.BACPHLIP_DIR/
            fi
        fi
        """

# 10. Create lifestyle consensus from BACPHLIP and Phabox2
rule lifestyle_consensus:
    input:
        bacphlip = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle.tsv",
        phabox = f"{config['output_dir']}/03_genomic_info/phabox_output/lifestyle.tsv"
    output:
        consensus = f"{config['output_dir']}/03_genomic_info/lifestyle_consensus.tsv"
    log:
        f"{config['output_dir']}/logs/lifestyle_consensus.log"
    conda:
        config["conda_envs"]["python"]
    params:
        threshold = 0.7  # Minimum confidence for BACPHLIP predictions
    shell:
        """
        echo "Creating lifestyle consensus..." > {log} 2>&1

        # Run the lifestyle consensus script
        python {workflow.basedir}/scripts/lifestyle_consensus.py \
            --bacphlip {input.bacphlip} \
            --phabox {input.phabox} \
            --output {output.consensus} \
            --threshold {params.threshold} \
            >> {log} 2>&1

        # Log results summary
        if [ -f {output.consensus} ]; then
            CONSENSUS_COUNT=$(tail -n +2 {output.consensus} | wc -l)
            echo "Successfully created lifestyle consensus for $CONSENSUS_COUNT contigs" >> {log} 2>&1

            # Count predictions by source
            echo "Predictions by source:" >> {log} 2>&1
            tail -n +2 {output.consensus} | cut -f4 | sort | uniq -c >> {log} 2>&1

            # Count predictions by lifestyle
            echo "Predictions by lifestyle:" >> {log} 2>&1
            tail -n +2 {output.consensus} | cut -f2 | sort | uniq -c >> {log} 2>&1
        else
            echo "ERROR: Lifestyle consensus file was not created" >> {log} 2>&1
            exit 1
        fi
        """

# 11. Create taxonomic consensus from multiple tools using R script with taxonomizr
rule taxonomic_consensus:
    input:
        mmseqs_raw = f"{config['output_dir']}/03_genomic_info/mmseqs_taxonomy.tsv",
        phabox_taxonomy = f"{config['output_dir']}/03_genomic_info/phabox_output/taxonomy.tsv",
        phabox_lifestyle = f"{config['output_dir']}/03_genomic_info/phabox_output/lifestyle.tsv",
        vcontact3_dir = f"{config['output_dir']}/03_genomic_info/vc3_output"
        # crassus_taxonomy would be added here when CrassUS integration is implemented
    output:
        consensus_taxonomy = f"{config['output_dir']}/03_genomic_info/consensus_taxonomy.tsv",
        consensus_summary = f"{config['output_dir']}/03_genomic_info/consensus_taxonomy_summary.json"
    log:
        f"{config['output_dir']}/logs/taxonomic_consensus.log"
    conda:
        config["conda_envs"]["r"]
    shell:
        """
        echo "Creating taxonomic consensus using R script with taxonomizr..." > {log} 2>&1
        
        # Check if taxonomizr database exists
        if [ ! -f "{config[databases][taxonomizr][db]}" ]; then
            echo "Error: Taxonomizr database not found: {config[databases][taxonomizr][db]}" >> {log} 2>&1
            echo "Please download the database first using R:" >> {log} 2>&1
            echo "  library(taxonomizr)" >> {log} 2>&1
            echo "  prepareDatabase('{config[databases][taxonomizr][db]}')" >> {log} 2>&1
            exit 1
        fi
        
        # Run the R script that replicates the original workflow with taxonomizr
        Rscript {workflow.basedir}/scripts/taxonomic_consensus.R \
            {input.mmseqs_raw} \
            {input.phabox_taxonomy} \
            {input.phabox_lifestyle} \
            {input.vcontact3_dir} \
            {config[databases][taxonomizr][db]} \
            {output.consensus_taxonomy} \
            {output.consensus_summary} \
            >> {log} 2>&1
        
        # Log results summary
        if [ -f {output.consensus_taxonomy} ]; then
            CONSENSUS_COUNT=$(tail -n +2 {output.consensus_taxonomy} | wc -l)
            echo "Successfully created consensus taxonomy for $CONSENSUS_COUNT contigs using taxonomizr" >> {log} 2>&1
        else
            echo "ERROR: Consensus taxonomy file was not created" >> {log} 2>&1
            exit 1
        fi
        """
