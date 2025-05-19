# Split protein files by phage ID for PHACTS array processing
rule split_proteins_by_phage:
    input:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
    output:
        split_dir = directory(f"{config['output_dir']}/03_split_proteins_by_phage"),
        split_list = f"{config['output_dir']}/03_split_proteins_by_phage/split_protein_list.txt"
    log:
        f"{config['output_dir']}/logs/split_proteins_by_phage.log"
    conda:
        config["conda_envs"]["phacts"]  # This env has biopython
    shell:
        """
        # Run the Python script to split proteins by phage ID
        python {workflow.basedir}/scripts/split_proteins_by_phage.py \
            --input {input.proteins} \
            --output-dir {output.split_dir} \
            --split-list {output.split_list} \
            --show-examples > {log} 2>&1
        
        # If no files were created (empty input), create an empty placeholder
        if [ ! -s "{output.split_list}" ]; then
            echo "Input was empty, creating placeholder file" >> {log} 2>&1
            touch "{output.split_dir}/empty.faa"
            echo "{output.split_dir}/empty.faa" > {output.split_list}
        fi
        """

# List all samples for phacts from split protein files (by phage ID)
def get_phacts_phage_samples():
    # After split_proteins_by_phage is run, this reads the split file list
    split_list = f"{config['output_dir']}/03_split_proteins_by_phage/split_protein_list.txt"
    
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
    
    # Fallback method 1: Use glob to find faa files directly
    try:
        split_dir = f"{config['output_dir']}/03_split_proteins_by_phage"
        if os.path.exists(split_dir) and os.path.isdir(split_dir):
            samples = [os.path.splitext(os.path.basename(f))[0] 
                    for f in glob.glob(f"{split_dir}/*.faa")]
            # Filter out any empty strings
            samples = [s for s in samples if s]
            if samples:
                return samples
    except Exception as e:
        print(f"Warning: Error using glob to find faa files: {e}")
    
    # Fallback method 2: Check tmp dir for existing predictions
    try:
        tmp_dir = f"{config['output_dir']}/03_phacts_results_by_phage/tmp"
        if os.path.exists(tmp_dir) and os.path.isdir(tmp_dir):
            # Look for directories in the tmp directory
            samples = [d for d in os.listdir(tmp_dir)
                      if os.path.isdir(os.path.join(tmp_dir, d))]
            if samples:
                return samples
    except Exception as e:
        print(f"Warning: Error checking tmp dir for predictions: {e}")
    
    # If all methods fail, return empty list
    return []

checkpoint wait_for_phacts_phage_splits:
    input:
        split_list = f"{config['output_dir']}/03_split_proteins_by_phage/split_protein_list.txt"
    output:
        touch(f"{config['output_dir']}/03_phacts_results_by_phage/.splits_ready")
    shell:
        "mkdir -p $(dirname {output})"

# Check if input files exist for PHACTS (phage-specific version)
rule check_phacts_phage_input_files:
    input:
        # Prerequisite - splits must be ready
        splits_ready = f"{config['output_dir']}/03_phacts_results_by_phage/.splits_ready",
        # Required protein predictions
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
        # Check that split files are created
        split_list = f"{config['output_dir']}/03_split_proteins_by_phage/split_protein_list.txt"
    output:
        touch(f"{config['output_dir']}/03_phacts_results_by_phage/.input_files_found")
    log:
        f"{config['output_dir']}/logs/check_phacts_phage_input_files.log"
    shell:
        """
        # Check if protein file exists
        if [ ! -f "{input.proteins}" ]; then
            echo "Error: Protein file does not exist: {input.proteins}" > {log} 2>&1
            exit 1
        fi
        
        # Check if split list exists and has content
        if [ ! -s "{input.split_list}" ]; then
            echo "Warning: Split protein file list is empty or doesn't exist: {input.split_list}" > {log} 2>&1
            echo "Creating a placeholder split list..." >> {log} 2>&1
            mkdir -p $(dirname {input.split_list})
            touch {input.split_list}
        fi
        
        # Everything checked out, all input files exist
        echo "All required input files for PHACTS (phage-specific) were found" > {log} 2>&1
        """

# Run PHACTS lifestyle prediction on phage-specific protein files
rule phacts_phage_prediction:
    input:
        checkpoint = f"{config['output_dir']}/03_phacts_results_by_phage/.splits_ready",
        input_check = f"{config['output_dir']}/03_phacts_results_by_phage/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins_by_phage/{{sample}}.faa"
    output:
        result_dir = directory(f"{config['output_dir']}/03_phacts_results_by_phage/tmp/{{sample}}"),
        result = f"{config['output_dir']}/03_phacts_results_by_phage/tmp/{{sample}}/{{sample}}.phacts.out"
    log:
        f"{config['output_dir']}/logs/phacts_phage_prediction/{{sample}}.log"
    threads: 4
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Set -e to stop script on first error, redirecting both stdout and stderr to the log
        set -e
        echo "Starting PHACTS rule for phage {wildcards.sample}" > {log} 2>&1

        # Create output directory
        mkdir -p {output.result_dir} >> {log} 2>&1

        # Define the path to phacts manually (change to match the updated value in 03_analysis.smk)
        phacts_path=/home/megan.j/PHACTS

        # Add phacts directory to PYTHONPATH
        export PYTHONPATH="${{phacts_path}}:${{PYTHONPATH:-}}"
        echo "PYTHONPATH is set to: $PYTHONPATH" >> {log} 2>&1

        # Run PHACTS 
        echo "Running PHACTS on phage: {wildcards.sample}" >> {log} 2>&1
        python /home/megan.j/PHACTS/phacts.py {input.protein_file} -o {output.result_dir}/prediction.txt >> {log} 2>&1

        # Debugging: List files in result_dir 
        echo "Listing files in {output.result_dir}:" >> {log} 2>&1
        ls -l {output.result_dir} >> {log} 2>&1

        # Check if output file exists and rename it to match the expected result name
        if [ -f "{output.result_dir}/prediction.txt" ]; then
            mv {output.result_dir}/prediction.txt {output.result}
            echo "Prediction completed successfully" >> {log} 2>&1
        else
            echo "PHACTS failed to produce output. Creating placeholder file." >> {log} 2>&1
            echo "No prediction was made for phage {wildcards.sample}" > {output.result}
        fi
        """

# Helper rule to force running all phacts predictions (phage-specific version)
rule run_all_phacts_phage_predictions:
    input:
        checkpoint = f"{config['output_dir']}/03_phacts_results_by_phage/.splits_ready",
        input_check = f"{config['output_dir']}/03_phacts_results_by_phage/.input_files_found",
        # For actual runs, get samples from the split files
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results_by_phage/tmp/{{sample}}/{{sample}}.phacts.out",
            sample=get_phacts_phage_samples()
        )
    output:
        touch(f"{config['output_dir']}/03_phacts_results_by_phage/.all_predictions_done")

# Aggregate PHACTS results (phage-specific version)
rule phacts_phage_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        all_done = f"{config['output_dir']}/03_phacts_results_by_phage/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results_by_phage/tmp/{{sample}}/{{sample}}.phacts.out",
            sample=get_phacts_phage_samples()
        )
    output:
        predictions = f"{config['output_dir']}/03_phacts_results_by_phage/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_phage_aggregate_results.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHACTS results (phage-specific)" > {log}
        echo -e "phage_id\tlifestyle\tprobability" > {output.predictions}

        if [ -n "$(ls -A $RESULTS_DIR/tmp 2>/dev/null)" ]; then
            # Process each phacts output file if it exists
            for file in $RESULTS_DIR/tmp/*/*.phacts.out 2>/dev/null; do
                if [ -f "$file" ]; then
                    # Check if the file contains useful predictions
                    if grep -q "Lifestyle:" "$file" && grep -q "Probability:" "$file"; then
                        # Extract phage ID from the filename (it's already phage-specific)
                        phage_id=$(basename "$(dirname "$file")")
                        lifestyle=$(grep "Lifestyle:" "$file" | awk '{{print $2}}')
                        probability=$(grep "Probability:" "$file" | awk '{{print $2}}')
                        
                        # Each file corresponds to a single phage
                        echo -e "$phage_id\t$lifestyle\t$probability" >> {output.predictions}
                    else
                        # Extract phage ID from filename
                        phage_id=$(basename "$(dirname "$file")")
                        echo "Warning: No valid prediction for phage $phage_id" >> {log}
                    fi
                fi
            done

            if [ "$(wc -l < {output.predictions})" -eq 1 ]; then
                echo "No valid PHACTS results were found, only header in output file" >> {log}
            else
                echo "Successfully compiled $(( $(wc -l < {output.predictions}) - 1 )) PHACTS results" >> {log}
            fi
        else
            echo "No PHACTS result files found, created empty result file" >> {log}
        fi

        # Clean up temporary directory after successful compilation
        if [ -d "$RESULTS_DIR/tmp" ]; then
            echo "Cleaning up temporary directory" >> {log}
            rm -rf "$RESULTS_DIR/tmp"
        fi
        """