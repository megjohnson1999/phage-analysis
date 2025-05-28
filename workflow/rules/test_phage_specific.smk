# Test rules for phage-specific PHACTS analysis
# These rules allow testing the phage-specific splitting and PHACTS analysis with arbitrary input files

# Rule to test phage-specific splitting with a custom input file
rule test_phage_splitting:
    input:
        protein_file = lambda wildcards: wildcards.input_file
    output:
        split_dir = directory("test_phage_specific/{input_file}/split_proteins"),
        split_list = "test_phage_specific/{input_file}/split_list.txt"
    log:
        "test_phage_specific/{input_file}/logs/test_phage_splitting.log"
    conda:
        config["conda_envs"]["phacts"]  # This env has biopython
    shell:
        """
        # Create output directory
        mkdir -p {output.split_dir}
        
        # Run the Python script to split proteins by phage ID
        python {workflow.basedir}/scripts/split_proteins_by_phage.py \
            --input {input.protein_file} \
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

# Run PHACTS on each split file in test mode
rule test_phage_prediction:
    input:
        protein_file = "test_phage_specific/{input_file}/split_proteins/{phage_id}.faa"
    output:
        result_dir = directory("test_phage_specific/{input_file}/phacts_results/tmp/{phage_id}"),
        result = "test_phage_specific/{input_file}/phacts_results/tmp/{phage_id}/{phage_id}.phacts.out"
    log:
        "test_phage_specific/{input_file}/logs/phacts_prediction/{phage_id}.log"
    threads: 1
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Set -e to stop script on first error, redirecting both stdout and stderr to the log
        set -e
        echo "Starting PHACTS rule for phage {wildcards.phage_id}" > {log} 2>&1

        # Create output directory
        mkdir -p {output.result_dir} >> {log} 2>&1

        # Define the path to phacts - use the same as in the main workflow
        phacts_path=/home/megan.j/PHACTS

        # Add phacts directory to PYTHONPATH
        export PYTHONPATH="${{phacts_path}}:${{PYTHONPATH:-}}"
        echo "PYTHONPATH is set to: $PYTHONPATH" >> {log} 2>&1

        # Run PHACTS 
        echo "Running PHACTS on phage: {wildcards.phage_id}" >> {log} 2>&1
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
            echo "No prediction was made for phage {wildcards.phage_id}" > {output.result}
        fi
        """

# Get the list of phage IDs from the split list
def get_test_phage_ids(wildcards):
    import os
    
    # Path to the split list file
    split_list = f"test_phage_specific/{wildcards.input_file}/split_list.txt"
    
    # Check if the split list exists
    if not os.path.exists(split_list):
        return []
    
    # Read the split list
    with open(split_list, "r") as f:
        files = [line.strip() for line in f if line.strip()]
    
    # Extract phage IDs from filenames
    phage_ids = [os.path.basename(file).replace(".faa", "") for file in files]
    
    # Filter out empty strings
    return [pid for pid in phage_ids if pid]

# Aggregate test PHACTS results
rule test_phage_aggregate_results:
    input:
        split_list = "test_phage_specific/{input_file}/split_list.txt",
        predictions = lambda wildcards: expand(
            "test_phage_specific/{input_file}/phacts_results/tmp/{phage_id}/{phage_id}.phacts.out",
            input_file=wildcards.input_file,
            phage_id=get_test_phage_ids(wildcards)
        )
    output:
        predictions = "test_phage_specific/{input_file}/phacts_results/phacts_predictions_compiled.tsv"
    log:
        "test_phage_specific/{input_file}/logs/phacts_aggregate_results.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHACTS results for test input" > {log}
        echo -e "phage_id\tlifestyle\tprobability" > {output.predictions}

        if [ -n "$(ls -A $RESULTS_DIR/tmp 2>/dev/null)" ]; then
            # Process each phacts output file if it exists
            for file in $RESULTS_DIR/tmp/*/*.phacts.out 2>/dev/null; do
                if [ -f "$file" ]; then
                    # Check if the file contains useful predictions
                    if grep -q "Lifestyle:" "$file" && grep -q "Probability:" "$file"; then
                        # Extract phage ID from the filename
                        phage_id=$(basename "$(dirname "$file")")
                        lifestyle=$(grep "Lifestyle:" "$file" | awk '{print $2}')
                        probability=$(grep "Probability:" "$file" | awk '{print $2}')
                        
                        # Each file corresponds to a single phage
                        echo -e "$phage_id\t$lifestyle\t$probability" >> {output.predictions}
                    else
                        # Extract phage ID from filename
                        phage_id=$(basename "$(dirname "$file")")
                        echo "Warning: No valid prediction for phage $phage_id" >> {log}
                    fi
                fi
            done

            LINECOUNT=$(wc -l < {output.predictions})
            if [ "$LINECOUNT" -eq 1 ]; then
                echo "No valid PHACTS results were found, only header in output file" >> {log}
            else
                RESULT_COUNT=$((LINECOUNT - 1))
                echo "Successfully compiled $RESULT_COUNT PHACTS results" >> {log}
            fi
        else
            echo "No PHACTS result files found, created empty result file" >> {log}
        fi
        """

# Complete test rule that runs the entire test pipeline for a given input file
rule test_phage_specific_phacts:
    input:
        protein_file = lambda wildcards: wildcards.input_file,
        split_list = "test_phage_specific/{input_file}/split_list.txt",
        predictions = "test_phage_specific/{input_file}/phacts_results/phacts_predictions_compiled.tsv"
    output:
        report = "test_phage_specific/{input_file}/test_report.txt"
    run:
        import os
        
        # Count the number of split files
        split_dir = f"test_phage_specific/{wildcards.input_file}/split_proteins"
        split_count = len([f for f in os.listdir(split_dir) if f.endswith(".faa")])
        
        # Check the number of predictions
        pred_file = f"test_phage_specific/{wildcards.input_file}/phacts_results/phacts_predictions_compiled.tsv"
        with open(pred_file, "r") as f:
            lines = f.readlines()
            pred_count = len(lines) - 1  # Subtract 1 for the header
        
        # Generate the report
        with open(output.report, "w") as f:
            f.write(f"=== Phage-specific PHACTS Analysis Test Report ===\n\n")
            f.write(f"Input file: {wildcards.input_file}\n")
            f.write(f"Number of phage-specific files created: {split_count}\n")
            f.write(f"Number of PHACTS predictions: {pred_count}\n\n")
            
            # Add the content of the predictions file
            f.write("PHACTS Predictions:\n")
            f.write("".join(lines))
            
            f.write("\n\nTest completed successfully!\n")
            f.write("To view detailed logs, check the logs directory in:\n")
            f.write(f"test_phage_specific/{wildcards.input_file}/logs/\n")