# Simple, robust PHACTS integration that relies on standalone scripts

# Make sure that we're also including the prerequisites from the original rule
# to maintain compatibility with existing workflows
rule phacts_simple_prediction:
    input:
        # Include all the prerequisites from the original rule for compatibility
        splits_ready = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        input_check = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa"
    output:
        result_dir = directory(f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}"),
        result = f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}/{{sample}}.phacts.out"
    log:
        f"{config['output_dir']}/logs/phacts_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        r"""
        # Create output directory
        mkdir -p {output.result_dir}
        
        # Use wrapper script for reliable execution
        WRAPPER_SCRIPT="$(dirname $(dirname {workflow.basedir}))/scripts/run_phacts.sh"
        
        # Check if wrapper script exists
        if [ ! -f "$WRAPPER_SCRIPT" ]; then
            echo "Error: PHACTS wrapper script not found at $WRAPPER_SCRIPT" > {log} 2>&1
            exit 1
        fi
        
        # Run PHACTS using the wrapper script
        "$WRAPPER_SCRIPT" {input.protein_file} -o {output.result_dir} --debug >> {log} 2>&1
        
        # Rename the output file to match expected format
        if [ -f "{output.result_dir}/prediction.txt" ]; then
            cp {output.result_dir}/prediction.txt {output.result}
            echo "Prediction completed successfully" >> {log} 2>&1
        else
            echo "Warning: PHACTS did not produce valid output. Creating placeholder file." >> {log} 2>&1
            echo "No prediction was made for this batch" > {output.result}
        fi
        """

# We'll reuse the original aggregation rule, but make sure our rule name is still the original one
# from the main workflow to avoid conflicts there too
rule phacts_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        all_done = f"{config['output_dir']}/03_phacts_results/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}/{{sample}}.phacts.out",
            sample=get_phacts_samples()
        )
    output:
        predictions = f"{config['output_dir']}/03_phacts_results/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_aggregate_results.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        r"""
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHACTS results" > {log}
        echo -e "phage_id\tlifestyle\tprobability" > {output.predictions}

        if [ -n "$(ls -A $RESULTS_DIR/tmp 2>/dev/null)" ]; then
            # Process each phacts output file if it exists
            for file in $RESULTS_DIR/tmp/*/*.phacts.out 2>/dev/null; do
                if [ -f "$file" ]; then
                    # Check if the file contains useful predictions
                    if grep -q "Lifestyle:" "$file" && grep -q "Probability:" "$file"; then
                        # Extract the contig information from the PHACTS output
                        # Look for lines containing >contig identifiers in the first part of the file
                        phage_ids=$(grep "^>" "$file" | head -n 10 | sed 's/^>//g' | cut -d "_" -f 1 | sort -u)
                        lifestyle=$(grep "Lifestyle:" "$file" | awk '{{print $2}}')
                        probability=$(grep "Probability:" "$file" | awk '{{print $2}}')
                        
                        # Add prediction for each phage in the batch
                        for phage_id in $phage_ids; do
                            echo -e "$phage_id\t$lifestyle\t$probability" >> {output.predictions}
                        done
                    else
                        # Extract batch number from filename
                        batch_id=$(basename "$file" .phacts.out)
                        echo "Warning: No valid prediction in $batch_id" >> {log}
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
        """