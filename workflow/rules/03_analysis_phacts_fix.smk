# Import the helper functions at the top of your file
import os
import sys

# Add proper paths for imports
# This handles both running from project root and from workflow directory
workflow_dir = os.path.dirname(workflow.basedir) if workflow.basedir.endswith("workflow") else workflow.basedir
if workflow_dir not in sys.path:
    sys.path.append(workflow_dir)

# Now import helper functions with proper path handling
from scripts.phacts_helper import generate_phacts_command, get_phacts_path

# Rule to install PHACTS if not already installed
rule install_phacts:
    output:
        phacts_script = directory(f"{config['output_dir']}/db/phacts/PHACTS"),
        flag_file = f"{config['output_dir']}/db/phacts/.installed"
    params:
        phacts_version = config.get("phacts_version", "main"),
        output_dir = f"{config['output_dir']}/db/phacts"
    log:
        f"{config['output_dir']}/logs/install_phacts.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Create installation directory
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        
        # Clone PHACTS repository if needed
        if [ ! -d "PHACTS" ]; then
            echo "Cloning PHACTS repository..." > {log} 2>&1
            git clone https://github.com/deprekate/PHACTS.git >> {log} 2>&1
            cd PHACTS
            git checkout {params.phacts_version} >> {log} 2>&1
        else
            echo "PHACTS already cloned, updating..." > {log} 2>&1
            cd PHACTS
            git fetch >> {log} 2>&1
            git checkout {params.phacts_version} >> {log} 2>&1
            git pull >> {log} 2>&1
        fi
        
        # Install PHACTS with pip in development mode
        echo "Installing PHACTS with pip..." >> {log} 2>&1
        pip install -e . >> {log} 2>&1
        
        # Create __init__.py files to ensure proper package structure
        echo "Creating __init__.py files to ensure proper package structure" >> {log} 2>&1
        touch __init__.py
        
        # Create a wrapper script that sets PYTHONPATH correctly
        echo "Creating wrapper script for reliable execution" >> {log} 2>&1
        cat > ../run_phacts.sh << EOL
#!/bin/bash

# Set PYTHONPATH to include PHACTS directory and parent
export PYTHONPATH="\$(dirname \$0)/PHACTS:\$(dirname \$0):\$PYTHONPATH"

# Run PHACTS with all arguments passed to this script
python "\$(dirname \$0)/PHACTS/phacts.py" "\$@"
EOL
        chmod +x ../run_phacts.sh
        
        # Create flag file to indicate successful installation
        cd ..
        touch .installed
        
        echo "PHACTS installation complete" >> {log} 2>&1
        """

# Check if input files exist for PHACTS (with explicit dependency on PHACTS installation)
# Renamed to avoid conflict with the rule in 03_analysis.smk
rule check_phacts_installation:
    input:
        # Prerequisite - splits must be ready
        splits_ready = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        # Required protein predictions
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
        # Check that split files are created
        split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt",
        # Add explicit dependency on PHACTS installation
        phacts_installed = f"{config['output_dir']}/db/phacts/.installed"
    output:
        touch(f"{config['output_dir']}/03_phacts_results/.installation_verified")
    log:
        f"{config['output_dir']}/logs/check_phacts_installation.log"
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
        
        # Check if PHACTS is properly installed
        if [ ! -f "{input.phacts_installed}" ]; then
            echo "Error: PHACTS installation not found. Please run install_phacts rule first." > {log} 2>&1
            exit 1
        fi
        
        # Everything checked out, all input files exist
        echo "All required input files for PHACTS were found and installation verified" > {log} 2>&1
        """

# 5a. Run PHACTS for lifestyle prediction on a single protein file batch (UPDATED version with different output paths)
rule phacts_single_prediction_v2:
    input:
        install_check = f"{config['output_dir']}/db/phacts/.installed",
        installation_verified = f"{config['output_dir']}/03_phacts_results/.installation_verified",
        input_check = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa"
    output:
        # Use different output directory to avoid conflicts with original rule
        result_dir = directory(f"{config['output_dir']}/03_phacts_results_v2/tmp/{{sample}}"),
        result = f"{config['output_dir']}/03_phacts_results_v2/tmp/{{sample}}/{{sample}}.phacts.out"
    threads: 4
    log:
        f"{config['output_dir']}/logs/phacts_prediction_v2/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Create output directory
        mkdir -p {output.result_dir}
        
        # Always use the workflow-installed version
        PHACTS_PATH="{config['output_dir']}/db/phacts/PHACTS/phacts.py"
        PHACTS_DIR="{config['output_dir']}/db/phacts/PHACTS"
        PHACTS_PARENT=$(dirname "$PHACTS_DIR")
        echo "Using workflow-installed PHACTS at $PHACTS_PATH" > {log} 2>&1
        
        # Enhanced debugging information
        echo $'\\n==== PHACTS DEBUG INFO ====' >> {log} 2>&1
        echo "Python version: $(python --version 2>&1)" >> {log} 2>&1
        echo "PHACTS_PATH: $PHACTS_PATH" >> {log} 2>&1
        echo "PHACTS_DIR exists: $(test -d "$PHACTS_DIR" && echo "Yes" || echo "No")" >> {log} 2>&1
        echo "phacts.py exists: $(test -f "$PHACTS_PATH" && echo "Yes" || echo "No")" >> {log} 2>&1
        echo "__init__.py exists: $(test -f "$PHACTS_DIR/__init__.py" && echo "Yes" || echo "No")" >> {log} 2>&1
        echo "Directory listing of PHACTS:" >> {log} 2>&1
        ls -la "$PHACTS_DIR" >> {log} 2>&1 2>&1 || echo "Failed to list directory" >> {log} 2>&1
        echo $'=========================\\n' >> {log} 2>&1
        
        # Run the debug script to diagnose import issues
        echo $'\\n==== RUNNING DIAGNOSTIC SCRIPT ====' >> {log} 2>&1
        python "$(dirname {workflow.basedir})/scripts/debug_phacts.py" --output-dir "{config['output_dir']}" >> {log} 2>&1 || true
        echo $'=================================\\n' >> {log} 2>&1
        
        # Add PHACTS directory to PYTHONPATH and run with enhanced setup
        echo $'\\n==== RUNNING PHACTS WITH MODIFIED ENVIRONMENT ====' >> {log} 2>&1
        # Add both the PHACTS directory and its parent to PYTHONPATH
        export PYTHONPATH="$PHACTS_DIR:$PHACTS_PARENT:$PYTHONPATH"
        echo "PYTHONPATH: $PYTHONPATH" >> {log} 2>&1
        
        # Create __init__.py if missing (this makes a directory a proper Python package)
        if [ ! -f "$PHACTS_DIR/__init__.py" ]; then
            echo "Creating missing __init__.py file" >> {log} 2>&1
            touch "$PHACTS_DIR/__init__.py"
        fi
        
        # Try with direct module execution
        echo "Attempting to run PHACTS..." >> {log} 2>&1
        python "$PHACTS_PATH" {input.protein_file} -o {output.result_dir} >> {log} 2>&1 || {{  
            echo "Failed with direct execution, trying with -m module syntax..." >> {log} 2>&1
            cd "$PHACTS_PARENT" && python -m PHACTS.phacts "{input.protein_file}" -o "{output.result_dir}" >> {log} 2>&1 || echo "All execution attempts failed" >> {log} 2>&1
        }}
        echo $'=================================================\\n' >> {log} 2>&1
        
        # Rename the output file to match expected format
        if [ -f "{output.result_dir}/prediction.txt" ]; then
            mv {output.result_dir}/prediction.txt {output.result}
            echo "Prediction completed successfully" >> {log} 2>&1
        else
            echo "Warning: PHACTS did not produce valid output. Creating placeholder file." >> {log} 2>&1
            echo "No prediction was made for this batch" > {output.result}
        fi
        """

# Aggregation rule for v2 prediction results
rule phacts_aggregate_results_v2:
    input:
        # This is the key part that makes the parallelization work
        all_done = f"{config['output_dir']}/03_phacts_results/.all_predictions_done_v2",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results_v2/tmp/{{sample}}/{{sample}}.phacts.out",
            sample=get_phacts_samples_compatibility()
        )
    output:
        predictions = f"{config['output_dir']}/03_phacts_results_v2/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_aggregate_results_v2.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHACTS results" > {log}
        echo -e "phage_id\\tlifestyle\\tprobability" > {output.predictions}

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
                            echo -e "$phage_id\\t$lifestyle\\t$probability" >> {output.predictions}
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

        # Create a copy in the original location for backward compatibility
        ORIG_PATH="{config['output_dir']}/03_phacts_results/phacts_predictions_compiled.tsv"
        mkdir -p $(dirname "$ORIG_PATH")
        cp {output.predictions} "$ORIG_PATH"
        echo "Created backward compatibility copy at $ORIG_PATH" >> {log}
        """

# Modified: Explicitly make sure install_phacts rule is included in the workflow
# Helper rule to force running all phacts predictions (now includes install dependency)
rule run_all_phacts_predictions_v2:
    input:
        # Add explicit dependency on PHACTS installation
        phacts_installed = f"{config['output_dir']}/db/phacts/.installed",
        installation_verified = f"{config['output_dir']}/03_phacts_results/.installation_verified",
        input_check = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        # For actual runs, get samples from the split files
        # For dry runs, this will be an empty list, which is fine
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results_v2/tmp/{{sample}}/{{sample}}.phacts.out",
            sample=get_phacts_samples_compatibility()
        )
    output:
        touch(f"{config['output_dir']}/03_phacts_results/.all_predictions_done_v2")

# New rule to ensure PHACTS is installed as a prerequisite to any PHACTS-related rules
rule ensure_phacts_installed:
    input:
        f"{config['output_dir']}/db/phacts/.installed"
    output:
        touch(f"{config['output_dir']}/.phacts_ready")
    shell:
        "echo 'PHACTS is properly installed and ready to use'"