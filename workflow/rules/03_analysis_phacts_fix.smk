# Import the helper functions at the top of your file
import os
from workflow.scripts.phacts_helper import generate_phacts_command, get_phacts_path

# 5a. Install PHACTS if needed
rule install_phacts:
    output:
        flag = f"{config['output_dir']}/03_phacts_results/.phacts_installed"
    log:
        f"{config['output_dir']}/logs/phacts_installation.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Try to locate PHACTS installation
        PHACTS_PATH=$(python -c "from workflow.scripts.phacts_helper import get_phacts_path; path, _ = get_phacts_path({}); print(path or '')")
        
        if [ -z "$PHACTS_PATH" ]; then
            echo "PHACTS not found. Installing..." > {log}
            # Run the installation script
            bash scripts/install_phacts.sh -e $(basename $CONDA_PREFIX) >> {log} 2>&1
        else
            echo "PHACTS found at $PHACTS_PATH" > {log}
        fi
        
        # Create flag file to indicate installation is complete
        touch {output.flag}
        """

# 5b. Run PHACTS for lifestyle prediction on a single protein file batch (UPDATED version)
rule phacts_single_prediction:
    input:
        checkpoint = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        input_check = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa",
        phacts_installed = f"{config['output_dir']}/03_phacts_results/.phacts_installed"
    output:
        result_dir = directory(f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}"),
        result = f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}/{{sample}}.phacts.out"
    threads: 4
    log:
        f"{config['output_dir']}/logs/phacts_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    script:
        """
        #!/bin/bash
        
        # Create output directory
        mkdir -p {output.result_dir}
        
        # Get the filename without extension
        NAME=$(basename {input.protein_file} .faa)
        
        # Run the PHACTS command with proper environment setup
        {generate_phacts_command(config, input.protein_file, output.result_dir)} > {log} 2>&1
        
        # Rename the output file to match expected format
        if [ -f "{output.result_dir}/prediction.txt" ]; then
            mv {output.result_dir}/prediction.txt {output.result}
        else
            echo "Warning: PHACTS did not produce valid output. Creating placeholder file." >> {log}
            echo "No prediction was made for this batch" > {output.result}
        fi
        """