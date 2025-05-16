# Import the helper functions at the top of your file
import os
from workflow.scripts.phacts_helper import generate_phacts_command, get_phacts_path

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
        
        # Create flag file to indicate successful installation
        cd ..
        touch .installed
        
        echo "PHACTS installation complete" >> {log} 2>&1
        """

# 5a. Run PHACTS for lifestyle prediction on a single protein file batch (UPDATED version)
rule phacts_single_prediction:
    input:
        install_check = f"{config['output_dir']}/db/phacts/.installed",
        checkpoint = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        input_check = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa"
    output:
        result_dir = directory(f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}"),
        result = f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}/{{sample}}.phacts.out"
    threads: 4
    log:
        f"{config['output_dir']}/logs/phacts_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Create output directory
        mkdir -p {output.result_dir}
        
        # Get PHACTS path from config or use the installed version
        if [ -n "{config.get('databases', {}).get('phacts', {}).get('path', '')}" ]; then
            PHACTS_PATH="{config['databases']['phacts']['path']}"
            PHACTS_DIR=$(dirname "$PHACTS_PATH")
            echo "Using configured PHACTS at $PHACTS_PATH" > {log} 2>&1
        else
            PHACTS_PATH="{config['output_dir']}/db/phacts/PHACTS/phacts.py"
            PHACTS_DIR="{config['output_dir']}/db/phacts/PHACTS"
            echo "Using workflow-installed PHACTS at $PHACTS_PATH" > {log} 2>&1
        fi
        
        # Add PHACTS directory to PYTHONPATH and run
        export PYTHONPATH="$PHACTS_DIR:$PYTHONPATH"
        python "$PHACTS_PATH" {input.protein_file} -o {output.result_dir} >> {log} 2>&1
        
        # Rename the output file to match expected format
        if [ -f "{output.result_dir}/prediction.txt" ]; then
            mv {output.result_dir}/prediction.txt {output.result}
            echo "Prediction completed successfully" >> {log} 2>&1
        else
            echo "Warning: PHACTS did not produce valid output. Creating placeholder file." >> {log} 2>&1
            echo "No prediction was made for this batch" > {output.result}
        fi
        """