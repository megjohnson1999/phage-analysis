# Import the helper function at the top of your file
from workflow.scripts.phacts_helper import generate_phacts_command

# 5a. Run PHACTS for lifestyle prediction on a single protein file batch (UPDATED version)
rule phacts_single_prediction:
    input:
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