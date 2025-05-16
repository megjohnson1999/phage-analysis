#!/usr/bin/env python
"""
PHACTS Helper functions for the phage-analysis workflow.
This script helps locate PHACTS and ensures it can be properly imported.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

def get_phacts_path(config):
    """
    Find the path to PHACTS executable based on configuration settings.
    
    Args:
        config: The Snakemake config dictionary
        
    Returns:
        tuple: (phacts_path, phacts_dir) where phacts_path is the path to the executable
               and phacts_dir is the directory containing PHACTS
    """
    # First check for workflow-installed PHACTS in the output directory
    if 'output_dir' in config:
        workflow_phacts = os.path.join(config['output_dir'], 'db/phacts/PHACTS/phacts.py')
        if os.path.isfile(workflow_phacts):
            return workflow_phacts, os.path.dirname(workflow_phacts)
    
    # If we didn't find the workflow installation, we'll return None
    # (the workflow should handle this by installing PHACTS)
    return None, None

def generate_phacts_command(config, input_file, output_dir):
    """
    Generate a proper command to run PHACTS with appropriate environment setup.
    
    Args:
        config: The Snakemake config dictionary
        input_file: Path to the input protein file
        output_dir: Path to the output directory
        
    Returns:
        str: A bash command that will properly run PHACTS
    """
    # Always use the workflow-installed version
    phacts_dir = os.path.join(config['output_dir'], 'db/phacts/PHACTS')
    phacts_path = os.path.join(phacts_dir, 'phacts.py')
    
    # Make phacts_dir an absolute path
    phacts_dir = os.path.abspath(phacts_dir)
    
    # Generate a bash command that sets up the environment properly
    cmd = f"""
# Set PYTHONPATH to include PHACTS directory
export PYTHONPATH="{phacts_dir}:$PYTHONPATH"

# Add parent directory to handle potential 'import phacts' cases
export PYTHONPATH="$(dirname {phacts_dir}):$PYTHONPATH"

# Debug info
echo "Using PHACTS path: {phacts_path}" >&2
echo "PYTHONPATH set to: $PYTHONPATH" >&2

# Run PHACTS
python "{phacts_path}" "{input_file}" -o "{output_dir}"
"""
    return cmd