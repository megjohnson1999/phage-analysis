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
    # First check if a path is provided in config
    if 'databases' in config and 'phacts' in config['databases'] and 'path' in config['databases']['phacts']:
        phacts_path = config['databases']['phacts']['path']
        if os.path.isfile(phacts_path):
            return phacts_path, os.path.dirname(phacts_path)
            
    # Next, try to find phacts.py in PATH
    phacts_in_path = shutil.which('phacts.py')
    if phacts_in_path:
        return phacts_in_path, os.path.dirname(phacts_in_path)
    
    # Check common installation locations
    default_locations = [
        "/home/luisalberto/Softwares/PHACTS/phacts.py",
        os.path.expanduser("~/Software/PHACTS/PHACTS/phacts.py"),
        os.path.expanduser("~/Softwares/PHACTS/phacts.py"),
        os.path.expanduser("~/software/PHACTS/phacts.py"),
        os.path.expanduser("~/phacts/phacts.py")
    ]
    
    for location in default_locations:
        if os.path.isfile(location):
            return location, os.path.dirname(location)
    
    # Could not find PHACTS
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
    phacts_path, phacts_dir = get_phacts_path(config)
    
    if not phacts_path:
        raise ValueError("PHACTS executable not found. Please install PHACTS or specify its path in config.yaml.")
    
    # Generate a bash command that sets up the environment properly
    cmd = f"""
# Set PYTHONPATH to include PHACTS directory
export PYTHONPATH="{phacts_dir}:$PYTHONPATH"

# Run PHACTS
python "{phacts_path}" "{input_file}" -o "{output_dir}"
"""
    return cmd