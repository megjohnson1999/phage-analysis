#!/bin/bash

set -e  # Exit immediately if a command exits with non-zero status

# Get the conda-prefix from slurm profile if it exists
SLURM_PROFILE="$(dirname $(pwd))/profile/slurm/config.yaml"
if [ -f "$SLURM_PROFILE" ]; then
    CONDA_PREFIX=$(grep "conda-prefix:" "$SLURM_PROFILE" | sed 's/conda-prefix: *//g' | tr -d '"')
    if [ -n "$CONDA_PREFIX" ]; then
        # Extract the last directory name from the path for the env name
        CONDA_ENV=$(basename "$CONDA_PREFIX")
        # Use the conda envs directory for installation
        INSTALL_DIR="$CONDA_PREFIX/PHACTS"
    else
        INSTALL_DIR="$HOME/Software/PHACTS"
        CONDA_ENV="phacts"
    fi
else
    INSTALL_DIR="$HOME/Software/PHACTS"
    CONDA_ENV="phacts"
fi

ACTIVATE_CONDA=true

# Display help message
function show_help {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -d DIR    Set installation directory (default: Auto-detected from slurm profile or $INSTALL_DIR)"
    echo "  -e ENV    Set conda environment name (default: Auto-detected from slurm profile or $CONDA_ENV)"
    echo "  -n        Do not activate conda environment"
    echo "  -h        Show this help message"
}

# Parse command-line options
while getopts "d:e:nh" opt; do
    case ${opt} in
        d)
            INSTALL_DIR=$OPTARG
            ;;
        e)
            CONDA_ENV=$OPTARG
            ;;
        n)
            ACTIVATE_CONDA=false
            ;;
        h)
            show_help
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            show_help
            exit 1
            ;;
    esac
done

# Create installation directory if it doesn't exist
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "=== Installing PHACTS to $INSTALL_DIR ==="

# Activate conda environment if requested
if [ "$ACTIVATE_CONDA" = true ]; then
    echo "=== Activating conda environment: $CONDA_ENV ==="
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV" || {
        echo "Error: Failed to activate conda environment $CONDA_ENV"
        echo "Please ensure the environment exists or create it first with:"
        echo "conda env create -f workflow/envs/phacts.yaml"
        exit 1
    }
fi

# Clone PHACTS from GitHub if needed
if [ ! -d "$INSTALL_DIR/PHACTS" ]; then
    echo "=== Cloning PHACTS repository ==="
    git clone https://github.com/deprekate/PHACTS.git
    cd PHACTS
else
    echo "=== PHACTS repository already exists, updating ==="
    cd PHACTS
    git pull
fi

# Install PHACTS using pip in development mode
echo "=== Installing PHACTS with pip ==="
pip install -e .

# Verify installation
echo "=== Verifying installation ==="
if python -c "import phacts" 2>/dev/null; then
    echo "=== PHACTS installation successful! ==="
else
    echo "=== Installation completed but module import failed. Adding to PYTHONPATH ==="
    
    # Create a wrapper script
    echo "Creating wrapper script for PHACTS..."
    cat > "$INSTALL_DIR/run_phacts.sh" << EOL
#!/bin/bash

# Set PYTHONPATH to include PHACTS directory
export PYTHONPATH="$INSTALL_DIR/PHACTS:\$PYTHONPATH"

# Run PHACTS with all arguments passed to this script
python "$INSTALL_DIR/PHACTS/phacts.py" "\$@"
EOL
    
    # Make wrapper script executable
    chmod +x "$INSTALL_DIR/run_phacts.sh"
    
    echo "=== Created wrapper script: $INSTALL_DIR/run_phacts.sh ==="
    echo "=== You can use this script to run PHACTS ==="
fi

# Update the config.yaml file if it exists
CONFIG_FILE="config/config.yaml"
if [ -f "$CONFIG_FILE" ]; then
    echo "=== Updating config.yaml with PHACTS path ==="
    
    # Check if phacts section exists
    if grep -q "phacts:" "$CONFIG_FILE"; then
        # Use sed to update the path
        sed -i.bak "s|path:.*\".*\"|path: \"$INSTALL_DIR/PHACTS/phacts.py\"|g" "$CONFIG_FILE"
    else
        echo "Warning: Could not update config.yaml automatically."
        echo "Please add the following to your config.yaml file:"
        echo ""
        echo "databases:"
        echo "  phacts:"
        echo "    path: \"$INSTALL_DIR/PHACTS/phacts.py\""
    fi
fi

echo "=== PHACTS Installation Complete ==="
echo "To use PHACTS, either:"
echo "1. Activate the conda environment with: conda activate $CONDA_ENV"
echo "2. Set PYTHONPATH with: export PYTHONPATH=\"$INSTALL_DIR/PHACTS:\$PYTHONPATH\""
echo "3. Use the wrapper script: $INSTALL_DIR/run_phacts.sh"
echo ""
echo "Installation directory: $INSTALL_DIR/PHACTS"
echo "Main executable: $INSTALL_DIR/PHACTS/phacts.py"