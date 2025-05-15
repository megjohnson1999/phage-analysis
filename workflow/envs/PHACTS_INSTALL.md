# PHACTS Installation Guide

PHACTS (PHAge Classification Tool Set) is not available through conda channels and needs to be installed manually. This guide explains how to properly install it for use with the pipeline.

## Prerequisites

The conda environment already includes these dependencies, but they are listed here for reference:
- Python 3.7+
- NumPy
- Pandas
- Biopython
- scikit-learn 0.24
- matplotlib
- libsvm

## Installation Options

### Option 1: Configure an Existing Installation

If PHACTS is already installed on your system (e.g., at `/home/luisalberto/Softwares/PHACTS/phacts.py`), you can simply specify this path in your `config.yaml`:

```yaml
databases:
  phacts:
    path: "/path/to/PHACTS/phacts.py"
```

### Option 2: Install PHACTS to Your PATH

1. First, activate the conda environment:
   ```bash
   conda activate phacts
   ```

2. Clone the PHACTS repository:
   ```bash
   git clone https://github.com/deprekate/PHACTS.git
   cd PHACTS
   ```

3. Install PHACTS so it's in your PATH:
   ```bash
   pip install -e .
   ```
   
   Alternatively, if the above doesn't work, you can install it manually:
   ```bash
   cp phacts.py $CONDA_PREFIX/bin/
   chmod +x $CONDA_PREFIX/bin/phacts.py
   ```

4. Verify installation:
   ```bash
   python -c "import phacts"
   which phacts.py
   ```

## Troubleshooting

If you encounter "phacts.py not found" errors:

1. Make sure the conda environment is activated before running the pipeline
2. Check if the path in your `config.yaml` is correct
3. If using PATH-based installation, check if phacts.py is in your PATH by running `which phacts.py`
4. If not found, manually install as described in Option 2, step 3

## How the Pipeline Finds PHACTS

The pipeline uses the following order to locate PHACTS:

1. First checks if a path is provided in `config.yaml`
2. If no path is provided or the file doesn't exist, it tries to find phacts.py in your PATH
3. If neither method works, the pipeline will fail with an error message

This approach provides flexibility while ensuring the pipeline can still work with existing installations.