# Phage Analysis Pipeline Guidelines

## Commands
- Run Python scripts: `python scripts/01_filterMmseqsLca.py [args]`
- Run R scripts: `Rscript scripts/01_phagePrediction.R [args]`
- Run Bash scripts: `bash scripts/02_clustering.sh [args]`
- SLURM submission: `sbatch scripts/03_iPhopArray.sh`

## Code Style
- Python: 4-space indentation, descriptive function names, argparse for CLI
- R: Use pipe operators (`%>%`), explicit package loading, validate arguments
- Bash: Use "${VARIABLE}" format, UPPERCASE for path variables, comment each step
- All scripts: Include descriptive comments, handle errors, validate inputs

## Naming Conventions
- Files: Numerical prefixes (01_, 02_, 03_) to indicate execution order
- Variables: Descriptive names that indicate purpose (e.g., FILTERING_WD)
- Functions: Verb-noun format that describes the action

## Notes
- Scripts use conda environments - ensure proper environment activation
- Pipeline integrates various bioinformatics tools (mmseqs2, checkV, genomad)
- PHACTS uses an existing installation at `/home/luisalberto/Softwares/PHACTS/phacts.py`
- The workflow sets PATH and PYTHONPATH variables for PHACTS access
- Follows bioinformatics best practices for sequence analysis