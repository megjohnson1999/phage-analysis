# Phage Analysis Pipeline

A Snakemake workflow for phage prediction, clustering, and characterization.

## Features

- Phage prediction using multiple tools (Jaeger, geNomad, Phold, CheckV)
- Viral contig clustering into vOTUs (optional)
- Host prediction using iPhop
- Lifestyle prediction using BACPHLIP and Phabox2
- Taxonomic classification using Phabox2 and vContact3

For a detailed overview of the entire workflow, see [WORKFLOW_SUMMARY.md](WORKFLOW_SUMMARY.md).

## Requirements

- Snakemake version 8+
- [Mamba](https://anaconda.org/conda-forge/mamba) for environment management
- [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for SLURM execution

## Installation

1. Clone this repository:
```
git clone https://github.com/megjohnson1999/phage-analysis.git
cd phage-analysis
```

2. Make sure you have Snakemake v8+, mamba, and the SLURM executor plugin installed in your environment:
```
conda install -c conda-forge -c bioconda snakemake=8 mamba
pip install snakemake-executor-plugin-slurm
```

3. Configure the pipeline in `config/config.yaml`:
   - Set input and output paths
   - Specify database locations
   - Configure workflow options

## Workflow Structure

The pipeline is organized into three main modules:

1. **Phage Prediction** (01_prediction.smk)
   - Filter contigs by length
   - Process assemblies (with optional Reneo graph-based processing)
   - Run mmseqs2 for taxonomy assignment
   - Execute multiple phage prediction tools (Jaeger, geNomad, Phold, CheckV)
   - Integrate results to identify phage contigs

2. **Clustering** (02_clustering.smk - optional)
   - Group similar phages into vOTUs using vclust
   - Select representative sequences

3. **Analysis** (03_analysis.smk)
   - Host prediction with iPhop
   - Lifestyle prediction with BACPHLIP and Phabox2
   - Taxonomic classification using Phabox2 and vContact3
   - Functional annotation

## Key Analysis Tools

### BACPHLIP
BACPHLIP (Bacteriophage lifestyle prediction tool) predicts phage lifestyle using a machine learning approach. It analyzes genomic features to classify phages as either temperate or virulent.

### Phabox2
Phabox2 provides comprehensive phage analysis including:
- Taxonomic classification at multiple levels
- Lifestyle prediction (lytic/lysogenic)
- Machine learning-based predictions from genomic features

### vContact3
vContact3 performs gene-content based taxonomic classification:
- Creates gene-sharing networks from protein sequences
- Assigns taxonomy based on network clustering
- Provides complementary classification to sequence-based methods

## Usage

### Configuration

1. Edit the configuration in `config/config.yaml` to set your input files and parameters:
   ```yaml
   # Key configuration options:
   output_dir: "/path/to/your/outputs"
   assembly_file: "/path/to/assembly.fasta"  # FASTA only workflow
   assembly_graph: "/path/to/graph.gfa"      # Optional for Reneo graph-based workflow
   reads_dir: "/path/to/reads/"              # Directory containing reads
   do_clustering: true                       # Set to false to skip clustering
   ```

2. Edit SLURM profile in `profile/slurm/config.yaml` if needed for your computing environment

### Running the Pipeline

Navigate to the repository directory and run:

```bash
# Standard execution with SLURM
snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"

# Run from workflow directory
cd workflow
snakemake --profile ../profile/slurm
```

### Input Options

The pipeline supports two mutually exclusive input modes:

- **FASTA-only workflow**: Provide `assembly_file` parameter (skips Reneo)
  ```
  snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

- **Reneo workflow**: Provide `assembly_graph` parameter (runs Reneo for graph-based binning)
  ```
  snakemake --profile profile/slurm --config assembly_graph="/path/to/graph.gfa" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```
  
  **Note**: When using GFA files, you don't need to explicitly set `assembly_file=""` - the pipeline automatically detects and clears placeholder paths from the default config.

**Important**: You must provide either `assembly_file` OR `assembly_graph`, not both. The pipeline will give an error if you try to provide both.

#### Reneo Integration

The pipeline includes a wrapper script (`workflow/scripts/run_reneo_wrapper.sh`) that handles Reneo execution. This wrapper:
- Manages Reneo's expected failures gracefully
- Creates empty output files if Reneo fails, allowing the pipeline to continue
- Logs detailed information about Reneo's execution status

**Important: Reneo requires a Gurobi license.** There are two ways to configure Reneo:

1. **Use an existing Reneo environment** (recommended if you already have Reneo with Gurobi set up):
   ```yaml
   # In config/config.yaml:
   conda_base_path: "/path/to/your/conda"  # e.g., "/ref/sahlab/software/miniforge3"
   conda_envs:
     reneo: "reneo"  # or full path: "/path/to/conda/envs/reneo"
   ```

2. **Let Snakemake create the environment** (requires Gurobi license setup):
   ```yaml
   # In config/config.yaml:
   conda_base_path: ""  # Leave empty
   conda_envs:
     reneo: "../envs/reneo.yaml"
   ```
   You'll need to set up the Gurobi license in the created environment.

### Additional Options

- Skip the clustering step:
  ```
  snakemake --profile profile/slurm --config do_clustering=false
  ```

- Override default SLURM settings:
  ```
  snakemake --profile profile/slurm --default-resources mem_mb=100000 runtime=2880
  ```

- Run a specific step of the pipeline:
  ```
  snakemake --profile profile/slurm <target_rule>
  ```

## SLURM Execution

The pipeline includes a SLURM profile configuration in `profile/slurm/config.yaml` with the following settings:

- Default resources: 50GB memory, 24 cores, 24-hour runtime
- SLURM account: sahlab
- Maximum concurrent jobs: 40
- Conda environment activation enabled
- Increased resources for memory-intensive tasks (mmseqs2, iPhop)


## Output

The pipeline produces the following key outputs:

- `01_phage_predictions/phageContigs.fasta`: Predicted phage contigs
- `02_clustering/vOTU_repSeqs.fasta`: Representative sequences for viral OTUs (if clustering enabled)
- `03_iphop_results/iphop_predictions_compiled.tsv`: Host predictions
- `03_bacphlip/bacphlip.predictions.tsv`: Lifestyle predictions from BACPHLIP
- `03_genomic_info/phabox_output/`: Phabox2 taxonomy and lifestyle predictions
- `03_genomic_info/vc3_output/`: vContact3 gene-content based taxonomy
- `03_genomic_info/`: Additional taxonomic and functional annotations

## Dependencies

The workflow automatically handles all dependencies through conda environments defined in the `workflow/envs/` directory.

## Recent Updates and Fixes

- **Replaced PHACTS with BACPHLIP**: Updated lifestyle prediction to use BACPHLIP for more accurate and efficient predictions
- **Improved Error Handling**: Enhanced PHOLD rule to handle prediction failures gracefully
- **Fixed Script Syntax**: Corrected shell script syntax in various pipeline rules
- **Added Comprehensive Documentation**: Created a detailed `WORKFLOW_SUMMARY.md` document explaining the entire pipeline
- **Enhanced Reneo Integration**: 
  - Added wrapper script to handle Reneo's expected failures gracefully
  - Fixed conda environment configuration for systems with existing Reneo/Gurobi setup
  - Resolved workflow.globals state issue during DAG construction
- **Fixed DAG Construction**: 
  - Resolved issues with dry-run when using GFA-only input
  - All rules now always defined with conditional execution logic
  - Input functions re-evaluate conditions during DAG construction
- **Automatic Placeholder Detection**: Pipeline now automatically detects and clears placeholder paths from default config
- **Flexible Conda Configuration**: Added `conda_base_path` option for using existing conda environments

## Troubleshooting

### Common Issues

1. **"Skipping Reneo - not configured for this run" despite providing GFA file**
   - **Cause**: workflow.globals state not available during DAG construction
   - **Solution**: Fixed in latest version - input functions now re-evaluate conditions locally

2. **Conda environment export error for Reneo**
   - **Cause**: Snakemake trying to export existing conda environment with full path as name
   - **Solution**: Configure `conda_base_path` in config.yaml for existing environments, or use YAML file for new environments

3. **MissingInputException with GFA files during dry-run**
   - **Cause**: DAG construction issues with conditional rules
   - **Solution**: All rules now always defined; conditional logic moved to input functions

4. **Placeholder paths causing errors**
   - **Cause**: Default config contains placeholder paths like `/path/to/assembly.fasta`
   - **Solution**: The pipeline now automatically detects and clears these. No manual intervention needed.

5. **Reneo fails but pipeline stops**
   - **Cause**: Reneo may fail on certain inputs but still produce partial outputs
   - **Solution**: The wrapper script now handles this by creating empty output files if needed

6. **PHOLD or iPhop failing on large datasets**
   - **Cause**: Memory or time limits exceeded
   - **Solution**: The pipeline chunks these analyses automatically. Adjust chunk sizes in the code if needed.

### Getting Help

- Check the logs in `<output_dir>/logs/` for detailed error messages
- Review `WORKFLOW_SUMMARY.md` for implementation details
- Submit issues to the GitHub repository with log files attached

## Documentation

- **README.md**: This file - basic usage and setup instructions
- **WORKFLOW_SUMMARY.md**: Comprehensive overview of the entire workflow and implementation details
- **workflow/TEST_INSTRUCTIONS.md**: General testing instructions for the workflow

## License

This project is licensed under the MIT License - see the LICENSE file for details.