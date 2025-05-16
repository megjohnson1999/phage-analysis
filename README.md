# Phage Analysis Pipeline

A Snakemake workflow for phage prediction, clustering, and characterization.

## Features

- Phage prediction using multiple tools (Jaeger, geNomad, Phold, CheckV)
- Viral contig clustering into vOTUs (optional)
- Host prediction using iPhop
- Lifestyle prediction using PHACTS
- Taxonomic classification using multiple approaches

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
   - Lifestyle prediction with PHACTS
   - Taxonomic classification and annotation

## PHACTS Integration

PHACTS (Phage Classification Tool Set) is used for predicting phage lifestyle. The pipeline provides two ways to use PHACTS:

### Option 1: Automatic Installation (Recommended)

The workflow now includes an automatic PHACTS installation:

- PHACTS is automatically installed during pipeline execution
- Installation happens in `output_dir/db/phacts/`
- Version control via `phacts_version` parameter in config.yaml

No manual installation is needed - simply run the pipeline and PHACTS will be installed and configured automatically.

### Option 2: Manual Installation

For users who prefer manual control:

1. Run the PHACTS installation script:
   ```
   bash scripts/install_phacts.sh
   ```

2. Or for custom installation path:
   ```
   bash scripts/install_phacts.sh -d /path/to/install/location
   ```

3. Update `config.yaml` to specify the path:
   ```yaml
   databases:
     phacts:
       path: "/path/to/phacts/phacts.py"
   ```

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
   
   # PHACTS configuration 
   phacts_version: "main"                    # Git branch or tag for PHACTS
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

The pipeline supports different input configurations:

- **FASTA only**: Provide `assembly_file` parameter (skips Reneo)
  ```
  snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

- **Assembly graph only**: Provide `assembly_graph` parameter (uses Reneo)
  ```
  snakemake --profile profile/slurm --config assembly_graph="/path/to/graph.gfa" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

- **Both files**: Provide both `assembly_file` and `assembly_graph` (uses Reneo and original FASTA)
  ```
  snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" assembly_graph="/path/to/graph.gfa" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

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
- `03_phacts_results/phacts_predictions_compiled.tsv`: Lifestyle predictions
- `03_genomic_info/`: Taxonomic and functional annotations

## Dependencies

The workflow automatically handles all dependencies through conda environments defined in the `workflow/envs/` directory, including PHACTS which is now installed as part of the workflow execution.

## License

This project is licensed under the MIT License - see the LICENSE file for details.