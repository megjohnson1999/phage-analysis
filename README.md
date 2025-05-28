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

**Note: The FASTA-only workflow is currently recommended.** The Reneo graph-based processing is undergoing troubleshooting and may have issues.

The pipeline supports different input configurations:

- **FASTA only (Recommended)**: Provide `assembly_file` parameter (skips Reneo)
  ```
  snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

- **Assembly graph only**: Provide `assembly_graph` parameter (uses Reneo) - *Currently being troubleshooted*
  ```
  snakemake --profile profile/slurm --config assembly_graph="/path/to/graph.gfa" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
  ```

- **Both files**: Provide both `assembly_file` and `assembly_graph` (uses Reneo and original FASTA) - *Currently being troubleshooted*
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

## Documentation

- **README.md**: This file - basic usage and setup instructions
- **WORKFLOW_SUMMARY.md**: Comprehensive overview of the entire workflow and implementation details
- **workflow/TEST_INSTRUCTIONS.md**: General testing instructions for the workflow

## License

This project is licensed under the MIT License - see the LICENSE file for details.