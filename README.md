# Phage Analysis Pipeline

A comprehensive Snakemake workflow for identifying, clustering, and characterizing phages from metagenomic assemblies.

## Features

- **Multi-tool phage prediction**: Jaeger, geNomad, PHOLD, CheckV
- **Optional graph-based assembly processing**: Reneo integration for improved contig binning
- **Viral clustering**: Group similar phages into vOTUs (optional)
- **Host prediction**: Identify bacterial hosts using iPhop
- **Lifestyle prediction**: Classify phages as temperate/virulent using BACPHLIP and Phabox2
- **Taxonomic classification**: Multiple approaches via MMseqs2, Phabox2, and vContact3
- **Taxonomic consensus**: Hierarchical integration of taxonomy predictions from multiple tools
- **Functional annotation**: Protein prediction and annotation
- **Progress tracking**: Automated summary collection and HTML report generation

### Taxonomic Consensus Integration

The pipeline includes a robust taxonomic consensus system that combines predictions from multiple tools:

- **MMseqs2**: Protein similarity-based taxonomy (highest priority)
- **Phabox2**: Machine learning-based phage-specific predictions  
- **vContact3**: Gene content-based clustering taxonomy

The consensus system:
- ✅ **Handles multiple input formats**: Automatically detects and processes LCA vs BLAST format outputs
- ✅ **Hierarchical validation**: Ensures taxonomic consistency across classification levels
- ✅ **Comprehensive coverage**: Integrates results from all available tools
- ✅ **Quality prioritization**: Uses tool-specific strengths (e.g., Phabox2 for phage-specific features)

For a detailed technical overview, see [WORKFLOW_SUMMARY.md](WORKFLOW_SUMMARY.md).

## Requirements

- Snakemake version 8+
- [Mamba](https://anaconda.org/conda-forge/mamba) for environment management
- [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for SLURM execution
- For Reneo workflow: Gurobi license (see Reneo section below)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/megjohnson1999/phage-analysis.git
cd phage-analysis
```

2. Install Snakemake and dependencies:
```bash
# Create a conda environment for running the pipeline
conda create -n phage-pipeline python=3.11
conda activate phage-pipeline

# Install Snakemake v8+, mamba, and SLURM executor
conda install -c conda-forge -c bioconda snakemake=8 mamba
pip install snakemake-executor-plugin-slurm
```

3. Download required databases (see Database Setup section below)

4. Configure the pipeline (see Configuration section below)

## Quick Start

**Tip**: Use `screen` or `tmux` to keep the pipeline running even if your SSH connection drops!

```bash
# Start a screen session
screen -S phage_run

# Navigate to the workflow directory
cd phage-analysis/workflow

# Option 1: Run with FASTA assembly
snakemake --profile ../profile/slurm \
  --config assembly_file="/path/to/assembly.fasta" \
           reads_dir="/path/to/reads/" \
           output_dir="/path/to/output/"

# Option 2: Run with GFA graph (Reneo workflow)
snakemake --profile ../profile/slurm \
  --config assembly_graph="/path/to/assembly.gfa" \
           reads_dir="/path/to/reads/" \
           output_dir="/path/to/output/"

# Option 3: Skip clustering step (works with FASTA or GFA input)
snakemake --profile ../profile/slurm \
  --config assembly_file="/path/to/assembly.fasta" \
           reads_dir="/path/to/reads/" \
           output_dir="/path/to/output/" \
           do_clustering=false

# Detach from screen: Ctrl+A, then D
# Check on it later: screen -r phage_run
```

## Configuration

### 1. Database Setup

The pipeline requires several databases. You can either:
- Use existing databases at your institution (update paths in config)
- Download databases yourself:

```bash
# Example database downloads (adjust paths as needed)
# CheckV
wget https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz

# geNomad
genomad download-database genomad_db/

# PHOLD
phold install -d phold_db/

# Other databases: see tool documentation
```

### 2. Edit Configuration File

Copy and modify the example configuration:
```bash
cp config/config.yaml config/my_config.yaml
```

Key settings in `config/my_config.yaml`:
```yaml
# Output location
output_dir: "/path/to/your/outputs"

# Input files (provide ONE of these):
assembly_file: "/path/to/assembly.fasta"  # For standard workflow
assembly_graph: ""                        # Leave empty if using FASTA

# OR for Reneo workflow:
assembly_file: ""                         # Leave empty if using GFA
assembly_graph: "/path/to/graph.gfa"      # For Reneo workflow

# Required for both workflows
reads_dir: "/path/to/reads/"              # Directory with FASTQ files

# Optional settings
do_clustering: true                        # Set false to skip vOTU clustering

# Database paths - update these!
databases:
  checkv:
    db: "/path/to/checkv-db-v1.5"
  genomad:
    db: "/path/to/genomad_db"
  # ... etc
```

### 3. SLURM Configuration (if using HPC)

Edit `profile/slurm/config.yaml` for your cluster:
```yaml
default-resources:
  - mem_mb=50000      # Default memory (50GB)
  - runtime=1440      # Default runtime (24 hours)
  - nodes=1
  - tasks=1
  - cpus_per_task=24

# Adjust for your SLURM account
slurm_account: "your_account_name"
```

## Running the Pipeline

### Basic Usage

**Important**: We recommend running the pipeline in a `screen` or `tmux` session to prevent interruption if your connection drops:

```bash
# Start a screen session
screen -S phage_pipeline

# Navigate to workflow directory
cd phage-analysis/workflow

# Run with your configuration file
snakemake --profile ../profile/slurm --configfile ../config/my_config.yaml

# Detach from screen with Ctrl+A then D
# Reattach later with: screen -r phage_pipeline
```

### Command Line Options

You can override config file settings via command line:

```bash
# Override specific parameters
snakemake --profile ../profile/slurm \
  --config assembly_file="/new/path/assembly.fasta" \
           output_dir="/new/output/path"

# Dry run to see what would be executed
snakemake -n --profile ../profile/slurm

# Run specific targets
snakemake --profile ../profile/slurm \
  results/01_phage_predictions/phageContigs.fasta

# Local execution (no SLURM)
snakemake --cores 8 --use-conda
```

### Input Modes

**Option 1: FASTA Workflow**
- Input: Assembled contigs in FASTA format
- Faster
- Example:
  ```bash
  snakemake --profile ../profile/slurm \
    --config assembly_file="assembly.fasta" \
             reads_dir="reads/" \
             output_dir="results/"
  ```

**Option 2: Reneo Graph-Based Workflow**
- Input: Assembly graph in GFA format
- Uses Reneo for enhanced contig binning
- Requires Gurobi license (see below)
- Example:
  ```bash
  snakemake --profile ../profile/slurm \
    --config assembly_graph="assembly.gfa" \
             reads_dir="reads/" \
             output_dir="results/"
  ```

### Reneo Configuration (for GFA input)

**Important: Reneo requires a Gurobi license.**

**Option 1: Use existing Reneo environment** (recommended):
```yaml
# In your config file:
conda_base_path: "/path/to/your/conda"  # e.g., "/home/user/miniforge3"
conda_envs:
  reneo: "reneo"  # name of existing environment
```

**Option 2: Let Snakemake create environment**:
```yaml
conda_base_path: ""  # Leave empty
conda_envs:
  reneo: "../envs/reneo.yaml"
```
Then set up Gurobi license in the created environment.

**Note**: If Reneo's virus detection fails (common with test data), the pipeline will continue using all processed sequences.

## Pipeline Outputs

The pipeline creates organized output directories:

```
output_dir/
├── Pipeline_Summary_Report.html         # 📊 Interactive HTML summary report
├── final_contig_summary.tsv            # ⭐ Comprehensive annotation table (all results in one file)
├── 01_checkv_output/                    # Initial quality assessment
│   └── quality_summary.tsv             # CheckV results (used for phage prediction)
├── 01_phage_predictions/
│   ├── phageContigs.fasta              # Final predicted phage sequences
│   └── phagePredictedContigs.tsv       # Detailed prediction scores
├── 02_clustering/                       # (if clustering enabled)
│   ├── vOTU_repSeqs.fasta              # Representative sequence per vOTU
│   └── clusters.tsv                     # Cluster membership information
├── 03_checkv_final/                     # Final quality assessment
│   └── quality_summary.tsv             # CheckV results matching analyzed sequences
├── 03_iphop_results/
│   └── iphop_predictions_compiled.tsv   # Bacterial host predictions
├── 03_genomic_info/
│   ├── consensus_taxonomy.tsv          # 🎯 Integrated taxonomy consensus
│   ├── consensus_taxonomy_summary.json # Consensus statistics and coverage
│   ├── lifestyle_consensus.tsv         # 🎯 Integrated lifestyle predictions (BACPHLIP + Phabox2)
│   ├── mmseqs_taxonomy.tsv             # Sequence-based taxonomy
│   ├── bacphlip_lifestyle.tsv          # BACPHLIP raw predictions
│   ├── bacphlip_lifestyle_with_completeness.tsv # BACPHLIP with CheckV quality
│   ├── phabox_output/                  # Phabox2 results
│   │   ├── final_prediction_summary.tsv # ML-based taxonomy and lifestyle
│   │   ├── taxonomy.tsv                # Legacy format (may be empty)
│   │   └── lifestyle.tsv               # Legacy format (may be empty)
│   └── vc3_output/                     # vContact3 gene-content taxonomy
│       └── exports/final_assignments.csv # Taxonomy predictions
├── pipeline_summaries/                 # JSON summary files for each step
└── logs/                               # Detailed logs for each step
```

### Key Output Files

**`final_contig_summary.tsv`** - The main results file consolidating all annotations:
- Contig ID and length
- Phage prediction rule
- CheckV quality metrics (completeness, contamination, provirus status)
- MMseqs2 LCA assignment from initial viral screening
- Consensus taxonomy (integrated from multiple tools)
- Lifestyle predictions (virulent/temperate/chronic) with confidence scores
- Host predictions from iPhop

This single file provides a complete overview of all analyzed sequences with their annotations.

**`lifestyle_consensus.tsv`** - Integrated lifestyle predictions combining:
- BACPHLIP predictions (prioritized when confidence ≥ 0.7)
- Phabox2 predictions (used as fallback)
- Source tracking and confidence scores

### Summary Report

The pipeline automatically generates a comprehensive HTML summary report (`Pipeline_Summary_Report.html`) that includes:

- **Configuration Overview**: Shows your input files and pipeline settings
- **Overall Statistics**: Key metrics like total input sequences, viral contigs found, and final phages
- **Progress Tracking**: Visual progress bar and completion status for each pipeline step
- **Detailed Results**: Step-by-step breakdown with sequence counts, quality metrics, and tool outputs

**To view the report**:
```bash
# Open in web browser (Linux/Mac)
open output_dir/Pipeline_Summary_Report.html

# Or copy to your local machine and open
scp user@cluster:path/to/output_dir/Pipeline_Summary_Report.html .
```

The report is generated automatically when the pipeline completes successfully.

## Troubleshooting

### Memory/Time Issues
- Increase resources for specific rules:
  ```bash
  snakemake --profile ../profile/slurm \
    --set-resources iphop_single_prediction:mem_mb=100000
  ```

### Reneo Issues
- Check Gurobi license: `python -c "import gurobipy"`
- Review wrapper log: `logs/reneo_binning.log`
- If virus detection fails, pipeline continues with all sequences

### Empty Results
- Check input file quality and size
- Review logs for tool-specific errors
- Ensure databases are properly downloaded

### Database Errors
- Verify database paths in config
- Check database versions match tool requirements
- Ensure sufficient disk space

## Example Commands

### Small Test Dataset
```bash
# Test with provided test data
cd workflow
snakemake --profile ../profile/slurm \
  --configfile ../test_data/test_config.yaml

# Test locally without SLURM
snakemake --use-conda --cores 4 \
  --configfile ../test_data/test_config_local.yaml
```

### Real Data Examples
```bash
# Always start in a screen session for long runs!
screen -S my_phage_analysis

# Metagenomic assembly from SPAdes
snakemake --profile ../profile/slurm \
  --config assembly_file="spades_output/scaffolds.fasta" \
           reads_dir="illumina_reads/" \
           output_dir="phage_results/"

# Assembly graph from metaFlye with Reneo
snakemake --profile ../profile/slurm \
  --config assembly_graph="flye_output/assembly_graph.gfa" \
           reads_dir="nanopore_reads/" \
           output_dir="phage_results_reneo/"

# Large dataset with custom resources
snakemake --profile ../profile/slurm \
  --config assembly_file="megahit_output/final.contigs.fa" \
           reads_dir="hiseq_reads/" \
           output_dir="large_dataset_results/" \
  --set-resources iphop_single_prediction:mem_mb=200000 \
                  iphop_single_prediction:runtime=2880

# Detach and let it run: Ctrl+A, D
```

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Support

- **Issues**: [GitHub Issues](https://github.com/megjohnson1999/phage-analysis/issues)
- **Documentation**: See [WORKFLOW_SUMMARY.md](WORKFLOW_SUMMARY.md) for technical details
- **Test Data**: Example files in `test_data/` directory

## Citation

If you use this pipeline, please cite:
- The individual tools used (see citations in tool documentation)
- This pipeline: [GitHub repository](https://github.com/megjohnson1999/phage-analysis)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
