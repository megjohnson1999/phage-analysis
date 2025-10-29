# Phage Analysis Pipeline

A comprehensive Snakemake workflow for identifying, clustering, and characterizing phages from metagenomic assemblies.

## Pipeline Overview

![Phage Analysis Workflow](phage-analysis-102925.drawio.svg)

## Features

- **Multi-tool phage prediction**: Jaeger, geNomad, PHOLD, CheckV
- **Optional graph-based assembly processing**: Reneo integration for improved contig binning
- **Viral clustering**: Group similar phages into vOTUs (optional)
- **Host prediction**: Identify bacterial hosts using iPhop
- **Lifestyle prediction**: Classify phages as temperate/virulent using BACPHLIP and Phabox2
- **Taxonomic classification**: Multiple approaches via MMseqs2, Phabox2, and vContact3
- **Taxonomic consensus**: Hierarchical integration of taxonomy predictions from multiple tools
- **Functional annotation**: Protein prediction and annotation
- **Per-sample abundance**: Calculate phage and ORF abundance across samples using CoverM (optional)
- **Progress tracking**: Automated summary collection and HTML report generation

The pipeline includes a robust taxonomic consensus system that combines predictions from MMseqs2 (protein similarity), Phabox2 (ML-based), and vContact3 (gene content), with automatic format detection and hierarchical validation. For detailed technical information, see [WORKFLOW_SUMMARY.md](WORKFLOW_SUMMARY.md).

## Requirements

- Snakemake version 8+
- [Mamba](https://anaconda.org/conda-forge/mamba) for environment management
- [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for SLURM execution
- For Reneo (GFA) workflow: Gurobi license

## Installation

```bash
# 1. Clone repository
git clone https://github.com/megjohnson1999/phage-analysis.git
cd phage-analysis

# 2. Create conda environment
conda create -n phage-pipeline python=3.11
conda activate phage-pipeline

# 3. Install Snakemake and dependencies
conda install -c conda-forge -c bioconda snakemake=8 mamba
pip install snakemake-executor-plugin-slurm

# 4. Download required databases (CheckV, geNomad, PHOLD, etc.)
# See tool documentation or use existing institutional databases
```

## Quick Start

### 1. Copy and Edit Configuration File

```bash
# Copy the example config
cp config/config.yaml config/my_config.yaml

# Edit with your paths (see Configuration section below for details)
nano config/my_config.yaml  # or use your preferred editor
```

**Key settings you MUST change:**
- `output_dir`: Where to save results
- `assembly_file` OR `assembly_graph`: Your input file (provide one, not both)
- `reads_dir`: Directory containing your FASTQ read files
- All database paths under `databases:` section

**Optional settings:**
- `do_clustering: false` - Skip clustering to speed up analysis
- SLURM resources in `profile/slurm/config.yaml`

### 2. Run the Pipeline

**Tip**: Use `screen` or `tmux` to keep the pipeline running if your SSH connection drops!

```bash
# Start a screen session
screen -S phage_run

# Navigate to workflow directory
cd phage-analysis/workflow

# Run with your config file
snakemake --profile ../profile/slurm --configfile ../config/my_config.yaml

# Detach from screen: Ctrl+A, then D
# Reattach later: screen -r phage_run
```

**Additional options:**
```bash
# Dry run to preview what will be executed
snakemake -n --profile ../profile/slurm --configfile ../config/my_config.yaml

# Local execution without SLURM
snakemake --use-conda --cores 8 --configfile ../config/my_config.yaml

# Override config values from command line (if needed)
snakemake --profile ../profile/slurm --configfile ../config/my_config.yaml \
  --config do_clustering=false
```

## Configuration

### Minimal Configuration

Edit `config/my_config.yaml` with your paths:

```yaml
# Output directory
output_dir: "/path/to/outputs"

# Input (provide ONE of these)
assembly_file: "/path/to/assembly.fasta"    # For FASTA workflow
assembly_graph: ""                          # Or GFA path for Reneo

# Required
reads_dir: "/path/to/reads/"

# Optional
do_clustering: true                         # Set false to skip clustering
calculate_abundance: false                  # Set true to calculate contig and ORF abundance
abundance_coverage_threshold: 0.75          # Min coverage to count abundance (0.0-1.0)

# Database paths (update these!)
databases:
  checkv:
    db: "/path/to/checkv-db"
  genomad:
    db: "/path/to/genomad_db"
  # ... see config/config.yaml for all databases
```

### SLURM Configuration

Edit `profile/slurm/config.yaml` for your HPC cluster:
- Adjust `default-resources` (memory, runtime, CPUs)
- Set `slurm_account` to your account name

See `config/config.yaml` for all available options and detailed comments.

## Key Output Files

The pipeline generates organized results in your specified `output_dir`:

**Main Results:**
- **`final_contig_summary.tsv`** - Comprehensive table with all annotations (contig info, phage predictions, CheckV quality, taxonomy consensus, lifestyle, host predictions)
- **`Pipeline_Summary_Report.html`** - Interactive HTML report with statistics, progress tracking, and visualizations

**Detailed Results:**
- `01_phage_predictions/phageContigs.fasta` - Predicted phage sequences
- `02_clustering/vOTU_repSeqs.fasta` - Representative sequences per vOTU (if clustering enabled)
- `03_checkv_final/quality_summary.tsv` - Quality metrics
- `03_genomic_info/consensus_taxonomy.tsv` - Integrated taxonomy from multiple tools
- `03_genomic_info/lifestyle_consensus.tsv` - Lifestyle predictions (BACPHLIP + Phabox2)
- `03_iphop_results/iphop_predictions_compiled.tsv` - Host predictions

**Abundance Outputs** (when `calculate_abundance: true`):
- `04_abundance/tpm_matrix.tsv` - Contig TPM abundance (contigs × samples)
- `04_abundance/count_matrix.tsv` - Contig read counts (contigs × samples)
- `04_abundance/orf_tpm_matrix.tsv` - ORF TPM abundance (ORFs × samples)
- `04_abundance/orf_count_matrix.tsv` - ORF read counts (ORFs × samples)
- `04_abundance/orf_annotations.tsv` - ORF metadata linking to contigs

**To view the HTML report:**
```bash
# On your local machine after copying from cluster
open Pipeline_Summary_Report.html

# Or copy from cluster
scp user@cluster:/path/to/output_dir/Pipeline_Summary_Report.html .
```

## Troubleshooting

**Memory/Time Issues**
```bash
# Increase resources for specific rules
snakemake --profile ../profile/slurm \
  --set-resources iphop_single_prediction:mem_mb=100000 \
                  iphop_single_prediction:runtime=2880
```

**Reneo Issues**
- Check Gurobi license: `python -c "import gurobipy"`
- Review `logs/reneo_binning.log`
- Pipeline continues if virus detection fails

**Empty Results**
- Verify input file quality and minimum contig lengths
- Check logs in `output_dir/logs/` for tool-specific errors
- Ensure databases are downloaded and paths are correct

**Database Errors**
- Verify all database paths in config file
- Check database versions match tool requirements
- Ensure sufficient disk space for databases and outputs

## Support

- **Issues**: [GitHub Issues](https://github.com/megjohnson1999/phage-analysis/issues)
- **Documentation**: See [WORKFLOW_SUMMARY.md](WORKFLOW_SUMMARY.md) for detailed technical information
- **Test Data**: Example configuration and data in `test_data/` directory

## Citation

If you use this pipeline, please cite:
- The individual tools used (see citations in tool documentation)
- This pipeline: [GitHub repository](https://github.com/megjohnson1999/phage-analysis)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
