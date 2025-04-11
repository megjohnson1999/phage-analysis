# Phage Analysis Pipeline

A Snakemake workflow for phage prediction, clustering and characterization from metagenomic assemblies.

## Features

- Phage prediction using multiple tools (Jaeger, geNomad, Phold, CheckV)
- Viral contig clustering into vOTUs (optional)
- Host prediction using iPhop
- Lifestyle prediction using PHACTS
- Taxonomic classification using multiple approaches

## Installation

1. Clone this repository:
```
git clone https://github.com/megjohnson1999/phage-analysis.git
cd phage-analysis
```

2. Update the `config/config.yaml` file with your input paths and parameters

## Usage

1. Edit the configuration in `config/config.yaml` to set your input files and parameters
2. Run the full pipeline:
```
snakemake --use-conda -j <cores>
```

3. Input Configuration Options:
   - **FASTA only**: Provide only `assembly_file` (skips Reneo)
   - **Assembly graph only**: Provide only `assembly_graph` (uses Reneo)
   - **Both files**: Provide both `assembly_file` and `assembly_graph` (uses Reneo and original FASTA)
```
# Example using only FASTA file (skips Reneo)
snakemake --use-conda -j <cores> --config assembly_file="/path/to/assembly.fasta" assembly_graph=""

# Example using only assembly graph (runs Reneo)
snakemake --use-conda -j <cores> --config assembly_file="" assembly_graph="/path/to/assembly.gfa"

# Example using both files
snakemake --use-conda -j <cores> --config assembly_file="/path/to/assembly.fasta" assembly_graph="/path/to/assembly.gfa"
```

4. Run the pipeline without clustering:
```
snakemake --use-conda -j <cores> --config do_clustering=false
```

### SLURM Execution

The pipeline includes a SLURM profile configuration in `profile/slurm/config.yaml` for running on a SLURM cluster:

1. Run the complete pipeline with the SLURM profile:
```
snakemake --profile profile/slurm
```

2. Run without clustering on SLURM:
```
snakemake --profile profile/slurm --config do_clustering=false
```

3. The SLURM profile includes the following default settings:
   - Default resources: 50GB memory, 24 cores, 24 hour runtime
   - SLURM account: sahlab
   - Maximum concurrent jobs: 40
   - Conda environment activation enabled
   - Increased resources for taxonomy-related tasks

4. To override default SLURM settings:
```
snakemake --profile profile/slurm --default-resources mem_mb=100000 runtime=2880
```

5. To run a specific step of the pipeline:
```
snakemake --profile profile/slurm <target_rule>
```

## Pipeline Steps

1. **Phage Prediction**: Identify phage contigs from metagenome assembly
   - Assembly processing (with or without Reneo)
   - mmseqs2 taxonomy assignment
   - Multiple phage prediction tools
   - Integration of prediction results

2. **Clustering** (Optional): Group similar phages into vOTUs
   - ANI-based clustering using vclust
   - Representative sequence selection

3. **Analysis**: Characterize phage genomes
   - Host prediction with iPhop
   - Lifestyle prediction with PHACTS
   - Taxonomic classification
   - Functional annotation

## Configuration

Key configuration options:

- `assembly_file`: Path to your input assembly FASTA (can be empty if only using assembly_graph)
- `assembly_graph`: Path to assembly graph file (leave empty to skip Reneo, required if no assembly_file)
- `do_clustering`: Set to false to skip the clustering step
- `reads_dir`: Directory containing sequencing reads
- Tool parameters and resource allocations

## Output

The pipeline produces the following key outputs:

- `01_phage_predictions/phageContigs.fasta`: Predicted phage contigs
- `02_clustering/vOTU_repSeqs.fasta`: Representative sequences for viral OTUs (if clustering enabled)
- `03_iphop_results/iphop_predictions_compiled.tsv`: Host predictions
- `03_phacts_results/phacts_predictions_compiled.tsv`: Lifestyle predictions
- `03_genomic_info/`: Taxonomic and functional annotations

## Dependencies

- Snakemake
- Python 3.7+
- R 4.0+
- Various bioinformatics tools (installed via conda environments)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
