# Phage Analysis Pipeline: Comprehensive Workflow Summary

## Overview

The phage analysis pipeline is a comprehensive, modular Snakemake workflow for identifying, clustering, and characterizing phages from metagenomic assemblies. It integrates multiple specialized bioinformatics tools to deliver robust predictions on phage taxonomy, host relationships, and lifestyle characteristics.

## Workflow Stages

The pipeline consists of three main stages:

1. **Prediction**: Identifying phage sequences from metagenomic assemblies
2. **Clustering** (Optional): Grouping similar phages into viral OTUs
3. **Analysis**: Characterizing phage properties (host range, lifestyle, taxonomy)

## Stage 1: Prediction

### Input Options
- **Assembly File**: Direct input of assembled contigs (FASTA)
- **Assembly Graph**: For graph-based assembly refinement with Reneo
- **Reads Directory**: Required for Reneo and coverage analysis

### Processing Steps

1. **Optional Assembly Refinement with Reneo**
   - Processes assembly graph files with raw reads for improved contig binning
   - Filters resulting contigs to ≥1KB length

2. **Initial Taxonomy Assignment (MMSeqs2)**
   - Performs sensitive taxonomic classification against reference databases
   - Parameters:
     ```
     --min-length 30 --e-value 1e-15 --search-type 2 --sensitivity 4.0
     --lca-mode 2 --add-taxa --tax-lineage 2
     ```

3. **Filtering for Viral Contigs**
   - Python script (`01_filterMmseqsLca.py`) analyzes taxonomy results
   - Selection criteria:
     - First taxID in lineage is 10239 (Viruses)
     - OR taxID is 1 (root)
     - OR taxID is 0 (unclassified)
   - Outputs filtered taxonomy table and viral contig IDs

4. **Multi-tool Phage Prediction**
   - **Jaeger**: Neural network-based prediction (score threshold 2.5)
   - **GeNomad**: Comprehensive viral detection (min-score 0.6)
   - **PHOLD**: Phage protein annotation and functional prediction
   - **CheckV**: Quality assessment and completeness estimation

5. **Integration of Predictions**
   - R script combines evidence from all prediction tools
   - Applies criteria for high-confidence phage contigs:
     - Functional diversity ≥3, OR
     - Special topology (DTR, ITR, Provirus), OR
     - Strong evidence from multiple tools

## Stage 2: Clustering (Optional)

### Purpose
Group similar phage sequences into viral operational taxonomic units (vOTUs)

### Processing Steps

1. **vClust-Based Clustering**
   - Length filtering (default: minimum 10KB)
   - Three-step process:
     - Initial prefiltering with k-mer similarity
     - Pairwise alignment of potential matches
     - Clustering with Leiden algorithm
   - Default parameters: 95% identity, 85% coverage

2. **Representative Sequence Selection**
   - Extracts one representative sequence per cluster
   - Creates collection of vOTU representatives for downstream analysis

## Stage 3: Analysis

### Purpose
Characterize phages by predicting host relationships, lifestyle, and taxonomy

### Processing Steps

1. **Parallel Processing Setup**
   - Splits phage sequences into manageable chunks (100 sequences/job)
   - Creates checkpoint system for workflow management

2. **Host Prediction (iPhop)**
   - Processed in parallel for better performance
   - Predicts bacterial hosts for each phage
   - Aggregates to single TSV with host predictions at genus level

3. **Protein Prediction (Prodigal)**
   - Identifies open reading frames in phage sequences
   - Uses metagenome mode for higher sensitivity
   - Outputs protein and gene sequences

4. **Lifestyle Prediction (PHACTS)**
   - **Phage-Specific Approach** (default):
     - Splits proteins by phage using `split_proteins_by_phage.py`
     - Each PHACTS prediction uses only proteins from one phage
     - Improves prediction accuracy significantly
     - Robust error handling for failed predictions
   - **Original Approach** (legacy):
     - Processes protein batches without phage-specific separation
     - Less accurate due to mixing proteins from different phages

5. **Taxonomic Classification**
   - **MMSeqs2**: Detailed sequence-based taxonomic classification
   - **Phabox2**: Phage-specific taxonomy and lifestyle analysis
   - **vContact3**: Gene content-based taxonomy using protein clustering

## Implementation Details

### Parallelization Strategy

The workflow implements efficient parallelization:

1. **Smart Chunking**
   - Divides large tasks into manageable pieces
   - Particularly important for PHOLD, iPhop, and PHACTS analyses

2. **Checkpoint System**
   - Uses Snakemake checkpoints for handling dynamic file generation
   - Three phases for parallel operations:
     ```
     .splits_ready → .input_files_found → .all_predictions_done
     ```

3. **Results Aggregation**
   - Merges individual results with proper header handling
   - Accounts for failed or empty results gracefully

### Error Handling

Comprehensive error handling throughout:

1. **Input Validation**
   - Checks for existence of required files
   - Validates configuration settings

2. **Robust Failure Recovery**
   - Creates empty outputs with proper headers when tools fail
   - Enables workflow to continue despite partial failures
   - Example from PHOLD rule:
     ```bash
     # Run with error trapping
     (phold run -i {input.contig_file} -o {output.results_dir} -d {config[databases][phold][db]} \
         -t {threads} --cpu --force > {log} 2>&1) || true
         
     # Create empty output if tool failed
     if [ ! -f "{output.predictions}" ]; then
         echo "WARNING: PHOLD failed to create output for {wildcards.sample}" >> {log}
         echo -e "contig_id\torf_id\tstart\tend\tstrand\taa_length\tcategory\tproduct\thit\tevalue\tidentity" > {output.predictions}
     fi
     ```

3. **Fallback Methods**
   - Multi-level sample identification:
     - Primary: Read from split file list
     - Fallback 1: Use glob to find files directly
     - Fallback 2: Check tmp directory for existing outputs

4. **Detailed Logging**
   - Each step generates comprehensive logs
   - Captures both stdout and stderr
   - Organized in logs/ directory for easy troubleshooting

### Environment Management

1. **Tool-Specific Conda Environments**
   - Separate environment for each major tool
   - Prevents version conflicts
   - Defined in workflow/envs/ directory

2. **PHACTS Configuration**
   - Special handling for PHACTS path
   - Sets PYTHONPATH to ensure modules are found:
     ```bash
     # Define the path to phacts manually
     phacts_path=/home/megan.j/PHACTS
     
     # Add phacts directory to PYTHONPATH
     export PYTHONPATH="${phacts_path}:${PYTHONPATH:-}"
     ```

## Key Features and Improvements

1. **Phage-Specific PHACTS Analysis**
   - Runs PHACTS on each phage separately for improved accuracy
   - Implemented as default analysis method

2. **Flexible Input Options**
   - Supports both direct assembly files and assembly graphs
   - Conditional workflow paths based on input types

3. **Robust Output Handling**
   - Properly manages file naming and paths
   - Fixes issues with mmseqs output naming:
     ```python
     # Output path that matches actual file structure
     lca_output = f"{config[output_dir]}/03_genomic_info/mmseqs_output_lca.tsv"
     ```

4. **Comprehensive Testing Framework**
   - Dedicated test rules for validating functionality
   - Special focus on phage-specific PHACTS analysis
   - Well-documented testing procedures in workflow/TEST_INSTRUCTIONS.md

## Configuration

Configure the pipeline through a YAML file with options for:

1. **Directory Paths**
   ```yaml
   output_dir: "results"
   assembly_file: "/path/to/assembly.fasta"  # Optional
   assembly_graph: "/path/to/assembly_graph.gfa"  # Optional
   reads_dir: "/path/to/reads"
   ```

2. **Database Locations**
   ```yaml
   databases:
     mmseqs2:
       db: "/path/to/mmseqs_db"
     genomad:
       db: "/path/to/genomad_db"
     checkv:
       db: "/path/to/checkv_db"
     # Additional database paths...
   ```

3. **Tool Parameters**
   ```yaml
   cluster:
     enabled: true
     min_length: 10000
     identity: 0.95
     coverage: 0.85
   ```

## Usage Examples

### Basic Execution
```bash
snakemake --use-conda --cores 24 --configfile config/config.yaml
```

### Running Specific Stages
```bash
# Just run phage prediction
snakemake --use-conda --cores 24 results/01_phage_predictions/phageContigs.fasta

# Run clustering
snakemake --use-conda --cores 24 results/02_clustering/vOTU_repSeqs.fasta

# Run host prediction
snakemake --use-conda --cores 24 results/03_iphop_results/iphop_predictions_compiled.tsv

# Run phage-specific PHACTS analysis
snakemake --use-conda --cores 24 results/03_phacts_results_by_phage/phacts_predictions_compiled.tsv
```

### Testing Phage-Specific PHACTS
```bash
# Navigate to workflow directory
cd workflow/

# Run test with a protein file
snakemake --configfile ../test_data/test_config.yaml --use-conda --cores 4 \
  test_phage_specific/results/test_report.txt
```