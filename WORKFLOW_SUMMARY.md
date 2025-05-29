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
- **Assembly File**: Direct input of assembled contigs (FASTA) - **Currently recommended**
- **Assembly Graph**: For graph-based assembly refinement with Reneo - *Currently being troubleshooted*
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

4. **Lifestyle Prediction**
   - **BACPHLIP**: Machine learning-based lifestyle prediction
     - Analyzes genomic features to classify phages as temperate or virulent
     - Runs directly on phage contigs
   - **Phabox2**: Comprehensive phage analysis including lifestyle prediction
     - Provides lytic/lysogenic classification
     - Integrates multiple prediction approaches

5. **Taxonomic Classification**
   - **MMSeqs2**: Detailed sequence-based taxonomic classification
   - **Phabox2**: Machine learning-based phage taxonomy at multiple levels
   - **vContact3**: Gene content-based taxonomy using network clustering
     - Creates gene-sharing networks from protein sequences
     - Assigns taxonomy based on network relationships

## Implementation Details

### Parallelization Strategy

The workflow implements efficient parallelization:

1. **Smart Chunking**
   - Divides large tasks into manageable pieces
   - Particularly important for PHOLD and iPhop analyses
   - Default chunk size: 100 sequences per batch

2. **Checkpoint System**
   - Uses Snakemake checkpoints for handling dynamic file generation
   - Three phases for parallel operations:
     ```
     .splits_ready → .input_files_found → .all_predictions_done
     ```
   - Particularly important for iPhop and PHOLD batch processing

3. **Results Aggregation**
   - Merges individual results with proper header handling
   - Accounts for failed or empty results gracefully

### Error Handling

Comprehensive error handling throughout:

1. **Input Validation**
   - Checks for existence of required files
   - Validates configuration settings
   - Automatically detects and clears placeholder paths from default config

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

3. **Reneo Wrapper Script**
   - Special handling for Reneo's expected failures
   - Located at `workflow/scripts/run_reneo_wrapper.sh`
   - Features:
     - Captures Reneo exit codes
     - Checks for expected output files
     - Creates empty outputs if Reneo fails, allowing pipeline continuation
     - Detailed logging of Reneo's execution status
   - Example usage:
     ```bash
     bash {workflow.basedir}/scripts/run_reneo_wrapper.sh \
         --input {input.assembly_graph} \
         --reads {input.reads_dir} \
         --minlength 1000 \
         --output {output_dir} \
         --threads {threads}
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

2. **Database Requirements**
   - Each tool requires specific databases
   - Configured in the config.yaml file
   - Includes paths for MMSeqs2, geNomad, CheckV, PHOLD, iPhop, Phabox2, and vContact3

## Key Features and Improvements

1. **Multiple Lifestyle Prediction Tools**
   - BACPHLIP for machine learning-based predictions
   - Phabox2 for integrated analysis
   - Provides complementary predictions for higher confidence

2. **Flexible Input Options**
   - Supports both direct assembly files and assembly graphs
   - Conditional workflow paths based on input types
   - Fixed DAG construction for GFA-only inputs:
     - All rules are always defined (no conditional rule definitions)
     - Uses input helper functions that return dummy inputs when rules shouldn't run
     - Rules check for dummy inputs and create empty outputs when skipped

3. **Robust Output Handling**
   - Properly manages file naming and paths
   - Fixes issues with mmseqs output naming:
     ```python
     # Output path that matches actual file structure
     lca_output = f"{config[output_dir]}/03_genomic_info/mmseqs_output_lca.tsv"
     ```

4. **Comprehensive Tool Integration**
   - Seamless integration of multiple prediction and analysis tools
   - Unified output format for easy downstream analysis
   - Parallel processing for computationally intensive steps

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

# Run BACPHLIP lifestyle prediction
snakemake --use-conda --cores 24 results/03_bacphlip/bacphlip.predictions.tsv

# Run Phabox2 analysis
snakemake --use-conda --cores 24 results/03_genomic_info/phabox_output/lifestyle.tsv
```

### Running Complete Analysis
```bash
# Navigate to workflow directory
cd workflow/

# Run full pipeline with FASTA input (recommended)
snakemake --profile ../profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
```