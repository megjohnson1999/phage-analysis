# Batch Processing in Phage Analysis Pipeline

This guide explains the batch processing approach implemented in the pipeline to efficiently handle large numbers of sequences when running computationally intensive tools.

## Overview

For tools that process individual sequences or have high computational requirements, we've implemented a batch processing framework that:

1. Splits input sequences into chunks
2. Processes each chunk in parallel
3. Aggregates results back into a single output file

This approach is particularly valuable for tools like PHOLD, iPhop, and PHACTS that can create thousands of jobs when processing large metagenomic datasets.

## Implementation Strategy

The batch processing framework consists of several components:

### 1. Sequence Splitting

Input sequences are split into manageable chunks rather than individual files:

```python
rule split_sequences_for_tool:
    input:
        sequences = "input.fasta"
    output:
        split_dir = directory("split_seqs_dir"),
        split_list = "split_file_list.txt"
    params:
        # Number of sequences per chunk - adjust based on expected sequence sizes
        chunk_size = 1000
    shell:
        """
        # Split FASTA file into chunked files with multiple sequences per file
        seqkit split {input.sequences} -O {output.split_dir} -p {params.chunk_size}
        
        # Create list of split files - use absolute paths for reliability
        find {output.split_dir} -name "*.fasta" -type f | sort > {output.split_list}
        """
```

### 2. Dynamic Sample Management

A helper function gets the list of chunk files to process:

```python
def get_samples():
    split_list = "split_file_list.txt"
    
    # Try reading from the split file list if it exists
    if os.path.exists(split_list):
        with open(split_list, "r") as f:
            files = [line.strip() for line in f if line.strip()]
        return [os.path.splitext(os.path.basename(file))[0] for file in files]
    
    # Fallback mechanisms for when the file list doesn't exist yet
    # ...
    
    # Return empty list if all methods fail
    return []
```

### 3. Checkpoints for Dynamic Workflow

Snakemake checkpoints enable dynamic dependency handling:

```python
checkpoint wait_for_splits:
    input:
        split_list = "split_file_list.txt"
    output:
        touch(".splits_ready")
```

### 4. Parallel Processing Rules

Each chunk is processed independently:

```python
rule process_single_chunk:
    input:
        checkpoint = ".splits_ready",
        chunk_file = "split_seqs_dir/{sample}.fasta"
    output:
        results = "tmp_results/{sample}/results.tsv"
    threads: 24  # Can be adjusted per rule
    shell:
        """
        # Run tool on chunk of sequences
        tool run -i {input.chunk_file} -o $(dirname {output.results}) -t {threads}
        """
```

### 5. Helper Rule for Dry Runs

A helper rule ensures all jobs will be scheduled:

```python
rule run_all_predictions:
    input:
        checkpoint = ".splits_ready",
        samples = lambda wildcards: expand(
            "tmp_results/{sample}/results.tsv",
            sample=get_samples()
        )
    output:
        touch(".all_predictions_done")
```

### 6. Results Aggregation

After all chunks are processed, results are aggregated:

```python
rule aggregate_results:
    input:
        all_done = ".all_predictions_done",
        predictions = lambda wildcards: expand(
            "tmp_results/{sample}/results.tsv",
            sample=get_samples()
        )
    output:
        compiled_results = "final_results.tsv"
    shell:
        """
        # Extract header from first file
        head -n 1 {input.predictions[0]} > {output.compiled_results}
        
        # Append all data (skipping headers) to output file
        for file in {input.predictions}; do
            awk 'NR>1' "$file" >> {output.compiled_results}
        done
        """
```

## Benefits of Batch Processing

- **Reduced job management overhead**: Fewer jobs to track and manage
- **Better resource utilization**: Each batch can use more resources effectively
- **Faster processing**: Eliminates per-file startup and initialization time
- **Reduced filesystem load**: Fewer small files created during processing
- **Simplified workflow**: Cleaner code and dependency structure

## Implemented Tools

Batch processing has been implemented for:

1. **PHOLD**: Phage protein function annotation
2. **iPhop**: Phage host prediction 
3. **PHACTS**: Phage lifestyle prediction

Each implementation follows the same pattern but with tool-specific adjustments for input/output formats and resource requirements.

## Adjusting Batch Size

The `chunk_size` parameter in the splitting rule controls how many sequences are included in each batch. This can be adjusted based on:

- Sequence length (longer sequences may need smaller batches)
- Tool memory requirements
- Available computational resources
- Total number of sequences

For most use cases, a batch size of 500-1000 sequences provides a good balance between parallelization and overhead reduction.