# Phage-Specific PHACTS Analysis

This is an enhanced version of the workflow that implements phage-specific grouping for PHACTS lifestyle prediction, addressing a significant limitation in the previous approach.

## Problem Being Solved

In the original workflow design:
1. Prodigal predicts proteins from all phages, creating a single `proteins.faa` file
2. This file is split into arbitrary batches based on sequence count, not phage identity
3. PHACTS makes a single lifestyle prediction for each batch
4. The output aggregation assigns this single prediction to all phage IDs found in the batch

This approach has several significant limitations:
- A single batch often contained proteins from multiple phages with potentially different lifestyles
- The arbitrary grouping broke the biological association between proteins from the same phage
- Each batch received a single lifestyle prediction that was inappropriately applied to all phages in that batch
- This reduced the accuracy and biological relevance of the predictions

## Solution: Phage-Specific Approach

This enhanced workflow addresses these issues by implementing a phage-specific approach:

1. **Phage ID Extraction**: The workflow extracts the phage identifier from each protein sequence header
2. **Per-phage Grouping**: Proteins are grouped by their source phage, ensuring all proteins from the same phage stay together
3. **Dedicated Files**: Separate files are created for each phage, named after the phage ID
4. **Independent Analysis**: PHACTS predictions are run on each phage-specific file separately
5. **Clear Mapping**: Results are aggregated with a direct 1:1 mapping between phage IDs and lifestyle predictions

## Technical Implementation

### Extraction of Phage IDs

Phage IDs are extracted from protein headers using a pattern matching approach:

```python
def extract_phage_id(header):
    """
    Extract phage ID from a FASTA header.
    
    Args:
        header: The FASTA header string (with or without leading '>')
        
    Returns:
        str: The extracted phage ID
    """
    # Extract the part before the last underscore and number
    # For headers like "read_13355_7", this will extract "read_13355"
    match = re.match(r'^>?([^_]+(?:_[^_]+)*?)_\d+', header)
    if match:
        return match.group(1).replace('>', '')
    
    # Fallback: use the entire header up to the first space/hash
    match = re.match(r'^>?([^ #]+)', header)
    if match:
        return match.group(1).replace('>', '')
        
    return "unknown"
```

This function handles common header formats and extracts the phage ID component:
- For headers like `>read_13355_7`, it extracts `read_13355` (the part before the last underscore and number)
- For other formats, it falls back to extracting everything up to the first space or hash character
- If no match is found, it returns "unknown" (these are filtered out later)

### File Organization

```python
def split_proteins(input_file, output_dir):
    """
    Split protein FASTA file by phage ID.
    
    Args:
        input_file: Path to input FASTA file
        output_dir: Directory to write split files
        
    Returns:
        list: Paths to created output files
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Group proteins by phage ID
    phage_proteins = {}
    
    # Parse the input file
    for record in SeqIO.parse(input_file, "fasta"):
        phage_id = extract_phage_id(record.id)
        
        if phage_id not in phage_proteins:
            phage_proteins[phage_id] = []
        
        phage_proteins[phage_id].append(record)
    
    # Write a file for each phage
    output_files = []
    for phage_id, records in sorted(phage_proteins.items()):
        if phage_id == "unknown":
            continue  # Skip records where we couldn't determine phage ID
            
        output_file = os.path.join(output_dir, f"{phage_id}.faa")
        SeqIO.write(records, output_file, "fasta")
        output_files.append(output_file)
    
    return output_files
```

This function:
1. Parses the original protein file using Biopython
2. Groups sequences by their phage ID in a dictionary
3. Writes separate FASTA files for each phage
4. Returns a list of the output files created

## Workflow Integration

The phage-specific approach is now the primary method in the workflow:

1. A rule `split_proteins_by_phage` creates phage-specific protein files
2. A function `get_phacts_phage_samples` lists all phage samples
3. Rules for processing each phage file run independently and in parallel
4. An aggregation step combines all phage-specific results

## Header Format Requirements

The script handles protein headers from Prodigal output, which typically follow formats such as:

```
>phage1_42 # 4873 # 5514 # 1 # ID=8902_7;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.352
```

The phage ID extraction is designed to be flexible:
1. First attempts to extract the part before the last underscore and number
2. If that fails, falls back to using everything before the first space or hash
3. The extraction is robust to headers with or without the leading '>' character

## Expected Improvements

This phage-specific approach yields several important benefits:

1. **Biological Relevance**: Predictions are made on biologically meaningful units (individual phages)
2. **Improved Accuracy**: By analyzing proteins from the same phage together, PHACTS can generate more coherent predictions
3. **Clear Attribution**: Each phage has a single, definitive lifestyle prediction
4. **Higher Confidence**: Using all proteins from a phage provides more comprehensive data for prediction
5. **Better Interpretability**: Results are easier to interpret and use in downstream analyses

## Output Format

The final output file `03_phacts_results_by_phage/phacts_predictions_compiled.tsv` contains:

```
phage_id    lifestyle    probability
phage1      lytic        0.85
phage2      temperate    0.76
phage3      lytic        0.92
...
```

Each row represents a single phage with:
- Its unique identifier (extracted from protein headers)
- The predicted lifestyle (typically "lytic" or "temperate")
- The probability/confidence score for that prediction (0-1)

## Usage

The workflow now uses this phage-specific approach by default. Simply run the Snakemake workflow as usual:

```bash
snakemake --profile profile/slurm --config assembly_file="/path/to/assembly.fasta" reads_dir="/path/to/reads/" output_dir="/path/to/output/"
```

The results will be available in the output directory:
```
{output_dir}/03_phacts_results_by_phage/phacts_predictions_compiled.tsv
```

## Testing

You can test the phage-specific PHACTS analysis with a specific protein file:

```bash
# Navigate to the workflow directory
cd workflow/

# Run the test with your protein file
snakemake --use-conda --cores 4 --config input_path=/path/to/your/protein/file.faa -- test_phage_specific/results/test_report.txt
```

## Future Enhancements

Potential improvements to this approach include:

1. More sophisticated phage ID extraction for complex header formats
2. Integration with taxonomic information to enhance prediction accuracy
3. Statistical analysis of prediction confidence based on protein count per phage
4. Addition of visualization tools for comparing phage lifestyles across samples

## Conclusion

The phage-specific PHACTS analysis represents a significant improvement in the workflow's ability to predict phage lifestyles. By respecting the biological origin of proteins and making predictions at the individual phage level, the results are more accurate, more meaningful, and more useful for downstream analyses.