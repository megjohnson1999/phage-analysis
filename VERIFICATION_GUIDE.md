# Stage 3 Output Verification Guide

This guide helps you verify that your phage analysis pipeline ran correctly when starting from clustered sequences (Stage 2 output → Stage 3 analysis).

## Quick Start

### 1. Upload the verification script to your HPC

```bash
# From your local machine
scp verify_stage3_outputs.sh your_username@hpc.university.edu:/path/to/your/project/
```

### 2. Run the verification script on HPC

```bash
# SSH to your HPC
ssh your_username@hpc.university.edu

# Navigate to your project directory
cd /path/to/your/project/

# Run the verification script
bash verify_stage3_outputs.sh /path/to/output/directory /path/to/input/vOTU_repSeqs.fasta
```

**Note:** The second argument (input FASTA file) is optional but recommended for comparing input vs output sequence counts.

### 3. Interpret the results

The script will check:
- ✓ All expected output files exist
- ✓ Files contain data (not empty)
- ✓ Sequence counts are consistent
- ✓ Data quality looks reasonable
- ⚠ Log files for errors
- ✓ Sequence IDs match across all outputs

## Expected Outputs

When starting from **clustered sequences** (Stage 2) and running **Stage 3 (Comprehensive Analysis)**, you should have:

### Quality Assessment
- `03_genomic_info/quality_summary.tsv` - CheckV completeness and quality metrics

### Host Prediction
- `03_iphop_results/iphop_predictions_compiled.tsv` - Bacterial host predictions

### Protein Prediction
- `03_genomic_info/proteins.faa` - Predicted protein sequences
- `03_genomic_info/genes.fna` - Predicted gene sequences

### Taxonomy (3 tools)
- `03_genomic_info/mmseqs_taxonomy.tsv` - MMseqs2 taxonomy
- `03_genomic_info/phabox_output/taxonomy.tsv` - Phabox2 taxonomy
- `03_genomic_info/vc3_output/` - vContact3 results directory
- `03_genomic_info/consensus_taxonomy.tsv` - **Consensus from all 3 tools**

### Lifestyle Prediction (2 tools)
- `03_genomic_info/phabox_output/lifestyle.tsv` - Phabox2 lifestyle
- `03_genomic_info/bacphlip_lifestyle.tsv` - BACPHLIP lifestyle
- `03_genomic_info/bacphlip_lifestyle_with_completeness.tsv` - BACPHLIP + CheckV completeness
- `03_genomic_info/lifestyle_consensus.tsv` - **Consensus from both tools**

### Final Integrated Results
- `final_contig_summary.tsv` - **Complete table with all characterization results**
- `Pipeline_Summary_Report.html` - Interactive HTML report (if enabled)

## What Success Looks Like

```
✓ ALL CHECKS PASSED!

Your Stage 3 outputs appear complete and consistent.
```

**Key indicators:**
- All 12+ expected files exist
- Sequence counts match across all outputs
- Final summary table has same number of sequences as input
- No error messages in log files
- Reasonable data distributions (not all "NA" or "No prediction")

## Troubleshooting Common Issues

### Missing Files
```
✗ MISSING: 03_genomic_info/consensus_taxonomy.tsv
```

**Solution:**
- Check if the pipeline completed: `tail -100 .snakemake/log/*.log`
- Look for failed Snakemake rules
- Re-run from the failed step if needed

### Empty Files
```
⚠ WARNING: final_contig_summary.tsv appears empty (1 lines)
```

**Solution:**
- Check tool-specific log files in `03_*/` directories
- Verify input sequences were valid FASTA format
- Check if tools failed due to database or dependency issues

### Sequence Count Mismatch
```
✗ Sequence count mismatch detected!
  CheckV: 150 sequences
  Summary: 145 sequences
```

**Solution:**
- Some sequences may have been filtered by individual tools
- Check which sequences are missing: the script shows first 5 IDs
- Review tool logs to understand why sequences were excluded

### Errors in Log Files
```
⚠ Found potential errors in: iphop_chunk_01.log
```

**Solution:**
- Review the specific log file mentioned
- Common issues:
  - Memory errors → increase memory allocation in SLURM/Snakemake config
  - Database errors → verify database paths in config.yaml
  - Timeout errors → increase time limits

## Manual Inspection

After the automated checks, you may want to manually inspect key outputs:

### 1. Check the final summary table
```bash
# View first few lines
head -20 /path/to/output/final_contig_summary.tsv | column -t -s $'\t'

# Count total sequences
tail -n +2 /path/to/output/final_contig_summary.tsv | wc -l
```

### 2. Review CheckV quality
```bash
# See completeness distribution
tail -n +2 /path/to/output/03_genomic_info/quality_summary.tsv | cut -f7 | sort | uniq -c
```

### 3. Check host prediction success rate
```bash
# Count sequences with host predictions
tail -n +2 /path/to/output/03_iphop_results/iphop_predictions_compiled.tsv | \
  awk -F'\t' '{if ($2 != "" && $2 != "NA") print $0}' | wc -l
```

### 4. Review taxonomy assignments
```bash
# See family-level taxonomy
tail -n +2 /path/to/output/03_genomic_info/consensus_taxonomy.tsv | cut -f5 | sort | uniq -c | sort -rn
```

### 5. Check lifestyle predictions
```bash
# Temperate vs virulent distribution
tail -n +2 /path/to/output/03_genomic_info/lifestyle_consensus.tsv | cut -f2 | sort | uniq -c
```

## Additional Resources

For more details on the pipeline workflow:
- See `WORKFLOW_EXECUTION_ORDER.md` for complete pipeline structure
- See `DIAGRAM_OUTLINES.md` for workflow diagrams
- See `README.md` for general pipeline documentation

## Questions or Issues?

If the verification script reports issues you can't resolve:
1. Check Snakemake log files: `.snakemake/log/`
2. Check individual tool logs in output subdirectories
3. Verify all required databases are accessible
4. Ensure sufficient compute resources (memory, time) were allocated

Common log locations:
- Main Snakemake log: `.snakemake/log/`
- Tool-specific logs: Often in the output subdirectories (e.g., `03_iphop_results/`)
- SLURM logs: Usually `slurm-*.out` files in your submission directory
