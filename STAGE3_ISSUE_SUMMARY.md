# Stage 3 Output Issue Summary

## The Problem

Your `test_bacphlip/` directory contains **mixed data from multiple pipeline runs**:

### Files with CORRECT sequence counts (28,375):
✓ `03_checkv_final/quality_summary.tsv` - Created Oct 21
✓ `03_orf_predictions/proteins.faa` - Created Oct 20
✓ `03_genomic_info/bacphlip_lifestyle.tsv` - Recent
✓ `final_contig_summary.tsv` - Created Oct 22 (most recent)

### Files with INCORRECT sequence counts (from old runs):
✗ `03_iphop_results/iphop_predictions_compiled.tsv` - **72,351 sequences** (should be 28,375)
✗ `03_genomic_info/mmseqs_taxonomy.tsv` - **190,747 sequences** (should be 28,375)

## Why This Happened

When you reused the `test_bacphlip/` directory, **Snakemake didn't regenerate all files** because:
1. Some files already existed from previous runs
2. Snakemake's logic determined they didn't need to be remade
3. The files have different timestamps showing they're from different runs

## How to Verify What Actually Worked

Run the updated verification script:

```bash
bash verify_stage3_outputs.sh test_bacphlip/ full_test_with_tax_consensus/02_clustering/vOTU_repSeqs.fasta
```

The script will now check the **correct file locations** and show:
- Which files have the right number of sequences (28,375)
- Which files are stale from old runs

## Key Questions to Answer

### 1. Check MMseqs2 sequence IDs
```bash
# See if MMseqs2 results match your actual input
head -20 test_bacphlip/03_genomic_info/mmseqs_taxonomy.tsv
```

Compare the sequence IDs with:
```bash
head -20 full_test_with_tax_consensus/02_clustering/vOTU_repSeqs.fasta | grep ">"
```

**Do they match?** If not, MMseqs2 file is from a previous run.

### 2. Check iPhop sequence IDs
```bash
# See if iPhop results match your actual input
head -20 test_bacphlip/03_iphop_results/iphop_predictions_compiled.tsv
```

**Do the sequence IDs match your input?** If not, it's from a previous run.

### 3. Check the final summary
```bash
# This should be correct since it was created most recently
head -20 test_bacphlip/final_contig_summary.tsv
```

## The Real Question

**Does the final_contig_summary.tsv contain the correct data for your 28,375 sequences?**

This is the most important file. Check:
1. Does it have 28,375 rows (+ 1 header)?
2. Do the sequence IDs match your input?
3. Does it have taxonomy and lifestyle predictions filled in?

```bash
# Count rows
wc -l test_bacphlip/final_contig_summary.tsv

# Check first sequence ID
head -2 test_bacphlip/final_contig_summary.tsv | tail -1 | cut -f1

# Compare with your input
head -1 full_test_with_tax_consensus/02_clustering/vOTU_repSeqs.fasta
```

## What to Do Next

### Option 1: If final_contig_summary.tsv looks correct
✓ Your pipeline DID work correctly
✓ The old iPhop/MMseqs2 files are just leftover artifacts
✓ The final summary integrates recent results from files that ran correctly

### Option 2: If final_contig_summary.tsv has wrong data
✗ Need to rerun with a **fresh output directory**:

```bash
# Use a completely new directory name
snakemake --config \
    start_from=clustering \
    input_clustered_seqs=full_test_with_tax_consensus/02_clustering/vOTU_repSeqs.fasta \
    outdir=test_bacphlip_CLEAN \
    ...
```

### Option 3: Force regenerate specific files
If you want to keep the directory but fix specific files:

```bash
# Remove stale files
rm test_bacphlip/03_iphop_results/iphop_predictions_compiled.tsv
rm test_bacphlip/03_genomic_info/mmseqs_taxonomy.tsv

# Rerun to regenerate just those files
snakemake --forcerun iphop_aggregate_results mmseqs_phage_taxonomy ...
```

## My Recommendation

**Check the final_contig_summary.tsv first.**

If the sequence IDs and data look correct for your 28,375 sequences, then your pipeline worked fine - you just have some leftover files from previous runs that don't matter.

The fact that:
- BACPHLIP has exactly 28,375 sequences ✓
- Final summary has exactly 28,375 sequences ✓
- CheckV has the right date (Oct 21) ✓
- Final summary has the most recent date (Oct 22) ✓

...suggests that the **actual analysis DID run correctly** and just integrated data from the correct recent files, not the stale ones.

Run this to check:
```bash
# Extract sequence IDs from final summary
tail -n +2 test_bacphlip/final_contig_summary.tsv | cut -f1 | head -10

# Compare with your input
grep ">" full_test_with_tax_consensus/02_clustering/vOTU_repSeqs.fasta | head -10 | sed 's/>//'
```

If these match, **you're good to go!** The stale files don't matter.
