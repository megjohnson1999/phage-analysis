# Entry Point Testing Guide

## Quick Reference: What Each Entry Point Skips

| Entry Point | Skips | Runs | Expected Jobs |
|-------------|-------|------|---------------|
| `assembly` | Nothing | Everything | ~50+ |
| `reneo_output` | Reneo binning | Filtering → Prediction → Clustering → Analysis | ~47 |
| `viral_contigs` | Reneo + MMseqs filtering | Prediction → Clustering → Analysis | ~40-45 |
| `predicted_phages` | Reneo + Filtering + Prediction | Clustering → Analysis | ~25-30 |
| `clustered_sequences` | Reneo + Filtering + Prediction + Clustering | Analysis only | ~20-25 |

## Testing Each Entry Point

### 1. Test `assembly` (default - everything)

```bash
cd workflow
snakemake -n --cores 1
```

**Key jobs to look for:**
- ✅ genomad_prediction
- ✅ jaeger_prediction
- ✅ filter_mmseqs_lca
- ✅ cluster_phages
- ✅ iphop_aggregate_results
- ✅ generate_summary_report

**Expected:** All ~50+ jobs

---

### 2. Test `reneo_output` (skip Reneo binning)

**What you need:** A FASTA file of contigs (min 1000bp)

```bash
cd workflow
snakemake -n --configfile ../config/examples/config_from_reneo_output.yaml
```

**Key jobs to look for:**
- ✅ filter_mmseqs_lca (still filters your contigs)
- ✅ genomad_prediction (still predicts phages)
- ✅ cluster_phages
- ❌ NO Reneo-specific binning rules

**Expected:** ~47 jobs (you already tested this ✓)

---

### 3. Test `viral_contigs` (skip filtering)

**What you need:** A FASTA file of viral contigs (min 1000bp)

**Setup:**
```bash
# Create test config with your viral contigs file
cp config/examples/config_from_predicted_phages.yaml config/test_viral_contigs.yaml
nano config/test_viral_contigs.yaml
# Change:
# - start_from: "viral_contigs"
# - viral_contigs_file: "/path/to/your/viral_contigs.fasta"
```

```bash
cd workflow
snakemake -n --configfile ../config/test_viral_contigs.yaml
```

**Key jobs to look for:**
- ❌ NO filter_mmseqs_lca (we're providing already-filtered viral contigs)
- ✅ genomad_prediction (still predicts phages from viral contigs)
- ✅ jaeger_prediction
- ✅ cluster_phages
- ✅ All analysis steps

**Expected:** ~40-45 jobs

---

### 4. Test `predicted_phages` (skip prediction)

**What you need:** A FASTA file of predicted phage contigs

```bash
cd workflow
snakemake -n --configfile ../config/examples/config_from_predicted_phages.yaml
```

**Key jobs to look for:**
- ❌ NO genomad_prediction (we're providing already-predicted phages)
- ❌ NO jaeger_prediction
- ❌ NO filter_mmseqs_lca
- ✅ cluster_phages (cluster your provided phages)
- ✅ All analysis steps (iphop, phabox, checkv, vcontact3)
- ✅ generate_summary_report

**Expected:** ~25-30 jobs

---

### 5. Test `clustered_sequences` (skip clustering - analysis only)

**What you need:** A FASTA file of vOTU representative sequences

```bash
cd workflow
snakemake -n --configfile ../config/examples/config_from_clustered_sequences.yaml
```

**Key jobs to look for:**
- ❌ NO genomad_prediction
- ❌ NO cluster_phages (we're providing pre-clustered vOTUs)
- ✅ iphop_aggregate_results (host prediction)
- ✅ phabox_prediction (taxonomy)
- ✅ checkv_final_assessment (quality)
- ✅ vcontact3_taxonomy
- ✅ generate_summary_report
- ✅ (Optional) abundance calculation if enabled

**Expected:** ~20-25 jobs (fewest jobs)

---

## Common Issues

### File validation errors
If you get "file does not exist" errors during dry-run, you need to:
1. Update the config with real file paths, OR
2. Create dummy files for testing

### Create dummy test files:
```bash
# Create test files for dry-run validation
mkdir -p /tmp/phage_test
touch /tmp/phage_test/reneo_contigs.fasta
touch /tmp/phage_test/viral_contigs.fasta
touch /tmp/phage_test/predicted_phages.fasta
touch /tmp/phage_test/vOTUs.fasta

# Then update your test configs to point to these
```

---

## Verification Checklist

For each entry point test:
- [ ] Dry-run completes without errors
- [ ] Job count matches expected range
- [ ] Skipped jobs are NOT in the list
- [ ] Required jobs ARE in the list
- [ ] Pipeline prints: "Pipeline entry point: [your_entry_point]"
