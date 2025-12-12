# Stage 3: Analysis & Characterization - Detailed Breakdown

## Overview
Stage 3 performs comprehensive analysis and characterization of phage sequences. This is the most complex stage with extensive parallel processing for host prediction, taxonomic classification, lifestyle prediction, and functional annotation.

---

## Purpose & Scope

### What Stage 3 Does:
1. **Quality Assessment** - Final CheckV validation
2. **Host Prediction** - Identify bacterial hosts (iPhop)
3. **Taxonomic Classification** - Multiple complementary methods
4. **Lifestyle Prediction** - Temperate vs lytic/virulent
5. **Functional Annotation** - ORF/gene prediction
6. **Consensus Building** - Integrate predictions from multiple tools
7. **Final Integration** - Create comprehensive summary table

### Input Flexibility:
Stage 3 works with **either**:
- vOTU representatives (if clustering was enabled)
- All phage contigs (if clustering was skipped)

This is handled automatically by the `get_phage_input()` function.

---

## Workflow Summary

```
Input Sequences
     ↓
┌────────────────────────────────────────────────────┐
│ Substage 3A: Quality Check                        │
│   CheckV Final Assessment                          │
└────────────────┬───────────────────────────────────┘
                 ↓
┌────────────────────────────────────────────────────┐
│ Substage 3B: Host Prediction (PARALLELIZED)       │
│   Split → iPhop Array → Aggregate                 │
└────────────────┬───────────────────────────────────┘
                 ↓
┌────────────────────────────────────────────────────┐
│ Substage 3C: Characterization (PARALLEL TOOLS)    │
│   • Prodigal (ORFs)                                │
│   • MMseqs2 (Taxonomy)                             │
│   • Phabox2 (Taxonomy + Lifestyle)                 │
│   • vContact3 (Gene Content Taxonomy)              │
│   • BACPHLIP (Lifestyle)                           │
└────────────────┬───────────────────────────────────┘
                 ↓
┌────────────────────────────────────────────────────┐
│ Substage 3D: Consensus Building                   │
│   • Lifestyle Consensus (BACPHLIP + Phabox2)      │
│   • Taxonomic Consensus (MMseqs2 + Phabox2 + VC3) │
└────────────────┬───────────────────────────────────┘
                 ↓
┌────────────────────────────────────────────────────┐
│ Substage 3E: Final Integration                    │
│   Create final_contig_summary.tsv                 │
└────────────────────────────────────────────────────┘
```

---

## Input Selection: `get_phage_input()` Function

Located in `workflow/rules/02_clustering.smk`

### Logic Tree:

```python
def get_phage_input(wildcards):
    start_from = config.get("start_from", "raw_contigs")

    # Priority 1: Starting from clustering
    if start_from == "clustering":
        return input_clustered_seqs  # User-provided file

    # Priority 2: Starting from phage_contigs without clustering
    if start_from == "phage_contigs" and not do_clustering:
        return input_phage_contigs  # User-provided file

    # Priority 3: Standard workflow
    if do_clustering:
        return "02_clustering/vOTU_repSeqs.fasta"  # From Stage 2
    else:
        return "01_phage_predictions/phageContigs.fasta"  # From Stage 1
```

### Possible Inputs:

| Scenario | Input File | Typical Size |
|----------|-----------|--------------|
| Clustered | `02_clustering/vOTU_repSeqs.fasta` | 500-2,000 seqs |
| Not clustered | `01_phage_predictions/phageContigs.fasta` | 2,000-10,000 seqs |
| User-provided (phage_contigs) | Custom file | Variable |
| User-provided (clustering) | Custom file | Variable |

---

## Substage 3A: Quality Assessment

### Purpose
Final quality check on the sequences being analyzed (post-clustering).

---

### **`checkv_final_assessment` Rule**

**Tool:** CheckV (same as Stage 1, but different input)

**Why run CheckV again?**
- Stage 1 CheckV ran on **all viral contigs**
- Stage 3 CheckV runs on **final analyzed sequences** (clustered or all phages)
- Different inputs may have different quality distributions
- This CheckV result matches the sequences being analyzed

**Input:** Output from `get_phage_input()` function

**Command:**
```bash
checkv end_to_end phage_seqs.fasta \
    03_checkv_final/ \
    -d /path/to/checkv/db \
    -t 24
```

**Output:** `03_checkv_final/quality_summary.tsv`

**Key Columns:**
- `contig_id`
- `contig_length`
- `provirus` - Yes/No
- `proviral_length`
- `gene_count`
- `viral_genes`
- `host_genes`
- `checkv_quality` - Complete, High-quality, Medium-quality, Low-quality, Not-determined
- `miuvig_quality` - Tier classification
- `completeness` - Percentage (or "Not-determined")
- `completeness_method` - How completeness was estimated
- `contamination` - Estimated contamination level
- `warnings` - Any QC warnings

**Processing Time:** 30-60 minutes for 1,500 sequences (24 threads)

**Used For:**
- Final summary table
- Quality filtering
- Understanding completeness of analyzed sequences

---

## Substage 3B: Host Prediction (PARALLELIZED)

### Purpose
Identify bacterial/archaeal hosts for each phage sequence using multiple prediction methods.

**Why parallelized?** iPhop is computationally intensive; splitting speeds up processing dramatically.

---

### Step 1: **`split_phage_sequences` Rule** - Chunking

**Tool:** seqkit

**What it does:** Splits sequences into fixed-size chunks for parallel processing

**Parameters:**
- `sequences_per_job: 100` (fixed chunk size)

**Strategy:**
```bash
seqkit split2 --by-size 100 phage_seqs.fasta -O 03_split_seqs/
```

**Output:**
- `03_split_seqs/*.fasta` - Multiple FASTA files (100 seqs each)
- `03_split_seqs/split_file_list.txt` - List of all split files

**Example:**
```
Input:  1,500 sequences
Output: 15 chunk files
        • phage_seqs.part_001.fasta (100 seqs)
        • phage_seqs.part_002.fasta (100 seqs)
        • ...
        • phage_seqs.part_015.fasta (100 seqs)
```

**Error Handling:**
- If input is empty, creates empty placeholder
- Continues pipeline gracefully

---

### Step 2: **Checkpoint** - Synchronization

**`wait_for_iphop_splits` checkpoint**

**Purpose:** Ensures all split files are created before launching parallel jobs

**How Snakemake Checkpoints Work:**
1. Workflow pauses at checkpoint
2. Evaluates what files were created
3. Dynamically creates jobs based on discovered files
4. Continues execution

**Why needed?** Number of chunks unknown until runtime (depends on input size)

---

### Step 3: **`check_iphop_input_files` Rule** - Validation

**What it does:** Verifies iPhop input files exist and are accessible

**Checks:**
- Split file list exists
- All split files are readable
- File paths are valid

**Output:** `.input_files_found` touchfile

**Purpose:** Early failure detection before launching expensive iPhop jobs

---

### Step 4: **`iphop_single_prediction` Rule** - Parallel Host Prediction

**Tool:** iPhop (Integrated Phage Host Predictor)

**What it does:**
- Predicts bacterial/archaeal hosts for phages
- Uses multiple complementary methods:
  - BLAST homology
  - k-mer similarity
  - CRISPR spacer matching
  - Transfer RNA matching
  - Prophage region matching

**Input:** Single chunk file (100 sequences)

**Command:**
```bash
iphop predict --fa_file chunk.fasta \
    --db_dir /path/to/iphop/db \
    --out_dir 03_iphop_results/tmp/{chunk_id}/ \
    --num_threads 24
```

**Resources per job:**
- Threads: 24
- Memory: ~50-100 GB (depends on database)
- Time: ~30-60 min per chunk

**Output:** `03_iphop_results/tmp/{chunk_id}/Host_prediction_to_genus_m90.csv`

**Output Columns (key ones):**
- `Virus` - Phage contig ID
- `Host genus` - Predicted host genus
- `Host taxonomy` - Full host lineage
- `Confidence score` - Overall prediction confidence
- `Method` - Which method(s) supported prediction
- `# of predicted hosts` - Number of possible hosts

**Error Handling:**
- If iPhop fails, creates empty output file
- Allows pipeline to continue
- Aggregation step handles empty files gracefully

---

### Step 5: **`run_all_iphop_predictions` Rule** - Coordination

**Purpose:** Aggregate trigger that waits for all parallel jobs

**Input:** All chunk predictions from `get_iphop_samples()`

**Output:** `.all_predictions_done` touchfile

**Why needed?** Ensures all predictions complete before aggregation

---

### Step 6: **`iphop_aggregate_results` Rule** - Compilation

**What it does:** Merges all chunk predictions into single file

**Process:**
1. Find all prediction CSV files in tmp directory
2. Extract header from first file
3. Concatenate all data (skip headers from subsequent files)
4. Create final compiled file

**Input:** All `Host_prediction_to_genus_m90.csv` files from chunks

**Command (simplified):**
```bash
# Get header
head -n 1 first_file.csv > compiled.csv

# Append all data (skip headers)
for file in */Host_prediction_to_genus*.csv; do
    tail -n +2 "$file" >> compiled.csv
done
```

**Output:** `03_iphop_results/iphop_predictions_compiled.tsv`

**Format:** TSV with columns:
- Virus ID
- Host predictions
- Confidence scores
- Supporting methods

**Processing:** Handles both CSV and TSV formats from iPhop

**Typical Results:**
- Input: 1,500 phages across 15 chunks
- Output: 1,500 rows in compiled file
- Prediction rate: 60-80% (varies by database and phage novelty)

---

## Substage 3C: Genomic Characterization (PARALLEL)

### Purpose
Characterize phages using multiple complementary approaches. All tools run **in parallel** as they're independent.

---

### Tool 1: **`prodigal_orf_prediction`** - Gene Calling

**What it does:** Predicts open reading frames (ORFs) and protein-coding genes

**Tool:** Prodigal (metagenomic mode)

**Input:** Phage sequences

**Command:**
```bash
prodigal -i phage_seqs.fasta \
    -o genes.gff \
    -a proteins.faa \
    -d nucleotides.fna \
    -p meta \
    -q
```

**Parameters:**
- `-p meta` - Metagenomic mode (for diverse/fragmented sequences)
- `-q` - Quiet mode
- `-a` - Output amino acid sequences
- `-d` - Output nucleotide sequences

**Outputs:**
- `03_orf_predictions/genes.gff` - Gene coordinates
- `03_orf_predictions/proteins.faa` - Protein sequences
- `03_orf_predictions/nucleotides.fna` - Gene nucleotide sequences

**Used For:**
- Functional annotation
- Gene count statistics
- Protein-based analyses

**Processing Time:** 5-15 minutes for 1,500 sequences

---

### Tool 2: **`mmseqs_phage_taxonomy`** - Protein-Based Taxonomy

**What it does:** Taxonomic classification using protein similarity against NCBI NR

**Tool:** MMseqs2 easy-taxonomy

**Input:** Phage sequences

**Command:**
```bash
mmseqs easy-taxonomy phage_seqs.fasta \
    /path/to/nr/db \
    03_genomic_info/mmseqs_taxonomy \
    /tmp \
    --lca-mode 2 \
    --tax-lineage 2 \
    --threads 24
```

**Parameters:**
- `--lca-mode 2` - 2bLCA algorithm (balanced sensitivity/specificity)
- `--tax-lineage 2` - Include full taxonomic lineage
- Search parameters: Similar to Stage 1 MMseqs2

**Output:** `03_genomic_info/mmseqs_taxonomy.tsv`

**Columns:**
- Query sequence ID
- Taxonomic assignment
- Lineage string
- Confidence metrics

**Difference from Stage 1 MMseqs2:**
- **Stage 1:** Filtered all contigs, looking for "Viruses"
- **Stage 3:** Classifies known phages to finer taxonomic levels

**Processing Time:** 1-3 hours for 1,500 sequences (24 threads)

---

### Tool 3: **`phabox_prediction`** - ML-Based Phage Characterization

**What it does:**
- Taxonomic classification using machine learning
- Lifestyle prediction (temperate vs lytic)
- Uses phage-specific features

**Tool:** Phabox2 (end-to-end mode)

**Input:** Phage sequences

**Command:**
```bash
phabox2 --task end_to_end \
    --dbdir /path/to/phabox/db \
    --contigs phage_seqs.fasta \
    --outpth 03_genomic_info/phabox_output \
    --len 1000 \
    --threads 24
```

**What Phabox2 Does:**
1. **Virus Identification** - Confirms sequences are viral (redundant here)
2. **Taxonomy Prediction** - ML-based classification
   - Uses protein/gene content features
   - Trained on phage-specific databases
   - Provides confidence scores
3. **Lifestyle Prediction** - Temperate vs Lytic
   - PhaTYP model
   - Considers integrase genes, lysogeny markers
4. **Host Prediction** - CHERRY model (not used in this pipeline)

**Outputs:**

**`03_genomic_info/phabox_output/taxonomy.tsv`**
```
contig_id    taxonomy_prediction          confidence
phage_001    Caudovirales;Siphoviridae   0.87
phage_002    Caudovirales;Myoviridae     0.92
phage_003    unclassified                0.23
```

**`03_genomic_info/phabox_output/lifestyle.tsv`**
```
contig_id    lifestyle_prediction    confidence
phage_001    temperate               0.78
phage_002    virulent                0.91
phage_003    unknown                 0.34
```

**Lifestyle Categories:**
- **Temperate** - Can integrate into host genome (lysogenic)
- **Virulent/Lytic** - Only lytic cycle
- **Unknown/Uncertain** - Cannot determine

**Format Handling:**
- Script handles both new format (final_prediction_summary.tsv)
- Falls back to legacy format (phamer_prediction.csv, cherry_prediction.csv)
- Always produces standardized TSV outputs

**Processing Time:** 1-2 hours for 1,500 sequences (24 threads)

---

### Tool 4: **`vcontact3_taxonomy`** - Gene Content Clustering

**What it does:**
- Clusters phages based on shared gene content
- Assigns taxonomy by viral clusters (VCs)
- Reference-based classification

**Tool:** vContact3

**Input:** Phage sequences

**Approach:**
1. Predict proteins (Prodigal)
2. Create protein clusters (similarity grouping)
3. Build network based on shared protein clusters
4. Identify viral clusters (communities)
5. Assign taxonomy based on reference genomes in clusters

**Command:**
```bash
vcontact3 run \
    --nucleotide phage_seqs.fasta \
    --db-domain 'prokaryotes' \
    --db-version 223 \
    --output-dir 03_genomic_info/vc3_output
```

**Parameters:**
- `--db-domain prokaryotes` - Use bacterial virus database
- `--db-version 223` - Database version
- `--nucleotide` - Input is nucleotide sequences

**Output Directory:** `03_genomic_info/vc3_output/`

**Key Output Files:**
- `viral_cluster_overview.tsv` - Cluster assignments
- `genome_by_genome_overview.csv` - Per-genome results
- Network files for visualization

**Viral Cluster (VC) Assignment:**
- Phages in same VC share significant gene content
- VCs often correspond to genera or species
- Provides complementary taxonomic signal to sequence similarity

**Processing Time:** 2-4 hours for 1,500 sequences

---

### Tool 5: **`bacphlip_lifestyle`** - Lifestyle Classification

**What it does:** Predicts phage lifestyle (temperate vs lytic) using protein features

**Tool:** BACPHLIP (Random Forest classifier)

**Input:** Phage sequences

**Method:**
1. Predicts proteins
2. Identifies lysogeny-associated proteins
3. Uses Random Forest model
4. Classifies as temperate or lytic

**Command:**
```bash
bacphlip -i phage_seqs.fasta \
    --multi_fasta \
    > bacphlip_lifestyle.tsv
```

**Output:** `03_genomic_info/bacphlip_lifestyle.tsv`

**Columns:**
- Phage ID
- Lifestyle prediction
- Confidence score

**Additional Processing:**
Script adds completeness flags from CheckV:
- Flags incomplete genomes (may affect lifestyle prediction accuracy)
- Incomplete genomes may lack lysogeny markers even if temperate

**Output:** `03_genomic_info/bacphlip_lifestyle_with_completeness.tsv`

**Added Columns:**
- `completeness` - From CheckV
- `checkv_quality` - Quality tier
- `completeness_flag` - Indicator for interpretation

**Processing Time:** 30-60 minutes for 1,500 sequences

---

## Substage 3D: Consensus Building

### Purpose
Integrate predictions from multiple tools to create more robust classifications.

---

### **`lifestyle_consensus` Rule** - Lifestyle Integration

**What it does:** Combines BACPHLIP and Phabox2 lifestyle predictions

**Tool:** Custom Python script (`scripts/lifestyle_consensus.py`)

**Inputs:**
- `03_genomic_info/bacphlip_lifestyle_with_completeness.tsv`
- `03_genomic_info/phabox_output/lifestyle.tsv`

**Consensus Logic:**

```python
def determine_consensus(bacphlip_pred, bacphlip_conf, phabox_pred, phabox_conf, completeness):
    threshold = 0.7  # Confidence threshold

    # Priority 1: Both tools agree with high confidence
    if bacphlip_pred == phabox_pred and bacphlip_conf >= threshold and phabox_conf >= threshold:
        return bacphlip_pred, "high_confidence"

    # Priority 2: One tool has high confidence, other is uncertain
    if bacphlip_conf >= threshold and phabox_conf < threshold:
        return bacphlip_pred, "bacphlip_only"
    if phabox_conf >= threshold and bacphlip_conf < threshold:
        return phabox_pred, "phabox_only"

    # Priority 3: Tools disagree
    if bacphlip_pred != phabox_pred:
        # Consider completeness
        if completeness < 50:
            return "uncertain_incomplete", "low_confidence"
        else:
            return "conflicting", "low_confidence"

    # Default: uncertain
    return "uncertain", "low_confidence"
```

**Output:** `03_genomic_info/lifestyle_consensus.tsv`

**Columns:**
- `contig_id`
- `consensus_lifestyle` - Final prediction
- `confidence_level` - High, medium, low
- `bacphlip_prediction`
- `bacphlip_confidence`
- `phabox_prediction`
- `phabox_confidence`
- `completeness`
- `agreement_flag` - Tools agree/disagree

**Interpretation:**
- **High confidence:** Both tools agree with >70% confidence
- **Medium confidence:** One tool confident, other uncertain
- **Low confidence:** Disagreement or both uncertain
- **Consider completeness:** Incomplete genomes may lack markers

**Processing Time:** <5 minutes

---

### **`taxonomic_consensus` Rule** - Taxonomy Integration (OPTIONAL)

**What it does:** Hierarchically integrates MMseqs2, Phabox2, and vContact3 taxonomy

**Tool:** Custom R script (`scripts/taxonomic_consensus.R`)

**Controlled By:** `run_consensus` config parameter (default: true)

**Inputs:**
- `03_genomic_info/mmseqs_taxonomy.tsv` (from Stage 3)
- `03_genomic_info/phabox_output/taxonomy.tsv`
- `03_genomic_info/vc3_output/` (vContact3 results)
- Optional: `01_filtered_mmseqs/filtered_lca.tsv` (from Stage 1, if available)

**Consensus Strategy:**

**Hierarchical Prioritization:**
1. **MMseqs2** - Highest priority (protein similarity to known sequences)
2. **Phabox2** - Medium priority (ML-based, phage-specific)
3. **vContact3** - Lower priority (gene content, clustering-based)

**Integration Process:**

```R
# For each taxonomic level (Kingdom → Species):
for each phage {
  # Start with MMseqs2
  if (MMseqs2 has assignment at this level AND confidence > threshold) {
    consensus = MMseqs2_assignment
  }
  # Fall back to Phabox2
  else if (Phabox2 has assignment AND confidence > threshold) {
    consensus = Phabox2_assignment
  }
  # Fall back to vContact3
  else if (vContact3 has viral cluster assignment) {
    consensus = VC_based_taxonomy
  }
  # Unclassified
  else {
    consensus = "unclassified"
  }

  # Validate hierarchical consistency
  validate_lineage(consensus)
}
```

**Hierarchical Validation:**
Ensures taxonomic consistency (e.g., if Family is assigned, Order must also be assigned)

**Output:** `03_genomic_info/consensus_taxonomy.tsv`

**Columns:**
- `contig_id`
- `kingdom`
- `phylum`
- `class`
- `order`
- `family`
- `genus`
- `species`
- `full_lineage`
- `confidence_level`
- `primary_source` - Which tool provided assignment (MMseqs2, Phabox2, or vContact3)
- `agreement_score` - How many tools agreed

**Why Optional?**
- Computationally intensive (requires taxonomizr database)
- May not be needed for all analyses
- Simple use cases can use individual tool outputs

**Processing Time:** 10-30 minutes for 1,500 sequences

---

## Substage 3E: Final Integration

### Purpose
Create a comprehensive summary table integrating all analysis results.

---

### **`create_final_contig_table` Rule** - Master Integration

**What it does:** Combines all analysis outputs into single summary table

**Tool:** Custom Python script (`scripts/create_final_contig_table.py`)

**Required Inputs:**
- `phage_seqs` - Sequence file (for sequence info)
- `consensus_taxonomy` - Integrated taxonomy
- `lifestyle_consensus` - Integrated lifestyle
- `iphop_predictions` - Host predictions

**Optional Inputs (auto-detected):**
- `phage_predictions` - Stage 1 prediction metadata (if available)
- `checkv_results` - Final CheckV quality assessment
- `mmseqs_lca` - Stage 1 viral filtering taxonomy (if available)

**Integration Logic:**

```python
# Start with sequence IDs
final_table = get_sequence_ids(phage_seqs)

# Add core annotations
final_table = final_table.merge(consensus_taxonomy, on='contig_id', how='left')
final_table = final_table.merge(lifestyle_consensus, on='contig_id', how='left')
final_table = final_table.merge(iphop_predictions, on='contig_id', how='left')

# Add optional annotations if available
if checkv_results exists:
    final_table = final_table.merge(checkv_results, on='contig_id', how='left')

if phage_predictions exists:
    final_table = final_table.merge(phage_predictions, on='contig_id', how='left')

# Add Stage 1 LCA if available
if mmseqs_lca exists:
    final_table = final_table.merge(mmseqs_lca, on='contig_id', how='left')

# Calculate summary statistics
final_table['completeness_category'] = categorize_completeness(completeness)
final_table['has_host_prediction'] = host_genus.notna()
final_table['taxonomy_level'] = get_deepest_taxonomy_level(lineage)

return final_table
```

**Output:** `final_contig_summary.tsv`

**Column Categories:**

**Basic Information:**
- `contig_id` - Sequence identifier
- `length` - Sequence length
- `gc_content` - GC percentage (if calculated)

**Taxonomy:**
- `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`
- `taxonomy_confidence`
- `taxonomy_source` (MMseqs2/Phabox2/vContact3)

**Quality:**
- `completeness` - Percentage from CheckV
- `checkv_quality` - Quality tier
- `contamination` - Estimated contamination
- `warnings` - QC warnings

**Lifestyle:**
- `lifestyle_consensus` - Temperate/Lytic/Unknown
- `lifestyle_confidence` - High/Medium/Low
- `bacphlip_lifestyle`
- `phabox_lifestyle`

**Host Prediction:**
- `host_genus` - Predicted host genus
- `host_family` - Predicted host family
- `host_prediction_confidence`
- `host_prediction_methods` - Which methods supported

**Stage 1 Metadata (if available):**
- `functional_diversity` - From PHOLD
- `prediction_rule` - Which rule classified as phage
- `genomad_topology` - DTR/ITR/Provirus/Linear
- `jaeger_score`
- `genomad_score`
- `initial_mmseqs_lineage` - From Stage 1 filtering

**Computed Fields:**
- `completeness_category` - High (>90%), Medium (50-90%), Low (<50%)
- `has_host_prediction` - Boolean
- `taxonomy_confidence_level` - Based on agreement
- `analysis_flags` - Any noteworthy features

**Example Row:**
```
contig_id: phage_00042
length: 48523
completeness: 87.3
checkv_quality: High-quality
order: Caudovirales
family: Siphoviridae
genus: unclassified
lifestyle_consensus: temperate
lifestyle_confidence: high
host_genus: Bacillus
host_prediction_confidence: 0.89
functional_diversity: 4
genomad_topology: Provirus
```

**Output Format:** Tab-separated values (TSV)

**File Size:** Typically 1-5 MB for 1,500 sequences

**Processing Time:** <5 minutes

---

## Stage 3 Outputs Summary

### Primary Output (Main Result)
**`final_contig_summary.tsv`**
- Comprehensive table with all analysis results
- One row per phage sequence
- All columns from all analyses integrated

### Quality Assessment
- `03_checkv_final/quality_summary.tsv` - Final quality metrics

### Host Prediction
- `03_iphop_results/iphop_predictions_compiled.tsv` - Host predictions

### Taxonomy (Multiple Sources)
- `03_genomic_info/mmseqs_taxonomy.tsv` - MMseqs2 taxonomy
- `03_genomic_info/phabox_output/taxonomy.tsv` - Phabox2 taxonomy
- `03_genomic_info/vc3_output/` - vContact3 results
- `03_genomic_info/consensus_taxonomy.tsv` - Integrated consensus

### Lifestyle
- `03_genomic_info/bacphlip_lifestyle_with_completeness.tsv` - BACPHLIP predictions
- `03_genomic_info/phabox_output/lifestyle.tsv` - Phabox2 lifestyle
- `03_genomic_info/lifestyle_consensus.tsv` - Integrated consensus

### Functional Annotation
- `03_orf_predictions/proteins.faa` - Predicted proteins
- `03_orf_predictions/genes.gff` - Gene coordinates

### Intermediate Files
- `03_split_seqs/*.fasta` - iPhop input chunks
- `03_iphop_results/tmp/` - Individual chunk predictions

---

## Processing Time Estimates

**For typical dataset (1,500 sequences, 24 threads):**

| Substage | Steps | Time | Bottleneck |
|----------|-------|------|------------|
| 3A: Quality | CheckV final | 30-60 min | Database search |
| 3B: Host Prediction | Split + iPhop (15 chunks) + Aggregate | 8-12 hours | iPhop array jobs |
| 3C: Characterization | All tools (parallel) | 3-5 hours | vContact3, MMseqs2 |
| 3D: Consensus | Python/R scripts | 15-30 min | I/O, data processing |
| 3E: Integration | Python script | <5 min | Fast |

**Total Stage 3:** ~12-18 hours (mostly parallel processing)

**Scaling:**
- Doubles with 3,000 sequences
- iPhop is the main bottleneck (can be sped up with more parallelization)
- Most other tools scale linearly

---

## Parallelization Summary

### What Runs in Parallel:

**Level 1: Across Substages**
- Substage 3A (CheckV) can run independently
- Substage 3B (iPhop) runs independently after splitting
- Substage 3C tools (Prodigal, MMseqs2, Phabox2, vContact3, BACPHLIP) all run in parallel

**Level 2: Within iPhop (Substage 3B)**
- 15 iPhop jobs run simultaneously (one per chunk)
- Each uses 24 threads
- Total: ~360 CPU cores at peak

**Level 3: Within Tools**
- Each tool uses 24 threads internally
- MMseqs2, CheckV, Phabox2, vContact3 all multi-threaded

**Total Parallelization:**
- Peak concurrent jobs: ~20 (5 characterization tools + 15 iPhop chunks)
- Peak CPU usage: ~500 cores (if unlimited resources)
- Typical HPC allocation: 100-200 cores (sequential chunks)

---

## Quality Control Checkpoints

### What to Check After Stage 3:

**1. CheckV Final Assessment:**
```bash
# Completeness distribution
cut -f10 03_checkv_final/quality_summary.tsv | sort | uniq -c

# Quality tier distribution
cut -f8 03_checkv_final/quality_summary.tsv | sort | uniq -c
```

Expected: Mix of quality tiers, some high-quality genomes

**2. iPhop Predictions:**
```bash
# Prediction rate
tail -n +2 03_iphop_results/iphop_predictions_compiled.tsv | wc -l
```

Expected: 60-80% of sequences have host predictions

**3. Taxonomy Coverage:**
```bash
# How many sequences got taxonomy?
awk -F'\t' '$2!="unclassified"' consensus_taxonomy.tsv | wc -l
```

Expected: >70% classified to at least Order level

**4. Lifestyle Predictions:**
```bash
# Lifestyle distribution
cut -f2 lifestyle_consensus.tsv | tail -n +2 | sort | uniq -c
```

Expected: Mix of temperate and virulent, some uncertain

**5. Final Table Completeness:**
```bash
# Check all sequences made it through
wc -l final_contig_summary.tsv

# Check for missing values in key columns
grep -c "NA" final_contig_summary.tsv
```

---

## Red Flags

⚠️ **All sequences missing host predictions**
- Check: iPhop database, logs for errors
- May indicate: Database issue, incompatible sequences

⚠️ **All sequences "unclassified" for taxonomy**
- Check: MMseqs2/Phabox2 outputs, database versions
- May indicate: Novel/divergent phages, database issues

⚠️ **All lifestyle predictions "uncertain"**
- Check: BACPHLIP/Phabox2 outputs
- May indicate: Incomplete genomes, novel phages

⚠️ **Very low completeness across all sequences**
- Check: CheckV database, input sequences
- May indicate: Fragmented assembly, short sequences

⚠️ **Final table has fewer rows than input**
- Check: Integration script logs
- Should never happen - investigate immediately

---

## Entry Points Impact on Stage 3

| Entry Point | What's Available | What's Missing |
|-------------|-----------------|----------------|
| `raw_contigs` | Everything | - |
| `reneo` | Everything | - |
| `viral_contigs` | Everything except Stage 1 LCA | Initial MMseqs2 taxonomy |
| `phage_contigs` | Stage 3 analyses only | Stage 1 predictions, CheckV from Stage 1 |
| `clustering` | Stage 3 analyses only | Earlier stages |

**Stage 1 LCA (`01_filtered_mmseqs/filtered_lca.tsv`):**
- Auto-detected and included if available
- Provides initial viral taxonomy assignment
- Shows taxonomic lineage before phage prediction
- Optional - pipeline works without it

---

## Tool Comparison Matrix

| Aspect | MMseqs2 | Phabox2 | vContact3 | BACPHLIP |
|--------|---------|---------|-----------|----------|
| **Purpose** | Taxonomy | Tax + Life | Taxonomy | Lifestyle |
| **Method** | Protein BLAST | ML (features) | Gene content | ML (proteins) |
| **Database** | NCBI NR | Phage-specific | Reference genomes | Training set |
| **Speed** | Slow | Medium | Slow | Fast |
| **Accuracy** | High (known) | High (phages) | Medium | High |
| **Novel phages** | Poor | Medium | Poor | Medium |
| **Output** | Full lineage | Order/Family | Viral cluster | Temp/Lytic |

**Complementary Strengths:**
- **MMseqs2:** Best for well-characterized phages with protein matches
- **Phabox2:** Handles novel phages better with ML features
- **vContact3:** Finds relatives through gene sharing
- **BACPHLIP + Phabox2:** Independent lifestyle predictions

---

## Files Generated (Complete List)

```
03_checkv_final/
└── quality_summary.tsv

03_split_seqs/
├── *.fasta (100 seqs/file)
└── split_file_list.txt

03_iphop_results/
├── tmp/
│   └── {chunk_id}/
│       └── Host_prediction_to_genus_m90.csv
├── .splits_ready
├── .input_files_found
├── .all_predictions_done
└── iphop_predictions_compiled.tsv

03_orf_predictions/
├── genes.gff
├── proteins.faa
└── nucleotides.fna

03_genomic_info/
├── mmseqs_taxonomy.tsv
├── phabox_output/
│   ├── taxonomy.tsv
│   └── lifestyle.tsv
├── vc3_output/
│   ├── viral_cluster_overview.tsv
│   ├── genome_by_genome_overview.csv
│   └── [network files]
├── bacphlip_lifestyle.tsv
├── bacphlip_lifestyle_with_completeness.tsv
├── lifestyle_consensus.tsv
└── consensus_taxonomy.tsv

final_contig_summary.tsv ← **MAIN OUTPUT**

logs/
├── checkv_final_assessment.log
├── split_phage_sequences.log
├── check_iphop_input_files.log
├── iphop_single_prediction/
│   └── {chunk_id}.log
├── iphop_aggregate_results.log
├── prodigal_orf_prediction.log
├── mmseqs_phage_taxonomy.log
├── phabox_prediction.log
├── vcontact3_taxonomy.log
├── bacphlip_lifestyle.log
├── lifestyle_consensus.log
├── taxonomic_consensus.log
└── create_final_contig_table.log
```

---

## Summary Diagram (Detailed)

```
                    INPUT (from Stage 2 or 1)
                 vOTUs (1.5K) or All Phages (5K)
                              ↓
                   ┌──────────────────┐
                   │ get_phage_input()│ Function selects
                   └────────┬─────────┘
                            ↓
    ┌───────────────────────────────────────────────┐
    │ SUBSTAGE 3A: Quality Assessment               │
    │  CheckV Final (30-60 min)                     │
    └────────────────┬──────────────────────────────┘
                     ↓
    ┌───────────────────────────────────────────────┐
    │ SUBSTAGE 3B: Host Prediction (8-12 hours)     │
    │                                               │
    │  Split (100 seq/chunk)                        │
    │        ↓                                       │
    │  Checkpoint (wait)                            │
    │        ↓                                       │
    │  ┌─────────────────────────────────┐          │
    │  │ iPhop Array (15 jobs parallel)  │ ∥∥∥∥∥    │
    │  │  - 24 threads each              │          │
    │  │  - 30-60 min per job            │          │
    │  └─────────────┬───────────────────┘          │
    │                ↓                               │
    │  Aggregate Results                            │
    └────────────────┬──────────────────────────────┘
                     ↓
    ┌───────────────────────────────────────────────┐
    │ SUBSTAGE 3C: Characterization (3-5 hours)     │
    │                                               │
    │  ┌──────────┬──────────┬──────────┬─────────┐│
    │  │ Prodigal │ MMseqs2  │ Phabox2  │vContact3││ Parallel
    │  │  ORFs    │Taxonomy  │Tax+Life  │Gene Cont││
    │  │ 15 min   │ 2 hours  │ 2 hours  │ 4 hours ││
    │  └────┬─────┴────┬─────┴────┬─────┴────┬────┘│
    │       └──────────┼──────────┼──────────┘     │
    │  ┌───────────────┴──────────┘                │
    │  │ BACPHLIP Lifestyle (60 min)               │
    │  └───────────────┬───────────────────────────┘
    └──────────────────┼────────────────────────────┘
                       ↓
    ┌───────────────────────────────────────────────┐
    │ SUBSTAGE 3D: Consensus (15-30 min)            │
    │                                               │
    │  ┌──────────────────────────────────────────┐ │
    │  │ Lifestyle Consensus                      │ │
    │  │  BACPHLIP + Phabox2 → consensus          │ │
    │  └──────────────┬───────────────────────────┘ │
    │                 ↓                              │
    │  ┌──────────────────────────────────────────┐ │
    │  │ Taxonomic Consensus (Optional)           │ │
    │  │  MMseqs2 + Phabox2 + vContact3 → lineage │ │
    │  └──────────────┬───────────────────────────┘ │
    └─────────────────┼────────────────────────────┘
                      ↓
    ┌───────────────────────────────────────────────┐
    │ SUBSTAGE 3E: Integration (<5 min)             │
    │                                               │
    │  Combine All Results                          │
    │        ↓                                       │
    │  final_contig_summary.tsv                     │
    │                                               │
    │  Columns:                                     │
    │   • ID, length, completeness                  │
    │   • Taxonomy (all levels)                     │
    │   • Lifestyle (consensus)                     │
    │   • Host prediction                           │
    │   • Quality metrics                           │
    │   • Functional diversity                      │
    └───────────────────────────────────────────────┘
```
