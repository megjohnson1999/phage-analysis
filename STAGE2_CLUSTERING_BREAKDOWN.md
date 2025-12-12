# Stage 2: Clustering - Detailed Breakdown

## Overview
Stage 2 groups similar phage sequences into viral Operational Taxonomic Units (vOTUs) to reduce redundancy and identify distinct phage populations. This stage is **optional** and controlled by the `do_clustering` configuration parameter.

---

## Purpose & Rationale

### Why Cluster Phages?

**Problem:** Metagenomic assemblies often contain:
- Multiple near-identical phages from the same population
- Variants of the same phage with minor mutations
- Fragmented assemblies of the same phage
- Strain-level diversity within a phage species

**Solution:** Cluster similar sequences into vOTUs (viral OTUs)

**Benefits:**
1. **Reduces redundancy** - Analyze one representative per phage population
2. **Improves computational efficiency** - Fewer sequences for downstream analysis
3. **Biological relevance** - vOTUs approximate phage species/populations
4. **Statistical power** - Better for abundance calculations and comparative analyses

**Trade-offs:**
- Loses strain-level diversity information
- May miss rare variants
- Depends on arbitrary similarity thresholds

---

## Decision Point: To Cluster or Not?

### `do_clustering` Configuration

```yaml
# In config.yaml or command line
do_clustering: true   # Default: Run clustering
do_clustering: false  # Skip clustering, analyze all phages
```

**Automatically set to FALSE when:**
- `start_from=clustering` (already have clustered sequences)

---

## Workflow Summary

```
Phage Contigs → Prefilter → Align → Cluster → Select Representatives → vOTU Sequences
   (Input)        (K-mer)    (ANI)  (Leiden)   (Longest)            (Output)
```

---

## Input Selection Logic

### Helper Function: `get_phage_contigs_input()`

Determines which file to use based on `start_from`:

| `start_from` | Input Source | File |
|--------------|--------------|------|
| `raw_contigs` or `reneo` or `viral_contigs` | Pipeline Stage 1 | `01_phage_predictions/phageContigs.fasta` |
| `phage_contigs` | User-provided | `input_phage_contigs` (from config) |
| `clustering` | Not applicable | Clustering is skipped |

---

## Detailed Steps

### Step 1: `cluster_phages` Rule - Multi-Stage Clustering

**Tool:** vClust (viral clustering tool)

**Algorithm:** Three-stage filtering and clustering

**Input:** Phage sequences from Stage 1 or user-provided

**Parameters (from config.yaml):**
```yaml
params:
  vclust:
    min_length: 10000      # Minimum sequence length (10 KB)
    identity: 0.95         # ANI threshold (95%)
    coverage: 0.85         # Coverage threshold (85%)
```

**Resources:**
- Threads: 24
- Memory: Scales with dataset size

---

#### Stage 1.1: Prefiltering (K-mer Based)

**Command:**
```bash
vclust prefilter -i phage_contigs.fasta \
    -o vclust_fltr.txt \
    --min-kmers 20 \
    --min-ident 0.95
```

**What it does:**
- Uses k-mer (25-mers) similarity for fast initial filtering
- Identifies pairs of sequences that might be similar
- Filters for pairs with ≥20 common k-mers and ≥95% k-mer identity

**Why this step?**
- **Speed:** K-mer comparison is much faster than alignment
- **Reduces search space:** Only align promising pairs
- **Scales well:** Can handle hundreds of thousands of sequences

**Output:** `vclust_fltr.txt`
- List of sequence pairs to align
- Format: `seq1 seq2` (one pair per line)

**Error Handling:**
- If prefilter produces no output (no similar sequences found):
  - Creates single-sequence clusters (each sequence is its own vOTU)
  - All sequences become representatives
  - Continues pipeline gracefully

**Typical result:** Filters out ~90-99% of possible pairwise comparisons

---

#### Stage 1.2: Alignment (ANI Calculation)

**Command:**
```bash
vclust align -i phage_contigs.fasta \
    -o vclust_ani.tsv \
    --filter vclust_fltr.txt
```

**What it does:**
- Aligns only the pre-filtered sequence pairs
- Calculates Average Nucleotide Identity (ANI) for each pair
- Calculates query coverage (qcov) and reference coverage (rcov)

**Output:** `vclust_ani.tsv`

**Columns:**
- `seq1` - First sequence ID
- `seq2` - Second sequence ID
- `ani` - Average Nucleotide Identity (0-1)
- `qcov` - Query coverage (fraction of seq1 aligned)
- `rcov` - Reference coverage (fraction of seq2 aligned)
- `aligned_length` - Number of aligned bases

**Example:**
```
seq1        seq2        ani    qcov   rcov   aligned_length
phage_001   phage_042   0.98   0.92   0.89   45231
phage_001   phage_153   0.96   0.87   0.91   43982
```

**Error Handling:**
- If alignment produces no output:
  - Creates single-sequence clusters
  - Continues pipeline

**Typical result:**
- Input: Thousands to millions of pairs (from prefilter)
- Output: Only pairs meeting similarity thresholds

---

#### Stage 1.3: Graph-Based Clustering

**Command:**
```bash
vclust cluster -i vclust_ani.tsv \
    -o vclust_clusters.tsv \
    --algorithm leiden \
    --metric ani \
    --ani 0.95 \
    --qcov 0.85 \
    --rcov 0.85
```

**What it does:**
- Builds a similarity graph from ANI results
- Applies Leiden community detection algorithm
- Groups sequences into clusters (vOTUs)

**Algorithm Details:**

**Leiden Algorithm:**
- Graph-based community detection
- Optimizes modularity to find natural groupings
- Better than traditional hierarchical clustering for large datasets
- Deterministic (same input → same output)

**Similarity Thresholds:**
- **ANI ≥ 95%** - Average nucleotide identity
- **Coverage ≥ 85%** - Both query and reference coverage

**How it works:**
1. Create graph: Nodes = sequences, Edges = ANI connections
2. Add edge only if: ANI ≥ 95% AND qcov ≥ 85% AND rcov ≥ 85%
3. Run Leiden algorithm to partition graph into communities
4. Each community = one vOTU

**Output:** `vclust_clusters.tsv`

**Format:**
```
sequence_id     cluster_id
phage_001       1
phage_042       1
phage_153       1
phage_002       2
phage_087       2
```

**Interpretation:**
- Each unique `cluster_id` is one vOTU
- Sequences with same `cluster_id` are similar enough to be grouped

**Error Handling:**
- If clustering fails or produces no clusters:
  - Each sequence becomes its own cluster
  - Preserves all sequences for downstream analysis

**Typical result:**
- Input: 5,000 phage sequences
- Output: 1,000-2,000 vOTUs (60-80% reduction)
- Reduction varies by:
  - Sample diversity
  - Presence of dominant phages
  - Assembly fragmentation

---

#### Stage 1.4: Representative Selection

**Strategy:** Select the **longest sequence** from each cluster

**Why longest?**
- More complete genomes
- Better for annotation and analysis
- More information for taxonomy/host prediction
- Likely less fragmented

**Process:**

1. **Get sequence lengths:**
```bash
seqkit fx2tab --length --name phage_contigs.fasta > seq_lengths.tsv
```

Output format:
```
phage_001    52341
phage_042    48923
phage_153    51287
```

2. **Join clusters with lengths:**
```bash
# Sort both files
sort -k1,1 vclust_clusters.tsv > vclust_sorted_clusters.tsv
sort -k1,1 seq_lengths.tsv > sorted_seq_lengths.tsv

# Join on sequence ID
join -1 1 -2 1 -t $'\t' vclust_sorted_clusters.tsv sorted_seq_lengths.tsv
```

Result:
```
sequence_id    cluster_id    length
phage_001      1             52341
phage_042      1             48923
phage_153      1             51287
phage_002      2             43128
phage_087      2             39821
```

3. **Select longest per cluster (AWK script):**
```bash
awk -F'\t' '{
    if (!($2 in max_len) || max_len[$2] < $3) {
        max_len[$2] = $3;
        max_seq[$2] = $1;
    }
} END {
    for (c in max_seq)
        print max_seq[c];
}'
```

**Logic:**
- For each cluster_id ($2), track the sequence_id ($1) with maximum length ($3)
- Output only the selected sequence IDs

**Output:** `vOTU_repSeqs.tsv`
```
phage_001
phage_002
phage_127
...
```

**Alternative Strategies (not implemented, but common in other tools):**
- Highest quality (from CheckV)
- Most complete
- Highest coverage/abundance
- Centroid (most similar to all others in cluster)

---

### Step 2: `extract_votu_representatives` Rule - Sequence Extraction

**Tool:** seqkit

**What it does:** Extracts representative sequences from the original FASTA file

**Input:**
- `phage_contigs.fasta` - Original phage sequences (all of them)
- `vOTU_repSeqs.tsv` - List of representative IDs (selected subset)

**Command:**
```bash
seqkit grep -f vOTU_repSeqs.tsv phage_contigs.fasta > vOTU_repSeqs.fasta
```

**Error Handling:**
- If representative list is empty:
  - Copies all phage contigs as representatives
  - Ensures pipeline continues

**Output:** `02_clustering/vOTU_repSeqs.fasta`
- FASTA file containing only representative sequences
- **This becomes the input for Stage 3 (Analysis)**

---

## Outputs

### Primary Output (Used by Stage 3)
**`02_clustering/vOTU_repSeqs.fasta`**
- Representative sequence for each vOTU
- Used as input for all downstream analyses

### Metadata Outputs
**`02_clustering/clusters.tsv`**
- Mapping of all sequences to their cluster IDs
- Format: `sequence_id \t cluster_id`
- Used for:
  - Abundance calculations
  - Understanding cluster composition
  - Tracking which sequences were grouped

**`02_clustering/vOTU_repSeqs.tsv`**
- Simple list of representative sequence IDs
- One ID per line
- Used for extracting sequences

### Intermediate Files (in clustering_results/)
```
02_clustering/
├── clustering_results/
│   ├── vclust_fltr.txt              # Prefilter pairs
│   ├── vclust_ani.tsv               # ANI alignments
│   ├── vclust_ani.ids.tsv           # ID mappings
│   ├── vclust_clusters.tsv          # Raw clusters
│   ├── seq_lengths.tsv              # Sequence lengths
│   ├── vclust_sorted_clusters.tsv   # Sorted clusters
│   ├── sorted_seq_lengths.tsv       # Sorted lengths
│   └── joined_Clusters_length.tsv   # Joined data
├── clusters.tsv                      # Final cluster mapping
├── vOTU_repSeqs.tsv                 # Representative IDs
└── vOTU_repSeqs.fasta               # Representative sequences ← MAIN OUTPUT
```

---

## Clustering Statistics & Metrics

### Typical Clustering Results

**Example 1: Highly diverse sample**
```
Input sequences:     5,000 phages
vOTUs formed:        3,800 clusters (76% remain)
Singleton clusters:  3,200 (84% of vOTUs)
Multi-seq clusters:  600 (16% of vOTUs)
Largest cluster:     23 sequences
Reduction:           24%
```

**Example 2: Sample with dominant phages**
```
Input sequences:     5,000 phages
vOTUs formed:        1,200 clusters (24% remain)
Singleton clusters:  800 (67% of vOTUs)
Multi-seq clusters:  400 (33% of vOTUs)
Largest cluster:     187 sequences
Reduction:           76%
```

**Example 3: Low diversity / high contamination**
```
Input sequences:     2,000 phages
vOTUs formed:        1,950 clusters (97.5% remain)
Singleton clusters:  1,920 (98% of vOTUs)
Multi-seq clusters:  30 (2% of vOTUs)
Largest cluster:     5 sequences
Reduction:           2.5%
```

### Cluster Size Distribution

**Typical pattern (diverse metagenome):**
```
Cluster Size    Count    Percentage
1 (singleton)   3,200    84%
2-5             450      12%
6-10            100      3%
11-50           45       1%
50+             5        0.1%
```

### Quality Indicators

**Good clustering results:**
✅ Reasonable reduction (20-80%, depends on sample)
✅ Mix of singletons and multi-sequence clusters
✅ Cluster sizes follow power-law distribution
✅ No excessively large clusters (>10% of total sequences)

**Potential issues:**
⚠️ >95% singletons → Too stringent thresholds or very diverse sample
⚠️ <5% singletons → Too permissive thresholds or low diversity
⚠️ One huge cluster → Possible contamination or misassembly
⚠️ No clusters formed → Parameter or data issues

---

## Configuration Options

### Clustering Parameters

**In `config.yaml`:**
```yaml
params:
  vclust:
    min_length: 10000    # Minimum sequence length (bp)
    identity: 0.95       # ANI threshold (0-1)
    coverage: 0.85       # Coverage threshold (0-1)
```

**Command-line override:**
```bash
snakemake --config do_clustering=true \
    params.vclust.identity=0.90 \
    params.vclust.coverage=0.80
```

### Parameter Guidelines

| Parameter | Recommended Range | Effect of Increasing | Effect of Decreasing |
|-----------|-------------------|----------------------|----------------------|
| `identity` | 0.90 - 0.98 | More stringent, more vOTUs | More permissive, fewer vOTUs |
| `coverage` | 0.70 - 0.90 | Requires more complete matches | Allows partial matches |
| `min_length` | 5000 - 15000 | Excludes shorter sequences | Includes more sequences |

**Community standards:**
- **Species-level:** 95% ANI, 85% coverage (default)
- **Genus-level:** 70-75% ANI
- **Strain-level:** 99% ANI, 90% coverage

**Current pipeline uses species-level thresholds (95/85)**

---

## Processing Time Estimates

**For typical datasets (24 threads):**

| Input Sequences | Prefilter | Align | Cluster | Extract | Total |
|----------------|-----------|-------|---------|---------|-------|
| 1,000 phages | <1 min | 2-5 min | <1 min | <1 min | ~5 min |
| 5,000 phages | 2-3 min | 10-20 min | 1-2 min | <1 min | ~15-25 min |
| 10,000 phages | 5-10 min | 30-60 min | 2-5 min | 1-2 min | ~40-75 min |
| 50,000 phages | 30-60 min | 3-6 hours | 10-30 min | 5 min | ~4-7 hours |

**Bottlenecks:**
- Alignment step scales O(n²) with number of pairs
- Prefiltering is crucial for large datasets
- Memory usage scales with number of sequences

**Optimization:**
- vClust is already highly optimized for viral sequences
- Parallelization across 24 threads provides good speedup
- For very large datasets (>100K sequences), consider:
  - Stricter prefilter thresholds
  - Higher identity/coverage thresholds
  - Pre-filtering by length

---

## Decision Flow

### When Clustering Runs

```
START
  │
  ├─→ [do_clustering = true?]
  │   │
  │   YES ─→ [start_from?]
  │          │
  │          ├─→ raw_contigs, reneo, viral_contigs → Use Stage 1 output → CLUSTER
  │          ├─→ phage_contigs → Use input_phage_contigs → CLUSTER
  │          └─→ clustering → SKIP CLUSTERING (use input_clustered_seqs)
  │
  │   NO ─→ SKIP CLUSTERING (use all phages)
  │
OUTPUT for Stage 3:
  - If clustered: vOTU_repSeqs.fasta
  - If not clustered: phageContigs.fasta
```

### Helper Function: `get_phage_input()`

Used by Stage 3 to determine which sequences to analyze:

```python
def get_phage_input(wildcards):
    start_from = config.get("start_from", "raw_contigs")

    # If starting from clustering, use user-provided file
    if start_from == "clustering":
        return input_clustered_seqs

    # If starting from phage contigs and clustering disabled
    if start_from == "phage_contigs" and not do_clustering:
        return input_phage_contigs

    # Otherwise, follow normal logic
    if do_clustering:
        return "02_clustering/vOTU_repSeqs.fasta"
    else:
        return "01_phage_predictions/phageContigs.fasta"
```

---

## Entry Points

### How Different Entry Points Affect Clustering

| Entry Point | Clustering Behavior |
|-------------|---------------------|
| `start_from=raw_contigs` | Uses `01_phage_predictions/phageContigs.fasta` |
| `start_from=reneo` | Uses `01_phage_predictions/phageContigs.fasta` |
| `start_from=viral_contigs` | Uses `01_phage_predictions/phageContigs.fasta` |
| `start_from=phage_contigs` | Uses `input_phage_contigs` (user-provided) |
| `start_from=clustering` | **SKIPS CLUSTERING** (auto-sets `do_clustering=false`)<br>Uses `input_clustered_seqs` directly |

**User-provided inputs:**
- `input_phage_contigs`: Pre-predicted phage sequences (for `start_from=phage_contigs`)
- `input_clustered_seqs`: Pre-clustered vOTUs (for `start_from=clustering`)

---

## Biological Interpretation

### What is a vOTU?

**vOTU = viral Operational Taxonomic Unit**

Analogous to bacterial OTUs in 16S rRNA analysis:
- Represents a **phage species** or **closely related population**
- 95% ANI ≈ species-level similarity
- Accounts for:
  - Strain variation within a species
  - Assembly artifacts
  - Sequencing errors

### Ecological Significance

**vOTUs enable:**
1. **Diversity estimates** - Number of distinct phage populations
2. **Abundance calculations** - Coverage/read depth per vOTU
3. **Comparative analysis** - vOTU presence/absence across samples
4. **Network analysis** - Phage-host interactions
5. **Biogeography** - Distribution of phage populations

**Trade-offs:**
- Loses fine-scale strain variation
- May merge slightly different phages
- Arbitrary threshold (95% is convention, not biological truth)

---

## Quality Control

### What to Check

**After Clustering:**

1. **Cluster Statistics:**
   ```bash
   # Number of vOTUs
   wc -l < vOTU_repSeqs.tsv

   # Cluster size distribution
   cut -f2 clusters.tsv | sort | uniq -c | sort -nr | head -20

   # Singleton percentage
   awk '{print $2}' clusters.tsv | sort | uniq -c | awk '$1==1' | wc -l
   ```

2. **Representative Lengths:**
   ```bash
   seqkit stats vOTU_repSeqs.fasta
   ```
   - Check min/max/mean length
   - Ensure representatives are reasonable sizes

3. **Reduction Rate:**
   ```
   Reduction = (Input - Output) / Input × 100%
   ```
   - Expected: 20-80% for most samples
   - Very low (<10%): Check if clustering is working
   - Very high (>90%): Check for dominant phages or parameter issues

### Red Flags

⚠️ **All sequences become singletons**
- Possible causes: Thresholds too stringent, very diverse sample, prefilter failure
- Check: vclust_fltr.txt and vclust_ani.tsv for content

⚠️ **One massive cluster**
- Possible causes: Contamination, misassembly, permissive thresholds
- Check: Cluster composition, inspect sequences in large cluster

⚠️ **No clusters.tsv output**
- Possible causes: vClust failure, empty input
- Check: Log file for errors, input sequence count

⚠️ **Zero representatives selected**
- Should never happen due to error handling
- Check: Log file, intermediate files

---

## Comparison: Clustering ON vs OFF

### Example Workflow Paths

**With Clustering (do_clustering=true):**
```
Stage 1 → 5,000 phages
    ↓
Stage 2 → Cluster → 1,500 vOTUs (70% reduction)
    ↓
Stage 3 → Analyze 1,500 sequences
    ↓
Result: Faster analysis, population-level view
```

**Without Clustering (do_clustering=false):**
```
Stage 1 → 5,000 phages
    ↓
Stage 2 → SKIP
    ↓
Stage 3 → Analyze 5,000 sequences
    ↓
Result: Strain-level resolution, more compute time
```

### When to Skip Clustering

**Skip clustering when:**
- Sample has very low diversity (few unique phages)
- You need strain-level resolution
- You're studying phage evolution/mutations
- Dataset is small (<1,000 phages)
- You'll do clustering separately with custom parameters

**Use clustering when:**
- Sample has high diversity
- You want population/species-level analysis
- Large dataset (>5,000 phages)
- Comparative analysis across multiple samples
- Need to reduce computational cost

---

## Logs & Debugging

### Key Log File: `logs/cluster_phages.log`

**What to look for:**

```bash
# Number of input sequences
Number of input sequences: 5000

# Prefilter results
Running vclust prefilter...
[vClust output showing pairs found]

# Alignment results
Running vclust align...
[vClust output showing ANI calculations]

# Clustering results
Running vclust cluster...
[Leiden algorithm output]

# Representative selection
Selecting longest sequence per cluster as representative...
Clustering complete. Found 1500 vOTU representatives.
```

**Warning messages:**
- `WARNING: vclust prefilter produced no output` → No similar sequences
- `WARNING: vclust cluster failed` → Algorithm issues
- `WARNING: No clusters formed` → Fallback to singletons

**Error messages:**
- `ERROR: ...` → Hard failure, check input files and parameters

---

## Files Generated (Complete List)

```
02_clustering/
├── clustering_results/              # Intermediate vClust files
│   ├── vclust_fltr.txt             # Prefiltered sequence pairs
│   ├── vclust_ani.tsv              # ANI alignment results
│   ├── vclust_ani.ids.tsv          # ID mappings
│   ├── vclust_clusters.tsv         # Raw cluster assignments
│   ├── seq_lengths.tsv             # Sequence length table
│   ├── vclust_sorted_clusters.tsv  # Sorted clusters
│   ├── sorted_seq_lengths.tsv      # Sorted lengths
│   └── joined_Clusters_length.tsv  # Joined cluster-length data
│
├── clusters.tsv                     # Final cluster mapping (all seqs)
├── vOTU_repSeqs.tsv                # Representative sequence IDs
└── vOTU_repSeqs.fasta              # Representative sequences ← MAIN OUTPUT FOR STAGE 3

logs/
├── cluster_phages.log              # Clustering process log
└── extract_votu_representatives.log # Extraction log
```

---

## Algorithm Details: vClust

### Why vClust?

**vClust advantages:**
1. **Virus-specific:** Optimized for viral genomes
2. **Fast:** K-mer prefiltering + efficient alignment
3. **Accurate:** ANI-based clustering (gold standard)
4. **Scalable:** Handles 100K+ sequences
5. **Graph-based:** Leiden algorithm finds natural communities

**Alternatives (not used in this pipeline):**
- CD-HIT: General-purpose, slower for large datasets
- MMseqs2 clustering: Protein-based, different metric
- VSEARCH: Similar to CD-HIT
- Manual BLAST + clustering: Very slow

### Leiden vs Other Clustering Algorithms

**Leiden Algorithm:**
- Community detection on similarity graph
- Optimizes modularity
- Handles complex cluster structures
- Better than greedy approaches

**vs Single-linkage:**
- Leiden: Groups by overall community structure
- Single-linkage: Chains sequences, can create huge clusters

**vs Complete-linkage:**
- Leiden: More flexible boundaries
- Complete-linkage: Very stringent, many singletons

**vs UPGMA:**
- Leiden: Non-hierarchical, natural communities
- UPGMA: Hierarchical tree, requires cutting at arbitrary threshold

---

## Summary Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                    STAGE 2: CLUSTERING                      │
│                      (Optional)                             │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  INPUT: phageContigs.fasta (5,000 sequences)               │
│                          ↓                                  │
│  ┌──────────────────────────────────────────────────────┐  │
│  │ Decision: do_clustering?                             │  │
│  └────┬─────────────────────────────────────────────┬───┘  │
│       │ YES                                      NO  │      │
│       ↓                                              ↓      │
│  ┌─────────────────────────┐              Use all phages   │
│  │ 1. Prefilter (K-mers)   │                   ↓           │
│  │    - 25-mers            │              Skip to Stage 3  │
│  │    - ≥20 common         │                               │
│  │    - 95% identity       │                               │
│  └──────────┬──────────────┘                               │
│             ↓                                               │
│  ┌─────────────────────────┐                               │
│  │ 2. Align (ANI)          │                               │
│  │    - Filtered pairs     │                               │
│  │    - Calculate ANI      │                               │
│  │    - Calculate coverage │                               │
│  └──────────┬──────────────┘                               │
│             ↓                                               │
│  ┌─────────────────────────┐                               │
│  │ 3. Cluster (Leiden)     │                               │
│  │    - Build graph        │                               │
│  │    - ANI ≥ 95%          │                               │
│  │    - Cov ≥ 85%          │                               │
│  │    - Community detect   │                               │
│  └──────────┬──────────────┘                               │
│             ↓                                               │
│  ┌─────────────────────────┐                               │
│  │ 4. Select Reps          │                               │
│  │    - Longest per vOTU   │                               │
│  │    - 1,500 reps         │                               │
│  └──────────┬──────────────┘                               │
│             ↓                                               │
│  OUTPUT: vOTU_repSeqs.fasta (1,500 sequences)              │
│          70% reduction                                      │
│                                                             │
│  Metadata: clusters.tsv (all seq → cluster mapping)        │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```
