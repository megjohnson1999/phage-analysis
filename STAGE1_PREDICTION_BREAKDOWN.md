# Stage 1: Phage Prediction - Detailed Breakdown

## Overview
Stage 1 identifies phage sequences from metagenomic assemblies through a multi-step filtering and prediction process. It uses multiple complementary tools to identify high-confidence phage contigs from viral sequences.

---

## Workflow Summary

```
Raw Assembly → Assembly Processing → Viral Filtering → Phage Classification → Integration → Phage Contigs
```

---

## Substage 1A: Assembly Processing

### Purpose
Process the input assembly and filter to reasonable contig sizes for downstream analysis.

### Decision Point: Input Type

**TWO PATHS:**

#### PATH 1: Reneo Workflow (Graph-based)
**When:** User provides `assembly_graph` (GFA file) + `reads_dir`

**Input:**
- `assembly_graph` - Assembly graph in GFA format (e.g., from metaSPAdes, MEGAHIT)
- `reads_dir` - Directory of sequencing reads (FASTQ files)

**Steps:**

1. **`reneo_binning` rule**
   - **Tool:** Reneo (requires Gurobi license)
   - **What it does:**
     - Uses assembly graph structure to resolve ambiguous paths
     - Bins sequences based on coverage and connectivity
     - Identifies virus-like edges in the graph
     - Produces enhanced assembly with better-resolved contigs
   - **Parameters:**
     - `--minlength 1000` (filter ≥1KB during Reneo)
     - `--threads 24`
   - **Output:** `01_reneo_output/genomes_and_unresolved_edges.fasta`
     - Contains: Resolved genomic paths + unresolved virus-like edges
   - **Known issue:** May fail at koverage_genomes step but still produce usable output
     - Wrapper script (`run_reneo_wrapper.sh`) handles this gracefully

2. **`contig_length_filter` rule**
   - **Tool:** seqkit
   - **What it does:** Additional filtering to ensure ≥1KB
   - **Command:** `seqkit seq --min-len 1000`
   - **Output:** `01_reneo_output/genomes_and_unresolved_edges_1KB.fasta`

**Alternative Entry:** If `start_from=reneo`, copies from `input_reneo_dir` instead of running Reneo

---

#### PATH 2: Direct FASTA Workflow
**When:** User provides `assembly_file` (FASTA)

**Input:**
- `assembly_file` - Assembled contigs in FASTA format

**Steps:**

1. **`direct_contig_filter` rule**
   - **Tool:** seqkit
   - **What it does:** Filter to ≥1KB contigs
   - **Command:** `seqkit seq --min-len 1000`
   - **Output:** `01_filtered_assembly/filtered_assembly_1KB.fasta`

---

### Substage 1A Output
Either:
- `01_reneo_output/genomes_and_unresolved_edges_1KB.fasta` (Reneo path)
- `01_filtered_assembly/filtered_assembly_1KB.fasta` (Direct path)

**Typical result:** Thousands to hundreds of thousands of contigs ≥1KB

---

## Substage 1B: Viral Contig Identification

### Purpose
Narrow down the assembly to only viral sequences using taxonomic classification.

### Steps:

#### 1. **`mmseqs_taxonomy` rule** - Taxonomic Assignment

**Tool:** MMseqs2 easy-taxonomy

**What it does:**
- Performs sensitive protein similarity search against NCBI NR database
- Assigns taxonomic lineage to each contig using LCA (Lowest Common Ancestor) algorithm
- Identifies which contigs have viral taxonomy

**Input:** Filtered assembly (≥1KB)

**Parameters:**
- Database: MMseqs2 NR (NCBI non-redundant protein database)
- `--min-length 30` (minimum alignment length)
- `-e 1e-15` (very stringent E-value)
- `--search-type 2` (translated search - DNA vs protein)
- `-s 4.0` (sensitivity)
- `--lca-mode 2` (2bLCA algorithm)
- `--tax-lineage 2` (full taxonomic lineage)
- `--threads 24`

**Output Columns:**
```
query, target, evalue, pident, fident, nident, mismatch,
qcov, tcov, qstart, qend, qlen, tstart, tend, tlen, alnlen,
bits, qheader, theader, taxid, taxname, taxlineage
```

**Output File:** `01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv`

**Typical result:** All contigs with taxonomic assignments (viral, bacterial, archaeal, eukaryotic, unclassified)

---

#### 2. **`filter_mmseqs_lca` rule** - Extract Viral Lineages

**Tool:** Custom Python script (`scripts/01_filterMmseqsLca.py`)

**What it does:**
- Parses MMseqs2 LCA output
- Filters for contigs with "Viruses" in taxonomic lineage
- Identifies contigs with no hits (potential novel viruses)
- Reports passing, filtered, and missing contigs

**Input:**
- `01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv`
- Filtered assembly FASTA

**Filtering Logic:**
- **PASS:** Lineage contains "Viruses"
- **PASS:** No taxonomic assignment (could be novel viral sequences)
- **FAIL:** Non-viral taxonomy (bacteria, archaea, eukaryotes)

**Outputs:**
- `01_filtered_mmseqs/filtered_lca.tsv` - Viral-only taxonomy table
- `01_filtered_mmseqs/passing_contig_ids.txt` - List of viral contig IDs
- `01_filtered_mmseqs/missing_contig_ids.txt` - Contigs with no hits

**Typical result:** ~10-30% of input contigs identified as viral (varies by sample)

---

#### 3. **`extract_viral_contigs` rule** - Extract Sequences

**Tool:** seqkit

**What it does:** Extracts viral sequences using ID list

**Command:** `seqkit grep -f passing_contig_ids.txt`

**Output:** `01_filtered_mmseqs/passing_Viralcontigs.fasta`

**Alternative Entry:** If `start_from=viral_contigs`, copies from `input_viral_contigs` instead

---

### Substage 1B Output
`01_filtered_mmseqs/passing_Viralcontigs.fasta` - All viral contigs (includes phages, other viruses, prophages)

---

## Substage 1C: Phage Classification (PARALLEL PROCESSING)

### Purpose
Use multiple complementary tools to predict which viral contigs are specifically phages and assess their quality.

### Why Multiple Tools?
Each tool has different strengths:
- **Jaeger:** Phage-specific, good sensitivity
- **geNomad:** Detects DTR/ITR topology, provirus detection
- **PHOLD:** Identifies functional phage proteins
- **CheckV:** Quality assessment, completeness estimation

### Parallel Rules (run simultaneously):

---

#### Tool 1: **`jaeger_prediction`** - Phage-Specific Prediction

**What it does:**
- Machine learning-based phage classifier
- Trained specifically on phage genomes
- Provides reliability score for predictions

**Input:** `01_filtered_mmseqs/passing_Viralcontigs.fasta`

**Parameters:**
- `-s 2.5` (score threshold for phage classification)
- `--fsize 1000` (fragment size)
- `--stride 1000` (sliding window stride)

**Output:** `01_jaeger_output/passing_Viralcontigs_default_jaeger.tsv`

**Output Columns:**
- `contig_id` - Sequence identifier
- `prediction` - "Phage" or "Non-phage"
- `realiability_score` - Confidence score (0-1)

**What gets used later:**
- Contigs with `prediction == "Phage"` AND `realiability_score >= 0.6`

---

#### Tool 2: **`genomad_prediction`** - Viral Structure Detection

**What it does:**
- Identifies viral genomes using gene content + genomic features
- Detects special genomic structures:
  - **DTR** (Direct Terminal Repeats) - circular phages
  - **ITR** (Inverted Terminal Repeats) - linear phages
  - **Provirus** - integrated prophages
- Provides virus score and false discovery rate

**Input:** `01_filtered_mmseqs/passing_Viralcontigs.fasta`

**Parameters:**
- `--min-score 0.6` (minimum virus score)
- `--cleanup` (remove intermediate files)
- `--threads 24`

**Command:** `genomad end-to-end`

**Output:** `01_genomad_output/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv`

**Key Columns Used Later:**
- `contig_id`
- `topology` - DTR, ITR, Provirus, or Linear
- `virus_score` - Confidence score
- `n_hallmarks` - Number of viral hallmark genes

**Special significance:** DTR/ITR/Provirus topologies are strong phage indicators

---

#### Tool 3: **`phold_single_prediction`** - Functional Annotation (PARALLELIZED)

**What it does:**
- Predicts protein functions using structural homology
- Identifies specific phage functional categories
- Uses PHROGs (Prokaryotic Virus Remote Homologous Groups) database

**Parallelization Strategy:**
1. **`split_viral_contigs_for_phold`** - Split into chunks of 1000 sequences
2. **Checkpoint:** Wait for splits to complete
3. **`phold_single_prediction`** - Run PHOLD on each chunk (parallel)
4. **`phold_aggregate_results`** - Merge all results

**Why parallelize?** PHOLD is computationally intensive; splitting speeds up processing

**Input per chunk:** Individual FASTA files from split

**Parameters:**
- `--cpu` (use CPU instead of GPU)
- `--force` (overwrite existing)
- `-t 24` (threads per job)

**Output per chunk:** `01_phold_output/tmp/{sample}/phold_per_cds_predictions.tsv`

**Aggregated Output:** `01_phold_output/phold_per_cds_predictions.tsv`

**Output Columns:**
- `contig_id` - Parent contig
- `orf_id` - Specific protein/ORF
- `category` - Functional category (see below)
- `product` - Protein function description
- `evalue`, `identity` - Match statistics

**Functional Categories (used for phage classification):**
- **Connector** - Structural proteins connecting head and tail
- **Head & Packaging** - Capsid and DNA packaging proteins
- **Integration & Excision** - Prophage integration machinery
- **Lysis** - Host cell lysis proteins
- **Tail** - Tail structural proteins
- Other: DNA/RNA metabolism, transcription regulation, etc.

**What gets counted:** Number of categories per contig = "Functional Diversity"

---

#### Tool 4: **`checkv_assessment`** - Quality Assessment

**What it does:**
- Estimates completeness of viral genomes
- Detects contamination (host sequences)
- Identifies proviral regions
- Provides quality tier classification

**Input:** `01_filtered_mmseqs/passing_Viralcontigs.fasta`

**Command:** `checkv end_to_end`

**Parameters:**
- Database: CheckV database v1.5
- `-t 24` (threads)

**Output:** `01_checkv_output/quality_summary.tsv`

**Key Columns Used Later:**
- `contig_id`
- `contig_length`
- `completeness` - Estimated completeness (%, or "Not determined")
- `checkv_quality` - Quality tier (Complete, High-quality, Medium-quality, Low-quality, Not-determined)

**Quality tiers:**
- **Complete:** >90% complete, circularized or has DTRs
- **High-quality:** >90% complete
- **Medium-quality:** 50-90% complete
- **Low-quality:** <50% complete

---

## Substage 1D: Phage Integration & Final Selection

### Purpose
Integrate predictions from all tools using evidence-based rules to identify high-confidence phage contigs.

---

### **`integrate_phage_predictions` rule** - Multi-Tool Integration

**Tool:** Custom R script (`scripts/01_phagePrediction.R`)

**Inputs (4 files):**
1. PHOLD predictions (functional annotations)
2. Jaeger predictions (phage classifier)
3. geNomad predictions (structure/topology)
4. CheckV results (quality/completeness)

**Processing Steps:**

#### Step 1: Calculate Functional Diversity from PHOLD

**Categories counted as evidence:**
- Connector
- Head & Packaging
- Integration & Excision
- Lysis
- Tail

**Filtered OUT (not counted):**
- "hypothetical protein"
- "other"
- "unknown function"
- "DNA, RNA and nucleotide metabolism"
- "moron, auxiliary metabolic gene and host takeover"
- "transcription regulation"

**Metric:** `Functional_diversity` = Count of categories present (0-5)

**Interpretation:**
- **3+ categories** → High confidence structural phage
- **1-2 categories** → Needs support from other tools
- **0 categories** → Needs strong support from both Jaeger + geNomad

---

#### Step 2: Apply Classification Rules

**Contigs are classified as PHAGE if ANY rule is met:**

| Rule # | Condition | Interpretation |
|--------|-----------|----------------|
| **Rule 1** | `Functional_diversity >= 3` | Strong functional evidence of complete phage structure |
| **Rule 2** | `Functional_diversity < 3` AND `topology in [DTR, ITR, Provirus]` | Special genomic structures indicate phage |
| **Rule 3** | `1 <= Functional_diversity < 3` AND (`Jaeger phage` OR `geNomad positive`) | Moderate functional evidence + tool support |
| **Rule 4** | `Functional_diversity < 1` AND `Jaeger phage` AND `geNomad positive` | Both specialized tools agree despite low functional annotation |

**Additional filters in classification:**
- Jaeger must have `realiability_score >= 0.6`

**Output labels:** Each contig gets a `prediction_rule` column indicating which rule classified it

---

#### Step 3: Create Final Output

**Outputs:**
1. **`01_phage_predictions/phagePredictedContigs.tsv`** - Full table with:
   - `contig_id`
   - `contig_length` (from CheckV)
   - `completeness` (from CheckV)
   - `Functional_diversity` (calculated from PHOLD)
   - `realiability_score` (from Jaeger, if available)
   - `genomad_topology` (from geNomad)
   - `genomad_virus_score` (from geNomad)
   - `prediction_rule` - Which rule classified this as phage

2. **`01_phage_predictions/contig_ids.txt`** - Simple list of phage contig IDs

**Typical results:**
- Input: Thousands of viral contigs
- Output: Hundreds to thousands of phage contigs (depends on sample)
- Reduction: ~50-80% of viral contigs classified as phages

---

### **`extract_phage_contigs` rule** - Extract Final Sequences

**Tool:** seqkit

**What it does:** Extract phage sequences using integrated ID list

**Input:**
- `01_filtered_mmseqs/passing_Viralcontigs.fasta` (all viral contigs)
- `01_phage_predictions/contig_ids.txt` (phage IDs only)

**Command:** `seqkit grep -f contig_ids.txt`

**Output:** `01_phage_predictions/phageContigs.fasta`

**Alternative Entry:** If `start_from=phage_contigs`, this would normally be skipped

---

## Stage 1 Final Outputs

### Primary Outputs (Used by Stage 2+)
1. **`01_phage_predictions/phageContigs.fasta`**
   - FASTA file of all predicted phage sequences
   - Used as input for clustering (Stage 2) or direct analysis (Stage 3)

2. **`01_phage_predictions/phagePredictedContigs.tsv`**
   - Metadata table with prediction details
   - Used in final summary reports

### Intermediate Outputs (May be used in summaries)
- `01_filtered_mmseqs/filtered_lca.tsv` - Viral taxonomy
- `01_checkv_output/quality_summary.tsv` - Quality metrics
- `01_phold_output/phold_per_cds_predictions.tsv` - Functional annotations

---

## Key Metrics & Statistics

### Typical Pipeline Progression (example numbers):

```
Input Assembly: 100,000 contigs
    ↓ [Filter ≥1KB]
Filtered Assembly: 50,000 contigs (50%)
    ↓ [MMseqs2 Viral Filter]
Viral Contigs: 10,000 contigs (20% of filtered)
    ↓ [Multi-tool Phage Prediction]
Phage Contigs: 5,000 contigs (50% of viral, 10% of filtered)
```

### Tool Agreement Statistics:

**High confidence (multiple rules):**
- ~30-40% pass on Functional_diversity alone (Rule 1)
- ~20-30% supported by topology (Rule 2)
- ~20-30% supported by multiple tools (Rules 3-4)

**Tool-specific contributions:**
- PHOLD: Provides functional evidence, most informative
- geNomad: Critical for detecting prophages and structural features
- Jaeger: Adds independent ML-based classification
- CheckV: Quality control, not used in classification but critical for QC

---

## Processing Time Estimates

**For a typical metagenomic dataset (50K contigs, 24 threads):**

| Step | Time | Bottleneck |
|------|------|------------|
| Reneo binning | 2-6 hours | Graph resolution, Gurobi optimization |
| MMseqs2 taxonomy | 1-3 hours | Database size, search sensitivity |
| Jaeger | 10-30 min | Relatively fast |
| geNomad | 30-60 min | Gene calling + classification |
| PHOLD (parallelized) | 1-2 hours | Protein structure search |
| CheckV | 30-60 min | Database search |
| Integration (R) | <5 min | Fast data processing |

**Total Stage 1:** ~4-12 hours (depending on dataset size and whether Reneo is used)

---

## Quality Control Checkpoints

### What to check in logs:

1. **After Assembly Processing:**
   - Number of contigs ≥1KB
   - Total assembly size

2. **After Viral Filtering:**
   - Percentage of viral contigs (should be reasonable for your sample type)
   - Check for known viral families in taxonomy

3. **After Tool Predictions:**
   - PHOLD: Check for empty results (can happen with novel phages)
   - Jaeger: Reliability score distribution
   - geNomad: Topology distribution (expect mostly "Linear")
   - CheckV: Completeness distribution

4. **After Integration:**
   - Number of phages per rule
   - Distribution of functional diversity scores
   - Completeness of final phage set

### Red Flags:
- ⚠️ Zero viral contigs → Check MMseqs2 database, taxonomy filtering
- ⚠️ Zero phage contigs → Check tool outputs, may need to adjust thresholds
- ⚠️ PHOLD entirely empty → Database issue or novel phages
- ⚠️ >90% of viral = phage → May be too permissive, check rules

---

## Entry Points Summary

| Entry Point | Skips | Required Input | Use Case |
|-------------|-------|----------------|----------|
| `raw_contigs` | None | `assembly_file` OR `assembly_graph`+`reads_dir` | Full pipeline from assembly |
| `reneo` | Reneo execution | `input_reneo_dir` | Reuse existing Reneo results |
| `viral_contigs` | Assembly processing, MMseqs2, viral filtering | `input_viral_contigs` | Starting from pre-identified viral sequences |

---

## Decision Tree Summary

```
START
  │
  ├─→ [Assembly Type?]
  │   ├─→ Graph (.gfa) → Reneo → Filter 1KB → VIRAL FILTERING
  │   └─→ FASTA → Filter 1KB → VIRAL FILTERING
  │
VIRAL FILTERING
  │
  ├─→ MMseqs2 Taxonomy → Filter Viral Lineages → Extract Viral Contigs
  │
PARALLEL PHAGE CLASSIFICATION
  │
  ├─→ Jaeger (phage prediction)
  ├─→ geNomad (topology + score)
  ├─→ PHOLD (functional annotation)
  └─→ CheckV (quality)
  │
INTEGRATION
  │
  └─→ Apply Rules:
      Rule 1: Functional_diversity ≥ 3 → PHAGE
      Rule 2: Special topology → PHAGE
      Rule 3: Moderate function + tool support → PHAGE
      Rule 4: Both tools agree → PHAGE
  │
OUTPUT: phageContigs.fasta
```

---

## Files Generated (Complete List)

```
01_reneo_output/
├── genomes_and_unresolved_edges.fasta
└── genomes_and_unresolved_edges_1KB.fasta

01_filtered_assembly/
└── filtered_assembly_1KB.fasta

01_mmseqs_output/
└── genomes_and_unresolved_edges_mmseqs_lca.tsv

01_filtered_mmseqs/
├── filtered_lca.tsv
├── passing_contig_ids.txt
├── missing_contig_ids.txt
└── passing_Viralcontigs.fasta

01_jaeger_output/
└── passing_Viralcontigs_default_jaeger.tsv

01_genomad_output/
└── passing_Viralcontigs_summary/
    └── passing_Viralcontigs_virus_summary.tsv

01_phold_split_seqs/
├── *.fasta (1000 seqs/file)
└── split_file_list.txt

01_phold_output/
├── tmp/
│   └── {sample}/
│       └── phold_per_cds_predictions.tsv
└── phold_per_cds_predictions.tsv (aggregated)

01_checkv_output/
└── quality_summary.tsv

01_phage_predictions/
├── phagePredictedContigs.tsv
├── contig_ids.txt
└── phageContigs.fasta ← **MAIN OUTPUT FOR STAGE 2**

logs/
├── reneo_binning.log
├── contig_length_filter.log
├── direct_contig_filter.log
├── mmseqs_taxonomy.log
├── filter_mmseqs_lca.log
├── extract_viral_contigs.log
├── jaeger_prediction.log
├── genomad_prediction.log
├── split_viral_contigs_for_phold.log
├── phold_prediction/
│   └── {sample}.log (per chunk)
├── phold_aggregate_results.log
├── checkv_assessment.log
├── integrate_phage_predictions.log
└── extract_phage_contigs.log
```
