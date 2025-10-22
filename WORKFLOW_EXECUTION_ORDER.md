# Phage Analysis Pipeline - Execution Order & Dependencies

## Overview

This document describes the execution order and dependencies for the phage analysis pipeline, with particular focus on Stage 3 (characterization) where multiple analyses run in parallel.

---

## Complete Pipeline Flow

### **Stage 1: Phage Identification & Prediction**

**Input Options:**
- **Option A**: FASTA file of assembled contigs (standard workflow)
- **Option B**: GFA assembly graph + raw sequencing reads (Reneo workflow for enhanced binning)

**Step-by-step process:**

1. **Optional Assembly Refinement (Reneo)**
   - If using GFA input, Reneo processes the assembly graph with reads for improved contig binning
   - Filters contigs to ≥1KB length
   - Otherwise, directly filters FASTA assembly to ≥1KB

2. **Initial Viral Screening (MMseqs2)**
   - Performs sensitive protein-based taxonomic classification against reference databases
   - Identifies sequences classified as:
     - Viruses (taxID 10239)
     - Unclassified (taxID 0)
     - Root (taxID 1)
   - **Output**: Filtered set of "candidate viral contigs"

3. **Multi-Tool Phage Prediction** (parallel execution on candidate virals)
   - **Jaeger**: Neural network-based phage detection (score ≥2.5)
   - **GeNomad**: Comprehensive viral detection (score ≥0.6)
   - **PHOLD**: Phage protein functional annotation (runs in parallel chunks of 1000 sequences)
   - **CheckV**: Quality assessment and genome completeness estimation

4. **Prediction Integration (R script)**
   - Combines evidence from all 4 tools
   - Applies high-confidence criteria:
     - Functional diversity ≥3 genes, OR
     - Special genomic features (DTR, ITR, Provirus), OR
     - Strong multi-tool support
   - **Output**: High-confidence phage contigs (phageContigs.fasta)

---

### **Stage 2: Clustering (Optional)**

**Purpose**: Group similar phages into viral Operational Taxonomic Units (vOTUs)

**Process:**
1. **Length filtering**: Minimum 10KB (configurable)
2. **vClust three-step clustering**:
   - **Prefilter**: k-mer similarity screening (≥20 common 25-mers, 95% identity)
   - **Align**: Pairwise alignment of filtered candidates
   - **Cluster**: Leiden algorithm clustering (default: 95% identity, 85% coverage)
3. **Representative selection**: Longest sequence per cluster
4. **Output**: Representative sequences for each vOTU

---

### **Stage 3: Comprehensive Characterization**

**Input**: Either clustered representatives OR all predicted phages (depending on clustering setting)

**Quality Assessment First:**
- **CheckV Final Assessment**: Runs on the actual sequences being analyzed (clustered reps OR all phages)
  - This is the "final checkV checking" from the original workflow
  - Results match the sequences in downstream analysis
  - Used by BACPHLIP for completeness-aware lifestyle prediction

**Parallel Processing Setup:**
- Sequences split into chunks of 100 for array-based processing
- Checkpoints ensure proper workflow coordination

---

## Stage 3 Analysis Steps: Execution Order & Dependencies

### **Quality Assessment: CheckV Final**
- **Runs first** - quality check on final analyzed set
- **Input**: Output from `get_phage_input()` (clustered reps OR all phages)
- **Output**: quality_summary.tsv matching analyzed sequences
- **Purpose**: Provides completeness/quality metrics for final outputs
- **Used by**: BACPHLIP (for lifestyle predictions), Final summary table

### **Parallel Track 1: iPhop (Host Prediction)**
- **Runs independently** - only needs phage sequences
- **Parallelization**: Split into chunks → Run in parallel → Aggregate
- **No dependencies** on other Stage 3 tools
- **Input**: Phage sequences (from Stage 1 or 2)
- **Output**: iphop_predictions_compiled.tsv (host predictions at genus level)

---

### **Parallel Track 2: Prodigal (Protein Prediction)**
- **Runs independently** - only needs phage sequences
- **Input**: Phage sequences (from Stage 1 or 2)
- **Output**:
  - proteins.faa
  - genes.fna
- **No dependencies** on other Stage 3 tools

---

### **Parallel Track 3a: MMseqs2 Taxonomy**
- **Runs independently** - only needs phage sequences
- **Input**: Phage sequences (from Stage 1 or 2)
- **Output**: mmseqs_taxonomy.tsv
- **No dependencies** on other Stage 3 tools

---

### **Parallel Track 3b: Phabox2**
- **Runs independently** - only needs phage sequences
- **Input**: Phage sequences (from Stage 1 or 2)
- **Output**:
  - phabox_output/taxonomy.tsv
  - phabox_output/lifestyle.tsv
  - **Note**: Both taxonomy AND lifestyle in one run
- **No dependencies** on other Stage 3 tools

---

### **Parallel Track 3c: vContact3 Taxonomy**
- **DEPENDS ON**: Prodigal
- **Input**:
  - proteins.faa (from Prodigal) ← **dependency**
  - Phage sequences (from Stage 1 or 2)
- **Output**: vc3_output/ directory
- **Must wait for**: Prodigal to complete

---

### **Parallel Track 4: BACPHLIP Lifestyle**
- **DEPENDS ON**: CheckV Final Assessment
- **Input**:
  - Phage sequences (from `get_phage_input()`)
  - quality_summary.tsv (from CheckV Final Assessment) ← **dependency**
- **Output**:
  - bacphlip_lifestyle.tsv
  - bacphlip_lifestyle_with_completeness.tsv
- **Must wait for**: CheckV Final Assessment to complete
- **Note**: Uses CheckV results that match the input sequences exactly

---

### **Sequential Step: Taxonomic Consensus**
- **DEPENDS ON**: MMseqs2, Phabox2, vContact3
- **Input**:
  - mmseqs_taxonomy.tsv (from MMseqs2) ← **dependency**
  - phabox_output/taxonomy.tsv (from Phabox2) ← **dependency**
  - phabox_output/lifestyle.tsv (from Phabox2) ← **dependency**
  - vc3_output/ (from vContact3) ← **dependency**
- **Output**:
  - consensus_taxonomy.tsv
  - consensus_taxonomy_summary.json
- **Must wait for**: All three taxonomy tools to complete
- **Process**:
  - Hierarchical integration using R script with taxonomizr
  - Prioritization: MMseqs2 > Phabox2 > vContact3
  - Validates taxonomic consistency across levels
  - Handles multiple input formats automatically

---

## Visual Execution Flow

```
START (with phage sequences from Stage 1 or 2)
    │
    ├─────────────────────────────────────────┐
    │                                         │
    ▼                                         ▼
┌─────────────────┐                    ┌──────────────┐
│ iPhop           │ (INDEPENDENT)      │ Prodigal     │ (INDEPENDENT)
│ (Host Pred)     │                    │ (Proteins)   │
└─────────────────┘                    └───────┬──────┘
                                               │
    ┌──────────────────────────────────────────┤
    │                                          │
    ▼                                          ▼
┌─────────────────┐                    ┌──────────────┐
│ MMseqs2         │ (INDEPENDENT)      │ vContact3    │ (DEPENDS: Prodigal)
│ (Taxonomy)      │                    │ (Taxonomy)   │
└────────┬────────┘                    └───────┬──────┘
         │                                     │
         │         ┌─────────────────┐         │
         │         │ Phabox2         │ (INDEPENDENT)
         │         │ (Tax + Lifestyle)│
         │         └────────┬────────┘
         │                  │
         └──────────────────┴─────────────────┘
                            │
                            ▼
                   ┌─────────────────────┐
                   │ Taxonomic Consensus │ (DEPENDS: MMseqs2, Phabox2, vContact3)
                   │ (R script)          │
                   └─────────────────────┘

    ┌──────────────────────────────┐
    │ BACPHLIP                     │ (INDEPENDENT - uses CheckV from Stage 1)
    │ (Lifestyle + Completeness)   │
    └──────────────────────────────┘
```

---

## Execution Summary

### **Run in PARALLEL (truly independent):**
1. **iPhop** - Host prediction
2. **Prodigal** - Protein prediction
3. **MMseqs2** - Taxonomy
4. **Phabox2** - Taxonomy + lifestyle (both in single run)
5. **BACPHLIP** - Lifestyle prediction

### **Run AFTER Prodigal completes:**
6. **vContact3** - Taxonomy (needs proteins from Prodigal)

### **Run AFTER all taxonomy tools complete:**
7. **Taxonomic Consensus** - Integration (needs MMseqs2, Phabox2, vContact3 outputs)

---

## Key Points

- **5 tools can start immediately in parallel** (iPhop, Prodigal, MMseqs2, Phabox2, BACPHLIP)
- **vContact3 waits for Prodigal** (needs protein sequences)
- **Taxonomic Consensus is the final step** (integrates all taxonomy results)
- **Phabox2 provides BOTH taxonomy and lifestyle** in a single run
- **BACPHLIP runs independently** but enriches results with CheckV completeness data from Stage 1
- **iPhop uses checkpoint-based parallelization** (splits sequences → parallel jobs → aggregate)

This design maximizes parallelization while respecting data dependencies!

---

## Stage 4: Summary & Reporting

**Automated Report Generation:**
- Collects step-by-step summaries into JSON files
- Generates interactive HTML report (Pipeline_Summary_Report.html) with:
  - Configuration overview
  - Overall statistics (sequences in/out at each step)
  - Progress tracking with visual indicators
  - Detailed results for each analysis step
  - Quality metrics and tool outputs

---

## Final Outputs

- **Phage sequences**: High-confidence phage genomes (optionally clustered)
- **Host predictions**: Bacterial hosts at genus/genome level
- **Lifestyle**: Temperate vs. virulent classification (with completeness)
- **Taxonomy**: Consensus classification from multiple tools
- **Functional annotation**: Predicted proteins and functions
- **Quality metrics**: Completeness and contamination estimates
- **HTML report**: Interactive summary of entire pipeline

---

## Technical Features

### **Robust Error Handling**
- Empty file handling with proper headers
- Graceful failures (creates placeholders to continue workflow)
- Comprehensive logging at each step

### **Efficient Parallelization**
- Smart chunking for computationally intensive steps (PHOLD, iPhop)
- Checkpoint system for dynamic file generation
- Aggregation with header management

### **Flexible Configuration**
- Optional steps (clustering, consensus)
- Configurable parameters (identity thresholds, chunk sizes)
- SLURM integration for HPC environments

### **Quality Control**
- CheckV completeness assessment integrated throughout
- Multi-tool validation for predictions
- Completeness flags on lifestyle predictions
