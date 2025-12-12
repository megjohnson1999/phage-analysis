# Detailed Flow Diagram Example (Tool → Output Style)

## Full Pipeline Data Flow

### Stage 1: Phage Prediction

```
┌──────────────────────┐
│ INPUT                │ 🟡 Yellow (Start)
│ assembly.gfa         │
│ + reads_dir/         │
└──────────┬───────────┘
           ↓
    ┌──────────────┐
    │    Reneo     │ 🔵 Blue (Tool)
    │   Binning    │ 24 threads
    └──────┬───────┘ ⏱️ 2-6 hours
           ↓
┌──────────────────────────────┐
│ genomes_and_unresolved_      │ ⚪ Gray (Output)
│ edges.fasta                  │
└──────────┬───────────────────┘ 50K contigs
           ↓
    ┌──────────────┐
    │   seqkit     │ 🔵 Blue (Tool)
    │  Filter 1KB  │
    └──────┬───────┘
           ↓
┌──────────────────────────────┐
│ genomes_and_unresolved_      │ ⚪ Gray (Output)
│ edges_1KB.fasta              │
└──────────┬───────────────────┘ 25K contigs (50% kept)
           ↓
    ┌──────────────┐
    │   MMseqs2    │ 🔵 Blue (Tool)
    │  Taxonomy    │ vs NR database
    │ easy-taxonomy│
    └──────┬───────┘ ⏱️ 1-3 hours
           ↓
┌──────────────────────────────┐
│ mmseqs_lca.tsv               │ ⚪ Gray (Output)
│ (All taxonomic assignments)  │
└──────────┬───────────────────┘ 25K rows
           ↓
    ┌──────────────┐
    │ Python       │ 🔵 Blue (Tool)
    │01_filterMm...│ Filter "Viruses"
    │seqsLca.py    │
    └──────┬───────┘
           ↓
┌──────────────────────────────┐
│ passing_Viralcontigs.fasta   │ ⚪ Gray (Output)
│ + passing_contig_ids.txt     │
│ + filtered_lca.tsv           │
└──────────┬───────────────────┘ 5K viral seqs (20%)
           ↓
           │
    ┌──────┴──────┬──────────┬──────────┐
    ↓             ↓          ↓          ↓
┌────────┐  ┌─────────┐ ┌────────┐ ┌────────┐
│ Jaeger │  │ geNomad │ │ PHOLD  │ │ CheckV │ 🔵 Tools (Parallel)
│        │  │         │ │(split) │ │        │
└───┬────┘  └────┬────┘ └───┬────┘ └───┬────┘ ⏱️ 1-2 hours each
    ↓            ↓          ↓          ↓
┌────────┐  ┌─────────┐ ┌────────┐ ┌────────┐
│jaeger  │  │genomad_ │ │phold_  │ │checkv_ │ ⚪ Outputs
│pred.tsv│  │summary  │ │pred.tsv│ │qual.tsv│
└───┬────┘  └────┬────┘ └───┬────┘ └───┬────┘
    └────────────┴──────────┴──────────┘
                    ↓
           ┌────────────────┐
           │   R Script     │ 🔵 Tool
           │01_phage...R    │ Integration logic
           │  Integration   │
           └────────┬───────┘
                    ↓
      ┌──────────────────────────┐
      │ phagePredictedContigs.tsv│ ⚪ Output
      │ + contig_ids.txt         │
      └──────────┬───────────────┘ 2.5K phages (50% of viral)
                 ↓
          ┌─────────────┐
          │   seqkit    │ 🔵 Tool
          │   Extract   │
          └──────┬──────┘
                 ↓
      ┌──────────────────────┐
      │ phageContigs.fasta   │ 🟢 Green (Stage output)
      └──────────┬───────────┘ STAGE 1 OUTPUT
                 ↓
```

---

### Stage 2: Clustering (Optional)

```
      ┌──────────────────────┐
      │ phageContigs.fasta   │ 🟡 Input (from Stage 1)
      └──────────┬───────────┘ 2.5K sequences
                 │
          ┌──────▼──────┐
          │ do_         │ 🟠 Orange (Decision)
          │ clustering? │ ◇
          └──┬───────┬──┘
             │       │
          YES│       │NO
             ↓       │
    ╔════════════════════╗    │
    ║ CLUSTERING         ║    │ ┄┄ Dashed border
    ║                    ║    │    (optional)
    ║  ┌──────────────┐  ║    │
    ║  │   vClust     │  ║    │
    ║  │  Prefilter   │  ║ 🔵 Tool
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │vclust_fltr   │  ║ ⚪ Output
    ║  │   .txt       │  ║    │
    ║  └──────┬───────┘  ║    │ K-mer pairs
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │   vClust     │  ║ 🔵 Tool
    ║  │    Align     │  ║    │
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │vclust_ani    │  ║ ⚪ Output
    ║  │   .tsv       │  ║    │ ANI values
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │   vClust     │  ║ 🔵 Tool
    ║  │   Cluster    │  ║    │ Leiden algorithm
    ║  │              │  ║    │ 95% ANI, 85% cov
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │vclust_       │  ║ ⚪ Output
    ║  │clusters.tsv  │  ║    │ Cluster IDs
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │   Select     │  ║ 🔵 Tool
    ║  │   Longest    │  ║    │ AWK script
    ║  └──────┬───────┘  ║    │
    ║         ↓          ║    │
    ║  ┌──────────────┐  ║    │
    ║  │vOTU_repSeqs  │  ║ ⚪ Output
    ║  │   .fasta     │  ║    │ 900 vOTUs (64% reduction)
    ║  └──────┬───────┘  ║    │
    ╚═════════╬══════════╝    │
              │               │
              └───────┬───────┘ Merge point
                      ↓
         ┌────────────────────────┐
         │ Analysis Input         │ Annotation box
         │ • IF clustered: 900    │
         │ • IF not: 2,500        │
         └────────┬───────────────┘
                  ↓
         ┌────────────────────┐
         │  Sequences for     │ 🟢 Green (Stage output)
         │  Analysis          │
         └────────────────────┘ STAGE 2 OUTPUT

```

---

### Stage 3: Analysis (Simplified - Full version would be very long)

```
         ┌────────────────────┐
         │ Analysis Sequences │ 🟡 Input (from Stage 2 or 1)
         │ (900 or 2,500)     │
         └────────┬───────────┘
                  ↓
          ┌───────────────┐
          │    CheckV     │ 🔵 Tool
          │     Final     │
          └───────┬───────┘ ⏱️ 30-60 min
                  ↓
         ┌────────────────────┐
         │ checkv_final/      │ ⚪ Output
         │ quality_summary.tsv│
         └────────┬───────────┘
                  ↓
          ┌───────────────┐
          │    seqkit     │ 🔵 Tool
          │   Split 100   │ 100 seqs/chunk
          └───────┬───────┘
                  ↓
         ┌────────────────────┐
         │ 03_split_seqs/     │ ⚪ Output
         │ *.fasta (9 files)  │
         └────────┬───────────┘
                  ↓
      ┌───────────┴───────────┐
      │                       │
      ├──→ ┌────────┐  ∥∥∥∥∥ │ Parallel array
      ├──→ │ iPhop  │  (9 jobs)
      ├──→ │  Host  │  24 threads each
      ├──→ │ Predict│  ⏱️ 8-12 hours total
      │    └───┬────┘         │
      │        ↓              │
      │    ┌────────┐         │
      │    │chunk_  │ ⚪ Outputs (per chunk)
      │    │pred.csv│         │
      │    └───┬────┘         │
      └────────┴───────────────┘
                  ↓
          ┌───────────────┐
          │   Aggregate   │ 🔵 Tool
          │    Results    │
          └───────┬───────┘
                  ↓
         ┌────────────────────────┐
         │ iphop_predictions_     │ ⚪ Output
         │ compiled.tsv           │
         └────────┬───────────────┘
                  │
      ┌───────────┼───────────┬──────────┬──────────┐
      ↓           ↓           ↓          ↓          ↓
  ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐
  │Prodigal│ │MMseqs2 │ │Phabox2 │ │vContact│ │BACPHLIP│ 🔵 Tools
  │  ORFs  │ │  Tax   │ │Tax+Life│ │   3    │ │Lifestyle│ (Parallel)
  └───┬────┘ └───┬────┘ └───┬────┘ └───┬────┘ └───┬────┘ ⏱️ 3-5 hours
      ↓          ↓          ↓          ↓          ↓
  ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐
  │proteins│ │mmseqs_ │ │phabox_ │ │vc3_    │ │bacphlip│ ⚪ Outputs
  │.faa    │ │tax.tsv │ │*.tsv   │ │output/ │ │life.tsv│
  └───┬────┘ └───┬────┘ └───┬────┘ └───┬────┘ └───┬────┘
      │          └────┬──────┴──────────┘          │
      │               ↓                             │
      │      ┌─────────────────┐                   │
      │      │   R Script      │ 🔵 Tool           │
      │      │ taxonomic_      │ (Optional)        │
      │      │ consensus.R     │                   │
      │      └────────┬────────┘                   │
      │               ↓                             │
      │      ┌─────────────────┐                   │
      │      │ consensus_      │ ⚪ Output          │
      │      │ taxonomy.tsv    │                   │
      │      └────────┬────────┘                   │
      │               │                             │
      │               └─────────┬───────────────────┘
      │                         ↓
      │              ┌─────────────────┐
      │              │   Python        │ 🔵 Tool
      │              │ lifestyle_      │
      │              │ consensus.py    │
      │              └────────┬────────┘
      │                       ↓
      │              ┌─────────────────┐
      │              │ lifestyle_      │ ⚪ Output
      │              │ consensus.tsv   │
      │              └────────┬────────┘
      │                       │
      └───────────────┬───────┴──────────────────────┐
                      ↓                              │
             ┌─────────────────┐                     │
             │   Python        │ 🔵 Tool             │
             │ create_final_   │                     │
             │ contig_table.py │                     │
             └────────┬────────┘                     │
                      ↓                              │
             ┌─────────────────────┐                 │
             │ final_contig_       │ 🟢 Green        │
             │ summary.tsv         │ (MAIN OUTPUT)   │
             └─────────────────────┘ STAGE 3 OUTPUT  │
                                                     │
```

---

## draw.io Layer Strategy

### Recommended Layers:

1. **Background** - Stage containers/swimlanes
2. **Data Flow** - Files (cylinders) and arrows
3. **Tools** - Process boxes
4. **Annotations** - Statistics, parameters, time estimates
5. **Highlights** - Entry points, decision points

**Why layers?**
- Toggle visibility for different detail levels
- Export simplified vs detailed versions
- Edit without affecting other elements

---

## Size Annotations Format

Use **small text boxes** next to cylinders:

```
┌──────────────────┐  ┌───────┐
│  filename.fasta  │  │ 2.5K  │ ← Annotation
└──────────────────┘  │ seqs  │
                      └───────┘
```

OR use **subscript inside**:

```
┌──────────────────────┐
│  filename.fasta      │
│  (2,500 sequences)   │ ← Inside cylinder
└──────────────────────┘
```

---

## Time Estimate Format

Use **clock emoji or small callout**:

```
┌──────────┐
│ MMseqs2  │  ⏱️ 2-3 hrs
└──────────┘

OR

┌──────────┐     ╭─────────╮
│ MMseqs2  │────▶│ 2-3 hrs │
└──────────┘     ╰─────────╯
```

---

## Parameters Format

Use **dashed callout box**:

```
┌─────────┐      ╔═ Parameters ═╗
│ vClust  │◀─────║ ANI: 95%     ║
│ Cluster │      ║ Cov: 85%     ║
└─────────┘      ║ Leiden       ║
                 ╚══════════════╝
```

---

## Complete Workflow at a Glance

For README, create a **zoomed-out version** showing:

```
INPUT → │ STAGE 1 │ → │ STAGE 2? │ → │ STAGE 3 │ → OUTPUT
        │ Predict │   │ Cluster  │   │ Analyze │
        │ (4 hrs) │   │ (1 hr)   │   │ (12 hrs)│
```

Then link to detailed versions for each stage.

---

This style makes it **crystal clear** what files exist at each step and what tools create them!
