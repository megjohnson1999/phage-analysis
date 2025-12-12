# Phage Analysis Workflow - draw.io Diagram Structure

## Overall Layout
- **Type**: Swimlane/Cross-functional flowchart
- **Orientation**: Top-to-bottom (recommended) or Left-to-right
- **Page Size**: A3 or custom wide format

---

## Components Breakdown

### HEADER SECTION
```
┌─────────────────────────────────────────────────────────────────┐
│              PHAGE ANALYSIS WORKFLOW                            │
│  Metagenomic Phage Identification, Clustering & Characterization│
└─────────────────────────────────────────────────────────────────┘
```
- Shape: Rectangle with bold border
- Color: Dark blue header

### ENTRY POINTS PANEL (Left sidebar or top bar)
```
┌─────────────────┐
│ START OPTIONS   │
├─────────────────┤
│ ● raw_contigs   │──────┐
│ ● reneo         │────┐ │
│ ● viral_contigs │──┐ │ │
│ ● phage_contigs │┐ │ │ │
│ ● clustering    ││ │ │ │
└─────────────────┘│ │ │ │
                   ↓ ↓ ↓ ↓
              (Arrows to respective stages)
```
- Shape: List container with bullet points
- Color: Yellow/gold background
- Arrows: Different colors for each entry point

---

## STAGE 1: PREDICTION

### Container
- Shape: Swimlane
- Color: Light blue (#E3F2FD)
- Label: "STAGE 1: PHAGE PREDICTION"

### Sub-section 1A: Assembly Processing
```
┌─ INPUT ─────────────────────────────────────────────────┐
│                                                          │
│  ┌────────────┐              ┌──────────────┐          │
│  │ Assembly   │              │   Reads      │          │
│  │  Graph     │              │  Directory   │          │
│  │  (.gfa)    │              │              │          │
│  └──────┬─────┘              └──────┬───────┘          │
│         │                           │                   │
│         └───────────┬───────────────┘                   │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │   Decision Point:     │                       │
│         │  Assembly Type?       │                       │
│         └─────────┬─────────────┘                       │
│              ┌────┴────┐                                │
│              ↓         ↓                                │
│    ┌──────────────┐  ┌────────────────┐               │
│    │    Reneo     │  │ Direct Filter  │               │
│    │   Binning    │  │   FASTA ≥1KB   │               │
│    │  (Gurobi)    │  │   (seqkit)     │               │
│    └──────┬───────┘  └────────┬───────┘               │
│           │                   │                         │
│           └─────────┬─────────┘                         │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │  Filter Contigs ≥1KB  │                       │
│         │      (seqkit)         │                       │
│         └───────────┬───────────┘                       │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │ Filtered Assembly     │ ← Entry: reneo        │
│         │      1KB.fasta        │                       │
│         └───────────────────────┘                       │
└──────────────────────────────────────────────────────────┘
```

### Sub-section 1B: Viral Identification
```
┌─ VIRAL FILTERING ───────────────────────────────────────┐
│                                                          │
│         ┌───────────────────────┐                       │
│         │  MMseqs2 Taxonomy     │                       │
│         │  (LCA Assignment)     │                       │
│         │   [NR Database]       │                       │
│         └───────────┬───────────┘                       │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │ Filter Viral Lineages │                       │
│         │   (Python script)     │                       │
│         └───────────┬───────────┘                       │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │ Extract Viral Contigs │                       │
│         │      (seqkit)         │                       │
│         └───────────┬───────────┘                       │
│                     ↓                                    │
│         ┌───────────────────────┐                       │
│         │ Viral Contigs FASTA   │ ← Entry: viral_contigs│
│         └───────────────────────┘                       │
└──────────────────────────────────────────────────────────┘
```

### Sub-section 1C: Phage Prediction (Parallel)
```
┌─ PHAGE CLASSIFICATION (Parallel Processing) ────────────┐
│                                                          │
│              ┌──────────────────┐                       │
│              │ Viral Contigs    │                       │
│              └────────┬─────────┘                       │
│                       │                                  │
│         ┌─────────────┼──────────────┬─────────┐       │
│         ↓             ↓              ↓         ↓       │
│   ┌─────────┐  ┌──────────┐  ┌──────────┐ ┌────────┐ │
│   │ Jaeger  │  │ geNomad  │  │  CheckV  │ │ PHOLD  │ │
│   │ (phage) │  │  (viral) │  │(quality) │ │(annot) │ │
│   │  Score  │  │  Score   │  │          │ │ ∥∥∥∥∥  │ │ ← Parallelized
│   │   2.5   │  │  ≥0.6    │  │          │ │1000/ch │ │   (checkpoint)
│   └─────┬───┘  └─────┬────┘  └────┬─────┘ └───┬────┘ │
│         │            │            │           │       │
│         └────────────┼────────────┼───────────┘       │
│                      ↓                                 │
│         ┌────────────────────────┐                    │
│         │  Integration (R)       │                    │
│         │ Combine Predictions    │                    │
│         └────────────┬───────────┘                    │
│                      ↓                                 │
│         ┌────────────────────────┐                    │
│         │  Phage Contigs FASTA   │ ← Entry: phage     │
│         │  + Predictions TSV     │    _contigs        │
│         └────────────────────────┘                    │
└──────────────────────────────────────────────────────┘
```

**Visual Elements:**
- Diamond: Decision point for assembly type
- Parallel boxes: Same height, connected from single point above, merge to single point below
- Dashed box around PHOLD with "∥∥∥∥∥" symbol to indicate parallelization
- Cylinder shapes for FASTA/TSV files
- Database icon (parallelogram) for NR Database

---

## STAGE 2: CLUSTERING

### Container
- Shape: Swimlane
- Color: Light green (#E8F5E9)
- Label: "STAGE 2: CLUSTERING (Optional)"

```
┌─ CLUSTERING (Optional) ──────────────────────────────────┐
│                                                           │
│              ┌──────────────────┐                        │
│              │ Phage Contigs    │                        │
│              └────────┬─────────┘                        │
│                       ↓                                   │
│              ┌──────────────────┐                        │
│              │  do_clustering?  │ ◇ Decision             │
│              └────┬────────┬────┘                        │
│                   │        │                              │
│              YES  │        │ NO                           │
│                   ↓        ↓                              │
│         ┌──────────────┐  │                              │
│         │    vClust    │  │                              │
│         │  Clustering  │  │                              │
│         ├──────────────┤  │                              │
│         │ • Prefilter  │  │                              │
│         │   95% ID     │  │                              │
│         │ • Align      │  │                              │
│         │ • Cluster    │  │                              │
│         │   Leiden     │  │                              │
│         │   85% cov    │  │                              │
│         └──────┬───────┘  │                              │
│                ↓           │                              │
│    ┌──────────────────┐   │                              │
│    │ Select Longest   │   │                              │
│    │  per Cluster     │   │                              │
│    └──────┬───────────┘   │                              │
│           │               │                               │
│           ↓               ↓                               │
│    ┌─────────────────────────┐                          │
│    │  vOTU Representatives   │ ← Entry: clustering       │
│    │    OR All Phages        │                           │
│    └─────────────────────────┘                          │
└───────────────────────────────────────────────────────────┘
```

**Visual Elements:**
- Diamond for decision point
- Dashed border for entire stage (optional)
- Two paths (YES/NO) with different arrow styles
- Callout box for vClust parameters

---

## STAGE 3: ANALYSIS

### Container
- Shape: Swimlane
- Color: Light orange (#FFF3E0)
- Label: "STAGE 3: ANALYSIS & CHARACTERIZATION"

```
┌─ ANALYSIS ────────────────────────────────────────────────┐
│                                                            │
│         ┌────────────────────────┐                        │
│         │ Phage Sequences        │                        │
│         │ (Clustered or All)     │                        │
│         └───────────┬────────────┘                        │
│                     ↓                                      │
│         ┌────────────────────────┐                        │
│         │  CheckV Final Quality  │                        │
│         └───────────┬────────────┘                        │
│                     ↓                                      │
│  ┌─ HOST PREDICTION (Parallelized) ──────────────┐       │
│  │         ┌────────────────┐                     │       │
│  │         │ Split Sequences│                     │       │
│  │         │  (100/chunk)   │                     │       │
│  │         └───────┬────────┘                     │       │
│  │                 ↓                               │       │
│  │         ┌────────────────┐                     │       │
│  │         │  iPhop Array   │  ∥∥∥∥∥              │       │
│  │         │  24 threads/job│  Parallel           │       │
│  │         └───────┬────────┘                     │       │
│  │                 ↓                               │       │
│  │         ┌────────────────┐                     │       │
│  │         │   Aggregate    │                     │       │
│  │         │   Predictions  │                     │       │
│  │         └───────┬────────┘                     │       │
│  └─────────────────┼──────────────────────────────┘       │
│                    │                                       │
│  ┌─ CHARACTERIZATION (Parallel) ─────────────────┐       │
│  │                 │                              │       │
│  │     ┌───────────┼──────────┬────────┐         │       │
│  │     ↓           ↓          ↓        ↓         │       │
│  │ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐  │       │
│  │ │MMseqs2 │ │Phabox2 │ │vContact│ │BACPHLIP│  │       │
│  │ │  Tax   │ │Tax+Life│ │   3    │ │Lifestyle│  │       │
│  │ └───┬────┘ └───┬────┘ └───┬────┘ └───┬────┘  │       │
│  │     │          │          │          │        │       │
│  │     └──────────┼──────────┼──────────┘        │       │
│  │                ↓                               │       │
│  │     ┌────────────────────┐                    │       │
│  │     │ Consensus Building │                    │       │
│  │     │ • Taxonomy         │                    │       │
│  │     │ • Lifestyle        │                    │       │
│  │     └────────┬───────────┘                    │       │
│  └──────────────┼────────────────────────────────┘       │
│                 ↓                                         │
│     ┌───────────────────────┐                            │
│     │  Final Summary Table  │                            │
│     │  (All Results + QC)   │                            │
│     └───────────────────────┘                            │
└────────────────────────────────────────────────────────────┘
```

**Visual Elements:**
- Nested containers for parallelized sections
- Parallel processing indicators (∥∥∥∥∥)
- Multiple parallel boxes at same height
- Different colored borders for sub-sections

---

## STAGE 4: REPORTING

### Container
- Shape: Swimlane
- Color: Light purple (#F3E5F5)
- Label: "STAGE 4: REPORTING (Optional)"
- Border: Dashed (optional stage)

```
┌─ REPORTING (Optional) ────────────────────────────────────┐
│                                                            │
│              ┌──────────────────┐                         │
│              │ All Results      │                         │
│              └────────┬─────────┘                         │
│                       ↓                                    │
│              ┌──────────────────┐                         │
│              │ generate_        │ ◇ Decision              │
│              │  summaries?      │                         │
│              └────┬────────┬────┘                         │
│                   │        │                               │
│              YES  │        │ NO (skip)                     │
│                   ↓        └─────────────────┐            │
│         ┌──────────────────┐                 │            │
│         │ Collect Step     │                 │            │
│         │  Summaries       │                 │            │
│         │  (JSON files)    │                 │            │
│         └────────┬─────────┘                 │            │
│                  ↓                            │            │
│         ┌──────────────────┐                 │            │
│         │  Generate HTML   │                 │            │
│         │  Summary Report  │                 │            │
│         └────────┬─────────┘                 │            │
│                  │                            │            │
│                  └────────────┬───────────────┘            │
│                               ↓                            │
│                      ┌────────────────┐                   │
│                      │  Pipeline Done │                   │
│                      └────────────────┘                   │
└────────────────────────────────────────────────────────────┘
```

---

## LEGEND / KEY (Bottom of diagram)

```
┌─ LEGEND ──────────────────────────────────────────────────┐
│                                                            │
│  Shapes:                                                   │
│  ┌────────┐  Process/Tool                                │
│  │        │                                                │
│  └────────┘                                                │
│                                                            │
│  ◇  Decision Point                                        │
│                                                            │
│  ┌────────┐  Input/Output File                           │
│  │ (cyl)  │                                                │
│  └────────┘                                                │
│                                                            │
│  ▱  Database                                              │
│                                                            │
│  ┌ ─ ─ ─ ┐  Optional Stage                               │
│  │       │                                                │
│  └ ─ ─ ─ ┘                                                │
│                                                            │
│  ∥∥∥∥∥  Parallelized Processing                          │
│                                                            │
│  Colors:                                                   │
│  🔵 Prediction  🟢 Clustering  🟠 Analysis  🟣 Reporting  │
│                                                            │
│  Entry Points:                                             │
│  ● raw_contigs  ● reneo  ● viral_contigs                 │
│  ● phage_contigs  ● clustering                            │
└────────────────────────────────────────────────────────────┘
```

---

## ADDITIONAL CALLOUT BOXES

Add these as floating text boxes with arrows pointing to relevant sections:

### Parallelization Details
```
┌─ PARALLEL PROCESSING ─┐
│ • PHOLD: 1000 seqs/job│
│   Checkpoint-based    │
│                       │
│ • iPhop: 100 seqs/job │
│   Checkpoint-based    │
│   24 threads each     │
└───────────────────────┘
```

### Key Parameters
```
┌─ CLUSTERING PARAMS ─┐
│ • ANI: 95%          │
│ • Coverage: 85%     │
│ • Algorithm: Leiden │
│ • Min length: 10KB  │
└─────────────────────┘
```

### Database Requirements
```
┌─ REQUIRED DATABASES ──┐
│ • MMseqs2 NR          │
│ • iPhop (Aug 2023)    │
│ • Phabox2 v2          │
│ • vContact3 v223      │
│ • geNomad             │
│ • CheckV v1.5         │
│ • PHOLD               │
│ • taxonomizr          │
└───────────────────────┘
```

---

## draw.io TIPS

### Recommended Settings:
1. **Grid**: Enable snap to grid (10px)
2. **Connectors**: Use orthogonal or curved connectors
3. **Shadows**: Subtle shadows on major boxes for depth
4. **Fonts**:
   - Headers: 14pt bold
   - Process boxes: 10pt
   - Labels: 8pt
5. **Arrow styles**:
   - Main flow: Thick solid arrows
   - Optional paths: Dashed arrows
   - Entry points: Colored arrows with distinct markers

### Layers to Use:
1. **Background**: Swimlanes and containers
2. **Processes**: All process boxes
3. **Files**: Input/output cylinders
4. **Decisions**: Diamond shapes
5. **Annotations**: Callout boxes and labels
6. **Entry Points**: Arrows showing where start_from options enter

### Export Settings:
- Format: PNG (for README) or SVG (for web/scaling)
- Resolution: 300 DPI minimum
- Transparent background: Optional
- Include border: 10px padding

---

## FILE ORGANIZATION

Save multiple versions:
1. `workflow_overview.drawio` - Simple 4-stage view
2. `workflow_detailed.drawio` - Full detailed version
3. `workflow_with_entrypoints.drawio` - Entry points highlighted

Export as:
- `workflow_overview.png` - For README
- `workflow_detailed.svg` - For detailed documentation
- `workflow_entrypoints.png` - For configuration docs

---

## RECOMMENDED APPROACH

**For README**: Create 2 diagrams
1. **Overview diagram**: Simple 4-stage flow with minimal detail
2. **Detailed diagram**: Full workflow with all components (link to separate file)

**For Documentation**: Create 3-4 diagrams
1. Entry points map
2. Stage 1 detailed
3. Stages 2-4 detailed
4. Parallelization focus

This gives users a quick understanding from README while providing detailed references when needed.
