# Modular Pipeline Documentation

## Overview

The phage-analysis pipeline now supports modular execution, allowing you to choose which analysis stages to run. This makes the pipeline more flexible and efficient for different use cases.

## Pipeline Modules

The pipeline consists of 5 modules:

1. **prediction** (required) - Phage prediction from assembly
2. **clustering** (optional) - Cluster predicted phages into vOTUs
3. **analysis** (optional) - Taxonomic, lifestyle, and host prediction
4. **summaries** (optional) - Generate HTML summary reports
5. **abundance** (optional) - Calculate per-sample abundance matrices

## Configuration Methods

### Method 1: Quick Presets (Recommended)

Use the `analysis_type` parameter for predefined workflows:

```yaml
# Prediction only - just identify phages, no downstream analysis
analysis_type: "prediction_only"

# Standard - complete analysis without abundance (recommended for most users)
analysis_type: "standard"

# Comprehensive - everything including abundance calculation
analysis_type: "comprehensive"
```

**Preset Details:**

| Preset | Prediction | Clustering | Analysis | Summaries | Abundance |
|--------|-----------|-----------|----------|-----------|-----------|
| `prediction_only` | ✓ | ✗ | ✗ | ✗ | ✗ |
| `standard` | ✓ | ✓ | ✓ | ✓ | ✗ |
| `comprehensive` | ✓ | ✓ | ✓ | ✓ | ✓ |

### Method 2: Granular Control (Advanced)

For fine-grained control, comment out `analysis_type` and use the `modules` section:

```yaml
# Comment out or remove analysis_type
# analysis_type: "standard"

modules:
  clustering: true      # Cluster predicted phages into vOTUs
  analysis: true        # Run taxonomic/lifestyle/host prediction
  summaries: true       # Generate summary reports
  abundance: false      # Calculate per-sample abundance
```

## Module Dependencies

Some modules depend on others:

- **abundance** requires **clustering** (needs vOTUs for abundance calculation)
- **analysis** requires **prediction** (needs predicted sequences)
- **summaries** works best with **analysis** (generates reports from analysis results)

The pipeline will validate these dependencies and show helpful error messages if they're not met.

## Examples

### Example 1: Quick Phage Identification

Just want to identify phages without any clustering or analysis:

```yaml
analysis_type: "prediction_only"
```

**Outputs:**
- Filtered assembly (1KB+ sequences)
- Phage predictions from geNomad

### Example 2: Standard Analysis (Recommended)

Complete analysis pipeline without abundance:

```yaml
analysis_type: "standard"
```

**Outputs:**
- Filtered assembly
- Phage predictions
- vOTU clustering
- Taxonomy (MMseqs2, Phabox2, vContact3)
- Lifestyle prediction (BACPHLIP, Phabox2)
- Host prediction (iPHoP)
- Quality assessment (CheckV)
- HTML summary report

### Example 3: Full Analysis with Abundance

Everything including per-sample abundance:

```yaml
analysis_type: "comprehensive"
```

**Outputs:**
- All outputs from "standard" plus:
- TPM abundance matrix
- Count abundance matrix
- ORF-level abundance matrices

### Example 4: Custom Configuration

Clustering and analysis but skip summaries and abundance:

```yaml
modules:
  clustering: true
  analysis: true
  summaries: false
  abundance: false
```

## Migration from Old Config

The old configuration flags are **deprecated** but still supported:

```yaml
# OLD (deprecated)
do_clustering: true
generate_summaries: true
calculate_abundance: false

# NEW (recommended)
analysis_type: "standard"
```

If you use the old flags, you'll see a deprecation warning but the pipeline will still work.

## Validation

The pipeline validates your configuration on startup:

1. Checks for invalid `analysis_type` values
2. Validates module dependencies
3. Prints a clear module status summary
4. Shows helpful error messages for configuration issues

Example output:
```
================================================================================
Pipeline Modules Configuration
================================================================================
  prediction     : ENABLED
  clustering     : ENABLED
  analysis       : ENABLED
  summaries      : ENABLED
  abundance      : DISABLED
================================================================================
```

## Testing Your Configuration

To test without actually running the pipeline, use a dry-run:

```bash
snakemake -n --cores 1
```

This will show you what rules would be executed without actually running them.

## See Also

- `config/config.yaml` - Main configuration file with examples
- `config/examples/` - Example configurations for different use cases
  - `config_prediction_only.yaml` - Prediction only example
  - `config_comprehensive.yaml` - Comprehensive analysis example
  - `config_granular.yaml` - Granular module control example
