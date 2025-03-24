# Parallelizing Array-Based Tasks in Snakemake

The current implementation processes files sequentially within a single rule, which differs from the original SLURM array-based approach. Here's how to modify the workflow to achieve better parallelization while using your existing Slurm profile.

## 1. Add Snakemake imports at the top of your workflow/Snakefile

```python
import os
import glob
from snakemake.utils import min_version

# Ensure minimum Snakemake version
min_version("7.0")
```

## 2. Replace the iPhop host prediction rule with these rules

Replace the current `iphop_host_prediction` rule in `workflow/rules/03_analysis.smk` with:

```python
# 2a. Run iPhop for host prediction on a single split file
rule iphop_single_prediction:
    input:
        phage_file = f"{config['output_dir']}/03_split_seqs/{{sample}}.fasta"
    output:
        prediction = f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv"
    log:
        f"{config['output_dir']}/logs/iphop_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["iphop"]
    resources:
        mem_mb = config["resources"]["iphop"]["mem_mb"],
        threads = config["resources"]["iphop"]["threads"],
        time = config["resources"]["iphop"]["time"]
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.prediction})
        
        # Run iPhop
        iphop predict --fa_file {input.phage_file} \
            --out_dir $(dirname {output.prediction}) \
            --num_threads {resources.threads} > {log} 2>&1
        """

# 2b. Aggregate iPhop results
rule iphop_aggregate_results:
    input:
        # Dynamic input based on the split files
        lambda wildcards: expand(
            f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv",
            sample=[os.path.splitext(os.path.basename(f))[0] 
                   for f in glob.glob(f"{config['output_dir']}/03_split_seqs/*.fasta")]
        )
    output:
        results_dir = directory(f"{config['output_dir']}/03_iphop_results"),
        predictions = f"{config['output_dir']}/03_iphop_results/iphop_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/iphop_aggregate_results.log"
    shell:
        """
        # Ensure results directory exists
        mkdir -p {output.results_dir}
        
        # Compile results
        echo "Compiling results" > {log}
        head -n 1 {input[0]} > {output.predictions}.tmp
        cat {input} | grep -v "query" >> {output.predictions}.tmp
        
        # Convert to TSV
        tr ',' '\\t' < {output.predictions}.tmp > {output.predictions}
        rm {output.predictions}.tmp
        """
```

## 3. Replace the PHACTS lifestyle prediction rule

Replace the current `phacts_lifestyle_prediction` rule with:

```python
# 5a. Run PHACTS for lifestyle prediction on a single protein file
rule phacts_single_prediction:
    input:
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa"
    output:
        result = f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}.phacts.out"
    log:
        f"{config['output_dir']}/logs/phacts_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    resources:
        mem_mb = config["resources"]["phacts"]["mem_mb"],
        threads = config["resources"]["phacts"]["threads"],
        time = config["resources"]["phacts"]["time"]
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.result})
        
        # Run PHACTS
        phacts.py {input.protein_file} {output.result} > {log} 2>&1
        """

# 5b. Aggregate PHACTS results
rule phacts_aggregate_results:
    input:
        # Dynamic input based on the split files
        lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}.phacts.out",
            sample=[os.path.splitext(os.path.basename(f))[0] 
                   for f in glob.glob(f"{config['output_dir']}/03_split_proteins/*.faa")]
        )
    output:
        results_dir = directory(f"{config['output_dir']}/03_phacts_results"),
        predictions = f"{config['output_dir']}/03_phacts_results/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_aggregate_results.log"
    shell:
        """
        # Ensure results directory exists
        mkdir -p {output.results_dir}
        
        # Compile results
        echo -e "phage_id\\tlifestyle\\tprobability" > {output.predictions}
        
        for file in {input}; do
            phage_id=$(basename "$file" .phacts.out)
            lifestyle=$(grep "Lifestyle:" "$file" | awk '{{print $2}}')
            probability=$(grep "Probability:" "$file" | awk '{{print $2}}')
            echo -e "$phage_id\\t$lifestyle\\t$probability" >> {output.predictions}
        done
        """
```

## 4. Update the all rule to reference aggregated results

The existing `all` rule already references the correct aggregated output files, so no changes are needed there.

## 5. How to run with your Slurm profile

To run the workflow with your existing Slurm profile, simply use:

```bash
snakemake --profile profile/slurm
```

This will use all the settings you've defined in your `profile/slurm/config.yaml`, including:
- Default resource allocation
- Job concurrency (40 jobs)
- Automatic conda environment activation
- Special resource settings for specific rules

## Benefits of this approach:

1. **True parallelization**: Each file gets processed as a separate Slurm job
2. **Efficient resource usage**: Resources are allocated per job, not shared across all files
3. **Fault tolerance**: If one job fails, only that one needs to be rerun
4. **Cluster integration**: Fully utilizes your Slurm profile
5. **Preserves workflow structure**: Maintains the same overall workflow structure

This approach brings the workflow closer to the original SLURM array design while leveraging Snakemake's dependency tracking and your Slurm profile configuration.

## Implementation Notes:

1. The `glob.glob()` function is used to dynamically find all split files
2. Resources for each job come from your existing configuration
3. The aggregation step only runs after all individual jobs complete
4. The pattern `{sample}` in the rule is a wildcard that Snakemake will match to file names