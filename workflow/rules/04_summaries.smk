"""
Rules for collecting pipeline progress summaries and generating reports.
These rules are lightweight and designed not to impact pipeline performance.
"""

import os

# Summary directory
SUMMARY_DIR = f"{config['output_dir']}/pipeline_summaries"

# Helper function to determine which summary files to collect based on start_from
def get_required_summaries():
    """
    Return list of summary files that should be collected based on start_from configuration.
    This prevents trying to collect summaries from steps that didn't run.
    """
    summaries = []
    start = config.get("start_from", "raw_contigs")

    # Prediction/filtering summaries (only if starting from raw_contigs)
    if start == "raw_contigs":
        summaries.extend([
            f"{SUMMARY_DIR}/input_stats.json",
            f"{SUMMARY_DIR}/reneo_stats.json",
            f"{SUMMARY_DIR}/filtering_stats.json",
            f"{SUMMARY_DIR}/jaeger_stats.json",
            f"{SUMMARY_DIR}/genomad_stats.json",
            f"{SUMMARY_DIR}/phold_stats.json",
            f"{SUMMARY_DIR}/checkv_stats.json",
            f"{SUMMARY_DIR}/integration_stats.json",
        ])

    # Analysis summaries (always collected)
    summaries.extend([
        f"{SUMMARY_DIR}/iphop_stats.json",
        f"{SUMMARY_DIR}/lifestyle_stats.json",
        f"{SUMMARY_DIR}/consensus_taxonomy.json",
        f"{SUMMARY_DIR}/final_phages.json"
    ])

    # Clustering summary (only if do_clustering is enabled)
    if config.get("do_clustering", True):
        summaries.append(f"{SUMMARY_DIR}/clustering_stats.json")

    return summaries

# Rule to collect input statistics
rule collect_input_stats:
    output:
        summary = f"{SUMMARY_DIR}/input_stats.json"
    params:
        assembly_file = config.get("assembly_file", ""),
        assembly_graph = config.get("assembly_graph", ""),
        reads_dir = config.get("reads_dir", ""),
        reneo_output = f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta"
    log:
        f"{config['output_dir']}/logs/summaries/collect_input_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Determine which assembly source to use
        ASSEMBLY_INPUT=""

        # Priority: 1) Reneo output (if assembly_graph was used), 2) assembly_file, 3) assembly_graph
        if [ -n "{params.assembly_graph}" ] && [ -f "{params.reneo_output}" ]; then
            ASSEMBLY_INPUT="{params.reneo_output}"
            echo "Using Reneo output as assembly input: $ASSEMBLY_INPUT" > {log}
        elif [ -n "{params.assembly_file}" ] && [ -f "{params.assembly_file}" ]; then
            ASSEMBLY_INPUT="{params.assembly_file}"
            echo "Using assembly file as input: $ASSEMBLY_INPUT" > {log}
        elif [ -n "{params.assembly_graph}" ]; then
            echo "Assembly graph specified but no Reneo output found. Using empty assembly." > {log}
            ASSEMBLY_INPUT=""
        else
            echo "No valid assembly input found" > {log}
            ASSEMBLY_INPUT=""
        fi

        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step input_stats \
            --output {output.summary} \
            --inputs assembly_fasta:"$ASSEMBLY_INPUT" reads_dir:{params.reads_dir} \
            >> {log} 2>&1
        """

# Rule to collect Reneo statistics (only if Reneo is used)
rule collect_reneo_stats:
    input:
        reneo_output = f"{config['output_dir']}/01_reneo_output/genomes_and_unresolved_edges_1KB.fasta" if workflow.globals.get("use_reneo", False) else "/dev/null"
    output:
        summary = f"{SUMMARY_DIR}/reneo_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_reneo_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        if [ "{input.reneo_output}" != "/dev/null" ] && [ -f "{input.reneo_output}" ]; then
            python {workflow.basedir}/scripts/collect_step_summary.py \
                --step reneo_stats \
                --output {output.summary} \
                --inputs reneo_contigs:{input.reneo_output} \
                > {log} 2>&1
        else
            echo "Reneo not used, creating empty summary" > {log}
            mkdir -p $(dirname {output.summary})
            echo '{{"step": "reneo_stats", "timestamp": "'$(date -Iseconds)'", "inputs": {{}}, "note": "Reneo not used in this run"}}' > {output.summary}
        fi
        """

# Rule to collect viral filtering statistics
rule collect_filtering_stats:
    input:
        lca_results = f"{config['output_dir']}/01_mmseqs_output/genomes_and_unresolved_edges_mmseqs_lca.tsv",
        viral_contigs = f"{config['output_dir']}/01_filtered_mmseqs/passing_Viralcontigs.fasta"
    output:
        summary = f"{SUMMARY_DIR}/filtering_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_filtering_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step filtering_stats \
            --output {output.summary} \
            --inputs mmseqs_results:{input.lca_results} viral_contigs:{input.viral_contigs} \
            > {log} 2>&1
        """

# Rule to collect Jaeger statistics
rule collect_jaeger_stats:
    input:
        jaeger_results = f"{config['output_dir']}/01_jaeger_output/passing_Viralcontigs_default_jaeger.tsv"
    output:
        summary = f"{SUMMARY_DIR}/jaeger_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_jaeger_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step jaeger_stats \
            --output {output.summary} \
            --inputs jaeger_results:{input.jaeger_results} \
            > {log} 2>&1
        """

# Rule to collect geNomad statistics
rule collect_genomad_stats:
    input:
        genomad_results = f"{config['output_dir']}/01_genomad_output/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv"
    output:
        summary = f"{SUMMARY_DIR}/genomad_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_genomad_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step genomad_stats \
            --output {output.summary} \
            --inputs genomad_results:{input.genomad_results} \
            > {log} 2>&1
        """

# Rule to collect PHOLD statistics
rule collect_phold_stats:
    input:
        phold_results = f"{config['output_dir']}/01_phold_output/phold_per_cds_predictions.tsv"
    output:
        summary = f"{SUMMARY_DIR}/phold_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_phold_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step phold_stats \
            --output {output.summary} \
            --inputs phold_results:{input.phold_results} \
            > {log} 2>&1
        """

# Rule to collect CheckV statistics
rule collect_checkv_stats:
    input:
        checkv_results = f"{config['output_dir']}/01_checkv_output/quality_summary.tsv"
    output:
        summary = f"{SUMMARY_DIR}/checkv_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_checkv_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step checkv_stats \
            --output {output.summary} \
            --inputs checkv_results:{input.checkv_results} \
            > {log} 2>&1
        """

# Rule to collect integration statistics
rule collect_integration_stats:
    input:
        phage_predictions = f"{config['output_dir']}/01_phage_predictions/phagePredictedContigs.tsv",
        phage_contigs = f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
    output:
        summary = f"{SUMMARY_DIR}/integration_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_integration_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step integration_stats \
            --output {output.summary} \
            --inputs phage_predictions:{input.phage_predictions} phage_contigs:{input.phage_contigs} \
            > {log} 2>&1
        """

# Rule to collect iPhop statistics
rule collect_iphop_stats:
    input:
        iphop_results = f"{config['output_dir']}/03_iphop_results/iphop_predictions_compiled.tsv"
    output:
        summary = f"{SUMMARY_DIR}/iphop_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_iphop_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step iphop_stats \
            --output {output.summary} \
            --inputs iphop_results:{input.iphop_results} \
            > {log} 2>&1
        """

# Rule to collect lifestyle prediction statistics
rule collect_lifestyle_stats:
    input:
        phabox_lifestyle = f"{config['output_dir']}/03_genomic_info/phabox_output/lifestyle.tsv"
        # bacphlip_results = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"  # Disabled for now
    output:
        summary = f"{SUMMARY_DIR}/lifestyle_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_lifestyle_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Count lifestyle predictions from Phabox2 (only source)
        PHABOX_COUNT=0
        if [ -f "{input.phabox_lifestyle}" ]; then
            # Count lines excluding header (subtract 1)
            PHABOX_COUNT=$(( $(wc -l < "{input.phabox_lifestyle}") - 1 ))
            echo "Found $PHABOX_COUNT Phabox2 lifestyle predictions" >> {log}
        fi
        
        # BACPHLIP disabled for now - uncomment below to re-enable
        # BACPHLIP_COUNT=0
        # BACPHLIP_FILE="{config[output_dir]}/03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"
        # if [ -f "$BACPHLIP_FILE" ]; then
        #     BACPHLIP_COUNT=$(( $(wc -l < "$BACPHLIP_FILE") - 1 ))
        #     echo "Found $BACPHLIP_COUNT BACPHLIP lifestyle predictions" >> {log}
        # fi
        
        # Use Phabox2 count (only source currently enabled)
        TOTAL_COUNT=$PHABOX_COUNT
        
        # Create summary JSON
        cat > {output.summary} << EOF
{{
  "step": "lifestyle_stats",
  "timestamp": "$(date -Iseconds)",
  "inputs": {{
    "lifestyle_results": {{
      "tool": "phabox2",
      "total_predictions": $TOTAL_COUNT,
      "phabox2_predictions": $PHABOX_COUNT
    }}
  }},
  "outputs": {{}},
  "statistics": {{}}
}}
EOF
        
        echo "Lifestyle stats summary created with $TOTAL_COUNT Phabox2 predictions" >> {log}
        """

# Rule to collect clustering statistics (only if clustering is enabled)
rule collect_clustering_stats:
    input:
        cluster_reps = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta" if config.get("do_clustering", True) else "/dev/null"
    output:
        summary = f"{SUMMARY_DIR}/clustering_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_clustering_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        if [ "{input.cluster_reps}" != "/dev/null" ] && [ -f "{input.cluster_reps}" ]; then
            python {workflow.basedir}/scripts/collect_step_summary.py \
                --step clustering_stats \
                --output {output.summary} \
                --inputs cluster_reps:{input.cluster_reps} \
                > {log} 2>&1
        else
            echo "Clustering not enabled, creating empty summary" > {log}
            mkdir -p $(dirname {output.summary})
            echo '{{"step": "clustering_stats", "timestamp": "'$(date -Iseconds)'", "inputs": {{}}, "note": "Clustering not enabled in this run"}}' > {output.summary}
        fi
        """

# Rule to collect consensus taxonomy statistics (only if consensus is enabled)
rule collect_consensus_stats:
    input:
        consensus_taxonomy = f"{config['output_dir']}/03_genomic_info/consensus_taxonomy.tsv" if config.get("run_consensus", True) else "/dev/null",
        consensus_summary = f"{config['output_dir']}/03_genomic_info/consensus_taxonomy_summary.json" if config.get("run_consensus", True) else "/dev/null"
    output:
        summary = f"{SUMMARY_DIR}/consensus_taxonomy.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_consensus_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        if [ "{input.consensus_taxonomy}" != "/dev/null" ] && [ -f "{input.consensus_taxonomy}" ]; then
            # Load the consensus summary JSON that was already created
            if [ -f "{input.consensus_summary}" ]; then
                # Copy and modify the existing summary to match our format
                python -c "
import json
import sys
from datetime import datetime

try:
    with open('{input.consensus_summary}', 'r') as f:
        data = json.load(f)
    
    # Reformat for our summary system
    summary = {{
        'step': 'consensus_taxonomy',
        'timestamp': datetime.now().isoformat(),
        'inputs': {{
            'consensus_taxonomy': data
        }}
    }}
    
    with open('{output.summary}', 'w') as f:
        json.dump(summary, f, indent=2)
    
    total_contigs = data.get('total_contigs', 0)
    print('Consensus taxonomy summary created: ' + str(total_contigs) + ' contigs')
    
except Exception as e:
    print('Error processing consensus summary: ' + str(e), file=sys.stderr)
    # Create minimal summary
    summary = {{
        'step': 'consensus_taxonomy',
        'timestamp': datetime.now().isoformat(),
        'inputs': {{}},
        'note': 'Error reading consensus taxonomy results'
    }}
    with open('{output.summary}', 'w') as f:
        json.dump(summary, f, indent=2)
" > {log} 2>&1
            else
                echo "Consensus summary file not found, creating basic summary" > {log}
                python {workflow.basedir}/scripts/collect_step_summary.py \
                    --step consensus_taxonomy \
                    --output {output.summary} \
                    --inputs consensus_file:{input.consensus_taxonomy} \
                    > {log} 2>&1
            fi
        else
            echo "Consensus taxonomy not enabled, creating empty summary" > {log}
            mkdir -p $(dirname {output.summary})
            echo '{{"step": "consensus_taxonomy", "timestamp": "'$(date -Iseconds)'", "inputs": {{}}, "note": "Consensus taxonomy not enabled in this run"}}' > {output.summary}
        fi
        """

# Helper function to get final sequences for summary
def get_final_sequences_for_summary(wildcards):
    """
    Get the final sequences file based on start_from and do_clustering settings.
    Handles all cases including starting from intermediate points.
    """
    start_from = config.get("start_from", "raw_contigs")

    # If starting from clustering with user-provided file
    if start_from == "clustering":
        clustered_file = config.get("input_clustered_seqs", "")
        if clustered_file:
            return clustered_file
        # Fallback to expected output location
        return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"

    # If starting from phage_contigs
    if start_from == "phage_contigs":
        if config.get("do_clustering", True):
            # Clustering was run on provided phage contigs
            return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
        else:
            # No clustering, use provided phage contigs directly
            return config.get("input_phage_contigs", "")

    # Default: starting from raw_contigs
    if config.get("do_clustering", True):
        return f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta"
    else:
        return f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"

# Rule to collect final statistics
rule collect_final_stats:
    input:
        # Use helper function to determine final sequences based on pipeline configuration
        final_seqs = get_final_sequences_for_summary
    output:
        summary = f"{SUMMARY_DIR}/final_phages.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_final_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step final_phages \
            --output {output.summary} \
            --inputs final_sequences:{input.final_seqs} \
            > {log} 2>&1
        """

# Rule to generate the final summary report
rule generate_summary_report:
    input:
        # Depend on all summaries being collected
        # The actual summaries collected depend on start_from configuration
        summaries_collected = f"{SUMMARY_DIR}/.all_summaries_collected"
    output:
        report = f"{config['output_dir']}/Pipeline_Summary_Report.html"
    params:
        assembly_file = config.get('assembly_file', ''),
        assembly_graph = config.get('assembly_graph', ''),
        reads_dir = config.get('reads_dir', ''),
        output_dir = config['output_dir'],
        do_clustering = config.get('do_clustering', True)
    log:
        f"{config['output_dir']}/logs/summaries/generate_summary_report.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Create config info for the report
        CONFIG_INFO="<tr><td>Output Directory</td><td>{params.output_dir}</td></tr>"
        
        # Only show the assembly input that was actually used
        if [ -n "{params.assembly_file}" ] && [ -f "{params.assembly_file}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly File</td><td>{params.assembly_file}</td></tr>"
        elif [ -n "{params.assembly_graph}" ] && [ -f "{params.assembly_graph}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly Graph</td><td>{params.assembly_graph}</td></tr>"
        elif [ -n "{params.assembly_file}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly File</td><td>{params.assembly_file}</td></tr>"
        elif [ -n "{params.assembly_graph}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly Graph</td><td>{params.assembly_graph}</td></tr>"
        fi
        
        CONFIG_INFO="$CONFIG_INFO<tr><td>Reads Directory</td><td>{params.reads_dir}</td></tr>"
        CONFIG_INFO="$CONFIG_INFO<tr><td>Clustering Enabled</td><td>{params.do_clustering}</td></tr>"
        
        python {workflow.basedir}/scripts/generate_summary_report.py \
            --summary-dir {SUMMARY_DIR} \
            --output {output.report} \
            --config "$CONFIG_INFO" \
            > {log} 2>&1
        
        echo "Summary report generated: {output.report}" >> {log}
        echo "You can open it in a web browser to view the results." >> {log}
        """

# Rule to run all summary collection (helper rule)
rule collect_all_summaries:
    input:
        get_required_summaries()
    output:
        flag = f"{SUMMARY_DIR}/.all_summaries_collected"
    shell:
        "touch {output.flag}"