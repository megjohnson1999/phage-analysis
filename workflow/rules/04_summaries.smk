"""
Rules for collecting pipeline progress summaries and generating reports.
These rules are lightweight and designed not to impact pipeline performance.
"""

import os

# Summary directory
SUMMARY_DIR = f"{config['output_dir']}/pipeline_summaries"

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

# Rule to collect lifestyle prediction statistics from consensus
rule collect_lifestyle_stats:
    input:
        lifestyle_consensus = f"{config['output_dir']}/03_genomic_info/lifestyle_consensus.tsv"
    output:
        summary = f"{SUMMARY_DIR}/lifestyle_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_lifestyle_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo "Collecting lifestyle consensus statistics..." > {log}

        python -c "
import pandas as pd
import json
from datetime import datetime

# Read lifestyle consensus
df = pd.read_csv('{input.lifestyle_consensus}', sep='\t')

# Count totals
total = len(df)

# Count by source
source_counts = df['source'].value_counts().to_dict() if 'source' in df.columns else {{}}

# Count by lifestyle
lifestyle_counts = df['lifestyle'].value_counts().to_dict() if 'lifestyle' in df.columns else {{}}

# Calculate mean confidence
mean_confidence = df['confidence'].mean() if 'confidence' in df.columns else 0

# Create summary
summary = {{
    'step': 'lifestyle_stats',
    'timestamp': datetime.now().isoformat(),
    'inputs': {{
        'lifestyle_consensus': {{
            'total_predictions': total,
            'by_source': source_counts,
            'by_lifestyle': lifestyle_counts,
            'mean_confidence': round(mean_confidence, 3)
        }}
    }}
}}

# Write to JSON
with open('{output.summary}', 'w') as f:
    json.dump(summary, f, indent=2)

print(f'Lifestyle consensus stats: {{total}} predictions')
print(f'Sources: {{source_counts}}')
print(f'Lifestyles: {{lifestyle_counts}}')
" >> {log} 2>&1

        echo "Lifestyle stats summary created from consensus" >> {log}
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

# Rule to collect final statistics
rule collect_final_stats:
    input:
        # Use clustering output if available, otherwise use phage predictions
        final_seqs = f"{config['output_dir']}/02_clustering/vOTU_repSeqs.fasta" if config.get("do_clustering", True) else f"{config['output_dir']}/01_phage_predictions/phageContigs.fasta"
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

# Rule to collect final contig summary table statistics
rule collect_final_summary_stats:
    input:
        final_summary = f"{config['output_dir']}/final_contig_summary.tsv"
    output:
        summary = f"{SUMMARY_DIR}/final_summary.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_final_summary_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        echo "Collecting final contig summary statistics..." > {log}

        python -c "
import pandas as pd
import json
from datetime import datetime

# Read final summary table
df = pd.read_csv('{input.final_summary}', sep='\t')

# Count totals
total_contigs = len(df)

# Count annotation coverage
has_taxonomy = df['phylum'].notna().sum() if 'phylum' in df.columns else 0
has_lifestyle = df['lifestyle'].notna().sum() if 'lifestyle' in df.columns else 0
has_host = df['host_genome'].notna().sum() if 'host_genome' in df.columns else 0
has_checkv = df['checkv_quality'].notna().sum() if 'checkv_quality' in df.columns else 0

# Count quality distribution
quality_dist = {{}}
if 'checkv_quality' in df.columns:
    quality_dist = df['checkv_quality'].value_counts().to_dict()

# Count lifestyle distribution
lifestyle_dist = {{}}
if 'lifestyle' in df.columns:
    lifestyle_dist = df['lifestyle'].value_counts().to_dict()

# Create summary
summary = {{
    'step': 'final_summary',
    'timestamp': datetime.now().isoformat(),
    'inputs': {{
        'final_contig_summary': {{
            'total_contigs': total_contigs,
            'annotation_coverage': {{
                'has_taxonomy': {{'count': has_taxonomy, 'percentage': round(has_taxonomy/total_contigs*100, 1) if total_contigs > 0 else 0}},
                'has_lifestyle': {{'count': has_lifestyle, 'percentage': round(has_lifestyle/total_contigs*100, 1) if total_contigs > 0 else 0}},
                'has_host_prediction': {{'count': has_host, 'percentage': round(has_host/total_contigs*100, 1) if total_contigs > 0 else 0}},
                'has_quality_assessment': {{'count': has_checkv, 'percentage': round(has_checkv/total_contigs*100, 1) if total_contigs > 0 else 0}}
            }},
            'quality_distribution': quality_dist,
            'lifestyle_distribution': lifestyle_dist
        }}
    }}
}}

# Write to JSON
with open('{output.summary}', 'w') as f:
    json.dump(summary, f, indent=2)

print(f'Final summary stats: {{total_contigs}} contigs')
print(f'Annotation coverage - Taxonomy: {{has_taxonomy}}, Lifestyle: {{has_lifestyle}}, Hosts: {{has_host}}')
" >> {log} 2>&1

        echo "Final summary table statistics collected" >> {log}
        """

# Rule to generate the final summary report
rule generate_summary_report:
    input:
        # Collect all available summaries
        summaries = [
            f"{SUMMARY_DIR}/input_stats.json",
            f"{SUMMARY_DIR}/filtering_stats.json",
            f"{SUMMARY_DIR}/jaeger_stats.json",
            f"{SUMMARY_DIR}/genomad_stats.json",
            f"{SUMMARY_DIR}/phold_stats.json",
            f"{SUMMARY_DIR}/checkv_stats.json",
            f"{SUMMARY_DIR}/integration_stats.json",
            f"{SUMMARY_DIR}/iphop_stats.json",
            f"{SUMMARY_DIR}/lifestyle_stats.json",
            f"{SUMMARY_DIR}/final_summary.json",
            f"{SUMMARY_DIR}/final_phages.json"
        ],
        # Optional summaries (may not exist)
        reneo_summary = f"{SUMMARY_DIR}/reneo_stats.json",
        consensus_summary = f"{SUMMARY_DIR}/consensus_taxonomy.json",
        clustering_summary = f"{SUMMARY_DIR}/clustering_stats.json"
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
        f"{SUMMARY_DIR}/input_stats.json",
        f"{SUMMARY_DIR}/reneo_stats.json",
        f"{SUMMARY_DIR}/filtering_stats.json",
        f"{SUMMARY_DIR}/jaeger_stats.json",
        f"{SUMMARY_DIR}/genomad_stats.json",
        f"{SUMMARY_DIR}/phold_stats.json",
        f"{SUMMARY_DIR}/checkv_stats.json",
        f"{SUMMARY_DIR}/integration_stats.json",
        f"{SUMMARY_DIR}/iphop_stats.json",
        f"{SUMMARY_DIR}/lifestyle_stats.json",
        f"{SUMMARY_DIR}/consensus_taxonomy.json",
        f"{SUMMARY_DIR}/clustering_stats.json",
        f"{SUMMARY_DIR}/final_summary.json",
        f"{SUMMARY_DIR}/final_phages.json"
    output:
        flag = f"{SUMMARY_DIR}/.all_summaries_collected"
    shell:
        "touch {output.flag}"