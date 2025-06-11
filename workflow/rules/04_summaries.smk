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
        assembly = config.get("assembly_file", ""),
        reads_dir = config.get("reads_dir", "")
    log:
        f"{config['output_dir']}/logs/summaries/collect_input_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step input_stats \
            --output {output.summary} \
            --inputs assembly_fasta:{params.assembly} reads_dir:{params.reads_dir} \
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
        """

# Rule to collect lifestyle prediction statistics
rule collect_lifestyle_stats:
    input:
        bacphlip_results = f"{config['output_dir']}/03_genomic_info/bacphlip_lifestyle_with_completeness.tsv"
    output:
        summary = f"{SUMMARY_DIR}/lifestyle_stats.json"
    log:
        f"{config['output_dir']}/logs/summaries/collect_lifestyle_stats.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Check if BACPHLIP results exist, otherwise try to find other lifestyle prediction results
        LIFESTYLE_FILE="{input.bacphlip_results}"
        if [ ! -f "$LIFESTYLE_FILE" ]; then
            # Look for alternative lifestyle prediction files
            for alt_file in "{config[output_dir]}/03_phacts_results/phacts_predictions_compiled.tsv" "{config[output_dir]}/03_genomic_info/phabox_output/lifestyle.tsv"; do
                if [ -f "$alt_file" ]; then
                    LIFESTYLE_FILE="$alt_file"
                    break
                fi
            done
        fi
        
        python {workflow.basedir}/scripts/collect_step_summary.py \
            --step lifestyle_stats \
            --output {output.summary} \
            --inputs lifestyle_results:"$LIFESTYLE_FILE" \
            > {log} 2>&1 || true
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
            > {log} 2>&1 || true
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
            f"{SUMMARY_DIR}/final_phages.json"
        ],
        # Optional summaries (may not exist)
        reneo_summary = f"{SUMMARY_DIR}/reneo_stats.json",
        clustering_summary = f"{SUMMARY_DIR}/clustering_stats.json"
    output:
        report = f"{config['output_dir']}/Pipeline_Summary_Report.html"
    log:
        f"{config['output_dir']}/logs/summaries/generate_summary_report.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Create config info for the report
        CONFIG_INFO="<tr><td>Output Directory</td><td>{config[output_dir]}</td></tr>"
        if [ -n "{config.get('assembly_file', '')}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly File</td><td>{config[assembly_file]}</td></tr>"
        fi
        if [ -n "{config.get('assembly_graph', '')}" ]; then
            CONFIG_INFO="$CONFIG_INFO<tr><td>Assembly Graph</td><td>{config[assembly_graph]}</td></tr>"
        fi
        CONFIG_INFO="$CONFIG_INFO<tr><td>Reads Directory</td><td>{config[reads_dir]}</td></tr>"
        CONFIG_INFO="$CONFIG_INFO<tr><td>Clustering Enabled</td><td>{config.get('do_clustering', True)}</td></tr>"
        
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
        f"{SUMMARY_DIR}/clustering_stats.json",
        f"{SUMMARY_DIR}/final_phages.json"
    output:
        flag = f"{SUMMARY_DIR}/.all_summaries_collected"
    shell:
        "touch {output.flag}"