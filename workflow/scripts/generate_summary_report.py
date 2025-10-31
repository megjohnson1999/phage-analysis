#!/usr/bin/env python3
"""
Generate a comprehensive HTML summary report from pipeline step summaries.
"""

import os
import json
import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd

def load_summaries(summary_dir):
    """Load all summary JSON files from the summary directory."""
    summaries = {}
    
    if not os.path.exists(summary_dir):
        return summaries
    
    for json_file in Path(summary_dir).glob("*.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                step_name = data.get('step', json_file.stem)
                summaries[step_name] = data
        except Exception as e:
            print(f"Warning: Could not load {json_file}: {e}")
    
    return summaries

def format_number(num):
    """Format numbers with commas for readability."""
    if isinstance(num, (int, float)):
        return f"{num:,}"
    elif isinstance(num, list) and len(num) == 1:
        # Handle single-element lists (from consensus taxonomy data)
        return format_number(num[0])
    elif isinstance(num, list):
        # Handle multi-element lists
        return ", ".join([format_number(n) for n in num])
    return str(num)

def generate_html_report(summaries, output_file, config_info=None):
    """Generate an HTML summary report."""
    
    # HTML template
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Phage Analysis Pipeline Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
            .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
            h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; margin-top: 30px; }}
            h3 {{ color: #2980b9; }}
            .summary-box {{ background-color: #ecf0f1; padding: 20px; border-radius: 5px; margin: 15px 0; }}
            .metric {{ display: inline-block; margin: 10px 20px 10px 0; }}
            .metric-value {{ font-size: 1.5em; font-weight: bold; color: #2980b9; display: block; }}
            .metric-label {{ font-size: 0.9em; color: #7f8c8d; }}
            table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
            th, td {{ padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #f8f9fa; font-weight: bold; }}
            .success {{ color: #27ae60; }}
            .warning {{ color: #f39c12; }}
            .error {{ color: #e74c3c; }}
            .step-section {{ margin: 20px 0; padding: 15px; border-left: 4px solid #3498db; background-color: #f8f9fa; }}
            .timestamp {{ color: #7f8c8d; font-size: 0.9em; }}
            .progress-bar {{ width: 100%; background-color: #ecf0f1; border-radius: 10px; overflow: hidden; }}
            .progress-fill {{ height: 20px; background-color: #3498db; text-align: center; line-height: 20px; color: white; font-size: 0.8em; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Phage Analysis Pipeline Report</h1>
            <div class="timestamp">Generated: {timestamp}</div>
            
            {config_section}
            
            <h2>Overall Summary</h2>
            <div class="summary-box">
                {overall_metrics}
            </div>
            
            <h2>Pipeline Progress</h2>
            {pipeline_steps}
            
            <h2>Detailed Results</h2>
            {detailed_sections}
            
            <h2>Notes</h2>
            <ul>
                <li>This report summarizes the key statistics from your phage analysis pipeline run.</li>
                <li>All sequence counts and metrics are based on the final outputs of each tool.</li>
                <li>Quality metrics from CheckV help assess the completeness of identified phage genomes.</li>
                <li>Host predictions and lifestyle classifications provide biological context for the phages.</li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    # Generate timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Generate config section
    config_section = ""
    if config_info:
        config_section = f"""
        <h2>Configuration</h2>
        <table>
            {config_info}
        </table>
        """
    
    # Calculate overall metrics
    overall_metrics = ""
    
    # Try to extract key numbers
    input_seqs = 0
    final_phages = 0
    viral_contigs = 0
    
    if "input_stats" in summaries:
        input_data = summaries["input_stats"].get("inputs", {})
        for key, value in input_data.items():
            if isinstance(value, dict) and "total_sequences" in value:
                input_seqs = max(input_seqs, value["total_sequences"])
    
    if "final_phages" in summaries:
        final_data = summaries["final_phages"].get("inputs", {})
        for key, value in final_data.items():
            if isinstance(value, dict) and "total_sequences" in value:
                final_phages = value["total_sequences"]
    
    if "filtering_stats" in summaries:
        filter_data = summaries["filtering_stats"].get("inputs", {})
        for key, value in filter_data.items():
            if isinstance(value, dict) and "total_sequences" in value:
                viral_contigs = value["total_sequences"]
    
    overall_metrics = f"""
        <div class="metric">
            <span class="metric-value">{format_number(input_seqs)}</span>
            <span class="metric-label">Input Sequences</span>
        </div>
        <div class="metric">
            <span class="metric-value">{format_number(viral_contigs)}</span>
            <span class="metric-label">Viral Contigs</span>
        </div>
        <div class="metric">
            <span class="metric-value">{format_number(final_phages)}</span>
            <span class="metric-label">Final Phages</span>
        </div>
    """
    
    # Generate pipeline steps section with key metrics
    def get_step_metric(step_name, summaries):
        """Extract key metric for each pipeline step."""
        if step_name not in summaries:
            return "Not completed"

        data = summaries[step_name].get("inputs", {})

        if step_name == "input_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} input sequences"

        elif step_name == "reneo_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} contigs"

        elif step_name == "filtering_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} viral contigs"

        elif step_name == "jaeger_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "phage_predictions" in value:
                    return f"{format_number(value['phage_predictions'])} phage predictions"

        elif step_name == "genomad_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_predictions" in value:
                    return f"{format_number(value['total_predictions'])} predictions"

        elif step_name == "phold_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_annotations" in value:
                    return f"{format_number(value['total_annotations'])} annotations"

        elif step_name == "checkv_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} sequences analyzed"

        elif step_name == "integration_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value and value["total_sequences"] > 0:
                    return f"{format_number(value['total_sequences'])} viral contigs"

        elif step_name == "iphop_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_predictions" in value:
                    return f"{format_number(value['total_predictions'])} host predictions"

        elif step_name == "lifestyle_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_predictions" in value:
                    total = value.get('total_predictions', 0)
                    # Try to get source breakdown (BACPHLIP vs Phabox2)
                    if "by_source" in value:
                        by_source = value["by_source"]
                        bacphlip = by_source.get("BACPHLIP", 0)
                        phabox2 = by_source.get("Phabox2", 0)
                        return f"{format_number(total)} lifestyle predictions (BACPHLIP: {bacphlip}, Phabox2: {phabox2})"
                    return f"{format_number(total)} lifestyle predictions"

        elif step_name == "consensus_taxonomy":
            for key, value in data.items():
                if isinstance(value, dict) and "total_contigs" in value:
                    return f"{format_number(value['total_contigs'])} contigs classified"

        elif step_name == "clustering_stats":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} clusters"

        elif step_name == "final_summary":
            for key, value in data.items():
                if isinstance(value, dict) and "total_contigs" in value:
                    total = value.get('total_contigs', 0)
                    # Try to get annotation coverage
                    if "annotation_coverage" in value:
                        coverage = value["annotation_coverage"]
                        # Get taxonomy percentage from nested structure
                        tax_pct = coverage.get("has_taxonomy", {}).get("percentage", 0)
                        return f"{format_number(total)} contigs ({tax_pct:.1f}% with taxonomy)"
                    return f"{format_number(total)} contigs in final summary"

        elif step_name == "final_phages":
            for key, value in data.items():
                if isinstance(value, dict) and "total_sequences" in value:
                    return f"{format_number(value['total_sequences'])} final phages"

        return "Completed"

    step_order = [
        "input_stats", "reneo_stats", "filtering_stats", "jaeger_stats",
        "genomad_stats", "phold_stats", "checkv_stats", "integration_stats",
        "iphop_stats", "lifestyle_stats", "consensus_taxonomy", "clustering_stats",
        "final_summary", "final_phages"
    ]

    pipeline_steps = ""
    completed_steps = 0
    total_steps = len([s for s in step_order if s in summaries])

    for step in step_order:
        if step in summaries:
            completed_steps += 1
            metric = get_step_metric(step, summaries)
            step_display = step.replace("_", " ").title()
            pipeline_steps += f'<div class="step-section">[DONE] {step_display}: {metric}</div>'
        else:
            step_display = step.replace("_", " ").title()
            pipeline_steps += f'<div class="step-section">- {step_display}: Not completed</div>'

    if total_steps > 0:
        progress_percent = int((completed_steps / total_steps) * 100)
        pipeline_steps = f"""
        <div class="progress-bar">
            <div class="progress-fill" style="width: {progress_percent}%">{progress_percent}% Complete</div>
        </div>
        {pipeline_steps}
        """
    
    # Generate detailed sections
    detailed_sections = ""
    
    for step_name, data in summaries.items():
        detailed_sections += f"""
        <div class="step-section">
            <h3>{step_name.replace('_', ' ').title()}</h3>
            <div class="timestamp">Completed: {data.get('timestamp', 'Unknown')}</div>
        """
        
        # Add step-specific details
        inputs = data.get("inputs", {})
        for input_name, input_data in inputs.items():
            # Skip reads directory information - not useful for reporting
            if input_name == 'reads_dir':
                continue

            if isinstance(input_data, dict):
                # For integration stats, skip provirus predictions (always 0 in Reneo workflow) and only show viral contigs
                if step_name == "integration_stats" and input_name == "phage_predictions":
                    continue  # Skip provirus predictions - not relevant for viral-only datasets

                # Use clearer labels
                display_name = input_name.replace('_', ' ').title()
                if step_name == "integration_stats" and input_name == "phage_contigs":
                    display_name = "Viral Contigs"

                detailed_sections += f"<h4>{display_name}</h4>"
                detailed_sections += "<table>"

                for key, value in input_data.items():
                    # Skip read-related fields if they somehow get through
                    if key in ['total_reads', 'read_files', 'reads_directory', 'method', 'note'] and 'read' in key.lower():
                        continue

                    if isinstance(value, dict):
                        # Handle nested dictionaries (like completeness_distribution, taxonomy_coverage)
                        detailed_sections += f"<tr><td>{key.replace('_', ' ').title()}</td><td>"
                        for sub_key, sub_value in value.items():
                            # Special handling for taxonomy coverage format
                            if isinstance(sub_value, dict) and 'count' in sub_value and 'percentage' in sub_value:
                                count = format_number(sub_value['count'])
                                percentage = format_number(sub_value['percentage'])
                                detailed_sections += f"{sub_key}: {count} ({percentage}%)<br>"
                            else:
                                detailed_sections += f"{sub_key}: {format_number(sub_value)}<br>"
                        detailed_sections += "</td></tr>"
                    else:
                        detailed_sections += f"<tr><td>{key.replace('_', ' ').title()}</td><td>{format_number(value)}</td></tr>"

                detailed_sections += "</table>"
        
        detailed_sections += "</div>"
    
    # Fill in the template
    html_content = html_template.format(
        timestamp=timestamp,
        config_section=config_section,
        overall_metrics=overall_metrics,
        pipeline_steps=pipeline_steps,
        detailed_sections=detailed_sections
    )
    
    # Write the HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate HTML summary report from pipeline step summaries')
    parser.add_argument('--summary-dir', required=True, help='Directory containing summary JSON files')
    parser.add_argument('--output', required=True, help='Output HTML file')
    parser.add_argument('--config', help='Configuration information to include in report')
    
    args = parser.parse_args()
    
    # Load all summaries
    summaries = load_summaries(args.summary_dir)
    
    if not summaries:
        print("Warning: No summary files found. Creating minimal report.")
        summaries = {}
    
    # Generate the report
    generate_html_report(summaries, args.output, args.config)

if __name__ == "__main__":
    main()