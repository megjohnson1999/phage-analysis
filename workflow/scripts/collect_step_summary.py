#!/usr/bin/env python3
"""
Collect statistics at each major pipeline step.
This script is designed to be lightweight and not impact pipeline performance.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from datetime import datetime

def count_sequences(fasta_file):
    """Count sequences and calculate basic statistics from a FASTA file."""
    if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
        return {
            "total_sequences": 0,
            "total_length": 0,
            "mean_length": 0,
            "min_length": 0,
            "max_length": 0,
            "n50": 0
        }
    
    lengths = []
    total_length = 0
    
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            length = len(record.seq)
            lengths.append(length)
            total_length += length
    
    if not lengths:
        return {
            "total_sequences": 0,
            "total_length": 0,
            "mean_length": 0,
            "min_length": 0,
            "max_length": 0,
            "n50": 0
        }
    
    # Calculate N50
    lengths.sort(reverse=True)
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= total_length * 0.5:
            n50 = length
            break
    
    return {
        "total_sequences": len(lengths),
        "total_length": total_length,
        "mean_length": int(total_length / len(lengths)),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "n50": n50
    }

def count_reads(reads_dir):
    """Count reads in a directory."""
    if not os.path.exists(reads_dir):
        return {"total_reads": 0, "read_files": 0}
    
    read_files = []
    total_reads = 0
    
    # Look for common read file extensions
    extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    
    for ext in extensions:
        read_files.extend(Path(reads_dir).glob(f"*{ext}"))
    
    # Count reads in each file (simplified - just count lines/4 for fastq)
    for read_file in read_files:
        try:
            if str(read_file).endswith('.gz'):
                import gzip
                with gzip.open(read_file, 'rt') as f:
                    lines = sum(1 for _ in f)
            else:
                with open(read_file, 'r') as f:
                    lines = sum(1 for _ in f)
            total_reads += lines // 4  # FASTQ has 4 lines per read
        except:
            continue  # Skip problematic files
    
    return {
        "total_reads": total_reads,
        "read_files": len(read_files)
    }

def parse_tool_results(tool, result_file):
    """Parse results from various prediction tools."""
    if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
        return {"total_predictions": 0, "tool": tool}
    
    try:
        if tool == "jaeger":
            df = pd.read_csv(result_file, sep='\t')
            phage_predictions = df[df['prediction'] == 'Phage'] if 'prediction' in df.columns else pd.DataFrame()
            return {
                "tool": tool,
                "total_predictions": len(df),
                "phage_predictions": len(phage_predictions),
                "mean_score": phage_predictions['realiability_score'].mean() if 'realiability_score' in phage_predictions.columns and len(phage_predictions) > 0 else 0
            }
        
        elif tool == "genomad":
            df = pd.read_csv(result_file, sep='\t')
            return {
                "tool": tool,
                "total_predictions": len(df),
                "mean_virus_score": df['virus_score'].mean() if 'virus_score' in df.columns else 0,
                "topology_counts": df['topology'].value_counts().to_dict() if 'topology' in df.columns else {}
            }
        
        elif tool == "checkv":
            df = pd.read_csv(result_file, sep='\t')
            return {
                "tool": tool,
                "total_sequences": len(df),
                "completeness_distribution": {
                    "complete": len(df[df['completeness'] == 'Complete']) if 'completeness' in df.columns else 0,
                    "high_quality": len(df[df['completeness'] == 'High-quality']) if 'completeness' in df.columns else 0,
                    "medium_quality": len(df[df['completeness'] == 'Medium-quality']) if 'completeness' in df.columns else 0,
                    "low_quality": len(df[df['completeness'] == 'Low-quality']) if 'completeness' in df.columns else 0
                }
            }
        
        elif tool == "phold":
            df = pd.read_csv(result_file, sep='\t')
            categories = df['category'].value_counts().to_dict() if 'category' in df.columns else {}
            return {
                "tool": tool,
                "total_annotations": len(df),
                "functional_categories": categories,
                "unique_contigs": df['contig_id'].nunique() if 'contig_id' in df.columns else 0
            }
        
        elif tool == "iphop":
            df = pd.read_csv(result_file, sep='\t')
            return {
                "tool": tool,
                "total_predictions": len(df),
                "unique_hosts": df['genus'].nunique() if 'genus' in df.columns else 0,
                "mean_confidence": df['score'].mean() if 'score' in df.columns else 0
            }
            
        elif tool == "mmseqs":
            df = pd.read_csv(result_file, sep='\t')
            return {
                "tool": tool,
                "total_assignments": len(df),
                "viral_assignments": len(df[df['taxlineage'].str.contains('Viruses', na=False)]) if 'taxlineage' in df.columns else 0
            }
    
    except Exception as e:
        return {
            "tool": tool,
            "total_predictions": 0,
            "error": str(e)
        }
    
    return {"tool": tool, "total_predictions": 0}

def collect_step_summary(step_name, inputs, output_file):
    """Collect summary statistics for a pipeline step."""
    
    timestamp = datetime.now().isoformat()
    summary = {
        "step": step_name,
        "timestamp": timestamp,
        "inputs": {},
        "outputs": {},
        "statistics": {}
    }
    
    # Process different types of inputs based on step
    for input_name, input_path in inputs.items():
        if input_name.endswith('_fasta') or input_name.endswith('_contigs'):
            summary["inputs"][input_name] = count_sequences(input_path)
        elif input_name == 'reads_dir':
            summary["inputs"][input_name] = count_reads(input_path)
        elif input_name.endswith('_results'):
            # Determine tool from step name or file path
            tool = step_name.split('_')[0] if '_' in step_name else "unknown"
            summary["inputs"][input_name] = parse_tool_results(tool, input_path)
    
    # Save summary
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary for step '{step_name}' saved to {output_file}")
    return summary

def main():
    parser = argparse.ArgumentParser(description='Collect pipeline step summary statistics')
    parser.add_argument('--step', required=True, help='Step name')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--inputs', nargs='*', help='Input files in format name:path')
    
    args = parser.parse_args()
    
    # Parse inputs
    inputs = {}
    if args.inputs:
        for inp in args.inputs:
            if ':' in inp:
                name, path = inp.split(':', 1)
                inputs[name] = path
    
    collect_step_summary(args.step, inputs, args.output)

if __name__ == "__main__":
    main()