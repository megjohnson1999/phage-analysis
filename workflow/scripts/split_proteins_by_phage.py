#!/usr/bin/env python
"""
Split protein FASTA files by phage ID.

This script parses protein FASTA files from Prodigal and splits them
into separate files based on the phage ID extracted from the sequence headers.

Usage:
    python split_proteins_by_phage.py --input proteins.faa --output-dir split_proteins/
"""

import os
import re
import argparse
from Bio import SeqIO

def extract_phage_id(header):
    """
    Extract phage ID from a FASTA header.
    
    Args:
        header: The FASTA header string (with or without leading '>')
        
    Returns:
        str: The extracted phage ID
    """
    # Remove leading '>' if present
    header = header.lstrip('>')
    
    # Get the part before any space or # (typical Prodigal header format)
    # Example: "contig_1_123 # start=456 end=789" -> "contig_1_123"
    base_header = header.split()[0]
    
    # Extract the contig ID by removing the last numeric suffix
    # This keeps the full contig name intact (e.g., "disjointig_1_123" -> "disjointig_1")
    match = re.match(r'^(.+)_\d+$', base_header)
    if match:
        return match.group(1)
    
    # Fallback: use the entire header up to the first space/hash
    return base_header

def split_proteins(input_file, output_dir):
    """
    Split protein FASTA file by phage ID.
    
    Args:
        input_file: Path to input FASTA file
        output_dir: Directory to write split files
        
    Returns:
        list: Paths to created output files
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Group proteins by phage ID
    phage_proteins = {}
    
    # Parse the input file
    total_records = 0
    for record in SeqIO.parse(input_file, "fasta"):
        phage_id = extract_phage_id(record.id)
        total_records += 1
        
        if phage_id not in phage_proteins:
            phage_proteins[phage_id] = []
        
        phage_proteins[phage_id].append(record)
    
    # Write a file for each phage
    output_files = []
    for phage_id, records in sorted(phage_proteins.items()):
        if phage_id == "unknown":
            continue  # Skip records where we couldn't determine phage ID
            
        output_file = os.path.join(output_dir, f"{phage_id}.faa")
        SeqIO.write(records, output_file, "fasta")
        output_files.append(output_file)
    
    # Print summary stats
    print(f"Processed {total_records} total protein sequences")
    print(f"Created {len(output_files)} phage-specific files")
    if 1 <= len(phage_proteins) <= 2:
        print(f"WARNING: Only {len(phage_proteins)} phage IDs detected: {', '.join(phage_proteins.keys())}")
        print("This may indicate an issue with contig naming or header parsing.")
    
    return output_files

def main():
    """Main function to parse arguments and run the script."""
    parser = argparse.ArgumentParser(description="Split protein FASTA files by phage ID")
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory for split files")
    parser.add_argument("--split-list", help="Optional file to write list of split files")
    parser.add_argument("--show-examples", action="store_true", help="Show example phage IDs extracted")
    
    args = parser.parse_args()
    
    # Process the file
    output_files = split_proteins(args.input, args.output_dir)
    
    # Print summary
    print(f"Split {args.input} into {len(output_files)} phage-specific files in {args.output_dir}")
    
    # Write the list of split files if requested
    if args.split_list:
        with open(args.split_list, 'w') as f:
            for output_file in sorted(output_files):
                f.write(f"{output_file}\n")
        print(f"Wrote list of split files to {args.split_list}")
    
    # Show examples if requested
    if args.show_examples and output_files:
        print("\nExample phage IDs extracted:")
        
        # Show the first few examples
        for i, output_file in enumerate(sorted(output_files)[:5]):
            phage_id = os.path.basename(output_file).replace('.faa', '')
            count = sum(1 for _ in SeqIO.parse(output_file, "fasta"))
            print(f"  {i+1}. {phage_id} - {count} proteins")
        
        if len(output_files) > 5:
            print(f"  ... and {len(output_files) - 5} more")

if __name__ == "__main__":
    main()