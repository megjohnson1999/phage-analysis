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

def extract_phage_id(header, use_full_contig_id=True):
    """
    Extract phage ID from a FASTA header.
    
    Args:
        header: The FASTA header string (with or without leading '>')
        use_full_contig_id: If True, uses complete contig ID as phage ID.
                           If False, tries to extract a more general phage ID.
        
    Returns:
        str: The extracted phage ID
    """
    # Remove leading '>' if present
    clean_header = header.lstrip('>')
    
    # Get first part before any spaces/tabs
    base_header = clean_header.split()[0]
    
    if use_full_contig_id:
        # For most assemblers, the full contig ID before any spaces is best
        # This preserves unique identities for each contig
        return base_header
    else:
        # Option to extract only part of the contig ID by various patterns
        # This is the original behavior and may be useful in special cases
        
        # Try to handle SPAdes format (NODE_X_length_Y_cov_Z)
        spades_match = re.match(r'NODE_(\d+)_', base_header)
        if spades_match:
            return f"NODE_{spades_match.group(1)}"
            
        # Try to handle MEGAHIT format (kX_Y)
        megahit_match = re.match(r'(k\d+_\d+)', base_header)
        if megahit_match:
            return megahit_match.group(1)
            
        # Original pattern: Extract the part before the last underscore and number
        # For headers like "read_13355_7", this will extract "read_13355"
        original_match = re.match(r'([^_]+(?:_[^_]+)*?)_\d+$', base_header)
        if original_match:
            return original_match.group(1)
        
        # Fallback: just use the entire header before any space/hash
        return base_header
    
    # If all else fails
    return "unknown"

def split_proteins(input_file, output_dir, use_full_contig_id=True):
    """
    Split protein FASTA file by phage ID.
    
    Args:
        input_file: Path to input FASTA file
        output_dir: Directory to write split files
        use_full_contig_id: Whether to use the complete contig ID (True) or 
                           try to extract a more general phage ID (False)
        
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
        phage_id = extract_phage_id(record.id, use_full_contig_id)
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
    if len(phage_proteins) == 1:
        only_key = next(iter(phage_proteins))
        print(f"WARNING: All proteins grouped under a single phage ID: {only_key}")
        print("This may indicate an issue with contig naming or header parsing.")
    elif len(phage_proteins) == 2:
        keys = list(phage_proteins.keys())
        print(f"WARNING: Only 2 phage IDs detected: {keys[0]} and {keys[1]}")
        print("This may indicate an issue with contig naming or header parsing.")
    
    return output_files

def main():
    """Main function to parse arguments and run the script."""
    parser = argparse.ArgumentParser(description="Split protein FASTA files by phage ID")
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory for split files")
    parser.add_argument("--split-list", help="Optional file to write list of split files")
    parser.add_argument("--show-examples", action="store_true", help="Show example phage IDs extracted")
    parser.add_argument("--use-full-contig-id", action="store_true", 
                        help="Use the full contig ID as the phage ID (recommended for most assemblies)")
    
    args = parser.parse_args()
    
    # Process the file
    output_files = split_proteins(args.input, args.output_dir, args.use_full_contig_id)
    
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