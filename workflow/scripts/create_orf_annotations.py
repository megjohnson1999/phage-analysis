#!/usr/bin/env python3
"""
Create ORF Annotations Table

Parses Prodigal output to create a table linking ORFs to their parent contigs
and providing positional information.

This enables linking ORF abundance data back to contig-level annotations.
"""

import argparse
import re
import sys
from pathlib import Path


def parse_prodigal_header(header):
    """
    Parse Prodigal FASTA header format.

    Example header:
    >contig_1_1 # 100 # 500 # 1 # ID=1_1;partial=00;start_type=ATG

    Returns:
        dict with: orf_id, contig_id, start, end, strand, partial, start_type
    """
    # Remove leading '>'
    header = header.lstrip('>')

    # Split by ' # ' delimiter
    parts = header.split(' # ')

    if len(parts) < 4:
        print(f"Warning: Could not parse header: {header}", file=sys.stderr)
        return None

    orf_id = parts[0].strip()
    start = parts[1].strip()
    end = parts[2].strip()
    strand = '+' if parts[3].strip() == '1' else '-'

    # Extract contig ID from ORF ID
    # Prodigal typically names ORFs as: contigName_orfNumber
    # Try to extract the contig name (everything before the last underscore and number)
    contig_id = '_'.join(orf_id.split('_')[:-1])
    if not contig_id:
        contig_id = orf_id  # Fallback to full ORF ID if parsing fails

    # Parse additional info if present
    partial = None
    start_type = None

    if len(parts) >= 5:
        info_str = parts[4]
        # Extract key=value pairs
        if 'partial=' in info_str:
            match = re.search(r'partial=(\d+)', info_str)
            if match:
                partial = match.group(1)

        if 'start_type=' in info_str:
            match = re.search(r'start_type=(\w+)', info_str)
            if match:
                start_type = match.group(1)

    return {
        'orf_id': orf_id,
        'contig_id': contig_id,
        'start': int(start),
        'end': int(end),
        'strand': strand,
        'partial': partial,
        'start_type': start_type
    }


def main():
    parser = argparse.ArgumentParser(
        description='Create ORF annotations from Prodigal output'
    )
    parser.add_argument('--genes', required=True,
                       help='Path to Prodigal genes.fna file')
    parser.add_argument('--proteins', required=True,
                       help='Path to Prodigal proteins.faa file')
    parser.add_argument('--output', required=True,
                       help='Output path for annotations TSV')

    args = parser.parse_args()

    print("=" * 70)
    print("Creating ORF Annotations")
    print("=" * 70)
    print()

    # Parse gene file
    print(f"Parsing ORF sequences from: {args.genes}")
    annotations = []

    if not Path(args.genes).exists():
        print(f"Error: Gene file not found: {args.genes}")
        sys.exit(1)

    with open(args.genes, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parsed = parse_prodigal_header(line.strip())
                if parsed:
                    # Calculate ORF length
                    orf_length = abs(parsed['end'] - parsed['start']) + 1
                    parsed['orf_length'] = orf_length
                    annotations.append(parsed)

    print(f"  Parsed {len(annotations)} ORFs")

    # Count contigs
    unique_contigs = set(ann['contig_id'] for ann in annotations)
    print(f"  From {len(unique_contigs)} contigs")
    print()

    # Write output
    print(f"Writing annotations to: {args.output}")

    with open(args.output, 'w') as f:
        # Write header
        f.write('\t'.join([
            'orf_id', 'contig_id', 'start', 'end', 'strand',
            'orf_length', 'partial', 'start_type'
        ]) + '\n')

        # Write data
        for ann in annotations:
            f.write('\t'.join([
                str(ann['orf_id']),
                str(ann['contig_id']),
                str(ann['start']),
                str(ann['end']),
                str(ann['strand']),
                str(ann['orf_length']),
                str(ann.get('partial', '')),
                str(ann.get('start_type', ''))
            ]) + '\n')

    print(f"  Wrote {len(annotations)} ORF annotations")
    print()

    # Print summary statistics
    print("Summary Statistics:")
    print(f"  Total ORFs: {len(annotations)}")
    print(f"  Total contigs: {len(unique_contigs)}")
    print(f"  ORFs per contig (avg): {len(annotations) / len(unique_contigs):.1f}")

    # Strand distribution
    plus_strand = sum(1 for ann in annotations if ann['strand'] == '+')
    minus_strand = sum(1 for ann in annotations if ann['strand'] == '-')
    print(f"  Plus strand: {plus_strand} ({100*plus_strand/len(annotations):.1f}%)")
    print(f"  Minus strand: {minus_strand} ({100*minus_strand/len(annotations):.1f}%)")

    # Length statistics
    lengths = [ann['orf_length'] for ann in annotations]
    print(f"  ORF length (mean): {sum(lengths) / len(lengths):.0f} bp")
    print(f"  ORF length (median): {sorted(lengths)[len(lengths)//2]:.0f} bp")
    print(f"  ORF length (min): {min(lengths)} bp")
    print(f"  ORF length (max): {max(lengths)} bp")

    print()
    print("=" * 70)
    print("ORF annotations created successfully!")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
