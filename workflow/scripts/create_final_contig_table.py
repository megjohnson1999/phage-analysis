#!/usr/bin/env python3
"""
Create Final Contig Summary Table

Consolidates key information about each contig into a single summary table.
Includes phage prediction details, CheckV quality, taxonomy, lifestyle, and host predictions.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def load_optional_table(file_path, description):
    """
    Load a table if it exists, otherwise return None.
    """
    if not file_path or not Path(file_path).exists():
        print(f"  {description}: Not available")
        return None

    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"  {description}: Loaded {len(df)} rows")
        return df
    except Exception as e:
        print(f"  Warning: Error loading {description}: {e}")
        return None


def get_contig_list(phage_seqs_file):
    """
    Extract list of contig IDs from the FASTA file.
    """
    contigs = []
    try:
        with open(phage_seqs_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract contig ID (everything after > until first space)
                    contig_id = line[1:].split()[0].strip()
                    contigs.append(contig_id)
        print(f"  Contig list: Found {len(contigs)} sequences")
        return contigs
    except Exception as e:
        print(f"  Error reading sequences file: {e}")
        return []


def main():
    parser = argparse.ArgumentParser(
        description='Create comprehensive contig summary table'
    )
    parser.add_argument('--phage-seqs', required=True,
                       help='Path to final phage sequences FASTA file')
    parser.add_argument('--phage-predictions', required=False,
                       help='Path to phage predictions TSV (optional)')
    parser.add_argument('--checkv-results', required=False,
                       help='Path to CheckV quality summary TSV (optional)')
    parser.add_argument('--consensus-taxonomy', required=False,
                       help='Path to consensus taxonomy TSV (optional)')
    parser.add_argument('--lifestyle-consensus', required=False,
                       help='Path to lifestyle consensus TSV (optional)')
    parser.add_argument('--iphop-predictions', required=False,
                       help='Path to iPhop host predictions TSV (optional)')
    parser.add_argument('--output', required=True,
                       help='Output path for summary table (TSV)')

    args = parser.parse_args()

    print("=" * 70)
    print("Creating Final Contig Summary Table")
    print("=" * 70)
    print()

    # Get list of all contigs from the sequences file
    print("Loading contig list...")
    contigs = get_contig_list(args.phage_seqs)

    if not contigs:
        print("Error: No contigs found in sequences file")
        sys.exit(1)

    # Initialize summary DataFrame with contig IDs
    summary = pd.DataFrame({'contig_id': contigs})
    print()

    # Load optional data files
    print("Loading annotation data...")

    # 1. Phage prediction rules (why was it called a phage?)
    phage_pred = load_optional_table(args.phage_predictions, "Phage predictions")
    if phage_pred is not None and 'contig_id' in phage_pred.columns:
        # Keep relevant columns: prediction_rule, genomad_score, jaeger_score, phold_category
        pred_cols = ['contig_id', 'prediction_rule']
        if 'genomad_score' in phage_pred.columns:
            pred_cols.append('genomad_score')
        if 'jaeger_score' in phage_pred.columns:
            pred_cols.append('jaeger_score')
        if 'phold_category' in phage_pred.columns:
            pred_cols.append('phold_category')

        phage_pred_subset = phage_pred[pred_cols]
        summary = summary.merge(phage_pred_subset, on='contig_id', how='left')

    # 2. CheckV quality assessment
    checkv = load_optional_table(args.checkv_results, "CheckV quality")
    if checkv is not None:
        # Rename first column to contig_id if needed
        if checkv.columns[0] != 'contig_id':
            checkv = checkv.rename(columns={checkv.columns[0]: 'contig_id'})

        # Keep key CheckV columns
        checkv_cols = ['contig_id']
        optional_cols = ['checkv_quality', 'completeness', 'contamination',
                        'provirus', 'contig_length', 'proviral_length']
        for col in optional_cols:
            if col in checkv.columns:
                checkv_cols.append(col)

        checkv_subset = checkv[checkv_cols]
        summary = summary.merge(checkv_subset, on='contig_id', how='left')

    # 3. Consensus taxonomy
    taxonomy = load_optional_table(args.consensus_taxonomy, "Consensus taxonomy")
    if taxonomy is not None and 'contig_id' in taxonomy.columns:
        # Keep taxonomy and confidence columns
        tax_cols = ['contig_id', 'consensus_taxonomy', 'confidence', 'agreement_level']
        tax_subset = taxonomy[[col for col in tax_cols if col in taxonomy.columns]]
        # Rename to avoid conflicts
        tax_subset = tax_subset.rename(columns={
            'consensus_taxonomy': 'taxonomy',
            'confidence': 'taxonomy_confidence',
            'agreement_level': 'taxonomy_agreement'
        })
        summary = summary.merge(tax_subset, on='contig_id', how='left')

    # 4. Lifestyle consensus
    lifestyle = load_optional_table(args.lifestyle_consensus, "Lifestyle consensus")
    if lifestyle is not None and 'contig_id' in lifestyle.columns:
        # Keep lifestyle, confidence, and source
        life_cols = ['contig_id', 'lifestyle', 'confidence', 'source']
        life_subset = lifestyle[[col for col in life_cols if col in lifestyle.columns]]
        # Rename to be more specific
        life_subset = life_subset.rename(columns={
            'confidence': 'lifestyle_confidence',
            'source': 'lifestyle_source'
        })
        summary = summary.merge(life_subset, on='contig_id', how='left')

    # 5. iPhop host predictions (take top hit)
    iphop = load_optional_table(args.iphop_predictions, "iPhop host predictions")
    if iphop is not None:
        # Rename query column to contig_id if needed
        if 'query' in iphop.columns:
            iphop = iphop.rename(columns={'query': 'contig_id'})
        elif 'Virus' in iphop.columns:
            iphop = iphop.rename(columns={'Virus': 'contig_id'})

        if 'contig_id' in iphop.columns:
            # Group by contig and take first (highest confidence) prediction
            iphop_cols = ['contig_id']
            optional_host_cols = ['host', 'score', 'genus', 'family', 'order', 'class', 'phylum']
            for col in optional_host_cols:
                if col in iphop.columns:
                    iphop_cols.append(col)

            # Take first prediction per contig (assumes sorted by confidence)
            iphop_subset = iphop[iphop_cols].groupby('contig_id').first().reset_index()

            # Rename columns to be more specific
            rename_dict = {col: f'host_{col}' for col in iphop_subset.columns if col != 'contig_id'}
            iphop_subset = iphop_subset.rename(columns=rename_dict)

            summary = summary.merge(iphop_subset, on='contig_id', how='left')

    print()

    # Save output
    summary.to_csv(args.output, sep='\t', index=False)
    print(f"Summary table saved: {args.output}")
    print(f"  Total contigs: {len(summary)}")
    print(f"  Total columns: {len(summary.columns)}")
    print()

    # Print column summary
    print("Columns in summary table:")
    for i, col in enumerate(summary.columns, 1):
        non_null = summary[col].notna().sum()
        pct = 100 * non_null / len(summary)
        print(f"  {i:2d}. {col:30s} ({non_null:5d} / {len(summary):5d} = {pct:5.1f}%)")

    print()
    print("=" * 70)
    print("Final contig summary table created successfully!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
