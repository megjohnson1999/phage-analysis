#!/usr/bin/env python3
"""
Lifestyle Consensus Script for Phage Analysis Pipeline

This script creates a consensus lifestyle prediction by combining BACPHLIP
and Phabox2 predictions with a hierarchical approach:

1. Use BACPHLIP if confidence >= 0.7 for either virulent or temperate
2. Otherwise, fall back to Phabox2 prediction

BACPHLIP predictions are preferred when confident because it's specifically
designed for lifestyle prediction using genomic features.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def load_bacphlip_results(file_path):
    """
    Load BACPHLIP lifestyle predictions.

    Expected format: Sequence\tVirulent\tTemperate
    No header in the raw BACPHLIP output.
    """
    if not Path(file_path).exists():
        print(f"Warning: BACPHLIP file not found: {file_path}")
        return pd.DataFrame()

    try:
        # BACPHLIP output has no header, so we need to add column names
        df = pd.read_csv(file_path, sep='\t', header=None,
                        names=['contig_id', 'virulent_conf', 'temperate_conf'])

        if df.empty:
            print("Warning: BACPHLIP file is empty")
            return pd.DataFrame()

        # Convert confidence columns to numeric
        df['virulent_conf'] = pd.to_numeric(df['virulent_conf'], errors='coerce')
        df['temperate_conf'] = pd.to_numeric(df['temperate_conf'], errors='coerce')

        # Remove any rows with invalid data
        df = df.dropna(subset=['virulent_conf', 'temperate_conf'])

        print(f"Loaded {len(df)} BACPHLIP predictions")
        return df

    except Exception as e:
        print(f"Error loading BACPHLIP results: {e}")
        return pd.DataFrame()


def load_phabox_lifestyle(file_path):
    """
    Load Phabox2 lifestyle predictions.

    Expected format: contig_id\tlifestyle_prediction\tconfidence
    """
    if not Path(file_path).exists():
        print(f"Warning: Phabox2 lifestyle file not found: {file_path}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(file_path, sep='\t')

        if df.empty:
            print("Warning: Phabox2 lifestyle file is empty")
            return pd.DataFrame()

        # Standardize column names if needed
        if 'contig_id' not in df.columns and df.columns[0] != 'contig_id':
            df = df.rename(columns={df.columns[0]: 'contig_id'})

        # Ensure confidence is numeric
        if 'confidence' in df.columns:
            df['confidence'] = pd.to_numeric(df['confidence'], errors='coerce')

        print(f"Loaded {len(df)} Phabox2 lifestyle predictions")
        return df

    except Exception as e:
        print(f"Error loading Phabox2 lifestyle results: {e}")
        return pd.DataFrame()


def create_lifestyle_consensus(bacphlip_df, phabox_df, confidence_threshold=0.7):
    """
    Create consensus lifestyle predictions.

    Logic:
    1. Use BACPHLIP if either virulent or temperate confidence >= threshold
    2. Otherwise, use Phabox2 prediction

    Args:
        bacphlip_df: DataFrame with BACPHLIP results
        phabox_df: DataFrame with Phabox2 results
        confidence_threshold: Minimum confidence for BACPHLIP (default 0.7)

    Returns:
        DataFrame with consensus predictions
    """
    # Get all unique contig IDs from both sources
    all_contigs = set()

    if not bacphlip_df.empty:
        all_contigs.update(bacphlip_df['contig_id'].tolist())

    if not phabox_df.empty:
        all_contigs.update(phabox_df['contig_id'].tolist())

    if not all_contigs:
        print("Warning: No contigs found in any lifestyle prediction results")
        return pd.DataFrame(columns=['contig_id', 'lifestyle', 'confidence', 'source'])

    print(f"Creating consensus for {len(all_contigs)} unique contigs")

    # Create results list
    results = []

    for contig_id in all_contigs:
        # Get BACPHLIP prediction for this contig
        bacphlip_row = bacphlip_df[bacphlip_df['contig_id'] == contig_id]

        # Get Phabox prediction for this contig
        phabox_row = phabox_df[phabox_df['contig_id'] == contig_id]

        lifestyle = None
        confidence = None
        source = None

        # Check if BACPHLIP has a confident prediction
        if not bacphlip_row.empty:
            virulent_conf = bacphlip_row.iloc[0]['virulent_conf']
            temperate_conf = bacphlip_row.iloc[0]['temperate_conf']

            # Use BACPHLIP if either confidence meets threshold
            if virulent_conf >= confidence_threshold or temperate_conf >= confidence_threshold:
                # Choose lifestyle with higher confidence
                if virulent_conf >= temperate_conf:
                    lifestyle = 'virulent'
                    confidence = virulent_conf
                else:
                    lifestyle = 'temperate'
                    confidence = temperate_conf
                source = 'BACPHLIP'

        # Fall back to Phabox if BACPHLIP didn't provide a confident prediction
        if lifestyle is None and not phabox_row.empty:
            lifestyle = phabox_row.iloc[0]['lifestyle_prediction']
            confidence = phabox_row.iloc[0].get('confidence', 0.0)
            source = 'Phabox2'

        # If still no prediction, mark as unknown
        if lifestyle is None:
            lifestyle = 'unknown'
            confidence = 0.0
            source = 'none'

        results.append({
            'contig_id': contig_id,
            'lifestyle': lifestyle,
            'confidence': confidence,
            'source': source
        })

    # Create results DataFrame
    consensus_df = pd.DataFrame(results)

    # Sort by contig_id for consistency
    consensus_df = consensus_df.sort_values('contig_id').reset_index(drop=True)

    return consensus_df


def generate_summary_stats(consensus_df):
    """Generate summary statistics for the consensus predictions."""
    summary = {}

    # Count by lifestyle
    lifestyle_counts = consensus_df['lifestyle'].value_counts().to_dict()
    summary['lifestyle_distribution'] = lifestyle_counts

    # Count by source
    source_counts = consensus_df['source'].value_counts().to_dict()
    summary['source_distribution'] = source_counts

    # Average confidence by lifestyle
    avg_conf_by_lifestyle = consensus_df.groupby('lifestyle')['confidence'].mean().to_dict()
    summary['avg_confidence_by_lifestyle'] = {k: round(v, 3) for k, v in avg_conf_by_lifestyle.items()}

    # Average confidence by source
    avg_conf_by_source = consensus_df.groupby('source')['confidence'].mean().to_dict()
    summary['avg_confidence_by_source'] = {k: round(v, 3) for k, v in avg_conf_by_source.items()}

    # Total contigs
    summary['total_contigs'] = len(consensus_df)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Create consensus lifestyle predictions from BACPHLIP and Phabox2'
    )
    parser.add_argument('--bacphlip', required=True,
                       help='Path to BACPHLIP lifestyle results (TSV)')
    parser.add_argument('--phabox', required=True,
                       help='Path to Phabox2 lifestyle results (TSV)')
    parser.add_argument('--output', required=True,
                       help='Output path for consensus lifestyle predictions (TSV)')
    parser.add_argument('--threshold', type=float, default=0.7,
                       help='Minimum confidence threshold for BACPHLIP (default: 0.7)')

    args = parser.parse_args()

    print("=" * 60)
    print("Lifestyle Consensus Analysis")
    print("=" * 60)
    print(f"BACPHLIP confidence threshold: {args.threshold}")
    print()

    # Load input files
    print("Loading lifestyle predictions...")
    bacphlip_df = load_bacphlip_results(args.bacphlip)
    phabox_df = load_phabox_lifestyle(args.phabox)
    print()

    # Create consensus
    print("Creating lifestyle consensus...")
    consensus_df = create_lifestyle_consensus(bacphlip_df, phabox_df, args.threshold)

    if consensus_df.empty:
        print("Warning: No consensus predictions could be created")
        print("Creating empty output file with proper headers...")
        consensus_df = pd.DataFrame(columns=['contig_id', 'lifestyle', 'confidence', 'source'])

    print(f"Created consensus for {len(consensus_df)} contigs")
    print()

    # Generate and print summary statistics
    if not consensus_df.empty:
        summary = generate_summary_stats(consensus_df)

        print("Summary Statistics:")
        print("-" * 60)
        print(f"Total contigs: {summary['total_contigs']}")
        print()

        print("Lifestyle distribution:")
        for lifestyle, count in summary['lifestyle_distribution'].items():
            pct = 100 * count / summary['total_contigs']
            print(f"  {lifestyle}: {count} ({pct:.1f}%)")
        print()

        print("Source distribution:")
        for source, count in summary['source_distribution'].items():
            pct = 100 * count / summary['total_contigs']
            print(f"  {source}: {count} ({pct:.1f}%)")
        print()

        print("Average confidence by lifestyle:")
        for lifestyle, conf in summary['avg_confidence_by_lifestyle'].items():
            print(f"  {lifestyle}: {conf:.3f}")
        print()

        print("Average confidence by source:")
        for source, conf in summary['avg_confidence_by_source'].items():
            print(f"  {source}: {conf:.3f}")
        print()

    # Save output
    consensus_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved consensus lifestyle predictions: {args.output}")
    print("=" * 60)
    print("Lifestyle consensus completed successfully!")


if __name__ == '__main__':
    main()
