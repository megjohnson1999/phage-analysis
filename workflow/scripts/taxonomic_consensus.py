#!/usr/bin/env python3
"""
Taxonomic Consensus Script for Phage Analysis Pipeline

This script implements a hierarchical consensus approach to combine taxonomic
predictions from multiple tools, following the logic from the R script by
Luis Chica.

Priority order (highest to lowest):
1. CrassUS (for Crassvirales, when available)
2. mmseqs2/BLAST (protein similarity-based)
3. Phabox2 (phage-specific ML predictions)  
4. vContact3 (gene content-based clustering)

The consensus ensures taxonomic consistency by only accepting lower-level
classifications when higher taxonomic levels match across tools.
"""

import argparse
import pandas as pd
import numpy as np
import sys
import json
from pathlib import Path


def load_mmseqs_taxonomy(file_path):
    """Load and process mmseqs2 taxonomy results."""
    if not Path(file_path).exists():
        print(f"Warning: mmseqs2 file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        # Load mmseqs2 results
        df = pd.read_csv(file_path, sep='\t')
        
        # Filter for viral hits with reasonable identity
        df_filtered = df[
            df['taxlineage'].str.contains('d_Viruses', case=False, na=False) &
            (df['pident'] >= 70)
        ]
        
        if df_filtered.empty:
            print("Warning: No viral hits found in mmseqs2 results")
            return pd.DataFrame()
        
        # Get best hit by bitscore for each query
        df_best = df_filtered.loc[df_filtered.groupby('query')['bits'].idxmax()]
        
        # Parse taxonomy from taxlineage
        def parse_taxonomy_line(taxline):
            """Parse NCBI-style taxonomy string."""
            tax_dict = {
                'blast_superkingdom': 'Viruses',  # Default for viral hits
                'blast_phylum': None,
                'blast_class': None,
                'blast_order': None,
                'blast_family': None,
                'blast_genus': None,
                'blast_species': None
            }
            
            if pd.isna(taxline) or taxline == 'root':
                return tax_dict
            
            # Split by semicolon and process each level
            parts = str(taxline).split(';')
            for part in parts:
                part = part.strip()
                if part.startswith('d_'):
                    tax_dict['blast_superkingdom'] = part[2:]
                elif part.startswith('p_'):
                    tax_dict['blast_phylum'] = part[2:]
                elif part.startswith('c_'):
                    tax_dict['blast_class'] = part[2:]
                elif part.startswith('o_'):
                    tax_dict['blast_order'] = part[2:]
                elif part.startswith('f_'):
                    tax_dict['blast_family'] = part[2:]
                elif part.startswith('g_'):
                    tax_dict['blast_genus'] = part[2:]
                elif part.startswith('s_'):
                    tax_dict['blast_species'] = part[2:]
            
            return tax_dict
        
        # Apply taxonomy parsing
        taxonomy_data = df_best['taxlineage'].apply(parse_taxonomy_line)
        taxonomy_df = pd.DataFrame(list(taxonomy_data))
        
        # Combine with query IDs
        result = pd.concat([df_best[['query']].reset_index(drop=True), 
                           taxonomy_df.reset_index(drop=True)], axis=1)
        result = result.rename(columns={'query': 'contigID'})
        
        return result
        
    except Exception as e:
        print(f"Error loading mmseqs2 results: {e}")
        return pd.DataFrame()


def load_phabox_taxonomy(taxonomy_file, lifestyle_file):
    """Load and process Phabox2 taxonomy and lifestyle results."""
    taxonomy_df = pd.DataFrame()
    lifestyle_df = pd.DataFrame()
    
    # Load taxonomy
    if Path(taxonomy_file).exists():
        try:
            taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
            # Rename columns to match expected format
            if 'contig_id' in taxonomy_df.columns:
                taxonomy_df = taxonomy_df.rename(columns={'contig_id': 'contigID'})
        except Exception as e:
            print(f"Warning: Error loading Phabox2 taxonomy: {e}")
    
    # Load lifestyle
    if Path(lifestyle_file).exists():
        try:
            lifestyle_df = pd.read_csv(lifestyle_file, sep='\t')
            if 'contig_id' in lifestyle_df.columns:
                lifestyle_df = lifestyle_df.rename(columns={'contig_id': 'contigID'})
        except Exception as e:
            print(f"Warning: Error loading Phabox2 lifestyle: {e}")
    
    # For now, we focus on taxonomy. Lifestyle will be handled separately.
    # Phabox2 provides taxonomy predictions, but the exact format depends on the tool's output
    # We'll create a basic structure for now
    
    if not taxonomy_df.empty:
        # Create standardized taxonomy columns
        phabox_tax = pd.DataFrame()
        phabox_tax['contigID'] = taxonomy_df['contigID']
        phabox_tax['phabox_superkingdom'] = 'Viruses'  # Phabox2 only predicts viral sequences
        
        # Map Phabox2 predictions to standard taxonomy levels
        # Note: This will depend on the actual Phabox2 output format
        for col in ['phabox_phylum', 'phabox_class', 'phabox_order', 
                   'phabox_family', 'phabox_genus']:
            phabox_tax[col] = None
            
        return phabox_tax
    
    return pd.DataFrame()


def load_vcontact3_taxonomy(results_dir):
    """Load and process vContact3 taxonomy results."""
    vc3_file = Path(results_dir) / "final_assignments.csv"
    
    if not vc3_file.exists():
        print(f"Warning: vContact3 results not found: {vc3_file}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(vc3_file)
        
        # Filter for relevant contigs
        df_filtered = df[
            (df['GenomeName'].str.contains('contig|edge|virus_comp', case=False, na=False)) &
            (df['realm (prediction)'].notna()) &
            (df['realm (prediction)'] != 'No Realm') &
            (df['realm (prediction)'] != 'No prediction')
        ]
        
        if df_filtered.empty:
            return pd.DataFrame()
        
        # Standardize column names
        df_filtered = df_filtered.rename(columns={
            'GenomeName': 'contigID',
            'realm (prediction)': 'vc3_realm',
            'phylum (prediction)': 'vc3_phylum',
            'class (prediction)': 'vc3_class',
            'order (prediction)': 'vc3_order',
            'family (prediction)': 'vc3_family',
            'subfamily (prediction)': 'vc3_subfamily',
            'genus (prediction)': 'vc3_genus'
        })
        
        # Clean up "novel" predictions
        taxonomy_cols = ['vc3_realm', 'vc3_phylum', 'vc3_class', 'vc3_order', 
                        'vc3_family', 'vc3_subfamily', 'vc3_genus']
        
        for col in taxonomy_cols:
            if col in df_filtered.columns:
                df_filtered[col] = df_filtered[col].apply(
                    lambda x: None if pd.isna(x) or 'novel' in str(x).lower() else x
                )
        
        # Add superkingdom
        df_filtered['vc3_superkingdom'] = 'Viruses'
        
        # Standardize specific taxonomy names
        if 'vc3_phylum' in df_filtered.columns:
            df_filtered['vc3_phylum'] = df_filtered['vc3_phylum'].apply(
                lambda x: 'Phixviricota' if pd.notna(x) and 'Phixviricota' in str(x) else x
            )
        
        if 'vc3_class' in df_filtered.columns:
            df_filtered['vc3_class'] = df_filtered['vc3_class'].apply(
                lambda x: 'Malgrandaviricetes' if pd.notna(x) and 'Malgrandaviricetes' in str(x) else x
            )
        
        # Select final columns
        result_cols = ['contigID', 'vc3_superkingdom', 'vc3_phylum', 'vc3_class', 
                      'vc3_order', 'vc3_family', 'vc3_subfamily', 'vc3_genus']
        
        return df_filtered[result_cols]
        
    except Exception as e:
        print(f"Error loading vContact3 results: {e}")
        return pd.DataFrame()


def load_crassus_taxonomy(file_path):
    """Load and process CrassUS taxonomy results (optional)."""
    if not Path(file_path).exists():
        print(f"CrassUS file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        # Filter valid CrassUS results
        df_filtered = df[
            (df['discard'] == False) & 
            (df['family'] != 'outgroup')
        ]
        
        if df_filtered.empty:
            return pd.DataFrame()
        
        # Create standardized taxonomy
        result = pd.DataFrame()
        result['contigID'] = df_filtered['contig']
        result['crassus_superkingdom'] = 'Viruses'
        result['crassus_phylum'] = 'Uroviricota'
        result['crassus_class'] = 'Caudoviricetes'
        result['crassus_order'] = 'Crassvirales'
        result['crassus_family'] = df_filtered['family']
        result['crassus_genus'] = df_filtered['genus']
        result['crassus_species'] = df_filtered['species']
        
        # Replace 'unknown' with None
        taxonomy_cols = ['crassus_family', 'crassus_genus', 'crassus_species']
        for col in taxonomy_cols:
            result[col] = result[col].apply(lambda x: None if x == 'unknown' else x)
        
        return result
        
    except Exception as e:
        print(f"Error loading CrassUS results: {e}")
        return pd.DataFrame()


def create_consensus_taxonomy(mmseqs_df, phabox_df, vc3_df, crassus_df=None):
    """
    Create consensus taxonomy using hierarchical priority system.
    
    Priority order: CrassUS > mmseqs2 > Phabox2 > vContact3
    """
    
    # Get all unique contig IDs
    all_contigs = set()
    for df in [mmseqs_df, phabox_df, vc3_df]:
        if not df.empty:
            all_contigs.update(df['contigID'].tolist())
    
    if crassus_df is not None and not crassus_df.empty:
        all_contigs.update(crassus_df['contigID'].tolist())
    
    if not all_contigs:
        print("Warning: No contigs found in any taxonomy results")
        return pd.DataFrame()
    
    # Create base dataframe
    consensus_df = pd.DataFrame({'contigID': list(all_contigs)})
    
    # Merge all taxonomy dataframes
    if not mmseqs_df.empty:
        consensus_df = consensus_df.merge(mmseqs_df, on='contigID', how='left')
    
    if not phabox_df.empty:
        consensus_df = consensus_df.merge(phabox_df, on='contigID', how='left')
    
    if not vc3_df.empty:
        consensus_df = consensus_df.merge(vc3_df, on='contigID', how='left')
    
    if crassus_df is not None and not crassus_df.empty:
        consensus_df = consensus_df.merge(crassus_df, on='contigID', how='left')
    
    # Helper function to check if two values match (including NA cases)
    def values_match(val1, val2):
        if pd.isna(val1) and pd.isna(val2):
            return True
        if pd.isna(val1) or pd.isna(val2):
            return True  # NA matches anything
        return val1 == val2
    
    # Apply consensus logic for each taxonomic level
    def get_consensus_value(row, level, sources):
        """Get consensus value for a taxonomic level."""
        for source in sources:
            col_name = f"{source}_{level}"
            if col_name in row and pd.notna(row[col_name]):
                # Check if this source's higher levels are consistent
                # with the already determined consensus levels
                valid = True
                
                # Check consistency with higher levels
                if level == 'phylum' and 'superkingdom' in consensus_df.columns:
                    higher_col = f"{source}_superkingdom"
                    if (higher_col in row and 
                        not values_match(row[higher_col], row.get('superkingdom'))):
                        valid = False
                
                elif level == 'class' and 'phylum' in consensus_df.columns:
                    for higher_level in ['superkingdom', 'phylum']:
                        higher_col = f"{source}_{higher_level}"
                        if (higher_col in row and 
                            not values_match(row[higher_col], row.get(higher_level))):
                            valid = False
                            break
                
                elif level == 'order':
                    for higher_level in ['superkingdom', 'phylum', 'class']:
                        higher_col = f"{source}_{higher_level}"
                        if (higher_col in row and 
                            not values_match(row[higher_col], row.get(higher_level))):
                            valid = False
                            break
                
                elif level == 'family':
                    for higher_level in ['superkingdom', 'phylum', 'class', 'order']:
                        higher_col = f"{source}_{higher_level}"
                        if (higher_col in row and 
                            not values_match(row[higher_col], row.get(higher_level))):
                            valid = False
                            break
                
                elif level == 'genus':
                    for higher_level in ['superkingdom', 'phylum', 'class', 'order', 'family']:
                        higher_col = f"{source}_{higher_level}"
                        if (higher_col in row and 
                            not values_match(row[higher_col], row.get(higher_level))):
                            valid = False
                            break
                
                elif level == 'species':
                    for higher_level in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']:
                        higher_col = f"{source}_{higher_level}"
                        if (higher_col in row and 
                            not values_match(row[higher_col], row.get(higher_level))):
                            valid = False
                            break
                
                if valid:
                    return row[col_name]
        
        return None
    
    # Define source priority (highest to lowest)
    sources = ['crassus', 'blast', 'phabox', 'vc3'] if crassus_df is not None else ['blast', 'phabox', 'vc3']
    
    # Apply consensus for each level
    taxonomy_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    for level in taxonomy_levels:
        consensus_df[level] = consensus_df.apply(
            lambda row: get_consensus_value(row, level, sources), axis=1
        )
    
    # Default superkingdom to Viruses if missing
    consensus_df['superkingdom'] = consensus_df['superkingdom'].fillna('Viruses')
    
    # Select final columns
    final_columns = ['contigID'] + taxonomy_levels
    result = consensus_df[final_columns].copy()
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Create consensus taxonomy from multiple tools')
    parser.add_argument('--mmseqs-taxonomy', required=True,
                       help='Path to mmseqs2 taxonomy results (TSV)')
    parser.add_argument('--phabox-taxonomy', required=True,
                       help='Path to Phabox2 taxonomy results (TSV)')
    parser.add_argument('--phabox-lifestyle', required=True,
                       help='Path to Phabox2 lifestyle results (TSV)')
    parser.add_argument('--vcontact3-dir', required=True,
                       help='Path to vContact3 results directory')
    parser.add_argument('--crassus-taxonomy', default=None,
                       help='Path to CrassUS taxonomy results (TSV, optional)')
    parser.add_argument('--output-taxonomy', required=True,
                       help='Output path for consensus taxonomy (TSV)')
    parser.add_argument('--output-summary', required=True,
                       help='Output path for consensus summary (JSON)')
    
    args = parser.parse_args()
    
    print("Loading taxonomy results from multiple tools...")
    
    # Load results from each tool
    mmseqs_df = load_mmseqs_taxonomy(args.mmseqs_taxonomy)
    phabox_df = load_phabox_taxonomy(args.phabox_taxonomy, args.phabox_lifestyle)
    vc3_df = load_vcontact3_taxonomy(args.vcontact3_dir)
    
    crassus_df = None
    if args.crassus_taxonomy:
        crassus_df = load_crassus_taxonomy(args.crassus_taxonomy)
    
    print(f"Loaded results:")
    print(f"  - mmseqs2: {len(mmseqs_df)} contigs")
    print(f"  - Phabox2: {len(phabox_df)} contigs")
    print(f"  - vContact3: {len(vc3_df)} contigs")
    if crassus_df is not None:
        print(f"  - CrassUS: {len(crassus_df)} contigs")
    
    # Create consensus taxonomy
    consensus_df = create_consensus_taxonomy(mmseqs_df, phabox_df, vc3_df, crassus_df)
    
    if consensus_df.empty:
        print("Error: No consensus taxonomy could be created")
        sys.exit(1)
    
    print(f"Created consensus taxonomy for {len(consensus_df)} contigs")
    
    # Save consensus taxonomy
    consensus_df.to_csv(args.output_taxonomy, sep='\t', index=False)
    print(f"Saved consensus taxonomy: {args.output_taxonomy}")
    
    # Create summary statistics
    summary = {
        'total_contigs': len(consensus_df),
        'tool_contributions': {
            'mmseqs2': len(mmseqs_df),
            'phabox2': len(phabox_df),
            'vcontact3': len(vc3_df)
        },
        'taxonomy_coverage': {}
    }
    
    if crassus_df is not None:
        summary['tool_contributions']['crassus'] = len(crassus_df)
    
    # Calculate coverage for each taxonomic level
    for level in ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        non_null = consensus_df[level].notna().sum()
        summary['taxonomy_coverage'][level] = {
            'count': int(non_null),
            'percentage': round(100 * non_null / len(consensus_df), 2)
        }
    
    # Save summary
    with open(args.output_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Saved consensus summary: {args.output_summary}")
    print("Taxonomic consensus completed successfully!")


if __name__ == '__main__':
    main()