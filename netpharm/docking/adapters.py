"""
docking/adapters.py

Adapters to harmonize enrichment results from different tools (e.g., DAVID)
into a standard format compatible with netpharm plotting utilities.

CONTEXT:
This module is primarily used for the CONFIRMATORY phase of analysis (Figure 4),
where a small, high-affinity gene subset (filtered by docking) is analyzed
using DAVID, as per the original study methodology.
"""

import pandas as pd
from pathlib import Path
import re

def parse_david_results(directory: Path, prefix: str = "") -> pd.DataFrame:
    """
    Parse DAVID enrichment results (CSV) into a standardized DataFrame.
    
    Args:
        directory: Path to results directory.
        prefix: Optional prefix for filenames (e.g., "mock_").
    
    Expected Input Files in directory:
    - {prefix}david_go_bp.csv
    - {prefix}david_kegg.csv
    - {prefix}david_reactome.csv
    
    Expected DAVID CSV Columns:
    - 'Term': e.g., "GO:0006954~inflammatory response"
    - 'Count': Gene count in term
    - '%': Percentage of input gene list
    - 'PValue': Unadjusted p-value
    - 'Benjamini': Adjusted p-value (FDR)
    
    Returns:
        pd.DataFrame: Standardized schema matching g:Profiler output.
        - Source (GO:BP, KEGG, REAC)
        - Term_Name (Cleaned name)
        - Term_ID (Cleaned ID)
        - Intersection_Size (Count)
        - Recall (Derived from %) -> Gene Ratio (0-1)
        - Adjusted_P_value (Benjamini)
    """
    directory = Path(directory)
    all_frames = []
    
    # Map filenames to Source identifiers
    file_map = {
        'david_go_bp.csv': 'GO:BP',
        'david_kegg.csv': 'KEGG',
        'david_reactome.csv': 'REAC'
    }
    
    for filename, source in file_map.items():
        filepath = directory / (prefix + filename)
        if not filepath.exists():
            continue
            
        try:
            # Load data
            df = pd.read_csv(filepath)
            
            # Validate essential columns
            required = ['Term', 'Count', 'Benjamini', '%']
            if not set(required).issubset(df.columns):
                print(f"Warning: {filename} missing columns. Expected {required}. Found {df.columns.tolist()}")
                continue
            
            # Parse 'Term' into ID and Name
            # Format is usually "ID~Name"
            def parse_term(t):
                if '~' in str(t):
                    return str(t).split('~', 1)
                return [str(t), str(t)] # Fallback
            
            parsed_terms = df['Term'].apply(parse_term)
            
            clean_df = pd.DataFrame()
            clean_df['Source'] = [source] * len(df)
            clean_df['Term_ID'] = [x[0] for x in parsed_terms]
            clean_df['Term_Name'] = [x[1] for x in parsed_terms]
            
            clean_df['Intersection_Size'] = df['Count']
            
            # Handle Recall / Gene Ratio
            # DAVID '%' column is "Percentage of involved genes / Total genes" * 100
            # We convert to 0-1 ratio for consistency with g:Profiler 'Recall'
            clean_df['Recall'] = df['%'] / 100.0
            
            clean_df['Adjusted_P_value'] = df['Benjamini']
            
            # Optional: Keep raw PValue if needed
            if 'PValue' in df.columns:
                clean_df['P_value'] = df['PValue']
            
            all_frames.append(clean_df)
            
        except Exception as e:
            print(f"Error parsing {filename}: {e}")
    
    if not all_frames:
        return pd.DataFrame()
        
    final_df = pd.concat(all_frames, ignore_index=True)
    return final_df.sort_values('Adjusted_P_value')
