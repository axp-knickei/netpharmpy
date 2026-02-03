"""
AlphaFold PDB Fetcher

Automates the retrieval of predicted protein structures from the AlphaFold Database.
"""

import os
import requests
import time
import argparse
import pandas as pd

def get_uniprot_id(gene_name, organism_id="9606"):
    """
    Retrieve UniProt ID for a given gene name using the UniProt API. 
    
    Args:
        gene_name (str): Gene symbol (e.g., "MAPK1")
        organism_id (str): NCBI Taxonomy ID (9606 for Human)
        
    Returns:
        str: UniProt Accession ID (e.g., "P28482") or None
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene_name} AND organism_id:{organism_id} AND reviewed:true",
        "format": "json",
        "fields": "accession",
        "size": 1
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if data['results']:
            return data['results'][0]['primaryAccession']
        else:
            print(f"  [WARN] No UniProt ID found for {gene_name}")
            return None
            
    except Exception as e:
        print(f"  [ERROR] UniProt API failed for {gene_name}: {e}")
        return None

def fetch_alphafold_structure(uniprot_id, output_dir):
    """
    Download the AlphaFold predicted structure (PDB format).
    
    Args:
        uniprot_id (str): UniProt Accession ID
        output_dir (str): Directory to save PDB
        
    Returns:
        str: Path to saved file or None
    """
    # AlphaFold naming convention (v4 is current standard)
    filename = f"AF-{uniprot_id}-F1-model_v4.pdb"
    url = f"https://alphafold.ebi.ac.uk/files/{filename}"
    
    save_path = os.path.join(output_dir, filename)
    
    if os.path.exists(save_path):
        print(f"  [SKIP] File already exists: {filename}")
        return save_path
        
    try:
        print(f"  Downloading {filename}...")
        response = requests.get(url)
        
        if response.status_code == 200:
            with open(save_path, "wb") as f:
                f.write(response.content)
            return save_path
        elif response.status_code == 404:
            print(f"  [WARN] Structure not found in AlphaFold DB (404) for {uniprot_id}")
            return None
        else:
            print(f"  [ERROR] Download failed: HTTP {response.status_code}")
            return None
            
    except Exception as e:
        print(f"  [ERROR] Connection error: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Fetch AlphaFold structures for a list of genes.")
    parser.add_argument("--genes", nargs="+", help="List of gene symbols (e.g. MAPK1 TNF)")
    parser.add_argument("--file", help="Path to a text file containing gene list (one per line) or CSV with 'gene_name' column")
    parser.add_argument("--top", type=int, help="Limit to the top N genes (requires CSV with 'Degree' column or pre-sorted list)")
    parser.add_argument("--out", default="./data/alphafold_pdbs", help="Output directory")
    
    args = parser.parse_args()
    
    os.makedirs(args.out, exist_ok=True)
    
    genes_to_fetch = []
    
    if args.genes:
        genes_to_fetch.extend(args.genes)
        
    if args.file:
        if args.file.endswith('.csv'):
            df = pd.read_csv(args.file)
            
            # Logic for Network Metrics (Sort by Degree if possible)
            if 'Degree' in df.columns and args.top:
                print(f"📊 Sorting by 'Degree' to find Top {args.top} Hubs...")
                df = df.sort_values('Degree', ascending=False)
            
            # Limit rows if requested
            if args.top:
                df = df.head(args.top)

            if 'gene_name' in df.columns:
                genes_to_fetch.extend(df['gene_name'].tolist())
            elif 'Protein' in df.columns: # Support network_metrics.csv
                genes_to_fetch.extend(df['Protein'].tolist())
        else:
            with open(args.file, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]
                if args.top:
                    lines = lines[:args.top]
                genes_to_fetch.extend(lines)
    
    # Deduplicate (but keep order if possible for priority?) 
    # Actually, set() destroys order. If 'top' logic was applied, we want those specific ones.
    # Let's simple deduplicate while preserving order
    seen = set()
    unique_genes = []
    for g in genes_to_fetch:
        if g not in seen:
            unique_genes.append(g)
            seen.add(g)
    genes_to_fetch = unique_genes
    
    print(f"🔍 Found {len(genes_to_fetch)} genes to process.")
    print(f"📂 Saving to: {args.out}\n")
    
    results = []
    
    for gene in genes_to_fetch:
        print(f"Processing {gene}...")
        
        # Step 1: Get UniProt ID
        uniprot_id = get_uniprot_id(gene)
        
        if uniprot_id:
            print(f"  UniProt ID: {uniprot_id}")
            # Step 2: Download PDB
            pdb_path = fetch_alphafold_structure(uniprot_id, args.out)
            
            status = "Success" if pdb_path else "Failed"
            results.append({
                "gene": gene,
                "uniprot": uniprot_id,
                "file": os.path.basename(pdb_path) if pdb_path else None,
                "status": status
            })
        else:
            results.append({
                "gene": gene,
                "uniprot": None,
                "file": None,
                "status": "UniProt ID Not Found"
            })
            
        # Be nice to APIs
        time.sleep(1)
        
    # Save log
    log_df = pd.DataFrame(results)
    log_path = os.path.join(args.out, "fetch_log.csv")
    log_df.to_csv(log_path, index=False)
    print(f"\n✅ Done! Log saved to {log_path}")

if __name__ == "__main__":
    main()
