"""
API wrapper functions for external services.
"""

import requests
import time
from typing import List, Dict, Optional


def query_pubchem(cid=None, smiles=None):
    """
    Query PubChem for compound information.
    
    Args:
        cid: PubChem Compound ID
        smiles: SMILES string
    
    Returns:
        dict: Compound information with standardized keys
    
    Raises:
        requests.RequestException: If query fails
    """
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # Request multiple SMILES representations to handle API variations
    properties = "CanonicalSMILES,IsomericSMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight,IUPACName"
    
    if cid:
        url = f"{base_url}/compound/cid/{cid}/property/{properties}/JSON"
    elif smiles:
        url = f"{base_url}/compound/smiles/{smiles}/property/{properties}/JSON"
    else:
        raise ValueError("Either CID or SMILES must be provided")
    
    response = requests.get(url, timeout=30)
    time.sleep(0.3)  # Respect PubChem rate limits
    
    if response.status_code == 200:
        data = response.json()
        result = data['PropertyTable']['Properties'][0]
        
        # Standardize: Ensure CanonicalSMILES exists
        if 'CanonicalSMILES' not in result and 'ConnectivitySMILES' in result:
            result['CanonicalSMILES'] = result['ConnectivitySMILES']
        
        return result
    else:
        raise requests.RequestException(f"PubChem query failed with status {response.status_code}")



def query_reactome(query_term, species="Homo sapiens"):
    """
    Search Reactome for pathways with robust error handling.
    
    Args:
        query_term: Search keyword or pathway ID
        species: Species name
    
    Returns:
        list: List of pathway dictionaries
    
    Raises:
        ConnectionError: If network issues occur
    """
    base_url = "https://reactome.org/ContentService/search/query"
    
    params = {
        'query': query_term,
        'species': species,
        'types': 'Pathway'
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        time.sleep(0.5)
        
        if response.status_code == 200:
            data = response.json()
            if 'results' in data and data['results']:
                return data['results']
            return []
        elif response.status_code == 503:
            raise ConnectionError("Reactome service temporarily unavailable (503)")
        else:
            raise ConnectionError(f"Reactome API returned status {response.status_code}")
            
    except requests.exceptions.ConnectionError as e:
        raise ConnectionError(f"Network connection failed: {str(e)}")
    except requests.exceptions.Timeout:
        raise ConnectionError("Request timed out - Reactome server may be slow")
    except requests.exceptions.RequestException as e:
        raise ConnectionError(f"Request failed: {str(e)}")



def get_pathway_proteins(pathway_id):
    """
    Get proteins in a specific Reactome pathway.
    
    Args:
        pathway_id: Reactome pathway ID
    
    Returns:
        list: List of protein dictionaries
    """
    # Use the 'participants' endpoint - it returns 200 OK for Pathways
    url = f"https://reactome.org/ContentService/data/participants/{pathway_id}"
    
    try:
        response = requests.get(url, timeout=30)
        time.sleep(0.5)
        
        proteins = []
        if response.status_code == 200:
            data = response.json()
            
            # Iterate through the physical entities (PEs)
            for pe in data:
                # UniProt info is nested inside 'refEntities'
                ref_entities = pe.get('refEntities', [])
                if not ref_entities:
                    continue
                    
                for ref in ref_entities:
                    # Robustly check if this is a UniProt entry
                    # The API response often omits 'databaseName', so we check multiple fields
                    is_uniprot = False
                    if ref.get('databaseName') == 'UniProt':
                        is_uniprot = True
                    elif ref.get('schemaClass') in ['ReferenceGeneProduct', 'ReferenceIsoform']:
                        is_uniprot = True
                    elif str(ref.get('stId', '')).startswith('uniprot:'):
                        is_uniprot = True
                    elif str(ref.get('displayName', '')).startswith('UniProt:'):
                        is_uniprot = True
                        
                    if is_uniprot:
                        # 1. Extract UniProt ID
                        uniprot_id = ref.get('identifier', 'N/A')
                        
                        # 2. Extract Gene Name (try multiple sources)
                        gene_name = None
                        
                        # Try explicit geneName list in refEntity
                        if 'geneName' in ref and ref['geneName']:
                            gene_name = ref['geneName'][0]
                        
                        # Try parsing refEntity displayName (e.g., "UniProt:P07333 CSF1R")
                        if not gene_name and 'displayName' in ref:
                            parts = ref['displayName'].split(' ')
                            if len(parts) > 1:
                                gene_name = parts[-1]
                        
                        # Try parsing PhysicalEntity displayName (e.g., "CSF1R [plasma membrane]")
                        if not gene_name and 'displayName' in pe:
                            # Take the part before the first bracket
                            gene_name = pe['displayName'].split(' [')[0]
                            
                        if not gene_name:
                            gene_name = 'Unknown'
                            
                        proteins.append({
                            'gene_name': gene_name,
                            'uniprot_id': uniprot_id,
                            'protein_name': pe.get('displayName', 'N/A')
                        })
        
        return proteins

    except requests.RequestException:
        # Return empty list on failure instead of crashing
        return []


def query_string(gene_list: List[str], species: int = 9606, required_score: int = 700):
    """
    Query STRING database for protein interactions.
    
    Args:
        gene_list: List of gene names
        species: NCBI taxonomy ID (9606 = Homo sapiens)
        required_score: Minimum interaction score (0-1000)
    
    Returns:
        dict: Interactions data
    """
    string_api_url = "https://string-db.org/api/tsv/network"
    
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": required_score,
        "network_type": "functional"
    }
    
    response = requests.post(string_api_url, data=params, timeout=60)
    time.sleep(1)
    
    if response.status_code == 200:
        return response.text
    else:
        raise requests.RequestException(f"STRING query failed with status {response.status_code}")


def query_gprofiler(gene_list: List[str], organism: str = "hsapiens"):
    """
    Query g:Profiler for functional enrichment.
    
    Args:
        gene_list: List of gene names
        organism: Organism code (hsapiens, mmusculus, etc.)
    
    Returns:
        dict: Enrichment results
    """
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    
    payload = {
        "organism": organism,
        "query": gene_list,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr"
    }
    
    response = requests.post(url, json=payload, timeout=60)
    time.sleep(1)
    
    if response.status_code == 200:
        return response.json()
    else:
        raise requests.RequestException(f"g:Profiler query failed with status {response.status_code}")
