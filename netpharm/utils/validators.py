"""
Input validation functions.
"""

import re


def validate_cid(cid):
    """
    Validate PubChem CID format.
    
    Args:
        cid: Compound ID (can be string or int)
    
    Returns:
        int: Valid CID
    
    Raises:
        ValueError: If CID is invalid
    """
    try:
        cid_int = int(cid)
        if cid_int <= 0:
            raise ValueError("CID must be a positive integer")
        return cid_int
    except (ValueError, TypeError):
        raise ValueError(f"Invalid CID format: {cid}. Must be a positive integer.")


def validate_smiles(smiles):
    """
    Basic SMILES format validation.
    
    Args:
        smiles: SMILES string
    
    Returns:
        str: Valid SMILES string
    
    Raises:
        ValueError: If SMILES is invalid
    """
    if not isinstance(smiles, str):
        raise ValueError("SMILES must be a string")
    
    if len(smiles) < 3:
        raise ValueError("SMILES string too short")
    
    # Basic character check
    allowed_chars = set('CNOPSFClBrI[]()=#@+-\\/0123456789cnops')
    if not all(c in allowed_chars for c in smiles):
        raise ValueError(f"SMILES contains invalid characters: {smiles}")
    
    return smiles.strip()


def validate_threshold(threshold, min_val=0.0, max_val=1.0):
    """
    Validate probability threshold.
    
    Args:
        threshold: Threshold value
        min_val: Minimum allowed value
        max_val: Maximum allowed value
    
    Returns:
        float: Valid threshold
    
    Raises:
        ValueError: If threshold is invalid
    """
    try:
        threshold_float = float(threshold)
        if not min_val <= threshold_float <= max_val:
            raise ValueError(f"Threshold must be between {min_val} and {max_val}")
        return threshold_float
    except (ValueError, TypeError):
        raise ValueError(f"Invalid threshold: {threshold}")


def validate_pathway_id(pathway_id):
    """
    Validate Reactome pathway ID format.
    
    Args:
        pathway_id: Reactome pathway ID
    
    Returns:
        str: Valid pathway ID
    
    Raises:
        ValueError: If pathway ID is invalid
    """
    # Reactome IDs format: R-HSA-XXXXXX or R-MMU-XXXXXX, etc.
    pattern = r'^R-[A-Z]{3}-\d+$'
    if not re.match(pattern, pathway_id):
        raise ValueError(f"Invalid Reactome pathway ID format: {pathway_id}")
    return pathway_id
