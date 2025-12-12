"""
Basic smoke tests for network pharmacology tool.
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from netpharm.utils.validators import validate_cid, validate_smiles, validate_threshold
from netpharm.utils.api_wrappers import query_pubchem


class TestValidators:
    """Test input validation functions."""
    
    def test_validate_cid_valid(self):
        """Test valid CID validation."""
        assert validate_cid(969516) == 969516
        assert validate_cid("969516") == 969516
    
    def test_validate_cid_invalid(self):
        """Test invalid CID validation."""
        with pytest.raises(ValueError):
            validate_cid(-1)
        with pytest.raises(ValueError):
            validate_cid("invalid")
        with pytest.raises(ValueError):
            validate_cid(0)
    
    def test_validate_smiles_valid(self):
        """Test valid SMILES validation."""
        smiles = "COC1=CC=CC=C1"
        assert validate_smiles(smiles) == smiles
    
    def test_validate_smiles_invalid(self):
        """Test invalid SMILES validation."""
        with pytest.raises(ValueError):
            validate_smiles("XY")  # Too short
        with pytest.raises(ValueError):
            validate_smiles(123)  # Not a string
    
    def test_validate_threshold(self):
        """Test threshold validation."""
        assert validate_threshold(0.5) == 0.5
        assert validate_threshold("0.7") == 0.7
        
        with pytest.raises(ValueError):
            validate_threshold(1.5)  # Too high
        with pytest.raises(ValueError):
            validate_threshold(-0.1)  # Too low


class TestAPIs:
    """Test API wrapper functions."""
    
    def test_query_pubchem_by_cid(self):
        """Test PubChem query by CID."""
        result = query_pubchem(cid=969516)  # Curcumin
        
        # PubChem API may return CanonicalSMILES or ConnectivitySMILES
        has_smiles = any(key in result for key in ['CanonicalSMILES', 'ConnectivitySMILES'])
        assert has_smiles, f"No SMILES found in result: {result.keys()}"
        
        assert 'MolecularFormula' in result
        assert result['MolecularFormula'] == 'C21H20O6'
        
        # Verify CID
        assert result.get('CID') == 969516

    
    def test_query_pubchem_invalid_cid(self):
        """Test PubChem query with invalid CID."""
        with pytest.raises(Exception):
            query_pubchem(cid=999999999999)


def test_imports():
    """Test that all modules can be imported."""
    from netpharm import NetworkPharmacology
    from netpharm.compound import CompoundRetriever
    from netpharm.targets import TargetPredictor
    from netpharm.pathways import PathwayAnalyzer
    from netpharm.network import NetworkAnalyzer
    from netpharm.enrichment import EnrichmentAnalyzer
    from netpharm.visualize import NetworkVisualizer
    
    assert NetworkPharmacology is not None
    assert CompoundRetriever is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
