"""
Compound information retrieval from PubChem.
"""

import pubchempy as pcp
import pandas as pd
from .utils.validators import validate_cid, validate_smiles
from .utils.api_wrappers import query_pubchem


class CompoundRetriever:
    """Handle compound information retrieval from PubChem."""
    
    def __init__(self, logger):
        """
        Initialize CompoundRetriever.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.compound_data = None
    
    def get_compound_info(self, cid=None, smiles=None):
        """
        Retrieve compound information from PubChem.
        
        Args:
            cid: PubChem Compound ID
            smiles: SMILES string
        
        Returns:
            dict: Compound information
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 1] RETRIEVING COMPOUND INFORMATION FROM PUBCHEM")
        self.logger.info("="*70)
        
        try:
            if cid:
                cid = validate_cid(cid)
                self.logger.info(f"Querying PubChem CID: {cid}")
                compound = pcp.Compound.from_cid(cid)
            elif smiles:
                smiles = validate_smiles(smiles)
                self.logger.info(f"Querying PubChem SMILES: {smiles}")
                compounds = pcp.get_compounds(smiles, 'smiles')
                if not compounds:
                    raise ValueError("No compound found for given SMILES")
                compound = compounds[0]
            else:
                raise ValueError("Either CID or SMILES must be provided")
            
            # Extract compound data
            # Extract compound data (compatible with all PubChemPy versions)
            def get_smiles_property(compound, preferred, fallback):
                """Get SMILES with fallback for different PubChemPy versions."""
                for attr in [preferred, fallback]:
                    if hasattr(compound, attr):
                        return getattr(compound, attr)
                return 'N/A'

            self.compound_data = {
                'cid': compound.cid,
                'canonical_smiles': get_smiles_property(compound, 'connectivity_smiles', 'canonical_smiles'),
                'isomeric_smiles': get_smiles_property(compound, 'smiles', 'isomeric_smiles'),
                'molecular_formula': compound.molecular_formula,
                'molecular_weight': compound.molecular_weight,
                'iupac_name': compound.iupac_name if hasattr(compound, 'iupac_name') else 'N/A'
            }

            
            self.logger.info("\n✓ Compound information retrieved successfully!")
            self.logger.info(f"  CID: {self.compound_data['cid']}")
            self.logger.info(f"  Formula: {self.compound_data['molecular_formula']}")
            self.logger.info(f"  Molecular Weight: {self.compound_data['molecular_weight']}")
            self.logger.info(f"  Canonical SMILES: {self.compound_data['canonical_smiles']}")
            
            return self.compound_data
            
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: Failed to retrieve compound information")
            self.logger.error(f"   Reason: {str(e)}")
            self.logger.error("\n💡 Suggestions:")
            self.logger.error("   - Verify CID at: https://pubchem.ncbi.nlm.nih.gov/")
            self.logger.error("   - Check SMILES format")
            raise SystemExit(1)
    
    def save_compound_info(self, output_path):
        """
        Save compound information to CSV.
        
        Args:
            output_path: Path to save CSV file
        """
        if self.compound_data is None:
            raise ValueError("No compound data to save")
        
        df = pd.DataFrame([self.compound_data])
        df.to_csv(output_path, index=False)
        
        self.logger.info(f"\n📄 Compound information saved to: {output_path}")
        self.logger.info("   This CSV contains: CID, SMILES, molecular formula, and weight")
