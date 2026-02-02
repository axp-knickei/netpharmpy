import pytest
import pandas as pd
import networkx as nx
from unittest.mock import MagicMock, patch, mock_open
import sys
import os

# Add parent directory to path to ensure imports work
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from netpharm.compound import CompoundRetriever
from netpharm.targets import TargetPredictor
from netpharm.pathways import PathwayAnalyzer
from netpharm.network import NetworkAnalyzer

# Mock Logger to avoid cluttering output during tests
@pytest.fixture
def mock_logger():
    return MagicMock()

# =============================================================================
# 1. CompoundRetriever Tests
# =============================================================================
class TestCompoundRetriever:
    
    @patch('netpharm.compound.pcp')
    def test_get_compound_info_cid(self, mock_pcp, mock_logger):
        """Test retrieving compound info by CID using mocked PubChemPy."""
        # Setup mock compound
        mock_compound = MagicMock()
        mock_compound.cid = 969516
        mock_compound.canonical_smiles = 'COC1=CC=CC=C1'
        mock_compound.isomeric_smiles = 'COC1=CC=CC=C1'
        mock_compound.molecular_formula = 'C21H20O6'
        mock_compound.molecular_weight = 368.38
        mock_compound.iupac_name = 'Test Compound'
        
        # Configure the mock to return our mock_compound
        mock_pcp.Compound.from_cid.return_value = mock_compound
        
        retriever = CompoundRetriever(mock_logger)
        result = retriever.get_compound_info(cid=969516)
        
        # Assertions
        assert result['cid'] == 969516
        assert result['molecular_formula'] == 'C21H20O6'
        mock_pcp.Compound.from_cid.assert_called_once_with(969516)

    @patch('netpharm.compound.pcp')
    def test_get_compound_info_smiles(self, mock_pcp, mock_logger):
        """Test retrieving compound info by SMILES."""
        mock_compound = MagicMock()
        mock_compound.cid = 123
        mock_compound.molecular_formula = 'C2H6O'
        mock_compound.molecular_weight = 46.07
        # Mock attributes accessed via getattr in the code
        mock_compound.connectivity_smiles = 'CCO'
        mock_compound.canonical_smiles = 'CCO'
        mock_compound.isomeric_smiles = 'CCO'
        
        mock_pcp.get_compounds.return_value = [mock_compound]
        
        retriever = CompoundRetriever(mock_logger)
        result = retriever.get_compound_info(smiles='CCO')
        
        assert result['cid'] == 123
        assert result['canonical_smiles'] == 'CCO'
        mock_pcp.get_compounds.assert_called_once_with('CCO', 'smiles')

# =============================================================================
# 2. TargetPredictor Tests
# =============================================================================
class TestTargetPredictor:
    
    def test_get_all_target_genes_deduplication(self, mock_logger):
        """Test that get_all_target_genes correctly deduplicates and cleans genes."""
        predictor = TargetPredictor(mock_logger)
        
        # Manually set the internal dataframes as if they were loaded
        predictor.swiss_targets = pd.DataFrame({
            'Target': ['GENE_A', 'gene_b', 'Gene_A ']
        })
        predictor.superpred_targets = pd.DataFrame({
            'Target': ['GENE_B', 'GENE_C']
        })
        
        genes = predictor.get_all_target_genes()
        
        # Expect unique, upper-case, stripped genes: GENE_A, GENE_B, GENE_C
        assert len(genes) == 3
        assert 'GENE_A' in genes
        assert 'GENE_B' in genes
        assert 'GENE_C' in genes

    @patch('netpharm.targets.os.path.exists')
    @patch('netpharm.targets.pd.read_csv')
    @patch('builtins.input') # Mock input so it doesn't wait
    def test_predict_targets_manual_parsing(self, mock_input, mock_read_csv, mock_exists, mock_logger):
        """Test the logic of filtering and parsing the 'manual' CSV files."""
        predictor = TargetPredictor(mock_logger)
        
        # Setup mocks
        mock_exists.return_value = True # Files exist
        
        # Create mock DataFrames for Swiss and SuperPred
        # Swiss: Needs a Probability column and Common name
        swiss_df = pd.DataFrame({
            'Common name': ['GENE_A', 'GENE_B'],
            'Probability': [0.1, 0.9] # 0.1 should be filtered out if threshold is default
        })
        
        # SuperPred: Needs Target Name
        superpred_df = pd.DataFrame({
            'Target Name': ['GENE_C'],
            'Probability': ['80%'] # Test percentage string handling
        })
        
        # side_effect to return different dfs for different calls is tricky because exact order matters
        # The code reads Swiss first, then SuperPred Known, then SuperPred Predicted.
        # We'll just return compatible structures for all calls or make the mock smarter if needed.
        # For simplicity, let's just test Swiss logic or rely on the fact that read_csv will be called multiple times.
        
        def read_csv_side_effect(filepath, *args, **kwargs):
            if 'swiss_results' in filepath:
                return swiss_df
            elif 'Targets' in filepath: # SuperPred
                return superpred_df
            return pd.DataFrame()
            
        mock_read_csv.side_effect = read_csv_side_effect
        
        # Run with threshold 0.5 for Swiss
        stats = predictor.predict_targets_manual(
            smiles="dummy", 
            swiss_threshold=0.5, 
            superpred_threshold=0.0,
            output_dir="./mock_data"
        )
        
        # Check Swiss filtering (0.1 should be gone, 0.9 remains)
        assert len(predictor.swiss_targets) == 1
        assert predictor.swiss_targets.iloc[0]['Target'] == 'GENE_B'
        
        # Check SuperPred parsing (GENE_C should be there)
        assert len(predictor.superpred_targets) > 0
        assert 'GENE_C' in predictor.superpred_targets['Target'].values

# =============================================================================
# 3. PathwayAnalyzer Tests
# =============================================================================
class TestPathwayAnalyzer:

    def test_find_overlapping_targets(self, mock_logger):
        """Test the intersection logic for finding overlapping targets."""
        analyzer = PathwayAnalyzer(mock_logger)
        
        # Mock pathway proteins extraction
        # We manually set the dataframe
        analyzer.pathway_proteins = pd.DataFrame({
            'gene_name': ['G1', 'G2', 'G3'],
            'pathway_id': ['P1', 'P1', 'P2']
        })
        
        target_genes = ['G2', 'G4', 'G1'] # G4 is not in pathways
        
        result = analyzer.find_overlapping_targets(target_genes)
        
        # Should contain G1 and G2
        assert len(result) == 2
        assert 'G1' in result['gene_name'].values
        assert 'G2' in result['gene_name'].values
        assert 'G4' not in result['gene_name'].values
        assert 'G3' not in result['gene_name'].values # G3 was in pathway but not in target list

    def test_find_overlapping_targets_no_overlap(self, mock_logger):
        """Test behavior when there is no overlap."""
        analyzer = PathwayAnalyzer(mock_logger)
        analyzer.pathway_proteins = pd.DataFrame({'gene_name': ['A', 'B']})
        
        with pytest.raises(SystemExit):
            analyzer.find_overlapping_targets(['X', 'Y'])

# =============================================================================
# 4. NetworkAnalyzer Tests
# =============================================================================
class TestNetworkAnalyzer:

    def test_build_network(self, mock_logger):
        """Test building a NetworkX graph from a dataframe."""
        analyzer = NetworkAnalyzer(mock_logger)
        
        # Mock interactions dataframe
        analyzer.interactions_df = pd.DataFrame({
            'preferredName_A': ['P1', 'P2'],
            'preferredName_B': ['P2', 'P3'],
            'score': [900, 700] # Scores from STRING are 0-1000
        })
        
        G = analyzer.build_network()
        
        assert isinstance(G, nx.Graph)
        assert G.number_of_nodes() == 3 # P1, P2, P3
        assert G.number_of_edges() == 2 # P1-P2, P2-P3
        
        # Check weights (should be normalized /1000)
        edge_p1_p2 = G.get_edge_data('P1', 'P2')
        assert edge_p1_p2['weight'] == 0.9

    def test_analyze_network_topology(self, mock_logger):
        """Test that metrics are calculated correctly."""
        analyzer = NetworkAnalyzer(mock_logger)
        
        # Create a simple star graph: Center connected to Leaf1, Leaf2
        G = nx.Graph()
        G.add_edge('Center', 'Leaf1')
        G.add_edge('Center', 'Leaf2')
        analyzer.network = G
        
        metrics = analyzer.analyze_network()
        
        # Verify Center has highest degree
        center_row = metrics[metrics['Protein'] == 'Center'].iloc[0]
        leaf_row = metrics[metrics['Protein'] == 'Leaf1'].iloc[0]
        
        assert center_row['Degree'] == 2
        assert leaf_row['Degree'] == 1
        assert center_row['Degree_Centrality'] > leaf_row['Degree_Centrality']

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
