"""
Main NetworkPharmacology class that orchestrates the entire pipeline.
"""

import os
from datetime import datetime
from .compound import CompoundRetriever
from .targets import TargetPredictor
from .pathways import PathwayAnalyzer
from .network import NetworkAnalyzer
from .enrichment import EnrichmentAnalyzer
from .visualize import NetworkVisualizer
from .utils.logger import setup_logger
from .utils.config_handler import get_config


class NetworkPharmacology:
    """
    Main class for network pharmacology analysis pipeline.
    
    This class orchestrates the complete workflow from compound input
    to enrichment analysis and network visualization.
    """
    
    def __init__(self, cid=None, smiles=None, config=None, output_base='./outputs'):
        """
        Initialize NetworkPharmacology pipeline.
        
        Args:
            cid: PubChem Compound ID
            smiles: SMILES string
            config: Configuration dictionary (optional)
            output_base: Base output directory
        """
        # Determine compound identifier
        if cid:
            self.compound_id = str(cid)
            self.cid = cid
            self.smiles = None
        elif smiles:
            self.compound_id = "compound_smiles"
            self.cid = None
            self.smiles = smiles
        else:
            raise ValueError("Either CID or SMILES must be provided")
        
        # Setup output directory
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.output_dir = os.path.join(output_base, f"compound_{self.compound_id}_{timestamp}")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Setup logger
        self.logger = setup_logger(self.output_dir, self.compound_id)
        
        # Load configuration
        self.config = config if config else {}
        
        # Initialize components
        self.compound_retriever = CompoundRetriever(self.logger)
        self.target_predictor = TargetPredictor(self.logger)
        self.pathway_analyzer = PathwayAnalyzer(self.logger)
        self.network_analyzer = NetworkAnalyzer(self.logger)
        self.enrichment_analyzer = EnrichmentAnalyzer(self.logger)
        self.visualizer = NetworkVisualizer(self.output_dir)
        
        # Data storage
        self.compound_data = None
        self.target_genes = None
        self.overlapping_targets = None
        self.network = None
        
        self.logger.info("\n" + "="*70)
        self.logger.info("NETWORK PHARMACOLOGY ANALYSIS PIPELINE v1.0.0")
        self.logger.info("="*70)
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Compound identifier: {self.compound_id}")
    
    def get_compound_info(self):
        """Step 1: Retrieve compound information from PubChem."""
        step_dir = os.path.join(self.output_dir, "step1_compound_info")
        os.makedirs(step_dir, exist_ok=True)
        
        self.compound_data = self.compound_retriever.get_compound_info(
            cid=self.cid,
            smiles=self.smiles
        )
        
        # Save results
        output_path = os.path.join(step_dir, "compound_info.csv")
        self.compound_retriever.save_compound_info(output_path)
        
        return self.compound_data
    
    def predict_targets(self, swiss_threshold=0.0, superpred_threshold=0.5):
        """Step 2: Predict molecular targets (manual workflow)."""
        if self.compound_data is None:
            raise ValueError("Must retrieve compound info first")
        
        step_dir = os.path.join(self.output_dir, "step2_targets")
        os.makedirs(step_dir, exist_ok=True)
        
        # Create data directory for manual downloads
        data_dir = os.path.join(self.output_dir, "data")
        os.makedirs(data_dir, exist_ok=True)
        
        smiles = self.compound_data['canonical_smiles']
        
        stats = self.target_predictor.predict_targets_manual(
            smiles=smiles,
            swiss_threshold=swiss_threshold,
            superpred_threshold=superpred_threshold,
            output_dir=data_dir
        )
        
        # Save processed targets
        self.target_predictor.save_targets(step_dir)
        
        # Get all target genes
        self.target_genes = self.target_predictor.get_all_target_genes()
        
        return stats
    
    def analyze_pathways(self, search_terms):
        """Step 3: Analyze Reactome pathways and find overlaps."""
        if self.target_genes is None:
            raise ValueError("Must predict targets first")
        
        step_dir = os.path.join(self.output_dir, "step3_pathways")
        os.makedirs(step_dir, exist_ok=True)
        
        # Search pathways
        pathways = self.pathway_analyzer.search_pathways(search_terms)
        
        # Extract proteins
        pathway_proteins = self.pathway_analyzer.extract_pathway_proteins()
        
        # Find overlaps
        self.overlapping_targets = self.pathway_analyzer.find_overlapping_targets(
            self.target_genes
        )
        
        # Save results
        self.pathway_analyzer.save_results(step_dir)
        
        return self.overlapping_targets
    
    def build_network(self, confidence=0.700, scope="all"):
            """
            Step 4: Build and analyze STRING network.
            
            Args:
                confidence: Minimum interaction confidence
                scope: 'all' for all predicted targets, 'overlap' for pathway overlaps only
            """
            step_dir = os.path.join(self.output_dir, "step4_network")
            os.makedirs(step_dir, exist_ok=True)
            
            # Select genes based on scope
            if scope == "all":
                if self.target_genes is None:
                    raise ValueError("Must predict targets first")
                network_genes = self.target_genes
                self.logger.info(f"\nUsing ALL predicted targets for network: {len(network_genes)} genes")
            else:
                if self.overlapping_targets is None:
                    raise ValueError("Must analyze pathways first")
                network_genes = self.overlapping_targets['gene_name'].unique().tolist()
                self.logger.info(f"\nUsing only pathway-overlapping targets: {len(network_genes)} genes")
            
            # Query STRING
            interactions = self.network_analyzer.query_string_network(
                network_genes, 
                confidence=confidence
            )


            # Build network FIRST
            self.network = self.network_analyzer.build_network()

            # --- TEMP DEBUG ---
            print("\n[DEBUG] Network object inspection")
            print("Type:", type(self.network))
            print("Is directed:", self.network.is_directed())
            print("Nodes:", self.network.number_of_nodes())
            print("Edges:", self.network.number_of_edges())

            u, v, data = next(iter(self.network.edges(data=True)))
            print("Example edge:", u, v, data)

            n = next(iter(self.network.nodes()))
            print("Example node:", n)
            print("[DEBUG END]\n")

            
            # Analyze network
            metrics = self.network_analyzer.analyze_network()
            
            # Save results
            self.network_analyzer.save_results(step_dir)
            
            # Create visualizations
            self.visualizer.create_all_visualizations(self.network, step_dir)
            
            return self.network
    
    def enrichment_analysis(self, method='gprofiler'):
        """Step 5: Functional enrichment analysis."""
        if self.overlapping_targets is None:
            raise ValueError("Must analyze pathways first")
        
        step_dir = os.path.join(self.output_dir, "step5_enrichment")
        os.makedirs(step_dir, exist_ok=True)
        
        # Get gene list
        gene_list = self.overlapping_targets['gene_name'].unique().tolist()
        
        if method == 'david':
            data_dir = os.path.join(self.output_dir, "data")
            stats = self.enrichment_analyzer.analyze_david_manual(
                gene_list=gene_list,
                output_dir=data_dir
            )
            return stats
        
        elif method == 'gprofiler':
            results = self.enrichment_analyzer.analyze_gprofiler(gene_list=gene_list)
            self.enrichment_analyzer.save_results(step_dir, method='gprofiler')
            return results
        
        else:
            raise ValueError(f"Unknown enrichment method: {method}")
    
    def run_full_pipeline(self, config_path=None):
        """
        Run the complete pipeline with configuration.
        
        Args:
            config_path: Optional path to configuration file
        """
        # Get configuration
        if config_path or not self.config:
            self.config = get_config(config_path)
        
        try:
            # Step 1: Compound info
            self.get_compound_info()
            
            # Step 2: Target prediction
            swiss_thresh = self.config.get('target_prediction', {}).get('swiss_threshold', 0.0)
            superpred_thresh = self.config.get('target_prediction', {}).get('superpred_threshold', 0.5)
            self.predict_targets(swiss_threshold=swiss_thresh, superpred_threshold=superpred_thresh)
            
            # Step 3: Pathway analysis
            search_terms = self.config.get('pathways', {}).get('search_terms', [])
            self.analyze_pathways(search_terms)
            
            # Step 4: Network analysis
            confidence = self.config.get('string', {}).get('confidence', 0.700)
            # CHANGE: Use scope='all' to generate the full complex network
            self.build_network(confidence=confidence, scope="all")
            
            # Step 5: Enrichment
            enrichment_method = self.config.get('enrichment', {}).get('method', 'gprofiler')
            self.enrichment_analysis(method=enrichment_method)
            
            # Final summary
            self.logger.info("\n" + "="*70)
            self.logger.info("✅ PIPELINE COMPLETED SUCCESSFULLY!")
            self.logger.info("="*70)
            self.logger.info(f"\nAll results saved to: {self.output_dir}")
            self.logger.info("\nNext steps:")
            self.logger.info("  - Review CSV files in each step directory")
            self.logger.info("  - Open network visualizations")
            self.logger.info("  - Analyze enrichment results")
            self.logger.info("  - Cite this tool in your publication! (MIT License)")
            
        except SystemExit:
            self.logger.error("\n❌ Pipeline stopped due to error")
            self.logger.error("Check log file for details")
            raise
        except Exception as e:
            self.logger.error(f"\n❌ Unexpected error: {str(e)}")
            raise
