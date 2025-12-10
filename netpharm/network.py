"""
STRING protein-protein interaction network analysis.
"""

import pandas as pd
import networkx as nx
import os
from .utils.api_wrappers import query_string


class NetworkAnalyzer:
    """Handle STRING network analysis."""
    
    def __init__(self, logger):
        """
        Initialize NetworkAnalyzer.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.interactions_df = None
        self.network = None
        self.network_metrics = None
    
    def query_string_network(self, gene_list, confidence=0.700):
        """
        Query STRING database for protein interactions.
        
        Args:
            gene_list: List of gene names
            confidence: Minimum confidence score (0.0-1.0)
        
        Returns:
            pd.DataFrame: Interaction data
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 4] STRING PROTEIN-PROTEIN INTERACTION NETWORK")
        self.logger.info("="*70)
        
        if not gene_list:
            self.logger.error("\n❌ ERROR: No genes provided for network analysis")
            raise SystemExit(1)
        
        self.logger.info(f"\nQuerying STRING for {len(gene_list)} proteins...")
        self.logger.info(f"Confidence threshold: {confidence} (score ≥ {int(confidence*1000)}/1000)")
        
        try:
            # Query STRING
            required_score = int(confidence * 1000)
            response_text = query_string(gene_list, species=9606, required_score=required_score)
            
            # Parse TSV response
            lines = response_text.strip().split('\n')
            data = [line.split('\t') for line in lines]
            
            self.interactions_df = pd.DataFrame(data[1:], columns=data[0])
            
            # Convert score to numeric
            self.interactions_df['score'] = pd.to_numeric(
                self.interactions_df['score'], 
                errors='coerce'
            )
            
            # Clean protein names (remove species prefix)
            self.interactions_df['preferredName_A'] = self.interactions_df['preferredName_A'].str.replace('9606.', '', regex=False)
            self.interactions_df['preferredName_B'] = self.interactions_df['preferredName_B'].str.replace('9606.', '', regex=False)
            
            self.logger.info(f"\n✓ Retrieved {len(self.interactions_df)} interactions")
            
            if len(self.interactions_df) == 0:
                self.logger.error("\n❌ ERROR: No interactions found at this confidence level")
                self.logger.error("💡 Suggestions:")
                self.logger.error("   - Lower the confidence threshold")
                self.logger.error("   - Check if gene names are correct")
                self.logger.error("   - Try with more genes")
                raise SystemExit(1)
            
            return self.interactions_df
            
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: STRING query failed")
            self.logger.error(f"   Reason: {str(e)}")
            raise SystemExit(1)
    
    def build_network(self):
        """
        Build NetworkX graph from interactions.
        
        Returns:
            nx.Graph: Network graph
        """
        if self.interactions_df is None:
            raise ValueError("Must query STRING first")
        
        self.logger.info("\nBuilding protein interaction network...")
        
        self.network = nx.Graph()
        
        for _, row in self.interactions_df.iterrows():
            protein1 = row['preferredName_A']
            protein2 = row['preferredName_B']
            score = float(row['score'])
            
            self.network.add_edge(protein1, protein2, weight=score/1000)
        
        n_nodes = self.network.number_of_nodes()
        n_edges = self.network.number_of_edges()
        
        self.logger.info(f"✓ Network built successfully")
        self.logger.info(f"  Nodes (proteins): {n_nodes}")
        self.logger.info(f"  Edges (interactions): {n_edges}")
        
        if n_nodes == 0:
            self.logger.error("\n❌ ERROR: Empty network")
            raise SystemExit(1)
        
        return self.network
    
    def analyze_network(self):
        """
        Analyze network topology and identify hub proteins.
        
        Returns:
            pd.DataFrame: Network metrics
        """
        if self.network is None:
            raise ValueError("Must build network first")
        
        self.logger.info("\nAnalyzing network topology...")
        
        # Calculate centrality measures
        degree_centrality = nx.degree_centrality(self.network)
        betweenness_centrality = nx.betweenness_centrality(self.network)
        closeness_centrality = nx.closeness_centrality(self.network)
        
        # Create metrics dataframe
        metrics_data = []
        for node in self.network.nodes():
            metrics_data.append({
                'Protein': node,
                'Degree': self.network.degree(node),
                'Degree_Centrality': degree_centrality[node],
                'Betweenness_Centrality': betweenness_centrality[node],
                'Closeness_Centrality': closeness_centrality[node]
            })
        
        self.network_metrics = pd.DataFrame(metrics_data)
        self.network_metrics = self.network_metrics.sort_values('Degree', ascending=False)
        
        # Identify hub proteins (top 10)
        top_hubs = self.network_metrics.head(10)
        
        self.logger.info("\n" + "="*70)
        self.logger.info("NETWORK HUB PROTEINS (Top 10):")
        self.logger.info("="*70)
        self.logger.info("Hub proteins are highly connected and may be key therapeutic targets.\n")
        
        for idx, row in top_hubs.iterrows():
            self.logger.info(f"  {row['Protein']}: {row['Degree']} connections "
                           f"(centrality: {row['Degree_Centrality']:.3f})")
        
        # Network statistics
        if self.network.number_of_nodes() > 1:
            avg_degree = sum(dict(self.network.degree()).values()) / self.network.number_of_nodes()
            
            self.logger.info("\nNetwork Statistics:")
            self.logger.info(f"  Average degree: {avg_degree:.2f}")
            self.logger.info(f"  Network density: {nx.density(self.network):.3f}")
            
            if nx.is_connected(self.network):
                self.logger.info(f"  Network diameter: {nx.diameter(self.network)}")
                self.logger.info(f"  Average shortest path: {nx.average_shortest_path_length(self.network):.2f}")
            else:
                n_components = nx.number_connected_components(self.network)
                self.logger.info(f"  Connected components: {n_components}")
        
        return self.network_metrics
    
    def save_results(self, output_dir):
        """
        Save network analysis results.
        
        Args:
            output_dir: Output directory
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Save interactions
        if self.interactions_df is not None:
            interactions_path = os.path.join(output_dir, "string_interactions.csv")
            self.interactions_df.to_csv(interactions_path, index=False)
            self.logger.info(f"\n📄 STRING interactions saved to: {interactions_path}")
            self.logger.info("   Contains: Protein pairs, interaction scores, evidence types")
        
        # Save network metrics
        if self.network_metrics is not None:
            metrics_path = os.path.join(output_dir, "network_metrics.csv")
            self.network_metrics.to_csv(metrics_path, index=False)
            self.logger.info(f"📄 Network metrics saved to: {metrics_path}")
            self.logger.info("   Contains: Hub analysis, centrality measures for each protein")
