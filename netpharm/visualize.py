"""
Network visualization functions.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import os
from pyvis.network import Network as PyVisNetwork


class NetworkVisualizer:
    """Handle network visualizations."""
    
    def __init__(self, logger):
        """
        Initialize NetworkVisualizer.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        sns.set_style("whitegrid")
    
    def create_basic_visualization(self, network, output_path):
        """
        Create basic NetworkX visualization.
        
        Args:
            network: NetworkX graph
            output_path: Path to save image
        """
        self.logger.info("\n--- Creating Basic Network Visualization ---")
        
        plt.figure(figsize=(14, 10))
        
        # Layout
        pos = nx.spring_layout(network, k=2, iterations=50, seed=42)
        
        # Node sizes by degree
        node_degrees = dict(network.degree())
        node_sizes = [v * 150 + 400 for v in node_degrees.values()]
        
        # Edge widths by weight
        edge_weights = [network[u][v]['weight'] * 3 for u, v in network.edges()]
        
        # Draw
        nx.draw_networkx_nodes(network, pos,
                              node_size=node_sizes,
                              node_color='lightblue',
                              edgecolors='darkblue',
                              linewidths=2,
                              alpha=0.9)
        
        nx.draw_networkx_edges(network, pos,
                              width=edge_weights,
                              alpha=0.5,
                              edge_color='gray')
        
        nx.draw_networkx_labels(network, pos,
                               font_size=10,
                               font_weight='bold')
        
        plt.title('Protein-Protein Interaction Network\n(STRING confidence ≥ 0.700)',
                 fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        self.logger.info(f"✓ Basic visualization saved: {output_path}")
        self.logger.info("  Pros: Simple, clean, publication-ready")
        self.logger.info("  Cons: Static, limited customization")
    
    def create_colored_visualization(self, network, output_path):
        """
        Create colored visualization by node degree.
        
        Args:
            network: NetworkX graph
            output_path: Path to save image
        """
        self.logger.info("\n--- Creating Colored Network Visualization ---")
        
        plt.figure(figsize=(16, 12))
        
        # Layout
        pos = nx.spring_layout(network, k=2.5, iterations=50, seed=42)
        
        # Node colors by degree
        node_degrees = dict(network.degree())
        node_colors = [node_degrees[node] for node in network.nodes()]
        
        # Node sizes
        node_sizes = [v * 150 + 400 for v in node_degrees.values()]
        
        # Edge widths
        edge_weights = [network[u][v]['weight'] * 3 for u, v in network.edges()]
        
        # Draw nodes
        nodes = nx.draw_networkx_nodes(network, pos,
                                      node_size=node_sizes,
                                      node_color=node_colors,
                                      cmap=plt.cm.YlOrRd,
                                      edgecolors='black',
                                      linewidths=2,
                                      alpha=0.9)
        
        # Draw edges
        nx.draw_networkx_edges(network, pos,
                              width=edge_weights,
                              alpha=0.5,
                              edge_color='gray')
        
        # Draw labels
        nx.draw_networkx_labels(network, pos,
                               font_size=9,
                               font_weight='bold')
        
        # Colorbar
        plt.colorbar(nodes, label='Node Degree (Number of Connections)')
        
        plt.title('Protein-Protein Interaction Network\nColored by Node Degree',
                 fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        self.logger.info(f"✓ Colored visualization saved: {output_path}")
        self.logger.info("  Pros: Shows hub proteins clearly, aesthetic")
        self.logger.info("  Cons: Static, color scheme may not suit all publications")
    
    def create_interactive_visualization(self, network, output_path):
        """
        Create interactive HTML visualization using PyVis.
        
        Args:
            network: NetworkX graph
            output_path: Path to save HTML file
        """
        self.logger.info("\n--- Creating Interactive Network Visualization ---")
        
        # Create PyVis network
        net = PyVisNetwork(height='750px', width='100%', bgcolor='#ffffff', font_color='black')
        
        # Convert NetworkX to PyVis
        net.from_nx(network)
        
        # Customize appearance
        net.set_options("""
        {
          "nodes": {
            "font": {"size": 14, "face": "Arial"},
            "scaling": {"min": 10, "max": 30}
          },
          "edges": {
            "color": {"inherit": true},
            "smooth": {"type": "continuous"}
          },
          "physics": {
            "barnesHut": {
              "gravitationalConstant": -30000,
              "centralGravity": 0.3,
              "springLength": 95
            },
            "minVelocity": 0.75
          }
        }
        """)
        
        # Save
        net.save_graph(output_path)
        
        self.logger.info(f"✓ Interactive visualization saved: {output_path}")
        self.logger.info("  Pros: Interactive, zoomable, can explore connections")
        self.logger.info("  Cons: Requires browser, not suitable for print publications")
        self.logger.info(f"  💡 Open in browser: file://{os.path.abspath(output_path)}")
    
    def create_all_visualizations(self, network, output_dir):
        """
        Create all three types of visualizations.
        
        Args:
            network: NetworkX graph
            output_dir: Output directory
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("CREATING NETWORK VISUALIZATIONS")
        self.logger.info("="*70)
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Basic
        basic_path = os.path.join(output_dir, "network_basic.png")
        self.create_basic_visualization(network, basic_path)
        
        # Colored
        colored_path = os.path.join(output_dir, "network_colored.png")
        self.create_colored_visualization(network, colored_path)
        
        # Interactive
        interactive_path = os.path.join(output_dir, "network_interactive.html")
        self.create_interactive_visualization(network, interactive_path)
        
        self.logger.info("\n✓ All visualizations created successfully!")
        self.logger.info("\n📊 Visualization Summary:")
        self.logger.info(f"  1. Basic (PNG): {basic_path}")
        self.logger.info(f"  2. Colored (PNG): {colored_path}")
        self.logger.info(f"  3. Interactive (HTML): {interactive_path}")
