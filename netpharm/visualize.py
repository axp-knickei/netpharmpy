"""
Network visualization functions.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
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

    def _get_top_hub_subgraph(self, network, top_n=25):
        """
        Return a subgraph containing the top N nodes by degree.
        Degree is computed on the full network.
        """
        degree_dict = dict(network.degree())

        top_nodes = sorted(
            degree_dict,
            key=degree_dict.get,
            reverse=True
        )[:top_n]

        core_net = network.subgraph(top_nodes).copy()
        return core_net, degree_dict

    
    def create_basic_visualization(self, network, output_path):
        """
        Create basic NetworkX visualization.
        
        Args:
            network: NetworkX graph
            output_path: Path to save image
        """
        self.logger.info("\n--- Creating Basic Network Visualization ---")
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Layout
        core_net, degree_dict = self._get_top_hub_subgraph(network, top_n=25)
        pos = nx.spring_layout(core_net, k=0.6, iterations=200, seed=42)

        
        # Node sizes by degree
        node_degrees = {n: degree_dict[n] for n in core_net.nodes()}
        ## Raw degree value (data, not visualization)
        raw_degrees = [degree_dict[n] for n in core_net.nodes()]
        
        min_node_size = 300
        max_node_size = 1800

        min_deg = min(raw_degrees)
        max_deg = max(raw_degrees)

        node_sizes = [
            min_node_size
            + (deg - min_deg) / (max_deg - min_deg) * (max_node_size - min_node_size)
            if max_deg > min_deg else min_node_size
            for deg in raw_degrees
        ]
        
        # Edge widths by weight
        edge_weights = [
            core_net[u][v].get('weight', 1) * 3
            for u, v in core_net.edges()
            ]
        
        # Draw nodes
        nodes = nx.draw_networkx_nodes(
            core_net,
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=node_colors,
            cmap=plt.cm.viridis,
            edgecolors='black'
        )    
        
        nx.draw_networkx_edges(
            core_net, 
            pos,
            ax = ax,
            width=edge_weights,
            alpha=0.5,
            edge_color='gray'
            )
        
        hub_labels = {
            n: n for n in node_degrees if node_degrees[n] >= 10
        }

        nx.draw_networkx_labels(
            core_net,
            pos,
            ax=ax,
            labels=hub_labels,
            font_size=9,
            font_weight='bold'
        )
        
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
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Layout
        core_net, degree_dict = self._get_top_hub_subgraph(network, top_n=25)
        pos = nx.spring_layout(core_net, k=0.6, iterations=200, seed=42)
        
        # Node colors by degree
        node_degrees = {n: degree_dict[n] for n in core_net.nodes()}
        
        # Node sizes
        raw_degrees = [degree_dict[n] for n in core_net.nodes()] # The duplication node_degrees and raw_degrees is intentional

        min_node_size = 300
        max_node_size = 1800

        min_deg = min(raw_degrees)
        max_deg = max(raw_degrees)

        node_sizes = [
            min_node_size
            + (deg - min_deg) / (max_deg - min_deg) * (max_node_size - min_node_size)
            if max_deg > min_deg else min_node_size
            for deg in raw_degrees
        ]
        
        # Edge widths
        edge_weights = [
            core_net[u][v].get('weight', 1) * 3 
            for u, v in core_net.edges()
            ]
        
        # Draw nodes

        nodes = nx.draw_networkx_nodes(
            core_net, 
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=node_colors,
            cmap=plt.cm.viridis,
            edgecolors='black',
            linewidths=2,
            alpha=0.9
        )

        
        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.viridis,
            norm = plt.Normalize(
                vmin=min(node_colors),
                vmax=max(node_colors)
            )
        )
        sm.set_array([])

        cbar = fig.colorbar(sm, ax = ax)

        nx.draw_networkx_edges(
            core_net, 
            pos,
            ax=ax,
            width=edge_weights,
            alpha=0.5,
            edge_color='gray'
            )


        
        # Draw labels
        hub_labels = {n: n for n in node_degrees if node_degrees[n] >= 10}

        nx.draw_networkx_labels(
            core_net, 
            pos,
            ax = ax,
            labels = hub_labels,
            font_size=9,
            font_weight='bold'
            )
        
        # Colorbar
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
