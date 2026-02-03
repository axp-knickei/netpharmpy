"""
Static network visualization for publication-quality figures.
"""

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import os
import numpy as np

class StaticVisualizer:
    """Generates high-resolution static network plots using Matplotlib."""

    def __init__(self, output_dir):
        """
        Initialize StaticVisualizer.

        Args:
            output_dir: Directory to save figures
        """
        self.output_dir = output_dir
        self.figures_dir = os.path.join(output_dir, "figures_publication")
        os.makedirs(self.figures_dir, exist_ok=True)

    def plot_network_hubs(self, G, hub_nodes=None, title="Target-Interaction Network"):
        """
        Plot network highlighting hub nodes.

        Args:
            G: NetworkX graph
            hub_nodes: List of node names identified as hubs
            title: Plot title
        """
        if hub_nodes is None:
            hub_nodes = []

        plt.figure(figsize=(12, 12))
        
        # Calculate layout
        # k parameter controls node spacing (optimum distance)
        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)

        # Node sizes based on degree
        degrees = dict(G.degree())
        # Scale: min size 100, factor 50 per degree
        node_sizes = [300 + (degrees[n] * 50) for n in G.nodes()]

        # Node colors
        # Hubs = Red (#E41A1C), Others = Blue (#377EB8)
        node_colors = []
        for n in G.nodes():
            if n in hub_nodes:
                node_colors.append('#E41A1C')  # Red for hubs
            else:
                node_colors.append('#377EB8')  # Blue for others

        # Draw edges first (so they are behind nodes)
        # Weight thickness by edge weight if available
        weights = [G[u][v].get('weight', 0.5) * 2 for u, v in G.edges()]
        
        nx.draw_networkx_edges(
            G, pos, 
            alpha=0.3, 
            edge_color='gray',
            width=weights
        )

        # Draw nodes
        nx.draw_networkx_nodes(
            G, pos,
            node_size=node_sizes,
            node_color=node_colors,
            alpha=0.9,
            linewidths=1,
            edgecolors='white'  # White border around nodes
        )

        # Draw Labels
        # Only label hubs to avoid clutter, or all if network is small (<50 nodes)
        if len(G.nodes()) < 50:
            labels_to_draw = {n: n for n in G.nodes()}
        else:
            labels_to_draw = {n: n for n in G.nodes() if n in hub_nodes}

        # Draw labels with a slight offset or box
        text_items = nx.draw_networkx_labels(
            G, pos,
            labels=labels_to_draw,
            font_size=10,
            font_weight='bold',
            font_color='black',
            font_family='sans-serif'
        )
        
        # Add a white halo to text for readability
        import matplotlib.patheffects as path_effects
        for t in text_items.values():
            t.set_path_effects([path_effects.withStroke(linewidth=3, foreground='white')])

        # Legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='Hub Targets',
                   markerfacecolor='#E41A1C', markersize=15),
            Line2D([0], [0], marker='o', color='w', label='Other Targets',
                   markerfacecolor='#377EB8', markersize=10)
        ]
        plt.legend(handles=legend_elements, loc='upper right', frameon=True)

        plt.title(title, fontsize=16, fontweight='bold')
        plt.axis('off')  # Hide axis
        
        # Save
        filename = "network_hubs_static.png"
        save_path = os.path.join(self.figures_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        # Also save PDF for vector graphics (journals often prefer this)
        pdf_path = save_path.replace('.png', '.pdf')
        plt.savefig(pdf_path, format='pdf', bbox_inches='tight')

        return save_path

    def plot_centrality_distribution(self, metrics_df):
        """
        Plot bar chart of top nodes by centrality.

        Args:
            metrics_df: DataFrame from NetworkAnalyzer.analyze_network()
        """
        top_n = metrics_df.head(15).copy()
        top_n = top_n.sort_values('Degree', ascending=True) # Sort for horiz bar chart

        plt.figure(figsize=(10, 8))
        
        # Create gradient colors based on values
        colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_n)))

        bars = plt.barh(
            top_n['Protein'], 
            top_n['Degree'],
            color=colors,
            alpha=0.8
        )

        plt.xlabel('Degree (Number of Interactions)', fontsize=12)
        plt.title('Top 15 Hub Proteins by Degree Centrality', fontsize=14, fontweight='bold')
        plt.grid(axis='x', linestyle='--', alpha=0.5)

        # Add value labels to bars
        for bar in bars:
            width = bar.get_width()
            plt.text(
                width + 0.5, 
                bar.get_y() + bar.get_height()/2, 
                f'{int(width)}',
                ha='left', va='center', fontsize=10
            )

        save_path = os.path.join(self.figures_dir, "top_hubs_bar_chart.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return save_path
