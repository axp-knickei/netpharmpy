"""
Network visualization utilities for netpharmpy.

This module is responsible ONLY for visualization.
All biological or statistical decisions are assumed
to be made upstream (core / analysis).
"""

import os
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D


class NetworkVisualizer:
    """
    Visualize STRING protein–protein interaction networks.

    Visualization philosophy (fixed):
    - Plot a biologically meaningful subgraph
    - Hubs + first neighbors only
    - Remove isolated or weak connectors
    - Node size reflects degree *within plotted subgraph*
    - Edge width reflects STRING confidence
    """

    def __init__(self, output_dir, logger=None):
        self.output_dir = output_dir
        self.logger = logger

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def create_all_visualizations(
        self,
        G,
        hub_nodes,
        top_n=15,
        label_connectors=False,
    ):
        """
        Entry point called by core.py
        """

        if self.logger:
            self.logger.info("\nCREATING NETWORK VISUALIZATIONS")

        self._plot_hubs_and_connectors(
            G,
            hub_nodes=hub_nodes[:top_n],
            label_connectors=label_connectors,
        )

    # ------------------------------------------------------------------
    # Core plotting logic
    # ------------------------------------------------------------------

    def _plot_hubs_and_connectors(
        self,
        G,
        hub_nodes,
        label_connectors=False,
    ):
        """
        Plot hubs + first neighbors (connectors),
        enforcing strict graph consistency.
        """

        # --------------------------------------------------------------
        # 1. Build node set: hubs + first neighbors
        # --------------------------------------------------------------
        nodes_to_plot = set(hub_nodes)
        for h in hub_nodes:
            nodes_to_plot.update(G.neighbors(h))

        # Induced subgraph
        G_plot = G.subgraph(nodes_to_plot).copy()

        # --------------------------------------------------------------
        # 2. Remove weak connectors (degree < 2 inside subgraph)
        # --------------------------------------------------------------
        for node in list(G_plot.nodes()):
            if node not in hub_nodes and G_plot.degree(node) < 2:
                G_plot.remove_node(node)

        if self.logger:
            self.logger.info(
                f"Hubs: {len(hub_nodes)}, "
                f"Connectors: {len(G_plot.nodes()) - len(hub_nodes)}"
            )

        # --------------------------------------------------------------
        # 3. Layout (ONLY on G_plot)
        # --------------------------------------------------------------
        pos = nx.spring_layout(G_plot, seed=42)

        # --------------------------------------------------------------
        # 4. Node sizes = degree in G_plot
        # --------------------------------------------------------------
        node_sizes = [
            300 + 80 * G_plot.degree(n)
            for n in G_plot.nodes()
        ]

        # --------------------------------------------------------------
        # 5. Node colors (hubs vs connectors)
        # --------------------------------------------------------------
        node_colors = [
            "#d62728" if n in hub_nodes else "#1f77b4"
            for n in G_plot.nodes()
        ]

        # --------------------------------------------------------------
        # 6. Edge widths = normalized STRING confidence
        # --------------------------------------------------------------
        weights = [
            G_plot[u][v].get("weight", 1.0)
            for u, v in G_plot.edges()
        ]

        if len(weights) > 0:
            w_min, w_max = min(weights), max(weights)
            edge_widths = [
                1.5 + 4.0 * (w - w_min) / (w_max - w_min + 1e-6)
                for w in weights
            ]
        else:
            edge_widths = []

        # --------------------------------------------------------------
        # 7. Plot
        # --------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(14, 14))

        nx.draw_networkx_edges(
            G_plot,
            pos,
            ax=ax,
            width=edge_widths,
            edge_color="gray",
            alpha=0.7,
        )

        nx.draw_networkx_nodes(
            G_plot,
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=node_colors,
            edgecolors="black",
            linewidths=1.2,
        )

        # --------------------------------------------------------------
        # 8. Labels
        # --------------------------------------------------------------
        labels = {n: n for n in hub_nodes}
        if label_connectors:
            for n in G_plot.nodes():
                if n not in hub_nodes:
                    labels[n] = n

        nx.draw_networkx_labels(
            G_plot,
            pos,
            labels=labels,
            font_size=10,
            font_weight="bold",
        )

        # --------------------------------------------------------------
        # 9. Legends & annotations
        # --------------------------------------------------------------
        legend_elements = [
            Line2D([0], [0], marker="o", color="w",
                   label="Hub proteins",
                   markerfacecolor="#d62728", markersize=12,
                   markeredgecolor="black"),
            Line2D([0], [0], marker="o", color="w",
                   label="Connector proteins",
                   markerfacecolor="#1f77b4", markersize=12,
                   markeredgecolor="black"),
            Line2D([0], [0], color="gray",
                   lw=2, label="Protein–protein interaction"),
        ]

        ax.legend(
            handles=legend_elements,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            title="Network elements",
        )

        ax.set_title(
            "Protein–Protein Interaction Network\n"
            "(Hubs + First Neighbors)",
            fontsize=16,
            pad=20,
        )

        ax.text(
            0.01,
            0.01,
            "Hubs = top-degree proteins in STRING network\n"
            "Connectors = first neighbors of hubs\n"
            "Node size ∝ degree in plotted subgraph\n"
            "Edge width ∝ STRING confidence score",
            transform=ax.transAxes,
            fontsize=9,
            va="bottom",
            ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9),
        )

        ax.axis("off")

        # --------------------------------------------------------------
        # 10. Save
        # --------------------------------------------------------------
        out_path = os.path.join(
            self.output_dir,
            "network_static_hubs_connectors.png",
        )
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()

        if self.logger:
            self.logger.info(f"✓ Saved network visualization: {out_path}")
