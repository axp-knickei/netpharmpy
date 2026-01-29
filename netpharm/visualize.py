# netpharm/visualize.py
import os
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
from pyvis.network import Network as PyvisNetwork
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


class NetworkVisualizer:
    def __init__(self, output_dir):
        if hasattr(output_dir, "info"):
            raise TypeError(
                "NetworkVisualizer received a Logger instead of a path. "
                "Pass output_dir, not logger."
            )
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------
    def _get_hubs_and_connectors(self, network: nx.Graph, top_n: int = 15):
        """
        Return an induced subgraph containing top_n hubs (by global degree)
        and their first neighbors, filtered so connectors shown have >1 edge
        inside the induced subgraph.

        Returns:
            subgraph (nx.Graph), global_degree_dict (dict), roles (dict)
        """
        # global degrees (used for ranking hubs)
        global_degree = dict(network.degree())
        hubs = sorted(global_degree, key=global_degree.get, reverse=True)[:top_n]

        # all first neighbors of hubs
        connectors = set()
        for h in hubs:
            connectors.update(network.neighbors(h))

        # initial candidate nodes (hubs + connectors)
        nodes_to_keep = set(hubs) | connectors

        # induced subgraph
        subG = network.subgraph(nodes_to_keep).copy()

        # Now filter connectors by degree *inside the subgraph*
        valid_connectors = {
            n for n in subG.nodes() if n not in hubs and subG.degree(n) >= 2
        }

        nodes_final = set(hubs) | valid_connectors
        subG = network.subgraph(nodes_final).copy()

        # remove isolates if any remain
        isolates = list(nx.isolates(subG))
        if isolates:
            subG.remove_nodes_from(isolates)

        # roles
        roles = {}
        for n in subG.nodes():
            roles[n] = "hub" if n in hubs else "connector"

        return subG, global_degree, roles

    # ------------------------------------------------------------------
    # Static visualization
    # ------------------------------------------------------------------
    def create_static_network(
        self,
        network: nx.Graph,
        filename: str = "network_static_hubs_connectors.png",
        top_n: int = 15,
        label_connectors: bool = False,
        edge_width_range: tuple = (0.5, 4.0),
        node_size_range: tuple = (300, 1400),
        figsize: tuple = (10, 10),
        layout_seed: int = 42,
    ):
        """
        Produce a publication-ready static PNG: hubs (red) + connectors (blue).
        - Only connectors with degree >= 2 inside the subgraph are kept.
        - Node sizes reflect degree *inside* the plotted subgraph.
        - Edge widths are normalized from edge 'weight' or 'score' attribute.
        """
        G, global_deg, roles = self._get_hubs_and_connectors(network, top_n=top_n)
        if G.number_of_nodes() == 0:
            raise ValueError("Subgraph is empty (no hubs/connectors found).")

        # layout (reproducible)
        pos = nx.spring_layout(G, seed=layout_seed)

        # subgraph degree (used for node size)
        subdeg = dict(G.degree())

        # node size scaling (map degrees -> node_size_range)
        degs = list(subdeg.values())
        dmin, dmax = (min(degs), max(degs)) if degs else (0, 1)
        min_size, max_size = node_size_range

        def scale_node(d):
            if dmax == dmin:
                return (min_size + max_size) / 2.0
            return min_size + (d - dmin) / (dmax - dmin) * (max_size - min_size)

        node_sizes = [scale_node(subdeg[n]) for n in G.nodes()]

        # node colors by role
        node_colors = ["#d62728" if roles[n] == "hub" else "#1f77b4" for n in G.nodes()]

        # edge widths: gather raw weight values
        raw_weights = []
        for u, v, data in G.edges(data=True):
            # try common names for strength
            w = data.get("weight", None)
            if w is None:
                w = data.get("score", None)
            if w is None:
                w = 1.0
            raw_weights.append(float(w))

        # normalize weights into a visible width range
        if raw_weights:
            wmin, wmax = min(raw_weights), max(raw_weights)
            if wmax == wmin:
                edge_widths = [ (edge_width_range[0] + edge_width_range[1]) / 2.0 ] * G.number_of_edges()
            else:
                edge_widths = [
                    edge_width_range[0] + (w - wmin) / (wmax - wmin) * (edge_width_range[1] - edge_width_range[0])
                    for w in raw_weights
                ]
        else:
            edge_widths = [1.0] * G.number_of_edges()

        # === Drawing: edges first (so nodes sit on top) ===
        fig, ax = plt.subplots(figsize=figsize)

        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            width=edge_widths,
            edge_color="gray",
            alpha=0.55,
            zorder=1,
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=node_colors,
            edgecolors="black",
            linewidths=1.2,
            zorder=2,
        )

        # labels: hubs always; connectors optionally if they have degree >= 2 in subgraph
        if label_connectors:
            labels = {
                n: n for n in G.nodes() if roles[n] == "hub" or (roles[n] == "connector" and subdeg[n] >= 2)
            }
        else:
            labels = {n: n for n in G.nodes() if roles[n] == "hub"}

        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_weight="bold", zorder=3)

        ax.set_title("Protein–Protein Interaction Network\n(Hubs + First Neighbors)")
        ax.axis("off")

        # Legend / annotation
        legend_elements = [
            Patch(facecolor="#d62728", edgecolor="black", label="Hub proteins"),
            Patch(facecolor="#1f77b4", edgecolor="black", label="Connector proteins"),
            Line2D([0], [0], color="gray", lw=max(edge_widths) if edge_widths else 2, label="Protein–protein interaction"),
        ]
        ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5))

        # textual annotation box
        annotation_text = (
            "Hubs = top-degree proteins in STRING network\n"
            "Connectors = first neighbors of hubs (shown only if ≥2 edges in subgraph)\n"
            "Node size ∝ degree (inside plotted subgraph)\n"
            "Edge width ∝ STRING confidence score (normalized)"
        )
        ax.text(
            0.01,
            0.01,
            annotation_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
        )

        # summary counts
        num_hubs = sum(1 for r in roles.values() if r == "hub")
        num_connectors = sum(1 for r in roles.values() if r == "connector")
        ax.text(
            0.99,
            0.01,
            f"Hubs: {num_hubs}\nConnectors: {num_connectors}",
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            horizontalalignment="right",
        )

        out = os.path.join(self.output_dir, filename)
        plt.tight_layout()
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return out

    # ------------------------------------------------------------------
    # Interactive visualization (pyvis)
    # ------------------------------------------------------------------
    def create_interactive_network(self, network: nx.Graph, filename: str = "network_interactive.html"):
        G, _, roles = self._get_hubs_and_connectors(network, top_n=15)

        net = PyvisNetwork(height="800px", width="100%", bgcolor="#ffffff", notebook=False)
        net.force_atlas_2based()

        # add nodes with sizes from subgraph degree
        subdeg = dict(G.degree())
        min_size, max_size = 10, 30
        dmin, dmax = (min(subdeg.values()), max(subdeg.values())) if subdeg else (0,1)

        def scale_vis_size(d):
            if dmax == dmin:
                return (min_size + max_size) / 2.0
            return min_size + (d - dmin) / (dmax - dmin) * (max_size - min_size)

        for n in G.nodes():
            size = scale_vis_size(subdeg.get(n, 1))
            color = "#d62728" if roles[n] == "hub" else "#1f77b4"
            net.add_node(n, label=n, size=size, color=color)

        # add edges (use 'weight' or 'score' if available)
        for u, v, d in G.edges(data=True):
            w = d.get("weight", d.get("score", 1.0))
            net.add_edge(u, v, value=float(w))

        out = os.path.join(self.output_dir, filename)
        net.write_html(out)
        return out

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------
    def create_all_visualizations(self, network: nx.Graph, top_n: int = 15):
        # static (PNG)
        self.create_static_network(network, filename="network_static_hubs_connectors.png", top_n=top_n)
        # interactive (HTML)
        self.create_interactive_network(network, filename="network_interactive.html")
