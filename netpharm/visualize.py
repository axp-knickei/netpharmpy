import os
import matplotlib.pyplot as plt
import networkx as nx
from pyvis.network import Network
from pathlib import Path

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
    # Core graph abstraction
    # ------------------------------------------------------------------
    def _get_hubs_and_connectors(self, network, top_n=15):
        """
        Build a hub-centered subgraph:
        - hubs = top_n nodes by global degree
        - connectors = all first neighbors of hubs
        """

        degree_dict = dict(network.degree())
        hubs = sorted(degree_dict, key=degree_dict.get, reverse=True)[:top_n]

        connectors = set()
        for h in hubs:
            connectors.update(network.neighbors(h))

        nodes_to_keep = set(hubs) | connectors
        subgraph = network.subgraph(nodes_to_keep).copy()

        node_roles = {}
        for n in subgraph.nodes():
            if n in hubs:
                node_roles[n] = "hub"
            else:
                node_roles[n] = "connector"

        return subgraph, degree_dict, node_roles

    # ------------------------------------------------------------------
    # Static visualization (summary)
    # ------------------------------------------------------------------
    def create_static_network(self, network, filename, top_n=15):
        LEGEND_STYLE = "outside"  # options: outside, inside, bottom

        G, degree_dict, roles = self._get_hubs_and_connectors(network, top_n)

        fig, ax = plt.subplots(figsize=(10, 10))
        pos = nx.spring_layout(G, seed=42)

        # ---- Node sizes (global degree, normalized) ----
        degrees = [degree_dict[n] for n in G.nodes()]
        min_size, max_size = 400, 2000
        dmin, dmax = min(degrees), max(degrees)

        def scale(d):
            if dmax == dmin:
                return min_size
            return min_size + (d - dmin) / (dmax - dmin) * (max_size - min_size)

        node_sizes = [scale(degree_dict[n]) for n in G.nodes()]

        # ---- Node colors (role-based) ----
        node_colors = []
        for n in G.nodes():
            if roles[n] == "hub":
                node_colors.append("#d62728")  # red
            else:
                node_colors.append("#1f77b4")  # blue

        # ---- Edges ----
        edge_widths = [G[u][v].get("weight", 1) * 2 for u, v in G.edges()]

        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            width=edge_widths,
            edge_color="gray",
            alpha=0.6,
        )

        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=node_colors,
            edgecolors="black",
            linewidths=1.5,
        )

        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor="#d62728", edgecolor="black", label="Hub proteins"),
            Patch(facecolor="#1f77b4", edgecolor="black", label="Connector proteins"),
            Line2D([0], [0], color="gray", lw=2, label="Protein–protein interaction"),
        ]

        if LEGEND_STYLE == "outside":
            ax.legend(
                handles=legend_elements,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
                title="Network elements",
                fontsize=9,
                title_fontsize=10
            )
        elif LEGEND_STYLE == "inside":
            ax.legend(
                handles=legend_elements,
                loc="upper left",
                frameon=True,
                fontsize=9,
                title_fontsize=10
            )
        elif LEGEND_STYLE == "bottom":
            ax.legend(
                handles=legend_elements,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.05),
                ncol=3,
                frameon=False,
                fontsize=9,
                title_fontsize=10
            )



        # ---- Labels (hubs only) ----
        hub_labels = {n: n for n in G.nodes() if roles[n] == "hub"}
        nx.draw_networkx_labels(
            G,
            pos,
            labels=hub_labels,
            font_size=10,
            font_weight="bold",
        )

        ax.set_title("Protein–Protein Interaction Network\n(Hubs + First Neighbors)")
        ax.axis("off")

        annotation_text = (
            "Hubs = top-degree proteins in STRING network\n"
            "Connectors = first neighbors of hubs\n"
            "Node size ∝ global degree\n"
            "Edge width ∝ STRING confidence score"
        )

        ax.text(
            0.01, 0.01,
            annotation_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
        )

        num_hubs = sum(1 for r in roles.values() if r == "hub")
        num_connectors = sum(1 for r in roles.values() if r == "connector")

        ax.text(
            0.99, 0.01,
            f"Hubs: {num_hubs}\nConnectors: {num_connectors}",
            transform=ax.transAxes,
            fontsize=9,
            ha="right",
            va="bottom"
        )

        out = os.path.join(self.output_dir, filename)
        plt.tight_layout()
        plt.savefig(out, dpi=300)
        plt.close()

        return out

    # ------------------------------------------------------------------
    # Interactive visualization (full network)
    # ------------------------------------------------------------------
    def create_interactive_network(self, network, filename="network_interactive.html"):
        net = Network(height="750px", width="100%", bgcolor="white")

        for n in network.nodes():
            net.add_node(n, label=n)

        for u, v, d in network.edges(data=True):
            w = d.get("weight", 1)
            net.add_edge(u, v, value=w)

        out = os.path.join(self.output_dir, filename)
        net.write_html(out)
        return out

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------
    def create_all_visualizations(self, network):
        self.create_static_network(
            network,
            filename="network_static_hubs_connectors.png",
            top_n=15,
        )
        self.create_interactive_network(network)
