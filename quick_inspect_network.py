# quick_inspect_network.py
from netpharm.network import build_string_network

# load any previously generated STRING interactions CSV
# adjust path if needed
G = build_string_network(
    "./outputs/compound_969516_*/step4_network/string_interactions.csv"
)

print(type(G))
print("Is directed:", G.is_directed())
print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())

# inspect one edge
u, v, data = next(iter(G.edges(data=True)))
print("Example edge:", u, v, data)

# inspect one node
n = next(iter(G.nodes()))
print("Example node:", n)
