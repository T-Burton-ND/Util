import networkx as nx
from networkx.algorithms.similarity import graph_edit_distance

G = nx.path_graph(5)
H = nx.path_graph(4)

res = graph_edit_distance(G, H, timeout=0.1)
try:
    next(res)  # advance once if it's a generator
except TypeError:
    pass
