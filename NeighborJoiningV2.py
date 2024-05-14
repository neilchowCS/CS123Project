import matplotlib.pyplot as plt
import networkx as nx
#import pygraphviz
import netgraph

edge_length = {
    (0, 1) : 0.3,
    (1, 2) : 0.4,
    (2, 0) : 0.5,
}

edges = list(edge_length.keys())

nodes = {0:1,1:2,2:3}

fig, ax = plt.subplots()
netgraph.Graph(edges, edge_labels=edge_length, node_layout='geometric',
      node_layout_kwargs=dict(edge_length=edge_length), ax=ax, node_labels=nodes)
ax.set_aspect('equal')
plt.show()