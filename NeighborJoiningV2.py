import matplotlib.pyplot as plt
import networkx as nx
#import pygraphviz
import netgraph
import numpy as np

import NeighborJoining

m= np.array([[0,5,4,7,6,8],[5,0,7,10,9,11], [4, 7,0,7,6,8], [7, 10, 7,0,5,9], [6, 9, 6, 5,0,8],  [8, 11, 8, 9, 8,0]])


def NeighborJoining2(matrix):
    n = len(matrix)
    clusters = [[i] for i in range(n)]

    G = {}
    nodes = {}
    for i in range(len(m)+1):
        G[(i, len(m)+2)] = 1
        nodes[i] = i

    while n > 2:
        div_matrix = NeighborJoining.divergence_matrix(matrix)
        i, j = NeighborJoining.smallestOTU(div_matrix)
        new_cluster = clusters[i] + clusters[j]
        clusters = [clusters[k] for k in range(len(clusters)) if k != i and k != j]
        clusters.append(new_cluster)
        new_matrix = NeighborJoining.newDistMatrix(matrix, i, j)
        matrix = new_matrix
        n -= 1

    edges = list(G.keys())
    fig, ax = plt.subplots()
    netgraph.Graph(edges, edge_labels=G, node_layout='geometric',
                   node_layout_kwargs=dict(edge_length=G), ax=ax, node_labels=nodes)
    ax.set_aspect('equal')
    plt.show()

    return clusters

NeighborJoining2(m)
#
# edge_length = {
#     (0, 1) : 0.3,
#     (1, 2) : 0.4,
#     (2, 0) : 0.5,
# }
#
# edges = list(edge_length.keys())
#
# nodes = {0:1,1:2,2:''}
#
# fig, ax = plt.subplots()
# netgraph.Graph(edges, edge_labels=edge_length, node_layout='geometric',
#       node_layout_kwargs=dict(edge_length=edge_length), ax=ax, node_labels=nodes)
# ax.set_aspect('equal')
# plt.show()