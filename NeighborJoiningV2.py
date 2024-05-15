import sys

import matplotlib.pyplot as plt
import networkx as nx
#import pygraphviz
import netgraph
import numpy as np

import NeighborJoining

m= np.array([[0,5,4,7,6,8],[5,0,7,10,9,11], [4, 7,0,7,6,8], [7, 10, 7,0,5,9], [6, 9, 6, 5,0,8],  [8, 11, 8, 9, 8,0]])

def max_index_neighbor(graph, node):
    neighbors = list(graph.neighbors(node))
    if neighbors:
        max_index_node = max(neighbors)
        return max_index_node
    else:
        return None

def NeighborJoining2(matrix):
    safe = len(matrix)
    n = len(matrix)
    clusters = [[i] for i in range(n)]

    G = nx.Graph()
    nodes = np.arange(len(matrix) + 1)
    counter = n

    for i in range(len(matrix)+1):
        G.add_edge(i, sys.maxsize, weight = 1)



    for i in range(5):
    #while n > 2:
        div_matrix = NeighborJoining.divergence_matrix(matrix)
        i, j = NeighborJoining.smallestOTU(div_matrix)
        new_cluster = clusters[i] + clusters[j]

        counter+=1
        print("combine " + str(i) + " and " + str(j) + " to " + str(counter))
        print("delete (" + str(nodes[i]) + ", " + str(max_index_neighbor(G, nodes[i])) + ")")
        print("delete (" + str(nodes[j]) + ", " + str(max_index_neighbor(G, nodes[j])) + ")")
        print(nodes)

        x = max_index_neighbor(G, nodes[i])
        G.remove_edge(nodes[i], x)
        G.remove_edge(nodes[j], max_index_neighbor(G, nodes[j]))

        #G.add_node(counter)
        G.add_edge(nodes[i], counter, weight = 1)
        G.add_edge(nodes[j], counter, weight=1)
        G.add_edge(counter, x, weight=1)

        nodes = np.delete(nodes,[i,j])
        print(nodes)
        nodes = np.append(nodes, counter)
        print(nodes)

        clusters = [clusters[k] for k in range(len(clusters)) if k != i and k != j]
        clusters.append(new_cluster)
        new_matrix = NeighborJoining.newDistMatrix(matrix, i, j)
        matrix = new_matrix
        n -= 1

    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.gca()
    ax.axis("off")
    # Calculate position of nodes in the plot
    pos = nx.kamada_kawai_layout(G)
    # Assign the gene names to the nodes that represent a reference index
    #node_labels = {i: name for i, name in enumerate(genes)}
    labels = {x: x for x in range(safe + 1)}

    nx.draw_networkx_edges(
        G, pos, ax=ax
    )
    nx.draw_networkx_labels(
        G, pos, ax=ax, labels= labels , font_size=10,
        # Draw a white background behind the labeled nodes
        # for better readability
        bbox=dict(pad=0, color="white")
    )
    fig.tight_layout()

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