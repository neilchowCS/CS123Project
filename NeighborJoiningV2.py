import sys

import matplotlib.pyplot as plt
import networkx as nx
#import pygraphviz
#import netgraph
import numpy as np

from mis import NeighborJoining


def max_index_neighbor(graph, node):
    neighbors = list(graph.neighbors(node))
    if neighbors:
        max_index_node = max(neighbors)
        return max_index_node
    else:
        return None

def divergence_matrix(matrix, r):
    n = len(matrix)
    div_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                div_matrix[i, j] = matrix[i, j] - (r[i] + r[j]) / (n - 2.0)
    return div_matrix

# d(CU) = (d(AC) + d(BC) - d(AB)) / 2 = 3
def newDistMatrix(matrix, i, j):
    n = len(matrix)
    newMatrix = np.zeros((n + 1, n + 1))
    index = 0
    for x in range(n):
            for y in range(n):
                newMatrix[x][y] = matrix[x][y]
    for x in range(len(matrix)):
        newMatrix[x][len(newMatrix) - 1] = (matrix[i][x] + matrix[j][x] - matrix[i][j]) / 2.0

    for x in range(len(matrix)):
        newMatrix[len(newMatrix) - 1][x] = (newMatrix[x][i] + newMatrix[x][j] - matrix[i][j]) / 2.0

    if j < i:
        temp = j
        j = i
        i = temp
    # Delete the column
    newMatrix = np.delete(newMatrix, i, axis=1)
    # Delete the row
    newMatrix = np.delete(newMatrix, i, axis=0)

    # Delete the column
    newMatrix = np.delete(newMatrix, j-1, axis=1)
    # Delete the row
    newMatrix = np.delete(newMatrix, j-1, axis=0)

    return newMatrix

def NeighborJoining2(matrix, digits_to_round):
    n = len(matrix)
    clusters = [[i] for i in range(n)]

    G = nx.Graph()
    nodes = np.arange(len(matrix))
    counter = n


    for i in range(len(matrix)):
        G.add_edge(i, sys.maxsize, minlen = 1)

    #for xx in range(3):
    while n > 2:
        r = [sum(matrix[c]) for c in range(len(matrix))]
        div_matrix = divergence_matrix(matrix, r)
        i, j = NeighborJoining.smallestOTU(div_matrix)
        new_cluster = clusters[i] + clusters[j]

        #print(div_matrix)

        #S(AU) =d(AB) / 2 + [r(A) - r(B)] / 2(N-2) = 1
        minlen1 = round((matrix[i][j]) / 2.0 + (r[i] - r[j]) / (2.0*(n-2)), digits_to_round)
        minlen2 = round((matrix[i][j]) - minlen1, digits_to_round)

        counter+=1
        #print("combine " + str(i) + " and " + str(j) + " to " + str(counter))
        #print("delete (" + str(nodes[i]) + ", " + str(max_index_neighbor(G, nodes[i])) + ")")
        #print("delete (" + str(nodes[j]) + ", " + str(max_index_neighbor(G, nodes[j])) + ")")

        x = max_index_neighbor(G, nodes[i])
        G.remove_edge(nodes[i], x)
        G.remove_edge(nodes[j], max_index_neighbor(G, nodes[j]))

        #G.add_node(counter)
        G.add_edge(nodes[i], counter, minlen = minlen1)
        G.add_edge(nodes[j], counter, minlen= minlen2)
        G.add_edge(counter, x, minlen=1)

        nodes = np.delete(nodes,[i,j])
        nodes = np.append(nodes, counter)

        clusters = [clusters[k] for k in range(len(clusters)) if k != i and k != j]
        clusters.append(new_cluster)

        new_matrix = newDistMatrix(matrix, i, j)
        matrix = new_matrix

        n -= 1

    if (n == 2):
        x = max_index_neighbor(G, nodes[0])
        G.remove_edge(nodes[0], x)
        G.remove_edge(nodes[1], max_index_neighbor(G, nodes[1]))

        G.add_edge(nodes[0], nodes[1], minlen= round(matrix[0,1], digits_to_round))

    return (clusters, G)

def displayGraph(G, names):

    nameLabels = {it :names[it] for it in range(len(names))}

    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.gca()
    ax.axis("off")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_edges(
        G, pos, ax=ax
    )
    minlens = nx.get_edge_attributes(G, 'minlen')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=minlens, font_size=5)
    nx.draw_networkx_labels(
        G, pos, ax=ax, labels=nameLabels, font_size=8,
        bbox=dict(pad=0, color="white")
    )
    fig.tight_layout()

    plt.show()

def saveGraph(G, names, file_name):
    nameLabels = {it: names[it] for it in range(len(names))}

    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.gca()
    ax.axis("off")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_edges(
        G, pos, ax=ax
    )
    minlens = nx.get_edge_attributes(G, 'minlen')
    nx.draw_networkx_edge_labels(G, pos, font_size=6, edge_labels=minlens)
    nx.draw_networkx_labels(
        G, pos, ax=ax, labels=nameLabels, font_size=6,
        bbox=dict(pad=0, color="white")
    )
    fig.tight_layout()

    plt.savefig(file_name)



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