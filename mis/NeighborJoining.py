# m= [[5], [4, 7], [7, 10, 7], [6, 9, 6, 5], [8, 11, 8, 9, 8]]
import numpy as np
"""Incorrect version of neighbor joining. Methods are reused in version 2"""
# r(A) = [1/(L-2)] * [d(AB) + d(AC) = d(AD)]
def divergence_matrix(matrix):
    n = len(matrix)
    div_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                div_matrix[i, j] = (sum(matrix[i]) + sum(matrix[j]) - matrix[i, j]) / (n - 2)
    return div_matrix

# D(ij) = d(ij) - [r(i)+ r(j)]
def newDistMatrix(matrix, i, j):
    n = len(matrix)
    new_matrix = np.zeros((n - 1, n - 1))
    index = 0
    for x in range(n):
        if x != i and x != j:
            for y in range(n):
                if y != i and y != j:
                    new_matrix[index, y - (1 if y > i else 0) - (1 if y > j else 0)] = matrix[x, y]
            index += 1
    return new_matrix

def smallestOTU(matrix):
    """Returns pair of indices with smallest value in matrix."""
    min_val = np.inf
    min_indices = (0, 1)
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i, j] < min_val:
                min_val = matrix[i, j]
                min_indices = (i, j)
    return min_indices

def NeighborJoining(matrix):
    n = len(matrix)
    clusters = [[i] for i in range(n)]
    while n > 2:
        div_matrix = divergence_matrix(matrix)
        i, j = smallestOTU(div_matrix)
        new_cluster = clusters[i] + clusters[j]
        clusters = [clusters[k] for k in range(len(clusters)) if k != i and k != j]
        clusters.append(new_cluster)
        new_matrix = newDistMatrix(matrix, i, j)
        matrix = new_matrix
        n -= 1
    return clusters
