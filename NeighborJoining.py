m= [[5], [4, 7], [7, 10, 7], [6, 9, 6, 5], [8, 11, 8, 9, 8]]
import numpy as np

#r(A) = [1/(L-2)] * [d(AB) + d(AC) = d(AD)]
def divergence_matrix(matrix):
    n = len(matrix)
    div_matrix = np.zeros((n,n))
    inner = 0
    outer = 0
    mat_index = 0
    for i in range(n):
        for j in range(i+1, n):
            div_matrix[i,j] = (matrix[i][0] + matrix[j][0] - matrix[i][j])/(n-2)
            div_matrix[j,i] = div_matrix[i,j]
    return div_matrix

# D(ij) = d(ij) - [r(i)+ r(j)]
def newDistMatrix(matrix):
    n = len(matrix)
    newMatrix = np.zeros((n-1,n-1))
    for i in range (n-1):
        for j in range (n-1):
            if i==j:
                newMatrix[i,j] = 0
            elif i < j:
                newMatrix[i,j] = matrix[i][j]
            else:
                newMatrix[j,i] = matrix[i,j]
    return newMatrix

def smallestOTU(matrix):
    small = np.min(matrix)
    index = np.where(matrix==small)
    return index[0][0], index[1][0]

def NeighborJoining(matrix):
    n = len(matrix)
    while n>2:
        divMatrix = divergence_matrix(matrix)
        smallestOTU(matrix)
        newDistMatrix(divMatrix)
