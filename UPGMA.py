import matplotlib
from Bio import Phylo

"""Contains methods to generate phylogenetic trees using the UPGMA algorithm."""

def UPGMA(matrix, names, digits_to_round):
    """Performs UPGMA and returns tree in Newick format string."""
    orig = names.copy()

    names1 = [str(i) for i in range(len(names))]
    nodes = {}
    nodesStr = {}
    newick = None

    m = matrix
    while len(m) > 1:
        minimum = findMinimum(m)

        #print("pair " + names[minimum[0]] + " and " + names[minimum[1]] + " distance " + str(m[minimum[0]][minimum[1]] / 2.0))

        newick = None
        if ' ' in names1[minimum[0]] and ' ' in names1[minimum[1]]:
            newick = "(" + nodesStr[names1[minimum[0]]][0] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names1[minimum[0]]][1], digits_to_round)) + "," + nodesStr[names1[minimum[1]]][0] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names1[minimum[1]]][1], digits_to_round)) + ")"
        elif ' ' in names1[minimum[1]]:
            newick = "(" + orig[int(names1[minimum[0]])] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0, digits_to_round)) + "," + nodesStr[names1[minimum[1]]][0] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names1[minimum[1]]][1], digits_to_round)) + ")"
        elif ' ' in names1[minimum[0]]:
            newick = "(" + nodesStr[names1[minimum[0]]] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names1[minimum[0]]][1], digits_to_round)) + "," + orig[int(names1[minimum[1]])] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0, digits_to_round)) + ")"
        else:
            newick = "(" + orig[int(names1[minimum[0]])] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0, digits_to_round)) + "," + orig[int(names1[minimum[1]])] + ":" + str(round(m[minimum[0]][minimum[1]] / 2.0, digits_to_round)) + ")"

        nodesStr[names1[minimum[0]] + " " + names1[minimum[1]]] = (newick, m[minimum[0]][minimum[1]] / 2.0)

        m = regenerateMatrix(m, minimum[0], minimum[1], names1)
        #print(m)

    #print(nodes)
    print(newick)
    return newick


def findMinimum(matrix):
    """Returns index of minimum in 2D array."""
    minimum = [0, 1]
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i == j:
                continue
            if matrix[i][j] < matrix[minimum[0]][minimum[1]]:
                minimum = [i, j]
    return minimum


def regenerateMatrix(matrix, i, j, names):
    """Recalculates distance matrix after combination and removal of nodes i and j (indices)."""
    newMatrix = [[0.0 for x in range(len(matrix) + 1)] for y in range(len(matrix) + 1)]
    names.append(str(names[i]) + " " + str(names[j]))

    #populates matrix with old distances
    for x in range(len(matrix)):
        for y in range(len(matrix)):
            newMatrix[x][y] = matrix[x][y]

    #calculates new distances for added node
    for x in range(len(matrix)):
        newMatrix[x][len(newMatrix) - 1] = (newMatrix[x][i] + newMatrix[x][j]) / 2.0

    for x in range(len(matrix)):
        newMatrix[len(newMatrix) - 1][x] = (newMatrix[i][x] + newMatrix[j][x]) / 2.0

    #removed nodes that were combined
    reduceMatrix(newMatrix, i, j, names)

    return newMatrix


def reduceMatrix(matrix, i, j, names):
    """Removes column/rows of index i and j from 2D array."""
    if j < i:
        temp = j
        j = i
        i = temp
    del matrix[i]
    del matrix[j - 1]
    for x in matrix:
        del x[i]
        del x[j - 1]
    del names[i]
    del names[j - 1]


# #testing
# n = ["A", "B", "C", "D", "E"]
# m = [[0,6,1,2,6], [6,0,6,6,4], [1,6,0,2,6], [1,6,2,0,6], [6,4,6,6,0]]
# print(m)
# string = UPGMA(m, n)

#%matplotlib inline
#parser = Phylo.NewickIO.Parser.from_string(string)
#tree = parser.parse()
#Phylo.draw(list(tree)[0], branch_labels=lambda c: c.branch_length)