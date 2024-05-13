def UPGMA(matrix, names1):

    orig = names1.copy()
    #newick format: (a:8, (b:2, d:9):4)

    names = [str(i) for i in range(len(names1))]
    nodes = {}
    nodesStr = {}
    newick = None
    counter = 0

    m = matrix
    while len(m) > 1:
        minimum = findMinimum(m)

        print("pair " + names[minimum[0]] + " " + names[minimum[1]] + " distance " + str(m[minimum[0]][minimum[1]] / 2.0))

        newick = None
        if ' ' in names[minimum[0]] and ' ' in names[minimum[1]]:
            newick = "(" + nodesStr[names[minimum[0]]][0] + ":" + str(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names[minimum[0]]][1]) + "," + nodesStr[names[minimum[1]]][0] + ":" + str(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names[minimum[1]]][1]) + ")"
        elif ' ' in names[minimum[1]]:
            newick = "(" + orig[int(names[minimum[0]])] + ":" + str(m[minimum[0]][minimum[1]] / 2.0) + "," + nodesStr[names[minimum[1]]][0] + ":" + str(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names[minimum[1]]][1]) + ")"
        elif ' ' in names[minimum[0]]:
            newick = "(" + nodesStr[names[minimum[0]]] + ":" + str(m[minimum[0]][minimum[1]] / 2.0 - nodesStr[names[minimum[0]]][1]) + "," + orig[int(names[minimum[1]])] + ":" + str(m[minimum[0]][minimum[1]] / 2.0) + ")"
        else:
            newick = "(" + orig[int(names[minimum[0]])] + ":" + str(m[minimum[0]][minimum[1]] / 2.0) + "," + orig[int(names[minimum[1]])] + ":" + str(m[minimum[0]][minimum[1]] / 2.0) + ")"

        nodesStr[names[minimum[0]] + " " + names[minimum[1]]] = (newick, m[minimum[0]][minimum[1]] / 2.0)

        m = regenerateMatrix(m, minimum[0], minimum[1], names)
        print(m)

    print(nodes)
    print(newick)
    return newick

#returns index of minimum in 2D array
def findMinimum(matrix):
    minimum = [0, 1]
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i == j:
                continue
            if matrix[i][j] < matrix[minimum[0]][minimum[1]]:
                minimum = [i, j]
    return minimum


#recalculates distance matrix after combination and removal of nodes i and j (indices)
def regenerateMatrix(matrix, i, j, names):

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


#removes column/rows of index i and j from 2D array
def reduceMatrix(matrix, i, j, names):
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


#testing
n = ["A", "B", "C", "D", "E", "F"]
m = [[0,6,8,1,2,6], [6,0,8,6,6,4], [8,8,0,8,8,8], [1,6,8,0,2,6], [1,6,8,2,0,6], [6,4,8,6,6,0]]
print(m)
UPGMA(m, n)
