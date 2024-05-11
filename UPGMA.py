
def UPGMA(matrix):

    m = matrix
    while len(m) > 1:
        minimum = findMinimum(m)

        print("pair " + str(minimum[0]) + " " + str(minimum[1]) + " distance " + str(m[minimum[0]][minimum[1]]/2.0))
        m = regenerateMatrix(m, minimum[0], minimum[1])
        print(m)

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
def regenerateMatrix(matrix, i, j):

    #create new matrix 1 size larger
    newMatrix = [[0.0 for x in range(len(matrix) + 1)] for y in range(len(matrix) + 1)]

    #populates matrix with old distances
    for x in range(len(matrix)):
        for y in range(len(matrix)):
            newMatrix[x][y] = matrix[x][y]

    #calculates new distances for added node
    for x in range(len(matrix)):
        newMatrix[x][len(newMatrix) - 1] = (newMatrix[x][i] + newMatrix[x][j])/2.0

    for x in range(len(matrix)):
        newMatrix[len(newMatrix) - 1][x] = (newMatrix[i][x] + newMatrix[j][x])/2.0

    #removed nodes that were combined
    reduceMatrix(newMatrix, i, j)

    return newMatrix

#removes distance information for nodes i and j from 2D array
def reduceMatrix(matrix, i, j):
    if j < i:
        temp = j
        j = i
        i = temp
    del matrix[i]
    del matrix[j - 1]
    for x in matrix:
        del x[i]
        del x[j - 1]

#testing
m = [[0,9,9,4,5],[9,0,16,7,9],[9,16,0,11,7],[4,7,11,0,13],[5,9,7,13,0]]
print(m)
UPGMA(m)
