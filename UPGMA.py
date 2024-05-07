
def UPGMA(matrix):
    min = [0,0]
    for i in range(len(m)):
        for j in range(len(m[0])):
            if m[i][j] < min:
                min = [i][j]

    regenerateMatrix(matrix, min[0], min[1])

def regenerateMatrix(matrix, i, j):
    newMatrix = [[0.0 for x in range(len(matrix) + 1)] for y in range(len(matrix) + 1)]

    for x in range(len(matrix)):
        for y in range(len(matrix)):
            newMatrix[x][y] = matrix[x][y]

    print(newMatrix)

    for x in range(len(matrix)):
        newMatrix[x][len(newMatrix) - 1] = (newMatrix[x][i] + newMatrix[x][j])/2.0

    for x in range(len(matrix)):
        newMatrix[len(newMatrix) - 1][x] = (newMatrix[i][x] + newMatrix[j][x])/2.0

    print(newMatrix)

def reduceMatrix(matrix, i, j):
    return 0

#testing
m = [[0,2,3],[2,0,4],[4,5,0]]
regenerateMatrix(m,0,1)
