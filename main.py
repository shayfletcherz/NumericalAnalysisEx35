# https://github.com/shayfletcherz/NumericalAnalysisEx35.git

# EX 35

from datetime import datetime

#Exchanging lines
def manualSwapRow(matrix, b, row1, row2):
    if row2 < len(matrix) and row1 < len(matrix):
        temp = matrix[row1]
        matrix[row1] = matrix[row2]
        matrix[row2] = temp
        if b is not None:
            temp = b[row1]
            b[row1] = b[row2]
            b[row2] = temp
    return matrix, b


#Exchanging cols
def manualSwapCol(matrix, col1, col2):
    if col2 < len(matrix) and col1 < len(matrix):
        for i in range(len(matrix)):
            temp = matrix[i][col1]
            matrix[i][col1] = matrix[i][col2]
            matrix[i][col2] = temp
    return matrix


#Copy matrix function
def copyMat(matrix):
    B = makeMatrics(len(matrix), len(matrix[0]))
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            B[i][j] = matrix[i][j]
    return B


def rowSum(line):
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


def createDominantDiagonal(matrix, b):
    max = 0
    maxIndex = 0
    sum = 0
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            sum += abs(matrix[i][j])
            if abs(matrix[i][j]) > max:
                max = abs(matrix[i][j])
                maxIndex = j
        if (sum - max) <= max:
            matrix = manualSwapCol(matrix, maxIndex, i)
        else:
            max = 0
            maxIndex = 0
            for j in range(len(matrix)):
                sum += abs(matrix[j][i])
                if abs(matrix[j][i]) > max:
                    max = abs(matrix[j][i])
                    maxIndex = j
            if rowSum(matrix[j]) - max <= max:
                matrix, b = manualSwapRow(matrix, b, i, maxIndex)
            else:
                return None, None
    return matrix, b



# A function will get coefficients of the polynomial and the X that we want to find and return an approximate value
def getCoefficientsCalcY(coefficients, X):
    sum = 0
    for x in range(len(coefficients)):
        sum += coefficients[x][0] * (X ** x)
    return sum

#Presentation matrix and vector result after pivot
def isolateVariables(coefficientMat, b):
    vectorB = makeMatrics(len(coefficientMat))
    matA = makeMatrics(len(coefficientMat), len(coefficientMat[0]))
    for i in range(len(coefficientMat)):
        vectorB[i][0] = b[i][0] / coefficientMat[i][i]
        j = 0
        while j < len(coefficientMat[0]):
            if j is i:
                matA[i][i] = 1
            else:
                matA[i][j] -= coefficientMat[i][j] / coefficientMat[i][i]
            j += 1
    return matA, vectorB

#Presentation of iterations using the Seidel method
def gaussSeidelIter(coefficientMat, b):
    epsilon = 0.00001
    iteration = 0
    flag = True
    prevX = makeMatrics(len(coefficientMat))  # start as zero vector
    currentX = makeMatrics(len(coefficientMat))
    matA, vectorB = isolateVariables(coefficientMat, b)
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        prevX[0][0] = currentX[0][0]
        if iteration >= 100:
            break
        for i in range(len(coefficientMat)):
            j = 0
            currentX[i][0] = vectorB[i][0]
            while j < len(coefficientMat[0]):
                if j is not i:
                    currentX[i][0] += matA[i][j] * currentX[j][0]
                j += 1
        iteration += 1
    return currentX

#Creating a zero matrix
def makeMatrics(row, col=1):
    c = []
    for i in range(row):
        c += [[0] * col]
    return c

#Creating a polynomial matrix
def makePolynomialMat(points):
    size = len(points)
    newMat = makeMatrics(size, size)
    newB = makeMatrics(size, size)
    for i in range(size):
        xi = points[i][0]
        for j in range(size):
            newMat[i][j] = xi ** j
        newB[i][0] = points[i][1]
    return newMat, newB

#*****************************************************************
#*****************************************************************

#Evaluation Point
XtoFind = 1.37

#Data table
table = [[1.2,3.5095],
       [1.3,3.6984],
       [1.4,3.9043],
       [1.5,4.1293],
       [1.6,4.3756]]

day = str(datetime.today().day)
hour = str(datetime.now().hour)
minute = str(datetime.now().minute)

#Neville's algorithm
def neville_interpolation(table,xf):
    n = len(table)
    x = 0
    y = 1
    print("*************** Neville interpolation ***************\n")
    print("The Formula: P m,n(X) = (X-Xm)*P m+1,n (X) - (X-Xn)*Pm,n-1(X)")
    print("                        ------------------------------------")
    print("                                    Xn-Xm \n")
    print("intermediate stages:")
    for i in range(1,n,+1):
        for j in range(n-1,i-1,-1):
            table[j][y] = ((xf-table[j-i][x])*table[j][y]-(xf-table[j][x])*table[j-1][y])/(table[j][x]-table[j-i][x])
            print(table[j][y])
    result = str(table[n-1][y])
    print("The Approximate value by Neville interpolation:")
    return result+"00000"+day+hour+minute

#Calculation of an estimated value by a polynomial method
def polynomial(table, X):
    a, b = makePolynomialMat(table)
    copyA = copyMat(a)
    copyB = copyMat(b)
    copyA, copyB = createDominantDiagonal(copyA, copyB)
    if (copyA is not None) and (copyB is not None):
        a = copyA
        b = copyB
    #
    coefficients = gaussSeidelIter(a, b)
    print("***************Polynomial interpolation***********\n")
    print("The Formula: y(f)= y1-y2        y2*x1-y1*x2")
    print("                   -----  * xf + -----------")
    print("                   x1-x2         x1-x2\n")
    print("The Approximate value by Polynomial interpolation:")
    x = str(getCoefficientsCalcY(coefficients, X))
    return x+"00000"+day+hour+minute


print(neville_interpolation(table,XtoFind))
print("-------------------------------------------------------------")
print(polynomial(table, XtoFind))
