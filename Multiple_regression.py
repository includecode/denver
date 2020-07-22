#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 02:19:55 2020

@author: pavel
"""

import math
import invert
class MultipleRegression:
    def __init__(self, dataList1):
        self.dataList1 = dataList1
        
    """Etape 2: Centrer les données""" 
    def centerData(self):
        print("\n######################################################################################")
        print("\n2##.  Centrer les données")
        print("\tA.  Moyenne de X1")
        #self.meanX1 = round((sum(self.dataList1) / len(self.dataList1)), 2)
        self.meanX1 = sum(self.dataList1) / len(self.dataList1)

        print("{} / {} = {} / {} = {}\n\n".format(self.dataList1, len(self.dataList1),sum(self.dataList1), len(self.dataList1),  self.meanX1))
        lengthX1 = len(self.dataList1)

        print("\n >> X <<                 |                  \n")
        print("\n--------------------------------------------------------------------------------")
        self.centeredListX1 = []
        for i in range(lengthX1):
                print("[X] = {} - {}/{} = {} - {} = {}                     |".format(self.dataList1[i], sum(self.dataList1), len(self.dataList1), self.dataList1[i], self.meanX1, self.dataList1[i] - self.meanX1))
                self.centeredListX1.append(self.dataList1[i] - self.meanX1)
        return self.centeredListX1, self.meanX1
           
        
    """Etape 3: Calcul de la matrice de covariance""" 
    def computeCovarianceMatrix(self, centeredListXi, centeredListY):
        self.centeredListX1 = centeredListXi
        self.centeredListX2 = centeredListY
        print("\n######################################################################################")
        print("3# Calcul de la matrice de covariance")
        print("V(X) = E(X - E(X))² = E[x - E(x)²] = {SOMMEDES i de 1 à N} [i - E(x)]² ")
        print("\n\nA.Calcul de V(x1) ///somme des carrés valeurs centrées plus haut.")
        tempX1 = []
        tempX2 = []
        for i in range(len(self.centeredListX1)):
            tempX1.append(pow(self.centeredListX1[i], 2));
            print("{} +".format(tempX1[i]), end='')
            if(i == (len(self.centeredListX1) -1)):
                print("0/ {}".format(len(self.centeredListX1) -1), end='')
        self.varianceX1 = sum(tempX1) / (len(self.centeredListX1) -1)
        print("= {} / {} = {}".format(sum(tempX1), (len(self.centeredListX1) -1), self.varianceX1))
        print("\n\nB. Calcul de V(x2).   ///somme des carrés valeurs centrées plus haut")
        for i in range(len(self.centeredListX2)):
            tempX2.append(pow(self.centeredListX2[i], 2));
            print("{} +".format(tempX2[i]), end='')
            if(i == (len(self.centeredListX2) -1)):
                print("0/ {}".format(len(self.centeredListX2) -1), end='')
        self.varianceX2 = sum(tempX2) / (len(self.centeredListX2) -1)
        print("= {} / {} = {}".format(sum(tempX2), (len(self.centeredListX2) -1), self.varianceX2))
        print("\n\nC.Calcul de la matrice de covariance.")
        print("COV(x1, x2) = 1 / (N-1){SOMMEDES i de 1 à N}(x1i - moy(x1)) * (x2i - moy(x2))")
       
        """Computing covariance matrix"""
        somme = 0.0
        for i in range(len(self.centeredListX1)):
            print("({} * {}) +".format(self.centeredListX1[i], self.centeredListX2[i]), end = '')
            if(i==(len(self.centeredListX1) -1)):
                print("0/ {}".format(len(self.centeredListX1) -1))
            somme = somme + self.centeredListX1[i] * self.centeredListX2[i]
        print("COV(x1, x2) = {} / {}".format(somme, len(self.centeredListX1)-1))
        self.covX1X2 = somme / (len(self.centeredListX1)-1)
        print("COV(x1, x2) = {}".format(self.covX1X2))
        
        """Storing covariance matrix"""
        self.covMatrix = []
        self.covMatrix.append(self.varianceX1)
        self.covMatrix.append(self.covX1X2)
        self.covMatrix.append(self.varianceX2)
        return self.covMatrix
"""
a = [16, 4, 3,  7, 2, 16, 18, 17, 13, 1, 13]
b = [7, 8, 13, 10, 7,  9,  7,  8,  8, 8, 3]
c = [6, 2, 5, 3, 6, 2, 5, 3, 7, 1, 4]
d = [6, 5, 1, 3, 5, 3, 3, 1, 6, 0, 0]"""
a1 = [5, 3, 1, 3,  1, 1, 3, 1, 5, 7, 3]
b1 = [6, 5, 3, 5,  5, 5, 5, 4, 5, 7, 5]
c1 = [7, 1, 7, 7,  7, 1, 1, 1, 7, 4, 1]
a1 = [2, 2, 1, 3, 1, 3]
b1 = [2, 2, 1, 3, 2, 2]
c1 = [2, 2, 2, 3, 7, 2]
d1 = [2, 2, 1, 2, 2, 3]
y1  = [3, 2, 1, 2, 1, 2]

x1 = MultipleRegression(a1)
x2 = MultipleRegression(b1)
x3 = MultipleRegression(c1)
x4 = MultipleRegression(d1)
y = MultipleRegression(y1)

centeredDataX1, meanX1 = x1.centerData()
centeredDataX2, meanX2 = x2.centerData()
centeredDataX3, meanX3 = x3.centerData()
centeredDataX4, meanX4 = x4.centerData()
centeredDataY, meanY =  y.centerData()


covMatrixX1Y = x1.computeCovarianceMatrix(centeredDataX1, centeredDataY)
covMatrixX2Y = x2.computeCovarianceMatrix(centeredDataX2, centeredDataY)
covMatrixX3Y = x3.computeCovarianceMatrix(centeredDataX3, centeredDataY)
covMatrixX4Y = x4.computeCovarianceMatrix(centeredDataX4, centeredDataY)

print("Le vecteur des moyennes est x barre = ")
print("( {} )".format(meanX1))
print("( {} )".format(meanX2))
print("( {} )".format(meanX3))
print("( {} )".format(meanX4))

print("Cov(x, y) = ///remonter pour voir les fractions)")
print("( {} )".format(covMatrixX1Y[1]))
print("( {} )".format(covMatrixX2Y[1]))
print("( {} )".format(covMatrixX3Y[1]))
print("( {} )".format(covMatrixX4Y[1]))


print("V(y) = {}".format(covMatrixX4Y[2]))

print("\n\n\t\t/---------------%/%/%/%/%/%/%//%//%/%/----------\n\nCov(x) = ")
print("Cov(a, a) Cov(a, b) Cov(a, c) Cov(a, d)")
print("Cov(b, a) Cov(b, b) Cov(b, c) Cov(b, d)")
print("Cov(c, a) Cov(c, b) Cov(c, c) Cov(c, d)")
print("Cov(d, a) Cov(d, b) Cov(d, c) Cov(d, d)")

print("| |")

covMatrixAa = x1.computeCovarianceMatrix(centeredDataX1, centeredDataX1)
covMatrixAb = x1.computeCovarianceMatrix(centeredDataX1, centeredDataX2)
covMatrixAc = x1.computeCovarianceMatrix(centeredDataX1, centeredDataX3)
covMatrixAd = x1.computeCovarianceMatrix(centeredDataX1, centeredDataX4)

covMatrixBa = x2.computeCovarianceMatrix(centeredDataX2, centeredDataX1)
covMatrixBb = x2.computeCovarianceMatrix(centeredDataX2, centeredDataX2)
covMatrixBc = x2.computeCovarianceMatrix(centeredDataX2, centeredDataX3)
covMatrixBd = x2.computeCovarianceMatrix(centeredDataX2, centeredDataX4)

covMatrixCa = x3.computeCovarianceMatrix(centeredDataX3, centeredDataX1)
covMatrixCb = x3.computeCovarianceMatrix(centeredDataX3, centeredDataX2)
covMatrixCc = x3.computeCovarianceMatrix(centeredDataX3, centeredDataX3)
covMatrixCd = x3.computeCovarianceMatrix(centeredDataX3, centeredDataX4)

covMatrixDa = x3.computeCovarianceMatrix(centeredDataX4, centeredDataX1)
covMatrixDb = x3.computeCovarianceMatrix(centeredDataX4, centeredDataX2)
covMatrixDc = x3.computeCovarianceMatrix(centeredDataX4, centeredDataX3)
covMatrixDd = x3.computeCovarianceMatrix(centeredDataX4, centeredDataX4)

print("(  ({})   ({})  ({})    ({})   )".format(covMatrixAa[1], covMatrixAb[1], covMatrixAc[1], covMatrixAd[1]))
print("(  ({})   ({})  ({})    ({})   )".format(covMatrixBa[1], covMatrixBb[1], covMatrixBc[1], covMatrixBd[1]))
print("(  ({})   ({})  ({})    ({})   )".format(covMatrixCa[1], covMatrixCb[1], covMatrixCc[1], covMatrixCd[1]))
print("(  ({})   ({})  ({})    ({})   )".format(covMatrixDa[1], covMatrixDb[1], covMatrixDc[1], covMatrixDd[1]))


"""Prepare matrices for inversion"""
M = []
Ma = []
Mb = []
Mc = []
Md = []
Ma.append(covMatrixAa[1])
Ma.append(covMatrixAb[1])
Ma.append(covMatrixAc[1])
Ma.append(covMatrixAd[1])
M.append(Ma)

Mb.append(covMatrixBa[1])
Mb.append(covMatrixBb[1])
Mb.append(covMatrixBc[1])
Mb.append(covMatrixBd[1])
M.append(Mb)

Mc.append(covMatrixCa[1])
Mc.append(covMatrixCb[1])
Mc.append(covMatrixCc[1])
Mc.append(covMatrixCd[1])
M.append(Mc)

Md.append(covMatrixDa[1])
Md.append(covMatrixDb[1])
Md.append(covMatrixDc[1])
Md.append(covMatrixDd[1])
M.append(Md)


print("\n\n# Calcul de l'inverse de cette matrice")
MI = invert.invert_matrix(M)
print("\n\n  =>")
#invert.print_matrix("L'inverse de cette matrice est ", MI)

print("\n\n# Calcul du vecteur des covariances empiriques de y avec x1, x2, x3 et x4")
print("Y = AX + b +  epsilon\n")
print("A^ = Cov(X)⁻¹ * Cov(X, Y) \n")
print("\t\t| |\n")
print("{}    {}    {}    {}".format(MI[0][0], MI[0][1], MI[0][2], MI[0][3]))
print("{}    {}    {}    {}".format(MI[1][0], MI[1][1], MI[1][2], MI[1][3]))
print("{}    {}    {}    {}".format(MI[2][0], MI[2][1], MI[2][2], MI[2][3]))
print("{}    {}    {}    {}".format(MI[3][0], MI[3][1], MI[3][2], MI[3][3]))

print("\n\t\t*")
print("( {} )".format(covMatrixX1Y[1]))
print("( {} )".format(covMatrixX2Y[1]))
print("( {} )".format(covMatrixX3Y[1]))
print("( {} )".format(covMatrixX4Y[1]))

print("\n\t\t| |\n(Vecteur)")
def multiplyOneVectorElt(elt, inv0, inv1, inv2, inv3):
    return (inv0 * elt) +  (inv1 * elt) +  (inv2 * elt) +  (inv3 * elt)

res = []
for j in range(len(MI)):
    line = []
    for d1InvertedMatrix in MI[j]:
        line.append(d1InvertedMatrix)
    prod = (line[0] * covMatrixX1Y[1]) + (line[1] * covMatrixX2Y[1]) + (line[2] * covMatrixX3Y[1]) + (line[3] * covMatrixX4Y[1])
    print("({} * {}) + ({} * {}) + ({} * {}) + ({} * {}) = # {} # ".format(line[0], covMatrixX1Y[1], line[1], covMatrixX2Y[1], line[2], covMatrixX3Y[1], line[3], covMatrixX4Y[1], prod))
    res.append(prod)
    
print("################################################################################")
print("b^ = yBARRE - transposée de a^ * xBarre")
print("\n\t | |")
print("\n{}".format(meanY))
print("\n\t - ")
print("\n{}".format(res))
print("\n\t\t*")
print("( {} )".format(meanX1))
print("( {} )".format(meanX2))
print("( {} )".format(meanX3))
print("( {} )".format(meanX4))
print("\n\t | |")
means = []
means.append(meanX1)
means.append(meanX2)
means.append(meanX3)
means.append(meanX4)
prods = []
for j in range(len(res)):
    prod = (res[j] * means[j])
    prods.append(prod)
    print("/// {} * {}  +".format(res[j], means[j]))

print("///=", sum(prods))
print("///moy(Y) - ", sum(prods), " = ", meanY, " - ", sum(prods))
bChapeau = meanY - sum(prods)
print("\n{}".format(bChapeau))

print("###################################")
print("#Y = Ax + b +ep ==> Y = {} + {}x1 + {}x2 + {}x3 + {}x4 +epsilon".format(bChapeau, res[0], res[1], res[2], res[3]))
print("#Y = Ax + b +ep ==> Y = {} + {}x1 + {}x2 + {}x3 + {}x4 ".format(bChapeau, res[0], res[1], res[2], res[3]))

print("\n\n###################################")
print("Déterminer une regression linéaire de y sur les variables explicatives")

print("6. Calcul du r2: proportion de variance expliquéé")
print(" r² = cov(x1, x2)² / (ecartTypeX1)² * (ecartTypeX2)²")
r2a1 = pow(covMatrixX1Y[1], 2) / ((math.sqrt(math.pow(covMatrixX1Y[0], 2))) * (math.sqrt(covMatrixX1Y[2]) * (math.sqrt(covMatrixX1Y[2]))))
print("R2_X1 = ({})² / (RACINEDE{})² * (RACINEDE{})² = {}".format(covMatrixX1Y[1], covMatrixX1Y[0], covMatrixX1Y[2], r2a1))
print(r2a1 * r2a1)

print("\n\n################################################################")
print("6. Calcul des résidus")
print("\t6.1 Calcul des valeurs observées à l'aide des coefficients")
observesList = []
for i in range(len(a1)):
        obs = (res[0] * a1[i]) + (res[1] * b1[i]) + (res[2] * c1[i]) + (res[3] * d1[i]) + bChapeau
        print("({} * {}) + ({} * {}) + ({} * {}) + ({} * {}) + {} = {}".format(res[0], a1[i], res[1], b1[i], res[2], c1[i], res[3], d1[i], bChapeau, obs))
        observesList.append(obs)
    
print("\n\nt6.1 Calcul des différences (residus)")
print("e(i) = y(i) – [ b x(i) + a ]")
print(observesList)
residus = []
for i in range(len(observesList)):
    residus.append(y1[i] - observesList[i])

print("\n\nLe vecteur des residus est: ///Ecrire en vecteur")
for resi in residus:
    print("({})".format(round(resi, 5)+0))

regressionRes = MultipleRegression(residus)
centeredDataRes, meanRes = regressionRes.centerData()
covMatrixRes = regressionRes.computeCovarianceMatrix(centeredDataRes, centeredDataY)
print(1 - (covMatrixRes[0] / (covMatrixRes[2])))


"""
round(x,3)+0


print("SCR = SOMMEDESi (x2i - x2^i)² = ", end='')
SCR = 0
for i in range(len(self.dataList1)):
    SCR = SCR + (pow(self.dataList2[i] - x2PrimeList[i], 2))

print(SCR)"""