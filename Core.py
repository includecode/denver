#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 13:04:06 2020

@author: pavel
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 20:38:19 2020

@author: pavel

 x.1299999999999999y = x.13
 1.99999999999999y   = 2
 EX:
-1.799999999999999y = -1.8

priviligier l'écriture de fraction au cas ou le resultat est flottant avec plusieurs chiffres VALIDES après la virgule
Si à une étape on a un résultat avec virgules, calculer  à la main pour ne laisser que la fraction
Suppose qu'on a toujours 2 solutions après résolutionde l'équation
"""
import math
class SimpleRegression:
    def __init__(self, dataList1, dataList2 = [0]):
        self.dataList1 = dataList1
        self.dataList2 = dataList2
        
    """Etape 2: Centrer les données""" 
    def centerData(self):
        print("\n######################################################################################")
        print("\n2##.  Centrer les données")
        print("\tA.  Moyenne de X1")
        #self.meanX1 = round((sum(self.dataList1) / len(self.dataList1)), 2)
        self.meanX1 = sum(self.dataList1) / len(self.dataList1)

        print("{} / {} = {} / {} = {}\n\n".format(self.dataList1, len(self.dataList1),sum(self.dataList1), len(self.dataList1),  self.meanX1))
        
        print("\tB.  Moyenne de X2")
        self.meanX2 = sum(self.dataList2) / len(self.dataList2)
        print("{} / {} = {} / {} = {}".format( self.dataList2, len(self.dataList2), sum(self.dataList2), len(self.dataList2), self.meanX2))
        lengthX1 = len(self.dataList1)
        lengthX2 = len(self.dataList2)
        maxLength = max(len(self.dataList1), len(self.dataList2))

        print("\n >> X1 <<                 |            >> X2 <<                     \n")
        print("\n--------------------------------------------------------------------------------")
        self.centeredListX1 = []
        self.centeredListX2 = []
        for i in range(maxLength):
            if(i < lengthX1):
                print("[X1] = {} - {}/{} = {} - {} = {}                     |".format(self.dataList1[i], sum(self.dataList1), len(self.dataList1), self.dataList1[i], self.meanX1, self.dataList1[i] - self.meanX1))
                self.centeredListX1.append(self.dataList1[i] - self.meanX1)
            if(i < lengthX2):
                print("                                             [X2] = {} - {}/{} = {} - {} = {}".format(self.dataList2[i], sum(self.dataList2), len(self.dataList2),  self.dataList2[i], self.meanX2, self.dataList2[i] - self.meanX2))
                self.centeredListX2.append(self.dataList2[i] - self.meanX2)
        
    """Etape 3: Calcul de la matrice de covariance""" 
    def computeCovarianceMatrix(self):
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
        self.covMatrix.append(self.covX1X2)
        self.covMatrix.append(self.varianceX2)
        print("\n\nCov(Xchapeau) = ")
        print("\n( {}          {} )\n( {}          {} )".format(self.covMatrix[0], self.covMatrix[1], self.covMatrix[2], self.covMatrix[3]))
        print("\nOn estime donc \"a\" par")
        print("â = Cov(y,x) / V(x)")
        aChapeau = self.covMatrix[1] / self.covMatrix[0]
        print("  = {} / {} = {}".format(self.covMatrix[1], self.covMatrix[0], aChapeau))
        
        print("\net b^ par")
        print("b^ = moy(y) - â * moy(x)")
        bChapeau = self.meanX2 - (aChapeau * self.meanX1)
        print("   = {} - (({}/{}) * {}) = {}".format(self.meanX2, self.covMatrix[1], self.covMatrix[0], self.meanX1, bChapeau))
        print("On approche donc notre phénomène par: y = âx + b^")
        print("                                        ={}x + {}".format(aChapeau, bChapeau))
        print(" /!\ (S'il faut représenter les données c'est avec cette equation (prendre quelques valeurs aléatoires de X et Y))")
        
        """Calcul du R²"""
        print("Calcul du R² (proportion de variance expliquée = R² * 100)")
        rCarree = pow(self.covX1X2, 2) / (self.varianceX1 * self.varianceX2)
        print("R² = Cov(x, y)² / V(X) * V(Y) = {}² / {} * {} = {}".format(self.covX1X2, self.varianceX1, self.varianceX2, rCarree))
        print("Proportion de la variance expliquée = {}%".format(rCarree * 100))
       
        print(" /!\ (S'il faut représenter les données c'est avec cette equation (prendre quelques valeurs aléatoires de X et Y))")

        
    """Etape 4: Diagonaliser la matrice""" 
    def diagonaliseMatrix(self):
        print("\n\n\n-######  ########----#############----############  ###########--########")
        print("-@@@-Suite demandée dans le contrôle et pas dans le TD-@@@-")
        print("\######  ########----#############----########### ###########--################")
        print("4# Diagonalisation de la matrice")
        print("1. Det(M - XI) = ")
        print("\t( {}-X          {} )\n\t( {}          {}-X )".format(self.covMatrix[0], self.covMatrix[1], self.covMatrix[2], self.covMatrix[3]))
        print("=({}-X)({}-X) - ({})²".format(self.covMatrix[0], self.covMatrix[3], self.covMatrix[1]))
        print("={} (+){}X (+){}X + X² -{}".format(self.covMatrix[0] * self.covMatrix[3], self.covMatrix[0]*-1, self.covMatrix[3]*-1, pow(self.covMatrix[1], 2)))
        
        #Disriminent
        a = 1
        b = (self.covMatrix[0]*-1) + (self.covMatrix[3]*-1) 
        c = (pow(self.covMatrix[1], 2)*-1) + (self.covMatrix[0]) * (self.covMatrix[3])
        print("= X² (+){}X (+){}".format(b, c))
        delta = (pow(b, 2)) - 4 * a * c
        print("\n\nDELTA = ({})² - 4({})({})".format(b, a, c), end='')
        print(" = {}".format(delta))
        solution = []
        if delta > 0:
            racineDeDelta=math.sqrt(delta)
            solution = [(-b+racineDeDelta)/(2*a), (-b-racineDeDelta)/(2*a)]
            print("Les solutions sont", solution, " = VALEURS PROPRES")
        elif delta < 0:
            print("Cette équation n'a pas de solution = PAS DE VALEURS PROPRES")
            solution = []  #liste vide
        else:
            solution = [-b/(2*a)] #liste d'un seul élément
            print("La solution unique est", solution, " = VALEUR PROPRE")
      
        def pgcd(a,b):
            """pgcd(a,b): calcul du 'Plus Grand Commun Diviseur' entre les 2 nombres entiers a et b"""
            while b!=0:
                r=a%b
                a,b=b,r
            return a
        
        
        """Résolution de l'équation MP = lP pour chaque l de la solution"""
        def resolveMPlP(l):
            print("\n\nRésolution de l'équation: M(alpha, beta) = {}(alpha, beta) ///(alpha, beta) en colonne".format(l))
            print("on a le système: <=>| {}alpha (+) {}beta = {}alpha ".format(self.covMatrix[0], self.covMatrix[1], l))
            print("                    | {}alpha (+) {}beta = {}beta  ".format(self.covMatrix[2], self.covMatrix[3], l))
            print("                      {}beta  = {}alpha".format(self.covMatrix[1], l - self.covMatrix[0]))
            print("                      {}alpha = {}beta \n".format(self.covMatrix[2], l - self.covMatrix[3]))
            
            """Find PGCDs"""
            pgcdEquation1 = pgcd(round((self.covMatrix[1]), 1)*10, round((l - (self.covMatrix[0])), 1)*10)
            pgcdEquation2 = pgcd(round((self.covMatrix[2]), 1)*10, round((l - (self.covMatrix[3])), 1)*10)
            pgcdEquation1 = abs(pgcdEquation1)
            pgcdEquation2 = abs(pgcdEquation2)
            """PGCD SHOULD BE POSITIVE"""
            print("(INFO)find positive PGCD switing equation 1 to int: PGCD({}[rounded to 1]*10 ET {}[rounded to 1]*10) = {}".format(round((self.covMatrix[1]), 1), round((l - (self.covMatrix[0])), 1), pgcdEquation1))
            print("(INFO)find positive PGCD switing equation 2 to int: PGCD({}[rounded to 1]*10 ET {}[rounded to 1]*10) = {}".format(round((self.covMatrix[2]), 1), round((l - (self.covMatrix[3])), 1), pgcdEquation2))

            coeffAlpha1 = (round((l - self.covMatrix[0]), 1) *10) / pgcdEquation1;
            coeffBeta1  = (round((self.covMatrix[1]), 1) * 10) / pgcdEquation1;
            coeffAlpha2 = (round((l - self.covMatrix[3]), 1) *10) / pgcdEquation2;
            coeffBeta2  = (round((self.covMatrix[2]), 1) * 10) / pgcdEquation2;
            print(pgcdEquation1)
            print("                  <=> {}beta  = {}alpha".format(coeffBeta1, coeffAlpha1))
            print("                  <=> {}alpha  = {}beta".format(coeffBeta2, coeffAlpha2))
            print("                   => beta  = {}/{}alpha".format(coeffAlpha1, coeffBeta1))
            print("                   => alpha  = {}/{}beta".format(coeffAlpha2, coeffBeta2))
            
            print("\n\nOn veut de plus que p soit orthogonale")
            print("alpha² + beta² = 1")
            newBeta = coeffAlpha2 / coeffBeta2
            print("or alpha = {}/{}beta = {}".format(coeffAlpha2, coeffBeta2, newBeta))
            print("=> ({}/{})beta² + beta² = 1".format(coeffAlpha2, coeffBeta2))
            print("=> ({}/{})beta² + beta² = 1".format(coeffAlpha2, coeffBeta2))
            print("=> {}beta² + {}beta² = {}".format(pow(coeffAlpha2, 2), pow(coeffBeta2, 2), pow(coeffBeta2, 2)))
            print("=> {}beta² = {}".format(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2), pow(coeffBeta2, 2)))

            print("beta² = {}/{}".format(pow(coeffBeta2, 2), pow(coeffAlpha2, 2) + pow(coeffBeta2, 2)))
            betaPositive = math.sqrt(pow(coeffBeta2, 2)) / math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2))
            print("\n\n=> beta = RACINE({})/RACINE({}) = {}/{} = {}".format(pow(coeffBeta2, 2), pow(coeffAlpha2, 2) + pow(coeffBeta2, 2), math.sqrt(pow(coeffBeta2, 2)), math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2)), betaPositive))
            betaNegative = math.sqrt(pow(coeffBeta2, 2)) *-1 / math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2))
            print("<<OU>>")
            print("   beta = (-)RACINE({})/RACINE({}) = {}/{} = {}".format(pow(coeffBeta2, 2), pow(coeffAlpha2, 2) + pow(coeffBeta2, 2), math.sqrt(pow(coeffBeta2, 2))*-1, math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2)), betaNegative))


            print("\nDONC alpha = {}/{} * {}/{} = {} ".format(coeffAlpha2, coeffBeta2, math.sqrt(pow(coeffBeta2, 2)), math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2)), (coeffAlpha2/ coeffBeta2) * betaPositive), end='')
            print("<<OU>> alpha = {}/{} * (-){}/{} = {}  /// A moi de choisir le bon (générlement le positif pour la premier solution et le négatif & positif pour la deuxieme solution)".format(coeffAlpha2, coeffBeta2, math.sqrt(pow(coeffBeta2, 2)), math.sqrt(pow(coeffAlpha2, 2) + pow(coeffBeta2, 2)), (coeffAlpha2/ coeffBeta2) * betaNegative))
            return ((coeffAlpha2/ coeffBeta2) * betaPositive), ((coeffAlpha2/ coeffBeta2) * betaNegative), self.covMatrix[0], self.covMatrix[1], self.covMatrix[2], self.covMatrix[3], self.centeredListX1, self.centeredListX2

        p1, p2, covMatrix0, covMatrix1, covMatrix2, covMatrix3, centeredListX1, centeredListX2 = resolveMPlP(solution[0])
        p3, p4, covMatrix0, covMatrix1, covMatrix2, covMatrix3, centeredListX1, centeredListX2 = resolveMPlP(solution[1])
        p4 = abs(p4)
        p2 = abs(p2)
        #permutIfNegative(p1, p2)
        #permutIfNegative(p3, p4)
        print("\n\nLa matrice orthogonale de passage <<p>> est donc la suivante")
        print("( {}             {} )".format(p1, p3))
        print("( {}             {} )".format(p4, p2))

        print("(Côv x~) = ")
        print("\n\t| |")
        print("\n       x1        x2 ")
        print("x1 ( {}        {} )".format(covMatrix0, covMatrix1))
        print("x2 ( {}        {} )".format(covMatrix2, covMatrix3))
        print("\n\t| |")
        print("(P)")
        print("\n       u1        u2 ")
        print("x1 ( {}        {} )".format(p1, p3))
        print("x2 ( {}        {} )\n".format(p4, p2))
        print("\t * \n")
        print("(D)")
        print("\n       u1        u2 ")
        print("x1 ({}       {}  )".format(solution[0], 0))
        print("x2 ( {}        {} )".format(0, solution[1]))
        print("\t * \n")
        print("(P⁻1 = tP)")
        print("\n       x1        x2 ")
        print("x1 ({}       {}  )".format(p1, p4))
        print("x2 ( {}        {} )".format(p3, p2))
        
        print("\n\n\n-######  ########----#############----############  ###########--########")
        print("\n\n5# Calculer les coordonnées des observations dans le plan factoriel")
        def multiplyCenteredByMatrixX1(xbar1, xbar2, coord0, coord2):
            return (xbar1 * coord0) + (xbar2 * coord2)
        def multiplyCenteredByMatrixX2(xbar1, xbar2, coord1, coord3):
            return (xbar1 * coord1) + (xbar2 * coord3)
        
        print("********** ********************************************************************* ")
        print("x~ |----------| x1~ |----------| x2~ |----------+ u1~ |----------| u2~ ")
        print("*********** ********************************************************************* ")
        correlationsX1 = []
        correlationsX2 = []
        for i in range(len(centeredListX1)):
            x1 = multiplyCenteredByMatrixX1(centeredListX1[i], centeredListX2[i], p1, p4)
            x2 = multiplyCenteredByMatrixX2(centeredListX1[i], centeredListX2[i], p3, p2)
            correlationsX1.append(x1)
            correlationsX2.append(x2)
            print("{} |________| {} |_________| {} |_________| {} |_________| {} ".format(i+1, centeredListX1[i], centeredListX2[i], x1, x2))
            
        print("\n\n********** ********************************************************************* ")
        print("5. Calcul de la correlation empirique. (ration entre la covariance et les ecarts types)")
        self.ecartTypeX1 = math.sqrt(self.varianceX1)
        self.ecartTypeX2 = math.sqrt(self.varianceX2)

        print("Ecart type de x1: ({})² = {}".format(self.varianceX1, self.ecartTypeX1))
        print("Ecart type de x2: ({})² = {}".format(self.varianceX2, self.ecartTypeX2))
        
        print("\n\nCorrelation empirique (x1,x2) = cov(x1, x2) / ecartypeX1 * ecartypeX2".format(self.covX1X2, self.ecartTypeX1, self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2)))
        print("Correlation empirique (x1,x2) = {} / {} * {} = {}".format(self.covX1X2, self.ecartTypeX1, self.ecartTypeX2, self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2)))
        print("Correlation empirique (x1) = {} / {} = {}".format(self.covX1X2, self.ecartTypeX1, self.covX1X2 / (self.ecartTypeX1)))
        print("Correlation empirique (x2) = {} / {} = {}".format(self.covX1X2, self.ecartTypeX2, self.covX1X2 / (self.ecartTypeX2)))
        
        
        print("\n\nUne correlation proche de 0 signifie qu'il n'existe pas de relation  linéaire entre X et Y")
        print("\n\nPartant des données xi et yi, on cherche à estimer des valeurs a^ et b^ reflettant le mieux possible les valeurs a et b des phénomènes sous-djacents\nainsi que la quantité de cette regression")
        print("\ncalcul de la somme des carrées de correlation de X1")
        print("Coefficient de correlation (le 'r') = cov(x1, x2) / ecartTypeX1 * ecartTypeX2")
        self.r  = self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2)
        self.r2 = self.covX1X2 * self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2 * self.ecartTypeX1 * self.ecartTypeX2)
        print(" r = {} / {} * {} = {}".format(self.covX1X2, self.ecartTypeX1, self.ecartTypeX2, self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2)))
        print(" r² = cov(x1, x2)² / (ecartTypeX1)² * (ecartTypeX2)²")
        print(" r²= {}² / {}² * {}² = {}".format(self.covX1X2, self.ecartTypeX1, self.ecartTypeX2, self.covX1X2 * self.covX1X2 / (self.ecartTypeX1 * self.ecartTypeX2 * self.ecartTypeX1 * self.ecartTypeX2)))
        
        print("\n\nLa matrice de correlation est donc la suivante")
        print("    |     X1           X2")
        print("_______________________________________")
        print("X1  | ( {}         {} )".format(1, self.r))
        print("    |")
        print("X2  | ( {}         {} )".format(self.r, 1))
        print(self.ecartTypeX1 * self.dataList1[0])
        print((self.ecartTypeX1 * self.ecartTypeX1) / (len(self.dataList1) -1) + (self.ecartTypeX2 * self.ecartTypeX2) / (len(self.dataList1)-1))
        print(pow((self.covX1X2 / (self.ecartTypeX1)), 2))
        lambda1 = 0
        for i in range(len(self.dataList1)):
            lambda1 = lambda1 + (self.r * self.r) * self.centeredListX1[i]
            print(lambda1)
        print(lambda1, "mmmmmmm")
        sumCarreesCorrX1 = 0
        sumCarreesCorrX2 = 0
        for i in range(len(correlationsX1)):
            sumCarreesCorrX1 = sumCarreesCorrX1 + math.pow(correlationsX1[i] - multiplyCenteredByMatrixX1(centeredListX1[i], centeredListX2[i], p1, p4), 2)
            sumCarreesCorrX2 = sumCarreesCorrX2 + math.pow(correlationsX2[i] - multiplyCenteredByMatrixX2(centeredListX1[i], centeredListX2[i], p3, p2), 2)
        print(sumCarreesCorrX1) 
        print(sumCarreesCorrX2)

a = [16, 4, 3,  7, 2, 16, 18, 17, 13, 1, 13]
b = [7, 8, 13, 10, 7,  9,  7,  8,  8, 8, 3]
c = [6, 2, 5, 3, 6, 2, 5, 3, 7, 1, 4]
d = [6, 5, 1, 3, 5, 3, 3, 1, 6, 0, 0]
test = SimpleRegression(c, d)
test.centerData()
test.computeCovarianceMatrix()
test.diagonaliseMatrix()
        
        
