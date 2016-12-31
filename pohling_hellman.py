# coding: utf-8

"""
Autor          : Milutin Cerovic
Implementacija : Pohlig-Hellman algoritma za nalazenje diskretnog logaritma
Predmet        : Kriptografija
Datum          : Januar xx, 2017
"""

from math import sqrt, ceil

def PrimeFactorization(n):
    """Funkcija za odredjivanje svih prostih delilaca broja n"""
    d, primeFactors = 2, []
    while d*d <= n:
        while (n % d) == 0:
            primeFactors.append(d)
            n //= d
        d += 1
    if n > 1:
       primeFactors.append(n)
    return primeFactors

def CountOccurences(primeFactors):
    """Funkcija za prebrojavanje broja puta pojavljivanja svakog prostog delioca"""
    return [[x, primeFactors.count(x)] for x in set(primeFactors)]

def ExtendedGCD(a, b):
    """Prosireni Euklidov algoritam za racunanje najveceg zajednickog delioca"""
    a2, a1 = 1, 0
    b2, b1 = 0, 1
    while b:
        q, r = divmod(a, b)
        a1, a2 = a2 - q * a1, a1
        b1, b2 = b2 - q * b1, b1
        a, b = b, r
    return a, a2, b2

def ChineseRemainder(pairs):
    """Kineska teorema o ostacima"""
    N, X = pairs[0][1], 0
    for ni in pairs[1:]:
        N *= ni[1]
    for (ai, ni) in pairs:
        mi = (N / ni)
        X += mi * ai * ExtendedGCD(mi, ni)[1]
    return X % N


def PrintFormated(arg1, arg2, arg3, arg4, arg5):
    print " {:1s} | {:1s} | {:13s} | {:13s} | {:45s}".format(str(arg1), str(arg2), str(arg3), str(arg4), str(arg5))
    print("-"*86)

def PohlingHellman(beta, alpha, p):
    """Glavna funkcija za pozivanje PohlingHellman algoritma"""
    a = CountOccurences(PrimeFactorization(p - 1))
    congruenceList = []

    print("\n")
    print("-"*86)
    print(" Solving %d ≡ %d^x (mod %d)" % (beta, alpha, p))
    print("-"*86)
    PrintFormated("q", "e", "g^((p-1)/q^e)", "h^((p-1)/q^e)", "Solve (g^((p-1)/q^e))^x = h^((p-1)/q^e) for x")

    for i in xrange(len(a)):
        congruenceList.append(getXModP(beta, alpha, p, a[i][0], a[i][1]))
        b = len(congruenceList)
        elem1 = (alpha ** ((p - 1) / (a[i][0] ** a[i][1]))) % p
        elem2 = (beta ** ((p - 1) / (a[i][0] ** a[i][1]))) % p
        prvi = congruenceList[b - 1][0] % congruenceList[b - 1][1]
        drugi = congruenceList[b - 1][1]
        PrintFormated(a[i][0], a[i][1], elem1, elem2, "x ≡ %2d (mod %2d)" % (prvi, drugi))

    # congruenceList = [getXModP(beta, alpha, p, q, r) for (q, r) in CountOccurences(PrimeFactorization(p - 1))]
    # print("\nGiven %d = %d^x (mod %d) => x = %d\n" % (beta, alpha, p, ChineseRemainder(congruenceList)))
    print(" Solution x = %d" % ChineseRemainder(congruenceList))
    print("-"*86)
    print("\n")

def DiscreteLogModP(a, b, p):
    """Vraca x tako da b = a**x (mod p)"""
    ax = 1
    b %= p
    for x in range(p-1):
        if ax == b: return x
        ax = ax * a % p
    return None

def ShanksAlgorithm(a, y, n):
    """Vraca x tako da y = a**x (mod n)"""
    s = int(ceil(sqrt(n)))
    A = [y * pow(a, r, n) % n for r in xrange(s)]
    for t in xrange(1,s+1):
        value = pow(a, t*s, n)
        if value in A:
            return (t * s - A.index(value)) % n

def getXModP(beta, alpha, p, q, r):
    """return (x, q**r) with (p-1) / q**r = k, 0 <= x < q**r, os beta^(x*k) = alpha^k mod p"""
    oDiv = (p-1) / q # first divided group order
    bCurrent = beta
    xFinal = 0  # returns x=x0+x1q+x2q^2+...+xiq^i with 0<=xi<q-1
    alphaRaisedModp = pow(alpha, oDiv, p)
    qPow = 1
    alphaInv = ExtendedGCD(alpha, p)[1]
    for i in range(0,r):
        betaRaisedModp = pow(bCurrent, oDiv, p)
        xCurrent = ShanksAlgorithm(alphaRaisedModp, betaRaisedModp, p)
        xFinal += xCurrent * qPow
        #now we calculate the next beta, power of q, order factor
        bCurrent = bCurrent * pow(alphaInv, xCurrent * qPow, p) % p
        qPow *= q
        oDiv /= q
    return (xFinal, qPow)

if __name__ == '__main__':

    print("\nPress CTRL + C to exit\n")
    print("="*86)
    print("Pohling-Hellman's algorithm for descrete logartihm")
    print("Formula : g^x ≡ h (mod p)")
    print("="*86)
    print("\n")

    # PohlingHellman(18, 2, 29)
    PohlingHellman(166, 7, 433)
    # PohlingHellman(12, 7, 41)

    """
    while True:
        g = int(raw_input("g: "))
        h = int(raw_input("h: "))
        p = int(raw_input("p: "))
        PohlingHellman(h, g, p)
    """
