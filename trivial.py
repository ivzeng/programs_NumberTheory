import math
import numpy as np

# returns n!
def factorial(n) :
    if n <= 1: return 1
    return n*factorial(n-1)

# gets the prime factorization of integer n stored in a list of 2-element lists [p, a], where p is the prime number and a is the order.
def primeFactorization(n):
    def primeOrder(i, n):
        a = 0
        while n%i == 0:
            n //= i
            a += 1
        return n, [[i, a]]
    p = []
    if n % 2 == 0:
        n, po = primeOrder(2, n)
        p += po
    i = 3
    while i**2 <= n:
        if n % i == 0:
            n, po = primeOrder(i, n)
            p += po
        i += 2
    if n != 1:
        p += [[n, 1]]
    return p

# prints the prime factorization
def printPrimeFactorization(n):
    print(n , '=', ''.join(str(x[0]) + '^' + str(x[1]) + ' + ' for x in primeFactorization(n))[:-3])

# computes the number of positive factors of n
def pFactorsCount(n):
    return np.prod([(p[1]+1) for p in primeFactorization(n)])

# return all positive factors of n:
def pFactors(n):
    res = []
    def getPF(cur, pos, res):
        if pos == l:
            res += [cur]
            return
        for i in range(pf[pos][1]+1):
            getPF(cur*pow(pf[pos][0],i), pos+1, res)
    pf = primeFactorization(n)
    l = len(pf)
    getPF(1,0, res)
    res.sort()
    return res

# gets the greatest common divisor of a and b using the Euclidean Algorithm
def gcd(a,b) :
    if a == 0: return b
    r = b%a
    if r < 0: r *= -1
    return gcd(r,a)

# finds x, y such that gcd(a, b) = ax + by
def EuclideanAlogorithm(a, b):
    if a < b:
        a, b = b, a
    proc = [a, b]
    q = []
    while proc[-1] != 0:
        q += [proc[-2]//proc[-1]]
        proc += [proc[-2]%proc[-1]]
        print(proc[-3], '=', proc[-2], '*', q[-1], '+', proc[-1])
    i = len(q)-1
    if i == 0:
        return [0,1] if a < b else [1,0]
    print('back substitution:')
    print(proc[i+1], '=', proc[i-1] , '-', proc[i], '*', q[i-1])
    xy = [1, -q[i-1]]
    while i > 1:
        print('   =', proc[i-1], '*', xy[0] , '+ (', proc[i-2],  '-' , proc[i-1], '*', q[i-2], ') *', xy[1])
        xy = [xy[1], xy[0] - q[i-2]*xy[1]]
        print('   =', proc[i-2], '*' ,xy[0], '+', proc[i-1], '*', xy[1])
        i -= 1
    return xy if a > b else [xy[1], xy[0]]

# returns the number of positive int k < n such that gcd(k,n) = 1
def EulerPhi(n):
    return int(n*np.prod([1-1/p[0] for p in primeFactorization(n)]))

def EulerPhiSilly(n):
    return sum(1 if gcd(i, n) == 1 else 0 for i in range(1, n))

# returns a list of positive integers < m that are invertible under multiplication modulo m
def inv_mod(m):
    feasible = []
    for i in range(0,m):
        if gcd(i, m) == 1:
            feasible += [[i, -1]]
    for i in range(len(feasible)):
        if feasible[i][1] == -1:
            for j in range(i, len(feasible)):
                if feasible[i][0] * feasible[j][0] % m == 1:
                    feasible[i][1] = feasible[j][0] 
                    feasible[j][1] = feasible[i][0]
                    break
    return feasible

# prints result by inv_mod
def printInv_mod(m):
    feasible = inv_mod(m)
    for f in feasible:
        print(f[0], "*", f[1], '=', f[0]*f[1], '=', m, "*", f[0]*f[1]//m, '+', f[0]*f[1]%m, "equiv", f[0]*f[1]%m, 'mod', m)


# returns the order l of a in alst mod m, stored in a list of 2-element list [a, l], l = -1 if such order doesn't exists (i.e. gcd(a, m) != 1)
def order_mod(alst, m):
    fpm = pFactors(EulerPhi(m))
    for i in range(len(alst)):
        a = alst[i]
        alst[i] = [a, -1]
        if gcd(a, m) == 1:
            for l in fpm:
                if a**l % m == 1:
                    alst[i][1] = l
                    break
    return alst

# return a list of primitive root mod m from alst:
def primitiveRoot(alst, m):
    pm = EulerPhi(m)
    return [i[0] for i in order_mod(alst, m) if i[1] == pm]

# prime factorization
# n = 28
# printPrimeFactorization(n)
# print(pFactorsCount(n))
# print(pFactors(n))

# gcd
#print(EuclideanAlogorithm(2022, 203))
#print(2022 * 76 + 203 * -757)
#print((2022-757), 203*(2022-757)%2022)

# phi function
#n = 672
#printPrimeFactorization(n)
#print(EulerPhi(n))
#print(EulerPhiSilly(n))

# invertible
printInv_mod(32)


# order mod m and primitive roots
# m = 32
# alst = [i for i in range(1, m)]
# om = order_mod(alst.copy(), m)
# print('order l of a mod', m, ':' ,om)
# print('a:', [o[0] for o in om])
# print("l:", [o[1] for o in om])
# pr = primitiveRoot(alst, m)
# print("number of primitive roots:", len(pr))
# print("primitive roots: ", pr)
