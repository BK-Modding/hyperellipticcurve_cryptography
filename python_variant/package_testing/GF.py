# https://stackoverflow.com/questions/48065360/interpolate-polynomial-over-a-finite-field/48067397#48067397

import itertools
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import (gf_irreducible_p, gf_add, \
                                     gf_sub, gf_mul, gf_rem, gf_gcdex)
from sympy.ntheory.primetest import isprime


class GF():
    def __init__(self, p, n=1):
        p, n = int(p), int(n)
        if not isprime(p):
            raise ValueError("p must be a prime number, not %s" % p)
        if n <= 0:
            raise ValueError("n must be a positive integer, not %s" % n)
        self.p = p
        self.n = n
        if n == 1:
            self.reducing = [1, 0]
        else:
            for c in itertools.product(range(p), repeat=n):
                poly = (1, *c)
                if gf_irreducible_p(poly, p, ZZ):
                    self.reducing = poly
                    break
    def __eq__(self, other):
        return self.p == other.p and self.n == other.n and self.reducing == other.reducing

    def add(self, x, y):
        return gf_add(x, y, self.p, ZZ)

    def sub(self, x, y):
        return gf_sub(x, y, self.p, ZZ)

    def mul(self, x, y):
        return gf_rem(gf_mul(x, y, self.p, ZZ), self.reducing, self.p, ZZ)

    def minus(self, x):
        return gf_sub([0], x, self.p, ZZ)

    def reduce(self, x):
        return gf_rem(gf_mul(x, [1], self.p, ZZ), self.reducing, self.p, ZZ)

#    def div(self, x, y):
#        return gf_div

    def inv(self, x):
        s, t, h = gf_gcdex(x, self.reducing, self.p, ZZ)
        return s

    def eval_poly(self, poly, point):
        val = []
        for c in poly:
            val = self.mul(val, point)
            val = self.add(val, c)
        return val


# Class PolyRing, polynomials over a field
class PolyRing():
    def __init__(self, field):
        self.K = field

    def add(self, p, q):
        s = [self.K.add(x, y) for x, y in \
             itertools.zip_longest(p[::-1], q[::-1], fillvalue=[])]
        return s[::-1]

    def sub(self, p, q):
        s = [self.K.sub(x, y) for x, y in \
             itertools.zip_longest(p[::-1], q[::-1], fillvalue=[])]
        return s[::-1]

    def mul(self, p, q):
        if len(p) < len(q):
            p, q = q, p
        s = [[]]
        for j, c in enumerate(q):
            s = self.add(s, [self.K.mul(b, c) for b in p] + \
                         [[]] * (len(q) - j - 1))
        return s

def poly_degree(poly):
    return len(poly) - 1

# to calculate invert elem of x modulo p - prime
def mod_inv(x, p):
    return pow(x, p - 2, p)

# apply Polynomial(defined by list of coeffs) to x
#The coefficients must be in ascending order (``x**0`` to ``x**n``).
def poly_calc(x, coeffs):
    if not coeffs :
        return 0*x
    n = len(coeffs)
    y = 0
    for i in range(0, n):
        y += coeffs[i]*x**i
    return y

#polynom derivation
#The coefficients must be in ascending order (``x**0`` to ``x**n``).
def deriv_poly(poly):
    return [poly[i] * i for i in range(1, len(poly))]
