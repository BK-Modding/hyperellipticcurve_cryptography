from sympy.polys.galoistools import gf_normal, gf_degree
from GF import *
import numpy as np
# import matplotlib.pyplot as plt

class HEC():
    field: GF
    h: list
    f: list
    genus: int
    _h: list
    _f: list
    def __init__(self, p, n, h, f):
        self.field = GF(p, n)
        self._h = h
        self._f = f
        self.h_canonic(h)
        self.f_canonic(h, f)
        self.simpify()
        self.genus = (poly_degree(self.f)-1) // 2 #gf_degree

    def h_canonic(self, h):
        if (self.field.p == 2):
            self.h = h
        else:
            self.h = []

    def f_canonic(self, h, f):
        if (self.field.p != 2):
            h_inv = [ i * mod_inv(4, self.field.p) for i in h]
        self.f = self.field.add(f, self.field.mul(h, h_inv))

    def is_canonic(self):
        return self.h == []

    def __eq__(self, other):
        return self.field == other.field and self.h == other.h and self.f == other.f

    def is_correct(self):
        u = np.arange(self.field.p)
        v = np.arange(self.field.p)
        h_reverse = self.h[::-1]
        f_reverse = self.f[::-1]
        h_deriv = deriv_poly(h_reverse)
        f_deriv = deriv_poly(f_reverse)
        res1 = pow(v, 2) + np.multiply(v, poly_calc(u, h_reverse)) - poly_calc(u, f_reverse)
        res2 = 2 * v + poly_calc(u, h_reverse)
        res3 = np.multiply(v, poly_calc(u, h_deriv)) - poly_calc(u, f_deriv)
        return sum((res1 == 0) & (res2 == 0) & (res3 == 0)) == 0

    def simpify(self):
        self.h = self.field.reduce(self.h)
        self.f = self.field.reduce(self.f)

    def _coef(self, real):
        if real:
            h = self._h[::-1]
            f = self._f[::-1]
        else :
            h = self.h[::-1]
            f = self.f[::-1]
        return h, f

    def draw_curve(self, real = True):
        y, x = np.ogrid[-10:10:300j, -10:10:300j]
        h, f = self._coef(real)
        plt.contour(x.ravel(), y.ravel(), pow(y, 2) + y * poly_calc(x, h) - poly_calc(x, f), [0])
        plt.grid()
        plt.show()

    def str(self, real = False):
        h, f = self._coef(real)
        h_str = [' {}u^{} '.format(h[i], i)  for i in range(len(h)) if h[i] != 0 ]
        f_str = [' {}u^{} '.format(f[i], i)  for i in range(len(f)) if f[i] != 0 ]
        if not f:
            f_str = '0'
        res = 'v^2 +' + '+'.join(h_str) + '= ' + '+'.join(f_str)
        res = res.replace('+= ', '= ')
        res = res.replace('+ -', '- ')
        return res
    def __str__(self):
        return self.str(False)

#    @classmethod

"""
#Example:
ff = GF(3, 4)
g = ff.mul([1], [1, 1, 1])
print(ff.reducing)
print(g)

crv = HEC(5, 12, [], [1, 0, -5, 0, 4, 0])
print(crv.h)
print(crv.f)
crv.draw_curve(True)
crv.draw_curve(False)
print(crv.is_correct())
print (crv.str(True))
"""