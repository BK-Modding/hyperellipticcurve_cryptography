from builtins import classmethod

from HEC import *
from GF import *
from ring import *

class PointHEC():
    curve: HEC
    field: GF
    x: int
    y: int
    is_inf: bool
    def __init__(self, x, y, curve):
        self.is_inf = False
        self.curve = curve
        self.field = curve.field
        self.x = x % self.field.p
        self.y = y % self.field.p

    def __str__(self):
        if self.is_inf:
            point_str = 'point = Infinity'
        else:
            point_str = 'x = ' + str(self.x) + ' y = ' + str(self.y)
        return  point_str + ' on ' + str(self.curve)

    def __eq__(self, other):
        if self.is_inf:
            return self.curve == other.curve and other.is_inf
        return self.curve == other.curve and self.x == other.x and self.y == other.y

    @classmethod
    def infinite_point(cls, curve):
        elem = cls(0, 0, curve)
        elem.is_inf = True
        return elem

    def is_finite(self):
        return not self.is_inf

    def opposite(self):
        if self.is_inf:
            return self
        x = self.x
        y = -self.y - poly_calc(x, self.curve.h[::-1])
        return PointHEC(x, y, self.curve)

    def belong_curve(self):
        if self.is_inf:
            return True
        x = self.x
        y = self.y
        return 0 == y^2 + y * poly_calc(x, self.curve.h[::-1]) - poly_calc(x, self.curve.f[::-1])

    def is_spacial(self):
        opposite = self.opposite()
        return self.x == opposite.x and self.y == opposite.y and self.is_inf == self.is_inf

    def is_ordinary(self):
        return not self.is_spacial()

class PointRing(PointHEC):
    def __init__(self, x, y, poly):
        super().__init__(x, y, poly.curve)
        self.poly = poly

    @classmethod
    def infinit_point(cls, poly):
        elem = cls(0, 0, poly)
        elem.is_inf = True
        return elem

    def __eq__(self, other):
        if self.is_inf:
            return self.poly == other.poly and other.is_inf
        return self.poly == other.poly and self.x == other.x and self.y == other.y

    def __str__(self):
        if self.is_inf:
            point_str = 'point = Infinity'
        else:
            point_str = 'x = ' + str(self.x) + ' y = ' + str(self.y)
        return  point_str + '\nof' + str(self.poly) + '\non ' + str(self.curve)

    def is_zero(self):
        if self.is_inf:
            return self.poly.get_value_inf() == 0
        return self.poly.get_value(self.x, self.y) == 0

    def is_pole(self):
        if self.is_inf:
            return self.poly.get_value_inf() == 'Inf'
        return self.poly.get_value(self.x, self.y) == 'Inf'

    def ord(self):
        if self.is_inf:
            return -max(2*poly_degree(self.poly.a), 2*self.curve.genus + 1 + 2*poly_degree(self.poly.b))
        a_deriv = self.poly.a[::-1]
        b_deriv = self.poly.b[::-1]
        times = 0
        while poly_degree(a_deriv) and poly_degree(b_deriv) and \
            poly_calc(self.x, a_deriv) == 0 and poly_calc(self.x, b_deriv) == 0:
            times = times + 1
            a_deriv = deriv_poly(a_deriv)
            b_deriv = deriv_poly(b_deriv)
        r = times
        s = 0
        if poly_calc(self.x, a_deriv) + self.y * poly_calc(self.x, b_deriv) == 0:
            norm_deriv = self.poly.norm()[::-1]
            while poly_degree(norm_deriv) and poly_calc(self.x, norm_deriv) == 0:
                s = s + 1
                norm_deriv = deriv_poly(norm_deriv)
        if self.is_ordinary():
            return r + s
        else:
            return 2*r + s

class PointFrac(PointHEC):
    def __init__(self, x, y, frac):
        super().__init__(x, y, frac.curve)
        self.frac = frac

    @classmethod
    def infinit_point(cls, frac):
        elem = cls(0, 0, frac)
        elem.is_inf = True
        return elem

    def __eq__(self, other):
        if self.is_inf:
            return self.frac == other.frac and other.is_inf
        return self.frac == other.frac and self.x == other.x and self.y == other.y

    def __str__(self):
        if self.is_inf:
            point_str = 'point = Infinity'
        else:
            point_str = 'x = ' + str(self.x) + ' y = ' + str(self.y)
        return  point_str + '\nof' + str(self.frac) + '\non ' + str(self.curve)

    def is_zero(self):
        if self.is_inf:
            return self.frac.get_value_inf() == 0
        return self.frac.get_value(self.x, self.y) == 0

    def is_pole(self):
        if self.is_inf:
            return self.frac.get_value_inf() == 'Inf'
        return self.frac.get_value(self.x, self.y) == 'Inf'

    def ord(self):
        return self.top.ord() - self.button.ord()




# curve = HEC(7, 4, [4, 1], [3, 3, 3, 3, 3,3, 3])
# pnt = PointHEC(32, -1200, curve)
# print("POIONT", pnt)
# print(pnt.is_spacial())
# point_inf = PointHEC.infinite_point(curve)
# print(point_inf)
# print(point_inf.is_spacial())

# elem = PolyRing([4, 0, 0], [1, 1, 0], curve)
# print (elem)
# point_ring = PointRing(0, 0, elem)
# point_ring_inf = PointRing.infinit_point(elem)
# print (point_ring)
# print(point_ring.ord())
# print(point_ring_inf.ord())

curve = HEC(7, 1, [0], [1, 0, 0, 2, 1, 3])
print(curve)
pnt = PointHEC(3, 1, curve)
print(pnt)

elem = PolyRing([4, 1], [1], curve)
print(elem)
point_ring = PointRing(0, 0, elem)
point_ring_inf = PointRing.infinit_point(elem)
print (point_ring)
print(point_ring.ord())
print(point_ring_inf.ord())

