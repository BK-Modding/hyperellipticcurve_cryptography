from HEC import *
from GF import *

class PolyRing():
    field: GF
    curve: HEC
    a: list
    b: list
    def __init__(self, a, b, curve):
        self.curve = curve
        self.field = curve.field
        self.a = a
        self.b = b
        self.simlify()
        self.deg = self.get_degree()

    def simlify(self):
        self.a = self.field.reduce(self.a)
        self.b = self.field.reduce(self.b)

    def conjugate(self):
        field = self.field
        b = field.minus(self.b)
        a = field.add(self.a, field.mul(b, self.curve.h))
        return PolyRing(a, b, self.curve)

    def norm(self):
        res = self * self.conjugate()
        if res.b:
            print("Ring.norm Warning: norm worked incorrect for some reasons")
        return res.a
        #first = self.field.mul(self.a, self.a)
        #second = self.field.mul(self.field.mul(self.a, self.b), self.curve.h)
        #thirt = self.field.mul(self.field.mul(self.b, self.b), self.curve.f)
        #return self.field.sub(self.field.add(first, second), thirt)
    #def get_value_inf(self):

    def get_value_inf(self):
        return 'Inf'

    def get_value(self, x, y):
        a = self.a[::-1]
        b = self.b[::-1]
        return poly_calc(x, a) + y * poly_calc(x, b)

    def get_degree(self):
        return max(2*poly_degree(self.a), 2*poly_degree(self.b) + 1 + 2*self.curve.genus)
        #return poly_degree(self.norm())

    def __eq__(self, other):
        return self.field == other.field and self.curve == other.curve and self.a == other.a and other.b == self.b

    def __mul__(self, other):
        if self.curve != other.curve or self.field != other.field:
            return None
        field = self.field
        curve = self.curve
        a = field.add(field.mul(self.a, other.a), field.mul(field.mul(self.b, other.b), curve.f))
        b = field.add(field.mul(self.a, other.b), field.mul(self.b, other.a))
        return PolyRing(a, field.sub(b, curve.h), curve)

    def __add__(self, other):
        if self.curve != other.curve or self.field != other.field:
            return None
        field = self.field
        a = field.add(self.a, other.a)
        b = field.add(self.b, other.b)
        return PolyRing(a, b, self.curve)

    def __sub__(self, other):
        if self.curve != other.curve or self.field != other.field:
            return None
        field = self.field
        a = field.sub(self.a, other.a)
        b = field.sub(self.b, other.b)
        return PolyRing(a, b, self.curve)

    def __str__(self):
        a = self.a[::-1]
        b = self.b[::-1]
        a_str = [' {}u^{} '.format(a[i], i) for i in range(len(a)) if a[i] != 0]
        b_str = [' {}u^{} '.format(b[i], i) for i in range(len(b)) if b[i] != 0]
        if not a_str and not b_str:
            res = '0'
        else:
            res = '+'.join(a_str) + '+ v(' + '+'.join(b_str) + ')'
        res = res.replace('= +', '= ')
        res = res.replace('+ -', '- ')
        res = res.replace('+ v()', '')
        return res


#Exanple
crv = HEC(5, 17, [3, 3, 3], [1, -5, 0, 4, 0])
elem = PolyRing([4, 0, 9], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], crv)
elem2 = PolyRing([1, 2, 34, -1, 4], [1, 34, -1, 1, 1, 1, 1, 1], crv)

# print(elem)
#print(elem.conjugate())
#print (elem.conjugate().norm())
#print(elem.norm())
# print (elem2.norm())
# print (elem2.get_degree())
#print (elem2.conjugate().get_degree())
#print ((elem*elem2).norm())
#print (elem.field.mul(elem.norm(), elem2.norm()))

#print(elem.get_degree())
#print (poly_degree(elem.norm()))
#print(elem.get_degree())
#print (elem2.get_degree())
#print (elem.get_degree() + elem2.get_degree())
#elem02 = elem * elem2
#print (elem02.get_degree())