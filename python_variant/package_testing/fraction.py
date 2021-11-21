from HEC import *
from GF import *
from ring import *

class PolyFract:
    field: GF
    curve: HEC
    top: PolyRing
    button: PolyRing

    def __new__(cls, top, button):
        top.simlify()
        button.simlify()
        if top.field ==  button.field and top.curve == button.curve and not (button.a == [] and button.b == []):
            #if a == [] then PolyFract with a , b notm
            #if b == [] then PolyRing
            return super(PolyFract, cls).__new__(cls)

    def __init__(self, top, button):
        self.field = top.field
        self.curve = top.curve
        self.top = top
        self.button = button
        self.simlify()

    def simlify(self):
         self.top.simlify()
         self.button.simlify()

    def __eq__(self, other):
        return self.top == other.top and self.button == other.button and self.field == other.field and self.curve == other.curve

    def __str__(self):
        top_str = str(self.top)
        button_str = str(self.button)
        len_sep = max(len(top_str), len(button_str))
        separator = "".join(['-' for i in range(len_sep)])
        return top_str + '\n' + separator + '\n' + button_str

    def get_value_inf(self):
        if self.top.degree > self.button.degree:
            return 0
        if self.top.degree < self.button.degree:
            return 'Inf'
        return self.top.b[0] / self.button.b[0] ## + ???

    def get_value(self, x, y):
        top = self.top.get_value(x, y)
        button = self.button.get_value(x, y)
        if button == 0 and top != 0:
            return 'Inf'
        if button != 0:
            return top / button
        # 0/0
        a_top_deriv = self.top.a[::-1]
        b_top_deriv = self.top.b[::-1]
        a_button_deriv = self.button.a[::-1]
        b_button_deriv = self.button.b[::-1]
        while poly_degree(a_top_deriv) and poly_degree(b_top_deriv) and \
            poly_degree(a_button_deriv) and poly_degree(b_button_deriv) and \
            poly_calc(x, a_top_deriv) == 0 and poly_calc(x, b_top_deriv) == 0 and \
            poly_calc(x, a_button_deriv) == 0 and poly_calc(x, b_button_deriv) == 0:

            a_top_deriv = deriv_poly(a_top_deriv)
            b_top_deriv = deriv_poly(b_top_deriv)
            a_button_deriv = deriv_poly(a_button_deriv)
            b_button_deriv = deriv_poly(b_button_deriv)
        top = poly_calc(x, a_top_deriv) + y*poly_calc(x, b_top_deriv)
        button = poly_calc(x, a_button_deriv) + y*poly_calc(x, b_button_deriv)
        if button != 0:
            return top / button
        if top != 0 and button == 0:
            return 'Inf'
        if top == 0 and button == 0:
            return 'Undefine'  #make better someday


#def __add__
#def __sub__
#def __mul__


curve = HEC(7, 4, [4, 1], [3, 3, 3, 3, 3,3, 3])
poly1 = PolyRing([1, 0, 0], [3, -5, 0], curve)
poly2 = PolyRing([2, 0, 0, 0], [4, 0, 3, 0, 0, 0], curve)
fract = PolyFract(poly1, poly2)
print (fract)
print (fract.get_value(0, 0))