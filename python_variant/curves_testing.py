from galois import Poly

import HyperellipticCurves

k = HyperellipticCurves.set_field(2, [1, 0, 1, 0, 0, 1], "reg")
h = Poly([1, 1, 0], field=k)
f = Poly([1, 1, 0, 1, 0, 1], field=k)
g = 2

print(HyperellipticCurves.find_points(k, h, f))
print(HyperellipticCurves.find_points_im(k, h, f))

P = (int(k.Elements()[10]), int(k.Elements()[14]))
a = Poly([1, 5, 13, 19, 18], field=k)
b = Poly([1, 19, 14, 9], field=k)

print(HyperellipticCurves.ord_ordinary(h, f, P, a, b))

P = (int(k.Elements()[0]), int(k.Elements()[1]))
a = Poly([1, 0, 0, 1, 0, 0, 0], field=k)
b = Poly([1, 0], field=k)

print(HyperellipticCurves.ord_special(h, f, P, a, b))

a = Poly([1, 5, 13, 19, 18], field=k)
b = Poly([1, 19, 14, 9], field=k)

print(HyperellipticCurves.ord_infinity(g, a, b))

a1 = Poly([1, 2, 4, 1], field=k)
b1 = Poly([5, 8, 4, 4, 2], field=k)
a2 = Poly([9, 3], field=k)
b2 = Poly([6, 1, 2, 4, 6, 8, 9], field=k)

print(HyperellipticCurves.add_div(a1, b1, a2, b2, h, f))
print(HyperellipticCurves.double_div(a1, b1, h, f))
print(HyperellipticCurves.red_div(a1, b1, g, h, f))

# print("Div to Pol")
# Div = [[1, k.Elements()[16], k.Elements()[31]], 
#         [2, k.Elements()[21], k.Elements()[5]], 
#         [1, k.Elements()[23], k.Elements()[12]], 
#         [1, k.Elements()[0], k.Elements()[1]]]

# print(div_to_pol(k, h, f, Div))

a1 = Poly([1, 21, 22], field=k)
b1 = Poly([4, 19], field=k)
a2 = Poly([1, 27, 0], field=k)
b2 = Poly([25, 1], field=k)

print(HyperellipticCurves.cantor(a1, b1, a2, b2, h, f, g))
print(HyperellipticCurves.cantor(a1, b1, a1, b1, h, f, g))