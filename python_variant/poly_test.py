import galois

def test_poly_gcd_2():
    GF = galois.GF(2**5, irreducible_poly=galois.Poly(list(reversed([1, 0, 1, 0, 0, 1]))))
    a = galois.Poly([1, 2, 4, 1], field=GF)
    b = galois.Poly([9, 3], field=GF)
    gcd, x, y = galois.poly_gcd(a, b)
    assert a*x + b*y == gcd

test_poly_gcd_2()