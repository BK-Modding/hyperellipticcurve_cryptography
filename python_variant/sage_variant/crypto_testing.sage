def ord_div(D):
    D_temp = D
    n = 0

    while True:
        D_temp = D_temp +  D
        n += 1
        print(n, D_temp)
        if (D_temp[0], D_temp[1]) == (1, 0):
            break
        
    return n

def phi(D):
    field_order = int(D.scheme().base_ring().order())
    pol_ord = int(D[0].degree())
    coeffs = D[0].coefficients()

    modcoeffsum = 0    
    for i in range(0, pol_ord + 1):
        modcoeffsum = modcoeffsum + (int(coeffs[i] % field_order) * field_order**i)

    return modcoeffsum

def jacobian_ord(F):
    q = 2
    n = 56
    E_f = lambda x, y: y^2 + (x^2 + x)*y + x^7 + x^5 + x^4 + x^3 + x^2 + x + 1

    # Step 1 
    # Note: + 1 in calcs. represents incluing the point at infinity.

    def count_points(F):
        return [(x, y) for x in F for y in F if E_f(x, y) == 0]
        
    ps1 = count_points(GF(2))
    m1 = len(ps1) + 1
    print('points in GF(2) =', ps1)
    print('M1 =', m1)

    ps2 = count_points(GF(4))
    m2 = len(ps2) + 1
    print('points in GF(4) =', ps2)
    print('M2 =', m2)

    # Step 2
    a1 = m1 - 1 - q
    a2 = (m2 - 1 - q^2 + a1^2) / 2
    print('a1 =', a1)
    print('a2 =', a2)

    # Step 3
    # X^2 + a1X + (a2 − 2q) = 0 => x^2 - x - 5 = 0 => zeros = (1 +- sqrt(21)) / 2
    var('X')
    gammas = list(map(lambda x: x.rhs(), solve([X^2 + a1*X + (a2 - 2 * 2) == 0], X)))
    print('gammas =', gammas)

    # Step 4
    # X^2 − gamma1X + q = 0
    alpha1 = list(map(lambda x: x.rhs(), solve([X^2 - gammas[0]*X + q == 0], X)))[0]
    alpha2 = list(map(lambda x: x.rhs(), solve([X^2 - gammas[1]*X + q == 0], X)))[0]
    print('alpha1 =', alpha1)
    print('alpha2 =', alpha2)

    # Step 5
    nj = int(abs(1-alpha1^n)^2 * abs(1-alpha2^n)^2)
    print('size of jacobian =', nj)
    return nj



p = next_prime(2^56)
F = GF(p)
x = F['x'].gen()
f = x^7 + x^5 + x^4 + x^3 + x^2 + x + 1
h = x^2 + x

H = HyperellipticCurve(f, h)
J = H.jacobian()(F)

# p1, p2 = randint(1, 2^56), randint(1, 2^56)
p1 = 26388590002444381
p2 = 58492723559085106

P1 = H.lift_x(F(p1))
P2 = H.lift_x(F(p2))

D1 = J(P1)
D2 = J(P2)   

D = D1 + D2

print(D)

order = int(jacobian_ord(F))

key = 3
D_public_key = key * D
int_hash = 325 # 1415821221623963719413415453263690387336440359920


def sign():
    w = int(int_hash % order)
    e = randint(1, order) # 5
    De = e * D
    F_e = int(phi(De) % order)
    print("Fe: ", F_e)
    r = int((w * F_e) % order)
    s = int((e + key * r) % order)

    print(r, s)
    return r, s

def verify(r, s):
    w2 = int(int_hash % order)
    print(w2)
    Ds = s * D
    Dr = r * D_public_key
    Dsum = Ds + Dr
    
    y = int(phi(Dsum) % order)
    v2 = w2 * y
    print(r, v2)
    print(r == v2)

r, s = sign()
verify(r, s)