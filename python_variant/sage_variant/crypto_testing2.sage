def ord_div(D):
    D_temp = D
    n = 1
    
    while True:
        D_temp = D_temp + D
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



p = next_prime(2^80)
F = GF(p)
x = F['x'].gen()
f = x^5 + x^3 + 1
H = HyperellipticCurve(f, 0)
J = H.jacobian()(F)
P = H.lift_x(F(3))
D = J(P)   

print(D)
print(ord_div(D))
'''


p = 2003 # 7
K = GF(p)
R.<x> = K[]
f = x^5 + 1184*x^3 + 1846*x^2 + 956*x + 560 # x^5 + 2 * x^2 + x + 3
h = 0

C = HyperellipticCurve( f, h )
J = C.jacobian()
X = J(K)

U = (x^2 + 180*x + 1562) # ( x + 4 )

V = (962*x + 938) # ( 1 )

D = X( [ U,V ] )

print(D[0].roots())

order = 2^168 # ord_div(D)


key = 3
D_public_key = key * D
int_hash = 325 # 1415821221623963719413415453263690387336440359920

def MMI(A, n, s=1, t=0, N=0):
    return (n < 2 and t%N or MMI(n, A%n, t, s-A//n*t, N or n),-1)[n<1]

def sign():
    # m_encode = message.encode("utf-8")
    # hashed = sha1(m_encode).digest()
    # int_hash = int.from_bytes(hashed, "big")
    message = int_hash

    k_rand = randint(1, order)
    while MMI(k_rand, order) == -1:
        k_rand = randint(1, order)
    
    Dq = k_rand * D
    print("k_rand, ", k_rand)
    s = (MMI(k_rand, order) * (int_hash + key * phi(Dq))) % order

    return (message, (Dq), s)

def verify(m, Dq, s):
    # m_encode = message.encode("utf-8")
    # hashed = sha1(m_encode).digest()
    # int_hash = int.from_bytes(hashed, "big")
    message = int_hash
    
    v1 = (MMI(s, order) * (int_hash)) % order
    v2 = (MMI(s, order) * phi(Dq)) % order

    Dv = v1 * D + v2 * D_public_key
    print(Dv)
    print(Dq)
    if Dv == Dq:
        print("Signature accepted")
    else:
        print("Signature rejected")

message = int_hash
m, Dq, s = sign()

verify(m, Dq, s)
'''