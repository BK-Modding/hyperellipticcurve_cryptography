from numpy.random import randint
from hashlib import sha1

import HyperellipticCurves

# Algorithm to compute the order of a reduced divisor in Mumford representation
def ord_div(ap, bp, k, h, f, g):
    '''
    input:
        ap: (galois.Poly) first polynomial of the reduced divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the reduced divisor in Mumford representation
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined

    returns:
        n: (int) order of the reduced divisor
    '''
    aq = ap
    bq = bp
    test = 0
    n = 1
    
    while test == 0:
        aq, bq = HyperellipticCurves.cantor(aq, bq, ap, bp, h, f, g)
        n += 1
        if (aq, bq) == (HyperellipticCurves.i_p(1, k), HyperellipticCurves.i_p(0, k)):
            test = 1

    return n

# multiply a divisor in Mumford representation with a scalar integer
def multi_div(a, b, multi, k, h, f, g):
    '''
    input:
        a: (galois.Poly) first polynomial of the divisor in Mumford representation
        b: (galois.Poly) second polynomial of the divisor in Mumford representation
        multi: (int) integer the divisor is to be multiplied with
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined

    returns:
        ad: (galois.Poly) first polynomial of the multiplied divisor in Mumford representation
        bd: (galois.Poly) second polynomial of the multiplied divisor in Mumford representation
    '''
    binary = list(map(int, list(bin(multi)[2:])))
    cc = [0 for i in range(1, len(binary) + 1)]
    
    li = []
    r = []
    for i in range(len(binary)):
        if binary[i] == 1:
            cc_cop = cc.copy()
            # cc[i] = 1
            cc_cop[i] = 1
            li.append(cc_cop)
            r.append(len(binary) - i - 1)
            # cc[i] = 0

    l = [int("".join(map(str, li[i])), 2) for i in range(len(li))]

    dDiv = []

    for i in range(len(r)):
        aa, bb = a, b
        for j in range(r[i]):
            aa, bb = HyperellipticCurves.double_div(aa, bb, h, f)
            aa, bb = HyperellipticCurves.red_div(aa, bb, 2, h, f)
        
        dDiv.append((aa, bb))

    ad, bd = HyperellipticCurves.i_p(1, field=k), HyperellipticCurves.i_p(0, field=k)

    for i in range(len(dDiv)):
        ad, bd = HyperellipticCurves.cantor(ad, bd, dDiv[i][0], dDiv[i][1], h, f, 2)

    return ad, bd

# ap, bp: polynomials of the principal divisor D1
# key: private key
# aq, bq: public key polynomials
# am, bm: m as the reduced divisor M of the curve's Jacobian
# scalar and divisor multiplication: multi_div
# addition of divisors: cantor's algorithm

# generate private and public key
def key_gen(ap, bp, order, k, h, f, g):
    '''
    input:
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        order: (int) order of the principal divisor
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined

    returns:
        key: (int) randomly generated secret key
        aq: (galois.Poly) first polynomial of the public key divisor in Mumford representation
        bq: (galois.Poly) second polynomial of the public key divisor in Mumford representation
    '''
    key = randint(1, order + 1)
    aq, bq = multi_div(ap, bp, key, k, h, f, g)
    return (key, (aq, bq))

# encrypt a message m represented as a divisor of the Jacobian in Mumford representation
def encrypt(am, bm, ap, bp, order, aq, bq, k, h, f, g):
    '''
    input:  
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        am: (galois.Poly) first polynomial of the message represented as a divisor of the Jacobian in Mumford representation
        bm: (galois.Poly) second polynomial of the message represented as a divisor of the Jacobian in Mumford representation
        aq: (galois.Poly) first polynomial of the public key divisor in Mumford representation
        bq: (galois.Poly) second polynomial of the public key divisor in Mumford representation
        order: (int) order of the principal divisor
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined

    returns:
        ac1: (galois.Poly) is the polynomial a of the divisor C1 in Mumford representation, which creates the ciphertext.
        bc1: (galois.Poly) is the polynomial b of the divisor C1 in Mumford representation, which creates the ciphertext.
        ac2: (galois.Poly) is the polynomial a of the divisor C2 in Mumford representation, which creates the ciphertext.
        bc2: (galois.Poly) is the polynomial b of the divisor C2 in Mumford representation, which creates the ciphertext.
    '''
    a = randint(1, order + 1)
    ac1, bc1 = multi_div(ap, bp, a, k, h, f, g)
    aac2, bbc2 = multi_div(aq, bq, a, k, h, f, g)
    ac2, bc2 = HyperellipticCurves.cantor(am, bm, aac2, bbc2, h, f, g)
    return ((ac1, bc1), (ac2, bc2))

# decrypt a cipthertext c represented as a reduced divisor of the Jacobian in Mumford representation
def decrypt(ac1, bc1, ac2, bc2, key, k, h, f, g):
    '''
    input:  
        ac1: (galois.Poly) is the polynomial a of the divisor C1 in Mumford representation, which creates the ciphertext.
        bc1: (galois.Poly) is the polynomial b of the divisor C1 in Mumford representation, which creates the ciphertext.
        ac2: (galois.Poly) is the polynomial a of the divisor C2 in Mumford representation, which creates the ciphertext.
        bc2: (galois.Poly) is the polynomial b of the divisor C2 in Mumford representation, which creates the ciphertext.
        key: (int) randomly generated secret key
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined

    returns:
        amFin: (galois.Poly) is the polynomial a of the decrypted divisor of the message represented in Mumford representation
        bmFin: (galois.Poly) is the polynomial b of the decrypted divisor of the message represented in Mumford representation
    '''
    ac1op, bc1op = (ac1, (-h-bc1) % (ac1))
    akc1, bkc1 = multi_div(ac1op, bc1op, key, k, h, f, g)
    amFin, bmFin = HyperellipticCurves.cantor(ac2, bc2, akc1, bkc1, h, f, g)
    return amFin, bmFin

# map a divisor of the Jacobian in Mumford representation to an element of the ring of Integers Z
def phi(ap, bp):
    '''
    input:  
        ap: (galois.Poly) is the polynomial a of the divisor in Mumford representation
        bp: (galois.Poly) is the polynomial b of the divisor in Mumford representation

    returns:
        modcoeffsum: Integer value after mapping
    '''
    pol_ord = ap.degree
    coeffs = list(map(int, list(ap.coeffs)))[::-1]

    modcoeffsum = 0
    for i in range(0, pol_ord):
        modcoeffsum += ((coeffs[i] % 2) * 2**i)

    return modcoeffsum

# utility function to compute modular multiplicative inverse
def MMI(A, n, s=1, t=0, N=0):
    '''
    input:
        A: (int) number to find the modular multiplicative inverse of
        n: (int) the modulus

    returns:
        A^-1: (int) the modular multiplicative inverse of A such that (A * A^-1) mod n = 1
    '''
    return (n < 2 and t%N or MMI(n, A%n, t, s-A//n*t, N or n),-1)[n<1]

# HECDSA algorithm to sign a message, derived from https://www.esat.kuleuven.be/cosic/publications/article-520.pdf
def sign(message, ap, bp, order, aq, bq, key, k, h, f, g, hashfunc=sha1):
    '''
    input:
        message: (str) The message to sign as a string
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        aq: (galois.Poly) first polynomial of the public key divisor in Mumford representation
        bq: (galois.Poly) second polynomial of the public key divisor in Mumford representation
        order: (int) order of the principal divisor
        key: (int) randomly generated secret key
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined
        hashfunc: (function of hashlib) the hash function to use, default is SHA1

    returns:
        message: (str) The message that was signed
        aqq: (galois.Poly) first polynomial of the signed message represented as a divisor in Mumford representation
        bqq: (galois.Poly) second polynomial of the signed message represented as a divisor in Mumford representation
        s: (int) The signature parameter s
    '''
    m_encode = message.encode("utf-8")
    hashed = sha1(m_encode).digest()
    int_hash = int.from_bytes(hashed, "big")

    k_rand = randint(1, order)
    while MMI(k_rand, order) == -1:
        k_rand = randint(1, order)
    
    aqq, bqq = multi_div(ap, bp, k_rand, k, h, f, g)
    print("k_rand, ", k_rand)
    s = (MMI(k_rand, order) * (int_hash + key * phi(aqq, bqq))) % order

    return (message, (aqq, bqq), s)

# HECDSA algorithm to verify a message, derived from https://www.esat.kuleuven.be/cosic/publications/article-520.pdf
def verify(message, ap, bp, order, aq, bq, aqq, bqq, s, k, h, f, g, hashfunc=sha1):
    '''
    input:
        message: (str) The message to sign as a string
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        aq: (galois.Poly) first polynomial of the public key divisor in Mumford representation
        bq: (galois.Poly) second polynomial of the public key divisor in Mumford representation
        order: (int) order of the principal divisor
        aqq: (galois.Poly) first polynomial of the signed message represented as a divisor in Mumford representation
        bqq: (galois.Poly) second polynomial of the signed message represented as a divisor in Mumford representation
        s: (int) The signature parameter s
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined
        hashfunc: (function of hashlib) the hash function to use, default is SHA1

    returns:
        Output indicating whether the signature was accepted or rejected
    '''
    m_encode = message.encode("utf-8")
    hashed = sha1(m_encode).digest()
    int_hash = int.from_bytes(hashed, "big")
    
    v1 = (MMI(s, order) * (int_hash)) % order
    v2 = (MMI(s, order) * phi(aqq, bqq)) % order

    av1, bv1 = multi_div(ap, bp, v1, k, h, f, g)
    av2, bv2 = multi_div(aq, bq, v2, k, h, f, g)
    av, bv = HyperellipticCurves.cantor(av1, bv1, av2, bv2, h, f, g)
    print(av, bv)
    print(aqq, bqq)
    if (av, bv) == (aqq, bqq):
        print("Signature accepted")
    else:
        print("Signature rejected")

# alternate HECDSA algorithm to sign a message, derived from https://doi.org/10.1109/TCSET.2006.4404442
def sign2(message, ap, bp, order, key, k, h, f, g, hashfunc=sha1):
    '''
    input:
        message: (str) The message to sign as a string
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        order: (int) order of the principal divisor
        key: (int) randomly generated secret key
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined
        hashfunc: (function of hashlib) the hash function to use, default is SHA1

    returns:
        r: (int) The signature parameter r
        s: (int) The signature parameter s
    '''
    m_encode = message.encode("utf-8")
    hashed = sha1(m_encode).digest()
    int_hash = int.from_bytes(hashed, "big")

    w = int_hash % order
    e = randint(1, order)
    u, v = multi_div(ap, bp, e, k, h, f, g)
    F_e = phi(u, v) % order
    r = (w * F_e) % order
    s = (e + key * r) % order

    return (r, s)

# alternate HECDSA algorithm to verify a message, derived from https://doi.org/10.1109/TCSET.2006.4404442
def verify2(message, ap, bp, order, aq, bq, r, s, k, h, f, g, hashfunc=sha1):
    '''
    input:
        message: (str) The message to sign as a string
        ap: (galois.Poly) first polynomial of the principal divisor in Mumford representation
        bp: (galois.Poly) second polynomial of the principal divisor in Mumford representation
        aq: (galois.Poly) first polynomial of the public key divisor in Mumford representation
        bq: (galois.Poly) second polynomial of the public key divisor in Mumford representation
        order: (int) order of the principal divisor
        r: (int) The signature parameter r
        s: (int) The signature parameter s
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve
        k: (galois.GF) Galois field over which the hyperelliptic curve is defined
        hashfunc: (function of hashlib) the hash function to use, default is SHA1

    returns:
        Output indicating whether the signature was accepted or rejected
    '''
    m_encode = message.encode("utf-8")
    hashed = sha1(m_encode).digest()
    int_hash = int.from_bytes(hashed, "big")

    w = int_hash % order
    ur1, vr1 = multi_div(ap, bp, s, k, h, f, g)
    ur2, vr2 = multi_div(aq, bq, r, k, h, f, g)
    ur, vr = HyperellipticCurves.cantor(ur1, vr1, ur2, vr2, h, f, g)
    y = phi(ur, vr) % order
    v = w * y
    print(r, v)
    if r == v:
        print("Signature accepted")
    else:
        print("Signature rejected")