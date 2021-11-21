import numpy as np
from galois import GF, Poly, poly_gcd # importing galois field operations from the galois package

# converts integer to the 0th degree polynomial over the specified Galois field
def i_p(i, field):
    '''
    input:
        i: (int) integer to convert to 0th degree polynomial
        field: (galois.GF) Galois field

    returns: 
        galois.Poly (0th degree polynomial over the specified Galois field)
    '''
    return Poly([i], field=field)

# sets the galois field with prime characteristic p, coefficients of the field's irreducible polynomial n,
def set_field(p, n):
    '''
    input:
        p: (int) prime characteristic
        n: (list) coefficients of the field's irreducible polynomial

    returns: 
        K: (galois.GF) Galois field of prime characteristic p, extension p^(len(n) - 1), and irreducible polynomial n
    '''
    K = GF(p**(len(n) - 1), irreducible_poly=Poly(list(reversed(n))))
    return K

# finds the points of a hyperelliptic curve defined over a galois field
def find_points(k, h, f):
    '''
    input:
        k: (galois.GF) Galois field object
        h: (galois.Poly) non-monic polynomial of the hyperelliptic curve such that y^2 + h * y = f
        f: (galois.Poly) monic polynomial of the hyperelliptic curve such that y^2 + h * y = f

    returns: 
        points: (set) sorted set of point pairs that satisfy the curve equation, where each point is an element of the galois field
    '''
    points = set()
    FieldSize = len(k.Elements())
    for i in range(0, FieldSize):
        for j in range(0, FieldSize):
            fa = f(i)
            ha = h(i)
            if fa == (k.Elements()[j]**2 + k.Elements()[j] * ha):
                points.add((i, j))

    return sorted(points)

# Improved algorithm to find the points of a hyperelliptic curve defined over a galois field
def find_points_im(k, h, f):
    '''
    input:
        k: (galois.GF) Galois field object
        h: (galois.Poly) non-monic polynomial of the hyperelliptic curve such that y^2 + h * y = f
        f: (galois.Poly) monic polynomial of the hyperelliptic curve such that y^2 + h * y = f

    returns: 
        points: (set) sorted set of point pairs that satisfy the curve equation, where each point is an element of the galois field
    '''
    points = set()
    FieldSize = len(k.Elements())
    for i in range(0, FieldSize):
        for j in range(0, FieldSize):
            fa = f(i)
            ha = h(i)
            if fa == (k.Elements()[j]**2 + k.Elements()[j] * ha):
                points.add((i, j))
                points.add((i, list(k.Elements()).index(-k.Elements()[j]-ha)))
                break
    
    return sorted(points)

# number of points on a hyperelliptic curve
def number_of_points(k, h, f):
    '''
    input:
        k: (galois.GF) Galois field object
        h: (galois.Poly) non-monic polynomial of the hyperelliptic curve such that y^2 + h * y = f
        f: (galois.Poly) monic polynomial of the hyperelliptic curve such that y^2 + h * y = f

    returns: 
        (int) number point pairs that satisfy the curve equation, where each point is an element of the galois field
    '''
    return 1 + len(find_points(k, h, f))

# Order of a rational polynomial function at an ordinary point
def ord_ordinary(h, f, P, a, b):
    '''
    input:
        a: (galois.Poly) polynomial a of the polynomial function G = a(x) − b(x)y.
        b: (galois.Poly) polynomial b of the polynomial function G = a(x) − b(x)y.
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        P: (list) ordinary point of the hyperelliptic curve, given as the list of its two coordinates

    returns: 
        OrdG: (int) order of the polynomial at the ordinary point
    '''
    k = h.field
    pra, prb, pr = i_p(0, k), i_p(0, k), i_p(0, k)
    ra, rb, s = 0, 0, 0

    i = 1
    while pra == i_p(0, k):
        pra = a % Poly([1, -P[0]], field=k) ** i
        ra += 1

        i += 1

    i = 1
    while prb == i_p(0, k):
        prb = b % Poly([1, -P[0]], field=k) ** i
        rb += 1
        
        i += 1

    if ra > rb:
        r = rb - 1
    else:
        r = ra - 1
    aa = a / Poly([1, -P[0]], field=k) ** r
    bb = b / Poly([1, -P[0]], field=k) ** r
    NormG0 = aa**2 + aa * bb * h - bb ** 2 * f
    
    i = 1
    while pr == i_p(0, k):
        pr = NormG0 % Poly([1, -P[0]], field=k) ** i
        s += 1

        i += 1

    OrdG = r + s - 1
    
    return OrdG

# Order of a rational polynomial function at a special point
def ord_special(h, f, P, a, b):
    '''
    input:
        a: (galois.Poly) polynomial a of the polynomial function G = a(x) − b(x)y.
        b: (galois.Poly) polynomial b of the polynomial function G = a(x) − b(x)y.
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        P: (list) special point of the hyperelliptic curve, given as the list of its two coordinates

    returns: 
        OrdG: (int) order of the polynomial at the special point
    '''
    k = h.field
    pra, prb, pr = i_p(0, k), i_p(0, k), i_p(0, k)
    ra, rb, s = 0, 0, 0

    i = 1
    while pra == i_p(0, k):
        pra = a % Poly([1, -P[0]], field=k) ** i
        ra += 1

        i += 1

    i = 1
    while prb == i_p(0, k):
        prb = b % Poly([1, -P[0]], field=k) ** i
        rb += 1
        
        i += 1

    if ra > rb:
        r = rb - 1
    else:
        r = ra - 1
    aa = a / Poly([1, -P[0]], field=k) ** r
    bb = b / Poly([1, -P[0]], field=k) ** r
    NormG0 = aa**2 + aa * bb * h - bb ** 2 * f
    
    i = 1
    while pr == i_p(0, k):
        pr = NormG0 % Poly([1, -P[0]], field=k) ** i
        s += 1

        i += 1

    OrdG = 2 * r + s - 1
    
    return OrdG

# Order of a rational polynomial function at the point at infinity
def ord_infinity(g, a, b):
    '''
    input:
        a: (galois.Poly) polynomial a of the polynomial function G = a(x) − b(x)y.
        b: (galois.Poly) polynomial b of the polynomial function G = a(x) − b(x)y.
        g: (int) genus of the hyperelliptic curve

    returns: 
        OrdG: (int) order of the polynomial at the point at infinity
    '''
    OrdG = -max(2 * a.degree, 2 * b.degree + 2 * g + 1)
    return OrdG

# Norm of a polynomial function
def norm_g(a, b, h, f):
    '''
    input:
        a: (galois.Poly) polynomial a of the polynomial function G = a(x) − b(x)y.
        b: (galois.Poly) polynomial b of the polynomial function G = a(x) − b(x)y.
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f

    returns: 
        NG: (galois.Poly) The norm of the polynomial function G
    '''
    aa = a
    bb = b
    NG = aa**2 + aa*bb*h - bb**2*f
    return NG

# find GCD of two polynomials defined over a Galois field via Extended Euclidean Algorithm
def poly_gcd_linear_combination(pol1, pol2):
    '''
    input:
        pol1: (galois.Poly) first polynomial
        pol2: (galois.Poly) second polynomial

    returns:
        d: (galois.Poly) GCD(pol1, pol2)
        amul: (galois.Poly) polynomial of a linear combination such that amul * pol1 + bmul * pol2 = d
        bmul: (galois.Poly) polynomial of a linear combination such that amul * pol1 + bmul * pol2 = d
    '''
    k = pol1.field
    d, a, b = poly_gcd(pol1, pol2)
    target_mul = None

    a = i_p(a, field=k) if type(a) == int else a
    b = i_p(b, field=k) if type(b) == int else b

    irred_lin_comb = a * pol1 + b * pol2
    for elem in range(0, k.order): # k.Elements()
        mul = i_p(int(elem), field=k)
        if mul * irred_lin_comb == d:
            target_mul = mul
            break

    amul, bmul = a * target_mul, b * target_mul
    return d, amul, bmul

# algorithm to add two reduced divisors in Mumford representation of a hyperelliptic curve
def add_div(a1, b1, a2, b2, h, f):
    '''
    input:
        a1: (galois.Poly) first polynomial of the first reduced divisor in Mumford representation
        b1: (galois.Poly) second polynomial of the first reduced divisor in Mumford representation
        a2: (galois.Poly) first polynomial of the second reduced divisor in Mumford representation
        b2: (galois.Poly) second polynomial of the second reduced divisor in Mumford representation
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f

    returns:
        a: (galois.Poly) first polynomial of the resultant divisor (non-reduced) in Mumford representation such that (a1, b1) + (a2, b2) = (a, b)
        b: (galois.Poly) second polynomial of the resultant divisor (non-reduced) in Mumford representation such that (a1, b1) + (a2, b2) = (a, b)
    '''
    d1, e1, e2 = poly_gcd_linear_combination(a1, a2)
    d, c1, c2 = poly_gcd_linear_combination(d1, b1 + b2 + h)
    s1 = c1 * e1
    s2 = c1 * e2
    s3 = c2
    a = (a1 * a2) / d**2
    bb = (s1 * a1 * b2 + s2 * a2 * b1 + s3 * (b1 * b2 + f)) / d
    b = bb % a
    
    return a, b

# algorithm to double a reduced divisor (add to itself) in Mumford representation of a hyperelliptic curve
def double_div(a1, b1, h, f):
    '''
    input:
        a1: (galois.Poly) first polynomial of the reduced divisor in Mumford representation
        b1: (galois.Poly) second polynomial of the reduced divisor in Mumford representation
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f

    returns:
        a: (galois.Poly) first polynomial of the resultant divisor (non-reduced) in Mumford representation such that (a1, b1) + (a1, b1) = (a, b)
        b: (galois.Poly) second polynomial of the resultant divisor (non-reduced) Mumford representation such that (a1, b1) + (a1, b1) = (a, b)
    '''
    d, s1, s3 = poly_gcd_linear_combination(a1, (b1 + b1 + h))
    a = (a1 * a1) / d ** 2
    b = ((s1 * a1 * b1 + s3 * (b1 * b1 + f)) / d) % a
    
    return a, b

# convert a non-reduced divisor of a hyperelliptic curve in Mumford representation to its reduced form
def red_div(a, b, g, h, f):
    '''
    input:
        a: (galois.Poly) first polynomial of the non-reduced divisor in Mumford representation
        b: (galois.Poly) second polynomial of the non-reduced divisor in Mumford representation
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: (int) genus of the hyperelliptic curve

    returns:
        aa: (galois.Poly) first polynomial of the reduced divisor in Mumford representation
        bb: (galois.Poly) second polynomial of the reduced divisor in Mumford representation
    '''
    aa = a
    bb = b
    while aa.degree > g:
        a1 = (f - bb*h - bb**2) / aa
        b1 = (-h-bb) % a1
        aa=a1
        bb=b1

    ex = aa.degree
    inv = aa.field(1) / list(aa.coeffs)[list(aa.degrees).index(ex)]
    aa = i_p(int(inv), field=aa.field) * aa

    return aa, bb

# Function to compute the order of the Hyperelliptic Curve's Jacobian Group
def jacobian_ord(k, M1, M2):
    '''
    input:
        k: (galois.GF) Galois field over which the Jacobian's Hyperelliptic curve is defined
        M1: (int) is the number of rational points of the Galois field GF(p) that satisfy the hyperelliptic curve's equation.
        M2: (int) is the number of rational points of the Galois field GF(p^2) that satisfy the hyperelliptic curve's equation.

    returns:
        JacobOrd: (int) Order of the Hyperelliptic Curve's Jacobian Group
    '''
    p = k.characteristic
    a1 = M1 - 1 - p
    a2 = (M2 - 1 - p**2 + a1**2) / 2
    gamma1 = (-a1 - np.lib.scimath.sqrt(a1**2 - 4 * (a2 - 2 * p))) / 2
    gamma2 = (-a1 + np.lib.scimath.sqrt(a1**2 - 4 * (a2 - 2 * p))) / 2
    alpha1 = (gamma1 - np.lib.scimath.sqrt(gamma1**2 - 4 * p)) / 2
    alpha2 = (gamma2 - np.lib.scimath.sqrt(gamma2**2 - 4 * p)) / 2

    n = k.degree
    JacobOrd = int(abs(1 - alpha1**n)**2 * abs(1 - alpha2**n)**2)

    return JacobOrd

# DivToPol[K_,k_,h_,f_,x_,Div_]:=(
# l=Length[Div];
# ad=Div[[All,2]];
# aa=k[1];
# For[i=1,i<Length[ad]+1,i++,aa=aa*(x-ad[[i]])^(Div[[i]][[1]])];
# Pola=Expand[aa];
# ap=Apart[1/(aa),x];
# app={};
# For[i=1,i<Length[ap]+1,i++,app=AppendTo[app,ap[[i]]]];
# pow=Div[[All,1]];
# m=Max[pow];
# For[q=1,q<m+1,q++,For[j=1,j<Length[app]+1,j++,For[i=1,i<Length[app]+1,i++,If[i!=j,If[PolynomialRemainder[Denominator[app[[j]]],Denominator[app[[i]]],x]==0,{tog=Together[app[[i]]+app[[j]]],app=ReplacePart[app,j->tog],app=Delete[app,i],i=i-1}]]]]];
# nDiv={};
# For[i=1,i<Length[Div]+1,i++,{{p,q}={x-Div[[i]][[2]],Div[[i]][[3]]},{pn,qn}={p,q},If[Div[[i]][[1]]>1,For[j=1,j<Div[[i]][[1]],j++,{pn,qn}=AddDiv[p,q,pn,qn,h,f,x]]],AppendTo[nDiv,{pn,qn}]}];
# den={};
# For[i=1,i<Length[app]+1,i++,AppendTo[den,Denominator[app[[i]]]]];
# den=Expand[den];
# a2={};
# For[i=1,i<Length[den]+1,i++,For[j=1,j<Length[nDiv]+1,j++,If[den[[i]]==nDiv[[j]][[1]],AppendTo[a2,nDiv[[j]][[2]]*app[[i]]]]]];
# ssuma={};
# For[i=1,i<Length[a2]+1,i++,AppendTo[ssuma,a2[[i]]*aa]];
# ssuma=Expand[ssuma];
# suma=0;
# For[i=1,i<Length[ssuma]+1,i++,suma=Expand[suma+ssuma[[i]]]];
# b=PolynomialRemainder[suma,Pola,x];
# {Pola,b}
# )

# function in progress - not necessary for cryptography
def div_to_pol(k, h, f, Div):
    l = len(Div)
    ad = [row[1] for row in Div]
    aa = i_p(k.Elements()[1], k)
    # For[i = 1, i < Length[ad] + 1, i++, aa = aa*(x - ad[[i]])^(Div[[i]][[1]])];
    for i in range(len(ad)):
        aa = aa * (Poly([1, -ad[i]], field=k))**(Div[i][0])
    # ...

# function in progress - not necessary for cryptography
def pol_to_div(k, a, b):
    pass

# Cantor's algorithm to add two reduced divisors in Mumford representation
def cantor(a1, b1, a2, b2, h, f, g):
    '''
    input:
        a1: (galois.Poly) first polynomial of the first reduced divisor in Mumford representation
        b1: (galois.Poly) second polynomial of the first reduced divisor in Mumford representation
        a2: (galois.Poly) first polynomial of the second reduced divisor in Mumford representation
        b2: (galois.Poly) second polynomial of the second reduced divisor in Mumford representation
        h: (galois.Poly) monic polynomial h of the hyperelliptic curve such that y^2 + h * y = f 
        f: (galois.Poly) non-monic polynomial f of the hyperelliptic curve such that y^2 + h * y = f
        g: genus of the Hyperelliptic curve

    returns:
        aa: (galois.Poly) first polynomial of the resultant reduced divisor in Mumford representation such that (a1, b1) + (a2, b2) = (a, b)
        bb: (galois.Poly) second polynomial of the resultant reduced divisor in Mumford representation such that (a1, b1) + (a2, b2) = (a, b)
    '''
    a, b = add_div(a1, b1, a2, b2, h, f)
    aa, bb = red_div(a, b, g, h, f)
    return aa, bb
