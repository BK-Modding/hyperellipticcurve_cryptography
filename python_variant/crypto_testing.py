from galois import Poly, GF

import HyperellipticCurves
import HyperellipticCryptography

k = HyperellipticCurves.set_field(2, [1, 0, 1, 0, 0, 1], "reg")
h = Poly([1, 1, 0], field=k)
f = Poly([1, 1, 0, 1, 0, 1], field=k)
g = 2
ap = Poly([1, 21], field=k)
bp = Poly([5], field=k)

order = 482 # HyperellipticCryptography.ord_div(ap, bp, k, h, f, g)

key, (aq, bq) = HyperellipticCryptography.key_gen(ap, bp, order, k, h, f, g)

# ap, bp: polynomials of the principal divisor
# aq, bq: key polynomials
# am, bm: m as the reduced divisor M of the curve's Jacobian

am = Poly([1, 27, 0], field=k)
bm = Poly([25, 1], field=k)

print("Message Divisor:\t", am, bm)

(ac1, bc1), (ac2, bc2) = HyperellipticCryptography.encrypt(am, bm, ap, bp, order, aq, bq, k, h, f, g)

print("Encrypted Message Divisor:\t", (ac1, bc1), (ac2, bc2))

amFin, bmFin = HyperellipticCryptography.decrypt(ac1, bc1, ac2, bc2, key, k, h, f, g)

print("Decrypted Message Divisor:\t", amFin, bmFin)

print("*" * 50)
print("Signing and Verification")

message = "Hello there"
m, (aqq, bqq), s = HyperellipticCryptography.sign(message, ap, bp, order, aq, bq, key, k, h, f, g)

HyperellipticCryptography.verify(message, ap, bp, order, aq, bq, aqq, bqq, s, k, h, f, g)

# r, s = HyperellipticCryptography.sign(message, ap, bp, order, key, k, h, f, g)

# HyperellipticCryptography.verify(message, ap, bp, order, aq, bq, r, s, k, h, f, g)