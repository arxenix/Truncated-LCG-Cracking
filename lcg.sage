# Regular LCG, a pRNG
class LCG:
    # a LCG has 3 parameters, m (modulus), a (multiplier), and c (increment)
    # the seed is the starting state
    def __init__(self, m, a, c, seed):
        self.m = m
        self.a = a
        self.c = c
        self.state = seed

    # x_{n+1} = a*x_n + c (mod m)
    def next(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state


# LCG with bitshift (only output leading bits)
class TruncatedLCG(LCG):
    def __init__(self, m, a, c, seed, shift):
        self.m = m
        self.a = a
        self.c = c
        self.state = seed
        self.shift = shift

    def next(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state >> self.shift


# Java's java.util.Random, POSIX [ln]rand48, glibc [ln]rand48[_r]
class JavaRNG(TruncatedLCG):
    m = 2 ** 48
    a = 25214903917
    c = 11
    shift = 16  # bits 47..16

    def __init__(self, seed):
        self.state = seed

import os
import binascii
class PlaidRNG(TruncatedLCG):
    shift = 96
    
    def __init__(self):
        
        self.m = int(binascii.hexlify(os.urandom(16)), 16)
        self.a = 0
        self.c = 0
        self.state = 0

        while not (1 <= self.a < self.m and gcd(self.a, self.m) == 1):
            self.a = int(binascii.hexlify(os.urandom(16)), 16)
        while not (1 <= self.c < self.m and gcd(self.c, self.m) == 1):
            self.c = int(binascii.hexlify(os.urandom(16)), 16)
        while not (1 <= self.state < self.m):
            self.state = int(binascii.hexlify(os.urandom(16)), 16)


# setup our LCG
M = 2 ** 48
A = 25214903917
C = 11

NUM_OUTPUTS = 10000
y = []

lcg = JavaRNG(0xDEADB00L)
#lcg = PlaidRNG()
for i in range(NUM_OUTPUTS):
    y.append(lcg.next())


# put the above into a function
def findPolynomial(y_chunk, t=3):
    vu = 48 # num bits of modulus
    # alpha = proportion of outbit bits = 0.6666
    # beta = proportion of bits truncated = 0.3333
    av = 32 # alpha*vu = number of outbit bits
    bv = 16 # beta*vu = number of bits truncated
    # x_i = 2^(beta*vu) * y_i + z_i
    # x_i is internal state, y_i is leading bits (lcg output), z_i is trailing bits
    
    V = []
    for i in range(len(y_chunk)-t):
        Vi = []
        for j in range(t):
            Vi.append(y_chunk[i+j+1]-y_chunk[i+j])
        Vi = vector(ZZ, Vi) # QQ for rationals, RR for reals
        V.append(Vi)
        # print Vi
    NUM_VECTORS = len(V)
    
    # apply techniques from sec 2.2 to find relation sum(lambda_i * V_i) = 0
    # IMPORTANT NOTE: Sage's LLL is row-based (the basis vectors for the lattice should be the rows of the matrix)

    # construct lattice from vectors
    K = 1 # approximate relations
    KV = [K*Vi for Vi in V] # does nothing because K=1,
    lattice = Matrix(ZZ, KV)
    assert lattice.ncols() == t
    
    II = matrix.identity(NUM_VECTORS)
    lattice = lattice.augment(II)
    
    lattice_red = lattice.LLL()
    # LLL result polynomial coefficients
    lambdas = lattice_red[0][t:]
    
    # check to make sure it works
    x = PolynomialRing((Integers(M)), 'x').gen()
    poly = 0
    for i in range(len(lambdas)):
        poly += lambdas[i] * x^i
    if poly(A) == 0:
        # we successfully found a polynomial such that P(A)=0 (mod M)
        return lambdas
    else:
        return None

        

# split output into chunks, make a polynomial for each
SIZE_PARAM = 10
T_PARAM = 3
VECTOR_LEN = SIZE_PARAM - T_PARAM

NUM_POLYS = 8

# total amount of data consumed is NUM_POLYS * SIZE_PARAM

#y_chunks = [y[i:i+SIZE_PARAM] for i in range(0, len(y), SIZE_PARAM)]
basis_zero_vector = vector(ZZ, [0]*VECTOR_LEN)

polys = []
for i in range(NUM_POLYS):
    y_chunk = y[i:i+SIZE_PARAM]
    poly = findPolynomial(y_chunk, T_PARAM)
    if poly is not None:
        polys.append(poly)

#Put the polynomials in a matrix, apply LLL
polyM = Matrix(ZZ, polys)
polyM_red = polyM.LLL()

# remove zero rows and form new lattice
new_basis = []
for row in polyM_red:
    if row!=basis_zero_vector:
        new_basis.append(row)

new_polyM = Matrix(ZZ, new_basis)

# find the determinant - it is not a square matrix so we utilize normal equations (A^T * A)
new_polyM_normal = new_polyM.transpose() * new_polyM
det = new_polyM_normal.determinant()
det = ZZ(det).nth_root(2, truncate_mode=1)[0] #integer sqrt

print "Computed determinant: ", det
print "actual M:", lcg.m
print "outputs consumed: ", SIZE_PARAM+NUM_POLYS



# Compute A
nrows = new_polyM.nrows()
ncols = new_polyM.ncols()

#K = 2**32
K = M * 2**(NUM_POLYS+3)
for i in range(nrows):
    for j in range(ncols):
        if j>=2:
            new_polyM[i,j] *= K


new_polyM_red = new_polyM.LLL()

sublattice = sublattice.matrix_from_rows([0,1])

a1 = sublattice[0,0]
a2 = sublattice[0,1]
b1 = sublattice[1,0]
b2 = sublattice[1,1]

# solve a1*x+b1*y = 1
d,u,v = xgcd(a2,b2)
assert d==1
assert u*a2+v*b2 == 1
coeff = u*a1+v*b1
print "computed a:", -coeff
print "actual a:", lcg.a



# TODO -- find C
