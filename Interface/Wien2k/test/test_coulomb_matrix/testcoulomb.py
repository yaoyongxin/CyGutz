import numpy as np
from itertools import product
from pyglib.mbody.coulomb_matrix import U_matrix



def dumptxt_v2e(fname, L, umat):
    with open(fname, 'w') as f:
        for i,j,k,l in product(range(2*L+1), repeat=4):
            if abs(umat[2*i,2*k,2*j,2*l])<1.e-8:
                continue
            f.write('{:3d}{:3d}{:3d}{:3d}{:14.8f} m,m2,m1,m3,Vee\n'.format(
                    i-L,j-L,k-L,l-L,umat[2*i,2*k,2*j,2*l]))


with open('param.inp', 'r') as f:
    L, U, J = f.readlines()[0].split()

L, U, J = int(L), float(U), float(J)
umat, _, _ = U_matrix('slater', L, U_int=U, J_hund=J)
assert np.all(np.abs(umat.imag) < 1.e-10), ' Error: existing imag part!'
umat = umat.real

dumptxt_v2e('g_v2e.txt', L, umat)

with open('param2.inp', 'r') as f:
    line = f.readlines()[0].split()
    L = int(line[0])
    F = map(float, line[1:])

umat2, _, _ = U_matrix('slater', L, radial_integrals=F)
assert np.all(np.abs(umat2.imag) < 1.e-10), ' Error: existing imag part!'
umat2 = umat2.real

dumptxt_v2e('g2_v2e.txt', L, umat2)

print('max diff = {}'.format(np.max(np.abs(umat-umat2))))
