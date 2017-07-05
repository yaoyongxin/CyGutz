from pyglib.symm.angular_momentum_1p import get_L_vector_CH


def dumptxt_matrix(f, a):
    for a1 in a:
        f.write(('{:12.6f}'*len(a1)+'\n').format(*a1.real))
    f.write('\n')
    for a1 in a:
        f.write(('{:12.6f}'*len(a1)+'\n').format(*a1.imag))
    f.write('\n')


with open('param.inp', 'r') as f:
    L = int(f.readlines()[0].split()[0])


Lx, Ly, Lz = get_L_vector_CH(L)

f = open('g_lxyz.txt', 'w')
f.write('Lx comp_sph_harm\n')
dumptxt_matrix(f, Lx)
f.write('Ly comp_sph_harm\n')
dumptxt_matrix(f, Ly)
f.write('Lz comp_sph_harm\n')
dumptxt_matrix(f, Lz)

# get complex to real sph_harm transformation
from pyglib.symm.angular_momentum_1p import get_complex_to_real_sph_harm
c2r = get_complex_to_real_sph_harm(L)
Lx = c2r.T.conj().dot(Lx).dot(c2r)
Ly = c2r.T.conj().dot(Ly).dot(c2r)
Lz = c2r.T.conj().dot(Lz).dot(c2r)

f.write('Lx real_sph_harm\n')
dumptxt_matrix(f, Lx)
f.write('Ly real_sph_harm\n')
dumptxt_matrix(f, Ly)
f.write('Lz real_sph_harm\n')
dumptxt_matrix(f, Lz)

f.close()
