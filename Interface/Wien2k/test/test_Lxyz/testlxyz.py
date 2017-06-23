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
Lp = Lx + 1.j*Ly
Ln = Lx - 1.j*Ly

f = open('g_lxyz.txt', 'w')
dumptxt_matrix(f, Lx)
dumptxt_matrix(f, Ly)
dumptxt_matrix(f, Lz)
dumptxt_matrix(f, Lp)
dumptxt_matrix(f, Ln)
