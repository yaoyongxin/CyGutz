import h5py
import sys
from pyglib.symm.angular_momentum_1p import get_L_vector_CH
from pyglib.symm.atom_symm import get_JR_representation

irot = int(sys.argv[1])
l = int(sys.argv[2])

with h5py.File('rotylm1.h5', 'r') as f:
    rotmat = f['/rotations'][irot]

with open('param.inp', 'w') as f:
    f.write('{}\n'.format(l))
    for row in rotmat:
        f.write('{} {} {}\n'.format(*row))

Lx, ly, Lz = get_L_vector_CH(l)
rotylm = get_JR_representation([rotmat], [Lx, ly, Lz])[0]

with open('g_rotylm.txt', 'w') as f:
    for row in rotylm:
        for elem in row:
            f.write('{:10.6f}{:10.6f}  '.format(elem.real, elem.imag))
        f.write('\n')
