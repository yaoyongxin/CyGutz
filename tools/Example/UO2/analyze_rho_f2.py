'''
Analyze the f2 block of the local many-body density matrix of the Mott
solution.
'''


import numpy as np
from scipy.sparse.linalg import eigsh
import h5py
from pyglib.io.h5io import get_csr_matrix
from pyglib.math.matrix_util import sym_csr_matrix
from pyglib.symm.atom_symm import get_characters_espace

f = h5py.File('glog.h5', 'r')
rho = get_csr_matrix(f, "/Impurity_1/RHO")
f.close()

# Take the f2 block.
rho = rho[:91, :91]

print(' Total weight =',np.sum(rho.diagonal()))

vals, vecs = eigsh(rho, k=12)
print(' Eigen-values before symmetrization:', vals)

f = h5py.File('hilbert_space_rotations.h5', 'r')

# Read rotations.
rotations = []
for i in range(f["Impurity_1/val_block=2/dim_rotations"][...]):
    rotations.append(get_csr_matrix(f, \
            "Impurity_1/val_block=2/rotation_{}".format(i)))

rho = sym_csr_matrix(rho, rotations)
vals, vecs = eigsh(rho, k=12)
print(' Eigen-values after symmetrization:', vals)

# calculate charater
chi, _ = get_characters_espace(rotations, vecs[:, :3].T)
print(' Character:')
print(''.join(" {:<4.1f}".format(x) for x in chi))

np.savetxt('gamma_5_vec.txt', vecs[:, :3].T, fmt='%6.3f')
