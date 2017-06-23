from __future__ import print_function
import h5py, numpy
from pyglib.symm.angular_momentum_1p import get_L_vector_CH
from pyglib.symm.atom_symm import get_representation, get_product_table
from pyglib.math.matrix_util import sym_array


# d complex spherical Harmonics
L = get_L_vector_CH(2)

# test convention
# utrans = numpy.diag([1, 1, 1, -1, 1.])
# L = [utrans.conj().T.dot(L1).dot(utrans) for L1 in L]

print('Lx_r\n',L[0].real)
print('Lx_i\n',L[0].imag)

print('Ly_r\n',L[1].real)
print('Ly_i\n',L[1].imag)

print('Lz_r\n',L[2].real)
print('Lz_i\n',L[2].imag)



# get Lie parameters and one-particle dm to test
with h5py.File('data_complex_spherical_harmonics.h5', 'r') as f:
    nks = f['/NKS'][()].T
    lie_params = f['/lie_even_params'][()]



# take one-spin component
nks = nks[::2, ::2]

# rotations
R_list = get_representation(L, lie_params)

# check product table
ptable = get_product_table(R_list)
print('ptable\n',ptable)

print('nks\n',nks.real)


nks_sym = sym_array(nks, R_list)
print('nks_sym\n',nks_sym.real)

print(numpy.max(numpy.abs(nks.imag)), numpy.max(numpy.abs(nks_sym.imag)))


print('max symm error = {}'.format(numpy.max(numpy.abs(nks-nks_sym))))
