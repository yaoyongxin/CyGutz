from __future__ import print_function
import h5py, numpy
from pyglib.symm.angular_momentum_1p import get_L_vector_CH
from pyglib.symm.atom_symm import get_representation, get_product_table
from pyglib.math.matrix_util import sym_array


# d complex spherical Harmonics
L = get_L_vector_CH(2)

# to complex sph_harmonics
from pyglib.symm.angular_momentum_1p import get_complex_to_real_sph_harm
c2r = get_complex_to_real_sph_harm(2)

for imp in range(4):
    # get Lie parameters and one-particle dm of FeO to test
    with h5py.File('data_complex_spherical_harmonics.h5', 'r') as f:
        lie_params = f['/IMPURITY_{}_lie_even_params'.format(imp+1)][()]

    # rotations
    R_list = get_representation(L, lie_params)

    # check product table
    ptable = get_product_table(R_list)
    print('ptable\n',ptable)

    # from vasp, in real sph_harm basis
    # convert to <Ylm | \rho | Ylm'>
    nks = numpy.loadtxt('dm_ldau_{}.txt'.format(imp+1)).T

    print('nks.real in real_sph_harm basis:')
    for row in nks:
        print(('{:10.4f}'*len(row)).format(*row.real))

    print('nks.imag in real_sph_harm basis:')
    for row in nks:
        print(('{:10.4f}'*len(row)).format(*row.imag))

    nks = c2r.dot(nks).dot(c2r.T.conj())

    print('nks.real in comp_sph_harm basis:')
    for row in nks:
        print(('{:10.4f}'*len(row)).format(*row.real))

    print('nks.imag in comp_sph_harm basis:')
    for row in nks:
        print(('{:10.4f}'*len(row)).format(*row.imag))

    for i,rot in enumerate(R_list):
        err = nks.dot(rot) -  rot.dot(nks)
        print('{} maxerr of commutation = {}'.format(i,
                numpy.max(numpy.abs(err))))

