from __future__ import print_function
# analyze the embeding hamiltonian
import os, h5py, numpy, shutil
from scipy.linalg import block_diag
from pyglib.math.matrix_util import \
        get_utrans_orbfast_supdn_to_spin_fast_supdn


def h5gen_embedh_spin_updn(imp=1, h5wrt=False):
    # embedding hamiltonian parameters
    # upper case path indicates fortran convention.
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
        if not os.path.isfile('V2E_{}.inp'.format(imp)) or h5wrt:
            v2e = f['/V2E'][()].T
        else:
            v2e = None

    # default basis to symmetry adapted basis transformation
    with h5py.File('GPARAM.h5', 'r') as f:
        u_db2sab = f['/IMPURITY_{}/DB_TO_SAB'.format(imp)][()].T

    u_sf2of = get_utrans_orbfast_supdn_to_spin_fast_supdn( \
            u_db2sab.shape[0]).T
    u_sf2sab = u_sf2of.dot(u_db2sab)
    u_sab2sf = u_sf2sab.T.conj()

    # basis transformation
    daalpha_ = u_sab2sf.T.dot(daalpha).dot(u_sab2sf.conj())
    lambdac_ = u_sab2sf.T.conj().dot(lambdac).dot(u_sab2sf)
    h1e_ = u_sab2sf.T.conj().dot(h1e).dot(u_sab2sf)

    print(h1e_[0,0]); quit()

    if v2e is not None:
        # tensordot, switch orders?
        v2e_ = numpy.tensordot(numpy.tensordot(numpy.tensordot( \
                numpy.tensordot(v2e, \
                u_sf2sab, axes=([0], [1])), \
                u_sab2sf, axes=([0], [0])), \
                u_sf2sab, axes=([0], [1])), \
                u_sab2sf, axes=([0], [0]))

        # double check
        #    v2e_ = numpy.einsum('ijkl,pi,jq,rk,ls->pqrs', v2e, u_sf2sab, \
        #            u_sab2sf, u_sf2sab, u_sab2sf)
    else:
        v2e_= None

    if h5wrt:
        shutil.copy('EMBED_HAMIL_{}.h5'.format(imp), \
                'EMBED_HAMIL_CSH_{}.h5'.format(imp))
        with h5py.File('EMBED_HAMIL_CSH_{}.h5'.format(imp), 'a') as f:
            f['/D'][()] = daalpha_.T
            f['/LAMBDA'][()] = lambdac_.T
            f['/H1E'][()] = h1e_.T
            f['/V2E'][()] = v2e_.T

    return h1e_, lambdac_, daalpha_, v2e_


def wrt_text_rembed(h1e, lambdac, daalpha, v2e, imp=1):
    '''write the real embedding hamiltonian parameters in text format.
    '''
    norb = h1e.shape[0]
    v1e = block_diag(h1e, -lambdac)
    v1e[:norb, norb:] = daalpha.T
    v1e[norb:, :norb] = daalpha.conj()
    max_imag = numpy.max(numpy.abs(v1e.imag))

    if v2e is not None:
        max_imag = max(max_imag, numpy.max(numpy.abs(v2e.imag)))
    if max_imag > 1.e-7:
        print(' maximal imaginary part = {}'.format(max_imag))
        raise AssertionError(' embedding hamiltonian not real!')

    # write one-body part
    with open('H1E_{}.inp'.format(imp), 'w') as f:
        for i, v1 in enumerate(v1e):
            for j, v12 in enumerate(v1):
                if numpy.abs(v12) > 1.e-7:
                    f.write('{:2d} {:2d} {:16.12f}\n'.format(\
                            i+1, j+1, v12.real))
    # write two-body part
    if v2e is not None:
        with open('V2E_{}.inp'.format(imp), 'w') as f:
            for i, v1 in enumerate(v2e):
                for j, v12 in enumerate(v1):
                    for k, v123 in enumerate(v12):
                        if i == k:
                            continue
                        for l, v1234 in enumerate(v123):
                            if j == l:
                                continue
                            if numpy.abs(v1234) > 1.e-8:
                                # order of cd_i cd_k c_l c_j
                                # absorbing 1/2 factor
                                f.write('{:2d} {:2d} {:2d} {:2d} {:16.12f}\n'\
                                        .format(i+1, k+1, l+1, j+1, \
                                        v1234.real/2))



if __name__ == '__main__':
    h1e, lambdac, daalpha, v2e = h5gen_embedh_spin_updn(h5wrt=True)
    wrt_text_rembed(h1e, lambdac, daalpha, v2e)
