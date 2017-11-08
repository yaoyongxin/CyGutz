from __future__ import print_function
# analyze the embeding hamiltonian
import os, h5py, numpy, shutil
from scipy.linalg import block_diag
from pyglib.math.matrix_util import \
        get_utrans_orbfast_supdn_to_spin_fast_supdn


def twrite_2d_array_real(a, fname, shreshold=1.e-7):
    '''write 2d array (real) in text form of i, j, a_ij.
    '''
    with open(fname, 'w') as f:
        for i, a1 in enumerate(a):
            for j, a12 in enumerate(a1):
                if numpy.abs(a12) > shreshold:
                    f.write('{:2d} {:2d} {:16.12f}\n'.format(\
                            i+1, j+1, a12.real))


def twrite_u_matrix_real(v2e, fname, shreshold=1.e-7):
    '''write u-matrix in text form of i, k, l, j, v2e_(ij)(kl).
    '''
    with open(fname, 'w') as f:
        for i, v1 in enumerate(v2e):
            for j, v12 in enumerate(v1):
                for k, v123 in enumerate(v12):
                    if i == k:
                        continue
                    for l, v1234 in enumerate(v123):
                        if j == l:
                            continue
                        if numpy.abs(v1234) > shreshold:
                            # order of cd_i cd_k c_l c_j
                            # absorbing 1/2 factor
                            f.write('{:2d} {:2d} {:2d} {:2d} {:16.12f}\n'\
                                    .format(i+1, k+1, l+1, j+1, \
                                    v1234.real/2))


def get_u_sab2cshupdn(imp=1):
    '''get unitary transformation from symmetry-adapted basis
    to complex spherical harmonics basis with
    faster spin index ordered as up and down.
    '''
    # default basis to symmetry adapted basis transformation
    with h5py.File('GPARAM.h5', 'r') as f:
        u_db2sab = f['/IMPURITY_{}/DB_TO_SAB'.format(imp)][()].T

    u_sf2of = get_utrans_orbfast_supdn_to_spin_fast_supdn( \
            u_db2sab.shape[0]).T
    u_sf2sab = u_sf2of.dot(u_db2sab)
    u_sab2sf = u_sf2sab.T.conj()
    return u_sab2sf


def h5gen_embedh_spin_updn(imp=1, h5wrt=False):
    # embedding hamiltonian parameters
    # upper case path indicates fortran convention.
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
        if not os.path.isfile('V2E_{}.INP'.format(imp)) or h5wrt:
            v2e = f['/V2E'][()].T
        else:
            v2e = None

    u_sab2sf = get_u_sab2cshupdn(imp=imp)
    # basis transformation
    daalpha_ = u_sab2sf.T.dot(daalpha).dot(u_sab2sf.conj())
    lambdac_ = u_sab2sf.T.conj().dot(lambdac).dot(u_sab2sf)
    h1e_ = u_sab2sf.T.conj().dot(h1e).dot(u_sab2sf)

    if v2e is not None:
        # tensordot, switch orders?
        v2e_ = numpy.tensordot(numpy.tensordot(numpy.tensordot( \
                numpy.tensordot(v2e, \
                u_sab2sf.T.conj(), axes=([0], [1])), \
                u_sab2sf, axes=([0], [0])), \
                u_sab2sf.T.conj(), axes=([0], [1])), \
                u_sab2sf, axes=([0], [0]))

        # double check
        #    v2e_ = numpy.einsum('ijkl,pi,jq,rk,ls->pqrs', v2e, \
        #            u_sab2sf.T.conj(), u_sab2sf, \
        #            u_sab2sf.T.conj(), u_sab2sf)
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


def get_whole_h1e(h1e, lambdac, daalpha):
    '''get the total one-body part of the embedding hamiltonian
    given the components of impurity, bath and hybridization.
    '''
    norb = h1e.shape[0]
    v1e = block_diag(h1e, -lambdac)
    v1e[:norb, norb:] = daalpha.T
    v1e[norb:, :norb] = daalpha.conj()
    return v1e


def wrt_text_rembed(h1e, lambdac, daalpha, v2e, imp=1):
    '''write the real embedding hamiltonian parameters in text format.
    '''
    v1e = get_whole_h1e(h1e, lambdac, daalpha)

    max_imag = numpy.max(numpy.abs(v1e.imag))

    # taking care of f_b fd_a = fd_a f_b + \delta_a,b.
    # v1e += numpy.eye(v1e.shape[0])*lambdac.trace()/lambdac.shape[0]
    # but here the convention is for a typical molecule.

    if v2e is not None:
        max_imag = max(max_imag, numpy.max(numpy.abs(v2e.imag)))
    if max_imag > 1.e-7:
        print(' maximal imaginary part = {}'.format(max_imag))

    # write one-body part
    twrite_2d_array_real(v1e, 'H1E_{}.INP'.format(imp), shreshold=1.e-7)

    # write two-body part
    if v2e is not None:
        twrite_u_matrix_real(v2e, 'V2E_{}.INP'.format(imp))


def get_dm_cshupdn(imp=1):
    '''get the density matrix in complex spherical harmonics basis
    with faster spin index ordered as up, down.
    '''
    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'r') as f:
        dm = f['/DM'][()].T
    u_sab2sf = get_u_sab2cshupdn(imp=imp)
    u2_sab2sf = block_diag(u_sab2sf, u_sab2sf)
    dm = u2_sab2sf.T.dot(dm).dot(u2_sab2sf.conj())
    if numpy.max(numpy.abs(dm.imag)) < 1.e-7:
        dm = dm.real
    numpy.savetxt('dm.dat', dm)


def h5wrt_dm_sab_rc(imp=1):
    '''transform the density matrix (real) into symmetry adapted basis
    and write in EMBEB_HAMIL_{IMP}.h5 file in complex format.
    '''
    # read in density matrix (real)
    fname = "GDMRG_{}.OUT".format(imp)
    dm = numpy.loadtxt(fname)

    # cast to complex
    dm = numpy.array(dm, dtype=complex)

    # get transformation matrix to symmetry adapted basis
    line = open(fname, "r").readline()
    e_mol = float(line.split()[1])
    u_sab2sf = get_u_sab2cshupdn(imp=imp)
    u2_sab2sf = block_diag(u_sab2sf, u_sab2sf)

    # unitary transformation
    dm = u2_sab2sf.conj().dot(dm).dot(u2_sab2sf.T)
    with h5py.File("EMBED_HAMIL_RES_{}.h5".format(imp), 'w') as f:
        f["/DM"] = dm.T
        f["/emol"] = [e_mol]


def chk_consistent_dm(imp=1):
    '''consistent check for dm with unitary transformations.
    '''
    dm = numpy.loadtxt('dm.dat')
    with h5py.File('EMBED_HAMIL_CSH_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
    v1e = get_whole_h1e(h1e, lambdac, daalpha)
    res1 = numpy.sum(dm*v1e)

    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'r') as f:
        dm = f['/DM'][()].T
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
    v1e = get_whole_h1e(h1e, lambdac, daalpha)
    res2 = numpy.sum(dm*v1e)
    print(res1, res2)



if __name__ == '__main__':
    h1e, lambdac, daalpha, v2e = h5gen_embedh_spin_updn(h5wrt=True)
    wrt_text_rembed(h1e, lambdac, daalpha, v2e)
