# analyze the embeding hamiltonian
import h5py, numpy
from pyglib.math.matrix_util import \
        get_utrans_orbfast_supdn_to_spin_fast_supdn


def h5gen_embedh_spin_updn(imp=1):

    # embedding hamiltonian parameters
    # upper case path indicates fortran convention.
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][()].T
        lambdac = f['/LAMBDA'][()].T
        h1e     = f['/H1E'][()].T
        v2e     = f['/V2E'][()].T

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

    # tensordot, switch orders?
    v2e_ = numpy.tensordot(numpy.tensordot(numpy.tensordot( \
            numpy.tensordot(v2e, \
            u_sf2sab, axes=([0], [1])), \
            u_sab2sf, axes=([0], [0])), \
            u_sf2sab, axes=([0], [1])), \
            u_sab2sf, axes=([0], [0]))

#    v2e_ = numpy.einsum('ijkl,pi,jq,rk,ls->pqrs', v2e, u_sf2sab, \
#            u_sab2sf, u_sf2sab, u_sab2sf)

    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'a') as f:
        f['/D'][()] = daalpha_.T
        f['/LAMBDA'][()] = lambdac_.T
        f['/H1E'][()] = h1e_.T
        f['/V2E'][()] = v2e_.T


if __name__ == '__main__':
    h5gen_embedh_spin_updn()
