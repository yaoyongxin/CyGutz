from __future__ import print_function
from builtins import range

import numpy as np
from scipy.linalg import sqrtm,inv
from itertools import product
import h5py
from dataproc import get_csr_matrix
from multiplets_analysis_lib import generate_symmetrized_rho_f


def get_tr_phid_cd_phi_f(phi, c_list):
    '''
    calculate Tr(\phi^\dagger c^\dagger_A phi f_a).
    phi and c in sparse matrix format.
    '''
    ecdf = []
    for c in c_list:
        shift = c.shape[0] - phi.shape[0]
        phid_cd_phi = phi.getH() * c[shift:,shift:].getH() * phi
        row = []
        for f in c_list:
            row.append(np.sum((phid_cd_phi * f[shift:,shift:]).diagonal()))
        ecdf.append(row)
    return np.asarray(ecdf)


def get_tr_phid_cd_phi_fdff(phi, c_list):
    '''
    calculate Tr(\phi^\dagger c^\dagger_A phi f^\dagger_a f_b f_c).
    phi and c in sparse matrix format.
    '''
    ecdfdff = []
    for c in c_list:
        shift = c.shape[0] - phi.shape[0]
        phid_cd_phi = phi.getH() * c[shift:,shift:].getH() * phi
        row1 = []
        for f1 in c_list:
            row2 = []
            for f2 in c_list:
                fdf = f1.getH() * f2
                row3 = []
                for f3 in c_list:
                    fdff = (fdf * f3)[shift:, shift:]
                    row3.append(np.sum((phid_cd_phi * fdff).diagonal()))
                row2.append(row3)
            row1.append(row2)
        ecdfdff.append(row1)
    return np.asarray(ecdfdff)


def get_ecdfdelta(ecdf, delta):
    '''
    r_Aa = ecdf_Aa' * (delta(1-delta))^{-1/2}_{a'a}.
    '''
    dd = delta.dot(np.identity(delta.shape[0]) - delta)
    dd = inv(sqrtm(dd))
    return ecdf.dot(dd)


def get_ecdfdffdelta(ecdfdff, ecdf, delta):
    '''
    r_Aabc = (ecdfdff_Aa'b'c' + ecdf_Ab' * delta_a'c' - ecdf_Ac' * delta_a'b')
            (delta(1-delta))^(-1/2)_aa' * (delta(1-delta))^(-1/2)_b'b *
            (delta(1-delta))^(-1/2)_c'c
    '''
    m_range = range(delta.shape[0])
    for j, a, b, c in product(m_range,m_range,m_range,m_range):
        ecdfdff[j,a,b,c] += ecdf[j,b] * delta[a,c] -\
                ecdf[j,c] * delta[a,b]
    dd = delta.dot(np.identity(delta.shape[0]) - delta)
    dd = inv(sqrtm(dd))
    return np.einsum('jdef,ad,eb,fc -> jabc', ecdfdff, dd, dd, dd)


def calc_save_ecdfdffdelta(imp, mode='overwrite'):
    f = h5py.File("glog.h5",'a')
    if 'overwrite' not in mode:
        if '/Impurity_' + str(imp) + '/ecdfdffdelta' in f:
            return
    else:
        if '/Impurity_' + str(imp) + '/ecdfdffdelta' in f:
            del f['/Impurity_' + str(imp) + '/ecdfdffdelta']

    r = f['/Impurity_' + str(imp) + '/GA_R'][...].T
    delta = f["/Impurity_" + str(imp) + "/GA_NKS"][...].T
    c_list = []
    for i in range(r.shape[0]):
        c_list.append(get_csr_matrix(f, "/Impurity_" + str(imp) + \
                "/annihi.op._" + str(i+1)))
    if "/Impurity_" + str(imp) + "/phi_SYM.data" not in f:
        f.close()
        print(" Generate symmetrized \phi.")
        generate_symmetrized_rho_f(imp, '/phi')
        f = h5py.File("glog.h5",'a')

    phi = get_csr_matrix(f, "/Impurity_" + str(imp) + "/phi_SYM")

    # normalize
    norm = np.sum((phi.getH()*phi).diagonal())
    phi = phi.multiply(1/np.sqrt(norm))

    ecdf = get_tr_phid_cd_phi_f(phi, c_list)
    ecdfdelta = get_ecdfdelta(ecdf, delta)
    assert np.allclose(r, ecdfdelta.T, atol=1.e-5), ' r_matrix error!'

    ecdfdff = get_tr_phid_cd_phi_fdff(phi, c_list)
    ecdfdffdelta = get_ecdfdffdelta(ecdfdff, ecdf, delta)

    f['/Impurity_' + str(imp) + '/ecdfdffdelta'] = ecdfdffdelta

    f.close()


if __name__=="__main__":
    calc_save_ecdfdffdelta(1)
