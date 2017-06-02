from __future__ import print_function

import h5py
import itertools as it
import numpy as np
import sys
import scipy.sparse as sp
import scipy.sparse as sps
from scipy.linalg import sqrtm
from scipy.sparse.linalg import expm

sys.path.append("/home/ykent/WIEN_GUTZ/bin/tools/Gutzwiller")
from matrix_util import get_func_matrix_eigh, get_matrix_trace


__times = 0


def get_s_xyz(c_list):
    '''
    Get spin s vector.
    '''
    from local_operator_factory import get_Sz_Sp_from_c_CH
    s_z, s_p = get_Sz_Sp_from_c_CH(c_list)
    s_x = (s_p.getH() + s_p) / 2.0
    s_y = (-s_p.getH() + s_p) / 2.0j
    return (s_x, s_y, s_z)


def orthonormalize_matrices(A_list):
    '''
    Orthonormalize list of Hermitian matrices.

    Parameters
    ----------
    A_list: list of Hermitian matrices

    Returns
    -------
    Ao_list: list of orthonormalized Hermitian matrices
    '''
    Ao_list = []
    for i, A in enumerate(A_list):
        tr = get_matrix_trace(A*A.getH())
        A /= np.sqrt(tr)
        for Ao in Ao_list:
            tr = get_matrix_trace(A*Ao.getH())
            A -= tr*Ao
        if i>0:
          tr = get_matrix_trace(A*A.getH())
          assert tr.real > 1.e-6, "Error: almost degenerate matrix added!"
          A /= np.sqrt(tr)
        Ao_list.append(A)
    return Ao_list


def check_orthonormality(A_list):
    for i, A in enumerate(A_list):
        for j, B in enumerate(A_list):
            tr = get_matrix_trace(A*B.getH())
            if i==j:
                assert np.abs(tr - 1.0) < 1.e-6, "Not normalized."
            else:
                assert np.abs(tr) < 1.e-6, "Not orthogonal."


def get_rho_1bk(A_list, lambda_list):
    for i, _A, _lambda in it.izip(it.count(), A_list, lambda_list):
        if i==0:
            A = -_A*_lambda
        else:
            A -= _A*_lambda
    res = get_func_matrix_eigh(A, func='exp')
    #res = expm(A)
    return res


def get_trace_bk(rho_bk):
    for i, rho in enumerate(rho_bk):
        if i==0:
            tr = get_matrix_trace(rho)
        else:
            tr += get_matrix_trace(rho)
    return tr


def get_rho_exp_f(lambda_list, A_bk_list):
    rhop_bk = [get_rho_1bk(A_list, lambda_list) for A_list in A_bk_list]
    # Calculate trace(exp(-F))
    tr = get_trace_bk(rhop_bk)
    # Scale by 1/trace(exp(-F))
    for i, rhop in enumerate(rhop_bk):
        rhop_bk[i] = rhop/tr
    return rhop_bk


def get_trace_dist_piece(rho, rhop, tau=1.e-30):
    '''
    Calculate the trace distance of one valence block of rho.

    Parameters
    ----------
    rho : (N,N) ndarray or csc_matrix
    rhop : (N,N) ndarray or csc_matrix
    tau : float
        Smearing factor used to approximate
        sqrt(x**2 + tau**2) ~ |x|

    Returns
    -------
    tr_dist : float
        Trace distance.
    '''
    rho_diff = rho - rhop

    # square-root approximation
    rho_diff2 = rho_diff * rho_diff.getH() + tau**2*np.identity(rho.shape[0])
    #abs_rho_diff = get_func_matrix_eigh(rho_diff2, func='sqrt')
    abs_rho_diff = sqrtm(rho_diff2)

    # Or absolute value difference.
    # abs_rho_diff = get_func_matrix_eigh(rho_diff, func='abs')

    tr_dist = get_matrix_trace(abs_rho_diff).real
    return tr_dist


def get_trace_dist_all(rho_bk, rhop_bk):
    for i, rho, rhop in it.izip(it.count(), rho_bk, rhop_bk):
        if i==0:
            tr_dist = get_trace_dist_piece(rho, rhop)
        else:
            tr_dist += get_trace_dist_piece(rho, rhop)
    return tr_dist


def get_rho_dist(lambda_list, A_bk_list, rho_bk):
    rhop_bk = get_rho_exp_f(lambda_list, A_bk_list)
    tr_dist = get_trace_dist_all(rho_bk, rhop_bk)
    global __times
    __times += 1
    print(" Iter = %10d, trace dist = %.2e"%(__times, tr_dist))
    return tr_dist


def rho_projection_coef(A_bk_list, rho_bk):
    coef_list = np.zeros(len(A_bk_list[0]), dtype=np.float)
    for A_list, rho in it.izip(A_bk_list, rho_bk):
        log_rho = get_func_matrix_eigh(rho, func='log')
        for i, A in enumerate(A_list):
            coef_list[i] -= get_matrix_trace(A*log_rho).real
    return coef_list


if __name__ == "__main__":
    '''
    Perform fitting.
    '''
    f = h5py.File("glog.h5", 'r')

    from dataproc import get_csr_matrix
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    ID_block = f["/Impurity_1/SEC_ID"][...] - 1

    from dataread import get_c_operators
    c_list = get_c_operators(1, f)
    f.close()

    print(' Max imag of rho =', np.max(np.abs(rho.data.imag)))

    # n operators
    n_op = []
    for c_dn, c_up in it.izip(
            it.islice(c_list, 0, None, 2),
            it.islice(c_list, 1, None, 2)):
        n_op.append(
                c_dn.getH() * c_dn +
                c_up.getH() * c_up )

    # S operators
    #s_i = get_s_xyz(c_list[:2])
    #s_0 = get_s_xyz(c_list[2:4])

    # pmee terms
    A_list = []

    # identity
    one = sp.identity(rho.shape[0]).tocsc()
    A_list.append(one.copy())

    # Quadratic terms
    A_list.extend(n_op)

    for i, c_0dn, c_0up in it.izip(it.count(),
            it.islice(c_list, 0, None, 2),
            it.islice(c_list, 1, None, 2)):
        for c_1dn, c_1up in it.izip(
                it.islice(c_list, i*2 + 2, None, 2),
                it.islice(c_list, i*2 + 3, None, 2)):
            A_list.append(
                    c_0dn.getH() * c_1dn +
                    c_1dn.getH() * c_0dn +
                    c_0up.getH() * c_1up +
                    c_1up.getH() * c_0up)

    # Quartic terms
#    for _n1 in n_op:
#        A_list.append((_n1 - one)**2)
#    for _n1, _n2 in it.izip(it.islice(n_op, 0, None),it.islice(n_op, 1, None)):
#        A_list.append((_n1 - one)*(_n2 - one))

    # clean some data
    del c_list, n_op

    # Orthonomalize matrix basis-set
    A_list = orthonormalize_matrices(A_list)
    print(" Orthonormalizing matrix basis done.")

    # Remove identity, which can be any scale factor
    del A_list[0]

    # Check matrix basis
    check_orthonormality(A_list)

    # Change to valence block-diagonal form
    A_bk_list = []
    rho_bk = []
    j_iter = it.islice(ID_block, 1, None)
    for i, j in it.izip(ID_block, j_iter):
        A_bk_list.append([a[i:j, i:j] for a in A_list])
        rho_bk.append(rho[i:j, i:j])
    print(" Trace(Rho) =", get_trace_bk(rho_bk))

    # Clean some data
    del A_list, rho

    coef_list = rho_projection_coef(A_bk_list, rho_bk)
    print(coef_list)
    print(" Initial trace distance =",
            get_rho_dist(coef_list, A_bk_list, rho_bk))

#    Fitting.
    from scipy.optimize import minimize
    res = minimize(get_rho_dist, coef_list, args=(A_bk_list, rho_bk),
            method='CG', jac=False)
    print(res.x, res.fun)
    print(" Final trace distance =",
            get_rho_dist(res.x, A_bk_list, rho_bk))
