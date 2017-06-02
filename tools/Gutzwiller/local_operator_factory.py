'''
Get lcoal operators given the {c} operators.
'''

import h5py
import itertools
import numpy as np
from scipy.sparse import csr_matrix


def linear_trans_operators(op_list, U):
    '''
    Linear transformation (U) of the operators op_list.
    '''
    opt_list = []
    for U1 in U:
        sum = csr_matrix(op_list[0].shape, dtype=np.complex128)
        for i, x in enumerate(U1):
            if np.abs(x) < 1.e-12:
                continue
            sum += x * op_list[i]
        opt_list.append(sum)
    return opt_list


def get_Sz_Sp_from_c_CH(c_list):
    '''
    Get S_z, S^\dagger operators given the c operators in complex Harmonics
    basis set (spin-up block + spin-down block).

    '''
    S_z = csr_matrix(c_list[0].shape, dtype=np.complex128)
    S_p = csr_matrix(S_z)
    num_orbitals = len(c_list) / 2
    c_up_iter = itertools.islice(c_list, 0, num_orbitals)
    c_dn_iter = itertools.islice(c_list, num_orbitals, None)
    for c_up, c_dn in itertools.izip(c_up_iter, c_dn_iter):
        S_z += 0.5 * (- c_dn.getH() * c_dn + c_up.getH() * c_up)
        S_p += c_up.getH() * c_dn
    return S_z, S_p


def get_Lz_Lp_from_c_CH(c_list):
    '''
    Get L_z, L^\dagger operators given the c operators in complex Harmonics
    basis set (spin-up block + spin-down block).

    '''
    L_z = csr_matrix(c_list[0].shape, dtype=np.complex128)
    L_p = csr_matrix(L_z)
    num_orbitals = len(c_list) / 2
    L = num_orbitals / 2
    for i in range(num_orbitals):
        M = i - L
        c_dn = c_list[i]
        c_up = c_list[i + num_orbitals]
        L_z += M * (c_dn.getH() * c_dn + c_up.getH() * c_up)
        if i == num_orbitals - 1:
            continue
        coef = np.sqrt(np.float64((L - M) * (L + M + 1)))
        L_p += coef * (c_list[i + 1].getH() * c_dn +
                       c_list[i + 1 + num_orbitals].getH() * c_up)
    return L_z, L_p


def get_N_from_c(c_list, op_list):
    '''
    Get N operators given the c operators in complex Harmonics basis set
    (spin-up block + spin-down block).

    '''
    res = {}
    if "N" in op_list:
        N_op = csr_matrix(c_list[0].shape, dtype=np.complex128)
        for c in c_list:
            N_op += c.getH() * c
        res["N"] = N_op
    if "Nab" in op_list:
        Nab_list = []
        for i in range(len(c_list)):
            Nab_row = []
            for j in range(i + 1):
                Nab_row.append(c_list[i].getH() * c_list[j])
            Nab_list.append(Nab_row)
        res["Nab"] = Nab_list
    return res


def get_AM_op_from_c_CH(c_list, op_list):
    '''
    Get angular momentum operators given the c operators in complex Harmonics
    basis set (spin-up block + spin-down block).

    '''
    if op_list is None:
        return None

    S_z, S_p = get_Sz_Sp_from_c_CH(c_list)
    L_z, L_p = get_Lz_Lp_from_c_CH(c_list)
    J_z, J_p = S_z + L_z, S_p + L_p
    res = {}
    if "Sx" in op_list:
        res["Sx"] = (S_p.getH() + S_p) / 2.0
    if "Sy" in op_list:
        res["Sy"] = (S_p - S_p.getH()) / 2.0j
    if "Sz" in op_list:
        res["Sz"] = S_z
    if "Lx" in op_list:
        res["Lx"] = (L_p.getH() + L_p) / 2.0
    if "Ly" in op_list:
        res["Ly"] = (L_p - L_p.getH()) / 2.0j
    if "Lz" in op_list:
        res["Lz"] = L_z
    if "Jx" in op_list:
        res["Jx"] = (J_p.getH() + J_p) / 2.0
    if "Jy" in op_list:
        res["Jy"] = (J_p - J_p.getH()) / 2.0j
    if "Jz" in op_list:
        res["Jz"] = J_z
    if "S2" in op_list:
        res["S2"] = S_p.getH() * S_p + S_z * S_z + S_z
    if "L2" in op_list:
        res["L2"] = L_p.getH() * L_p + L_z * L_z + L_z
    if "J2" in op_list:
        res["J2"] = J_p.getH() * J_p + J_z * J_z + J_z
    for key in res.keys():
        res[key] = res[key].tocsr()
    return res


def get_local_operators(imp, f, op_list):
    '''
    Get the local operators, like N, Nab, Sx, Sy, Sz, Lx, Ly, Lz, S^2, L^2
    and J^2 operators for impurity imp.
    '''
    from dataproc import get_csr_matrix

    # Read in one-body rotation matrix.
    U_CH_cur = f[
        "/Impurity_" + str(imp) + "/ComplexHarmonicsToCurrentBasis"][...].T

    # Read annihilation operators in Fock basis.
    from dataread import get_c_operators
    c_list = get_c_operators(imp, f)

    try:
        U_fock = get_csr_matrix(
            f, "/Impurity_" + str(imp) + "/FockToCurrentBasis")
    # Rotate to the same basis as rho.
        for i in range(len(c_list)):
            c_list[i] = np.conj(U_fock.T) * c_list[i] * U_fock
    except Exception:
        pass

    res = {}
    res.update(get_N_from_c(c_list, op_list))
    c_CH_list = linear_trans_operators(c_list, U_CH_cur)
    res.update(get_AM_op_from_c_CH(c_CH_list, op_list))
    return res


def get_label_list(J, vecs):
    '''
    Get list of expection values of J for vecs.
    '''
    if J is None:
        label = None
    else:
        label = []
        for vec in vecs:
            res = np.vdot(vec, J.dot(vec))
            assert np.abs(res.imag) < 1.e-8, " error in get_label_list!"
            label.append(res.real)
        label = np.array(label)
    return label
