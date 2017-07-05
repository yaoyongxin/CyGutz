# -*- coding: utf-8 -*-

'''
Get the matrix representation of angular momentum operator
in one-particle basis.
'''

import numpy as np
from scipy.linalg import block_diag
from pyglib.math.matrix_util import trans_orbital_fast_to_spin_fast


def get_J_vector(l_list, basis):
    '''
    Get L(J) vector from a l list with CH (JJ) basis.
    '''
    J = [np.empty([0, 0], dtype=complex), np.empty(
        [0, 0], dtype=complex), np.empty([0, 0], dtype=complex)]
    for l in l_list:
        if basis == 'CH':
            Jz = get_matrix_Lz_CH(l)
            Jp = get_matrix_Lp_CH(l)
            Jn, Jx, Jy = get_other_op(Jz, Jp)
        else:
            Jz = get_matrix_Jz_JJ(l)
            Jp = get_matrix_Jp_JJ(l)
            Jn, Jx, Jy = get_other_op(Jz, Jp)
        J[0] = block_diag(J[0], Jx)
        J[1] = block_diag(J[1], Jy)
        J[2] = block_diag(J[2], Jz)
    return J


def get_L_vector_CH(l):
    '''
    Get (L_x, L_y, L_z) with CH basis.
    '''
    Lz = get_matrix_Lz_CH(l)
    Lp = get_matrix_Lp_CH(l)
    Ln, Lx, Ly = get_other_op(Lz, Lp)
    return Lx, Ly, Lz


def get_J_vector_JJ(l):
    '''
    Get (J_x, J_y, J_z) with JJ basis.
    '''
    Jz = get_matrix_Jz_JJ(l)
    Jp = get_matrix_Jp_JJ(l)
    Jn, Jx, Jy = get_other_op(Jz, Jp)
    return Jx, Jy, Jz


def get_Lp_coef(l, m):
    '''
    L^\dagger |l, m> = Lp_coef |l, m+1>.
    '''
    return np.sqrt(l * (l + 1) - m * (m + 1))


def get_diag_Lz_CH(l):
    '''
    Get diagonal elements of Lz for Complex spherical Harmomnics basis.
    '''
    Lz = []
    for j in range(int(2 * l + 1)):
        Lz.append(-l + j)
    return Lz


def get_matrix_Lz_CH(l):
    '''
    Get matrix_{z}.
    '''
    return np.diag(get_diag_Lz_CH(l))


def get_matrix_Lp_CH(l, m_rising=True):
    '''
    Get matrix L_{\dagger}.
    '''
    lm = int(2 * l + 1.5)
    Lp = np.zeros([lm, lm])
    if m_rising:
        for i in range(int(2 * l + 0.5)):
            Lp[i + 1, i] = get_Lp_coef(l, -l + i)
    else:
        for i in range(int(2 * l + 0.5)):
            Lp[i, i + 1] = get_Lp_coef(l, l - i - 1)
    return Lp


def get_matrix_Sz_CH_orbital_fast(l_list):
    '''
    Get matrix Sz for Complex spherical Harmomnics basis
    with spin-down + up block.
    '''
    num_lm = np.sum(2 * np.array(l_list) + 1)
    # spin-up + spin_dn. Wien2k convention.
    Sz = np.diag([ 0.5 for i in range(num_lm)] + [-0.5 for i in range(num_lm)])
    return Sz


def get_matrix_Lzp_CH_orbital_fast(l_list):
    '''
    Get matrix Lz, Lp for Complex spherical Harmomnics basis
    with spin-up + dn block.
    '''
    for i, _l in enumerate(l_list):
        _Jz = get_matrix_Lz_CH(_l)
        _Jp = get_matrix_Lp_CH(_l)
        if i == 0:
            Jz = _Jz.copy()
            Jp = _Jp.copy()
        else:
            Jz = block_diag(Jz, _Jz)
            Jp = block_diag(Jp, _Jp)
    # Adding the spin-dn block.
    Jz = block_diag(Jz, Jz)
    Jp = block_diag(Jp, Jp)
    return Jz, Jp


def get_matrix_Sp_CH_spin_fast(l_list):
    '''
    Get matrix Sp for Complex spherical Harmomnics basis with spin-fast-index.
    '''
    num_lm = np.sum(2 * np.array(l_list) + 1)
    Sp_sub = get_matrix_Lp_CH(0.5, m_rising=False)
    Sp_list = [Sp_sub for i in range(num_lm)]
    from scipy.linalg import block_diag
    return block_diag(*Sp_list)


def get_S_vector_CH_orbital_fast(l_list):
    '''
    Get S-vector in CH basis with fast orbital index.
    '''
    Sz = get_matrix_Sz_CH_orbital_fast(l_list)
    Sp = get_matrix_Sp_CH_spin_fast(l_list)
    Sp = trans_orbital_fast_to_spin_fast(Sp, lback=True)
    _, Sx, Sy = get_other_op(Sz, Sp)
    return [Sx, Sy, Sz]


def get_L_vector_CH_orbital_fast(l_list):
    '''
    Get L-vector in CH basis with fast orbital index.
    '''
    Lz, Lp = get_matrix_Lzp_CH_orbital_fast(l_list)
    _, Lx, Ly = get_other_op(Lz, Lp)
    return [Lx, Ly, Lz]


def get_trans_JJ_to_CH_orbital_fast(l_list):
    '''
    Get the unitary transformation from JJ to CH_orbital_fast basis.
    '''
    from pyglib.math.matrix_util import trans_JJ_to_CH_sup_sdn
    U_list = []
    for l in l_list:
        # watch out the convention.
        U = trans_JJ_to_CH_sup_sdn(l)
        U_list.append(U)
    return block_diag(*U_list)


def get_S_vector_JJ(l_list):
    '''
    Get S-vector in JJ basis.
    '''
    S_vec = get_S_vector_CH_orbital_fast(l_list)
    U = get_trans_JJ_to_CH_orbital_fast(l_list)
    for i, _S in enumerate(S_vec):
        S_vec[i] = U.dot(_S.dot(U.conj().T))
    return S_vec


def get_L_vector_JJ(l_list):
    '''
    Get L-vector in JJ basis.
    '''
    L_vec = get_L_vector_CH_orbital_fast(l_list)
    U = get_trans_JJ_to_CH_orbital_fast(l_list)
    for i, _L in enumerate(L_vec):
        L_vec[i] = np.dot(U, np.dot(_L, np.conj(U.T)))
    return L_vec


def get_matrix_Jz_JJ(l):
    '''
    Get matrix J_z with relativisitc Harmonics basis
    for l = 0, the basis is assumed to be (s_down and s_up)
    '''
    if l == 0:
        return np.diag([-0.5, 0.5])
    else:
        j1 = l - 0.5
        J1z = np.diag(get_diag_Lz_CH(j1))
        j2 = l + 0.5
        J2z = np.diag(get_diag_Lz_CH(j2))
        return block_diag(J1z, J2z)


def get_matrix_Jp_JJ(l):
    '''
    Get matrix J_p with relativisitc Harmonics basis
    for l = 0, the basis is assumed to be (s_down and s_up)
    '''
    if l == 0:
        alpha = get_Lp_coef(0.5, -0.5)
        return(np.array([[0., 0.], [alpha, 0.]]))
    else:
        j1 = l - 0.5
        J1p = get_matrix_Lp_CH(j1)
        j2 = l + 0.5
        J2p = get_matrix_Lp_CH(j2)
        return block_diag(J1p, J2p)


def get_other_op(Lz, Lp):
    '''
    Get other related operators given L_z and L^\dagger.
    '''
    Ln = Lp.conj().T
    Lx = (Lp + Ln) / 2.
    Ly = (Lp - Ln) / (2.j)
    return Ln, Lx, Ly


def get_l_list_from_string(code):
    '''
    Gt l list from letters, e.g., "spd" -> [0, 1, 2].
    '''
    l_list = []
    for c in code:
        l = get_l_from_char(c)
        if l >= 0:
            l_list.append(l)
    return l_list


def get_l_from_char(c):
    return {'s': 0, 'p': 1, 'd': 2, 'f': 3}.get(c, -1)


def get_S_vector(l_list, iso):
    '''
    Dispatcher.
    '''
    if iso > 1:
        return get_S_vector_JJ(l_list)
    else:
        return get_S_vector_CH_orbital_fast(l_list)


def get_L_vector(l_list, iso):
    '''
    Dispatcher.
    '''
    if iso > 1:
        return get_L_vector_JJ(l_list)
    else:
        return get_L_vector_CH_orbital_fast(l_list)


def get_complex_to_real_sph_harm(l):
    '''get the unitary transformation from compex spherical harmonics
    (Condon–Shortley phase convention) to real harmonics
    (https://en.wikipedia.org/wiki/Spherical_harmonics,
    consistent with VASP).
    '''
    mdim = 2*l + 1
    c2r = np.zeros((mdim, mdim), dtype=np.complex)
    for m in range(mdim):
        m_ = m - l
        if m_ > 0:
            c2r[ m_+l, m] = (-1)**m_/np.sqrt(2.)
            c2r[-m_+l, m] = 1./np.sqrt(2.)
        elif m_ == 0:
            c2r[l, l] = 1.
        else:
            c2r[ m_+l, m] = 1.j/np.sqrt(2.)
            c2r[-m_+l, m] = -1.j*(-1)**m_/np.sqrt(2.)
    return c2r


if __name__ == "__main__":
    Lx, Ly, Lz = get_L_vector_CH(0)
    print Lx, Ly, Lz
