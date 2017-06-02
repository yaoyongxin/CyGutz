from __future__ import print_function

'''
Tools to analyse the local many-body density matrix (multiplet structure).
'''

import os
from builtins import range
import numpy as np
from scipy.sparse import csr_matrix
import h5py
from dataproc import get_csr_matrix, h5write_csr_matrix


def get_rho_histogram(rho, ID_block, N=None, S=None, L=None, J=None,
                      num_ev=0, Rpr_list=None):
    '''
    Get the histogram of the local density matrix with labels.

    This function diagonalizes the reduced local many-body density matrix rho,
    in the order of valence block specified by the ID_block. The label of the
    resultant eigen-space, such as valence N, (averaged) J, character chi of
    the irreducible representation and degeneracy, will be computed. The final
    results are ordered accordig to descending order of the eigen-values of
    the eigen-spaces.

    Parameters
    ----------
    rho : csr_matrix
        reduced local many-body density matrix in sparse matrix format.
    ID_block : 1d array
        Indices to specify the valence blocks. For instance, ID_block will be
        equal to [1, 2, 16, 107] for a f-shell with valence range from 0-2.
        Specifically, the valence block i corresponds to the states from
        ID_block[i] to ID_block[i + 1] - 1.
    N : csr_matrix
        Total valence electron number operator in sparse matrix format.
    J : csr_matrix
        Total angular momentum :math:`J^{2}` operator in sparse matrix format.
    num_ev : integer
        Number of significant eigen-vectors of rho to be calculated.
    Rpr_list : list.
        List of rotation operations in the local Hilbert space of difference
        valence block.

    Returns
    -------
    vals : array
        Eigen-values of the eigen-spaces of rho.
    n_label : array
        Labels of valence of the eigen-spaces of rho.
    j_label : array
        Averaged J values of the eigen-spaces of rho.
    chi_label : array
        Characters of the eigen-spaces of rho as the irreducible
        representation of the rotation group defined by Rpr_list.
    multiplet_degeneracies : array
        Degeneracies of the eigen-spaces of rho.

    '''
    print(" Get rho histogram.")
    if num_ev > rho.shape[0]:
        num_ev = 0
    if num_ev <= 0:
        from dataproc import block_diagonalize
        vals, vecs = block_diagonalize(rho, ID_block)
    else:
        from scipy.sparse.linalg import eigsh
        vals, vecs = eigsh(rho, num_ev, which='LM')
        vecs = vecs.T
        # Poor lanczos, eigen-values may not be properly ordered.
        idx = vals.argsort()[::-1]
        vals = vals[idx]
        vecs = vecs[idx, :]

    from local_operator_factory import get_label_list
    print(" Get labels.")
    n_label = get_label_list(N, vecs)
    s_label = get_label_list(S, vecs)
    l_label = get_label_list(L, vecs)
    j_label = get_label_list(J, vecs)
    if s_label is not None:
        s_label = np.sqrt(s_label + 0.25) - 0.5
    if l_label is not None:
        l_label = np.sqrt(l_label + 0.25) - 0.5
    if j_label is not None:
        j_label = np.sqrt(j_label + 0.25) - 0.5

    from matrix_util import set_eigen_space, shrink_label
    idx = set_eigen_space(vals)
    vals = shrink_label(vals, idx, method="sum")
    n_label = shrink_label(n_label, idx, method="average")
    s_label = shrink_label(s_label, idx, method="average")
    l_label = shrink_label(l_label, idx, method="average")
    j_label = shrink_label(j_label, idx, method="average")
    multiplet_degeneracies = np.array(
        [idx[i + 1] - idx[i] for i in range(len(idx) - 1)])

    if Rpr_list is not None:
        # check commutation
        from AtomSymmetry import check_commute_G
        for ni in range(min(3,len(ID_block)-1)):
            print(' Check commutation for ni = ', ni)
            check_commute_G(rho[ID_block[ni]:ID_block[ni + 1], \
                    ID_block[ni]:ID_block[ni + 1]].todense(), Rpr_list[ni])

        from AtomSymmetry import get_characters_espace, check_sum_chi2_1
        chi_label = []

        for i, n in enumerate(n_label):
            ni = int(n + 0.5 - N[0,0].real)  # Relative due to possible shift.
            chi, repr = get_characters_espace(
                    Rpr_list[ni], vecs[idx[i] : idx[i + 1], \
                    ID_block[ni] : ID_block[ni + 1]])
            chi_label.append(chi)
            check = check_sum_chi2_1(chi)
            if check != "OK":
                print(" Warning: ", check, idx[i], idx[i + 1])
        chi_label = np.array(chi_label)
    else:
        chi_label = None

    idx = vals.argsort()[::-1]
    vals = vals[idx]
    if s_label is not None:
        s_label = s_label[idx]
    if l_label is not None:
        l_label = l_label[idx]
    if j_label is not None:
        j_label = j_label[idx]
    if n_label is not None:
        n_label = n_label[idx]
    if chi_label is not None:
        chi_label = chi_label[idx]
    if multiplet_degeneracies is not None:
        multiplet_degeneracies = multiplet_degeneracies[idx]
    return vals, n_label, s_label, l_label, j_label, chi_label, \
            multiplet_degeneracies


def generate_symmetrized_rho(rho, ID_block, VAL_N_block, Rpr_list):
    '''
    Explicitly symmetrize rho.
    '''
    from scipy.sparse import block_diag
    from matrix_util import sym_csr_matrix
    rho_list = []
    num_block = len(ID_block) - 1
    for i in range(num_block):
        print("   calculating block ", i)
        if ID_block[i + 1] <= ID_block[i]:
            continue
        rho_list.append(sym_csr_matrix(
            rho[ID_block[i]:ID_block[i + 1], ID_block[i]:ID_block[i + 1]], Rpr_list[i]))
    return block_diag(rho_list).tocsr()


def update_csr_matrix(f, path, rho_sym):
    try:
        del f[path + "_SYM.base"], \
            f[path + "_SYM.data"], \
            f[path + "_SYM.indices"], \
            f[path + "_SYM.indptr"], \
            f[path + "_SYM.ncol"], \
            f[path + "_SYM.nrow"]
    except:
        pass

    h5write_csr_matrix(f, path, rho_sym)


def generate_symmetrized_rho_f(imp,rpath):
    '''
    Generate and store the symmetrized rho-like arrays, e.g., rho, phi.
    '''
    f = h5py.File("hilbert_space_rotations.h5", 'r')
    from hs_rotation_representation import get_hs_rotations
    Rpr_list = get_hs_rotations(f, imp)
    f.close()
    f = h5py.File("glog.h5", 'r+')
    # one-based to zero-based
    ID_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1
    VAL_N_block = f["/Impurity_" + str(imp) + "/Sec_VAL"][...]
    rho = get_csr_matrix(f, "/Impurity_" + str(imp) + rpath)
    rho_sym = generate_symmetrized_rho(rho, ID_block, VAL_N_block, Rpr_list)

    update_csr_matrix(f, "/Impurity_" + str(imp) + rpath + "_SYM", rho_sym)
    f.close()


def generate_symmetrized_rho_by_hmexpand_f(imp=1,rpath='/RHO'):
    '''
    Generate and store the symmetrized rho-like arrays, e.g., rho, phi,
    by expanding the array with repect to a symmetry-adapted Hermitian
    matrix basis set stored in glog_sym.h5 file.
    '''
    f = h5py.File("glog.h5", 'r+')
    if "/Impurity_" + str(imp) + rpath + "_SYM.data" not in f:
        rho_fock = get_csr_matrix(f, "/Impurity_" + str(imp) + rpath)
        with h5py.File('glog_sym.h5', 'r') as fp:
            u_fock_sab = get_csr_matrix(fp, "/Impurity_" + str(imp)
                    + "/FockToCurrentBasis")
            nbase = u_fock_sab.shape[0] - rho_fock.shape[0]
            u_fock_sab = u_fock_sab[nbase:, nbase:]

            # To symmetry adapted basis
            rho_sab = u_fock_sab.getH().multiply(rho_fock.multiply(u_fock_sab))
            num_hm_basis = fp["/Impurity_" + str(imp) + "/NUM_HMAT_BASIS"][0]
            for i in range(num_hm_basis):
                hm_basis = get_csr_matrix(fp, "/Impurity_" + str(imp) +
                        "/HMAT_BASIS_" + str(i+1))[nbase:, nbase:]
                coef = np.sum(hm_basis.getH().multiply(rho_sab).diagonal())
                if i == 0:
                    rho_sym_sab = hm_basis.multiply(coef)
                else:
                    rho_sym_sab += hm_basis.multiply(coef)

            # Back to Fock basis
            rho_sym_fock = u_fock_sab.multiply(rho_sym_sab.multiply( \
                    u_fock_sab.getH()))

        update_csr_matrix(f, "/Impurity_" + str(imp) + rpath + "_SYM",
                rho_sym_fock)
    f.close()


def get_local_histogram(imp, f, op_list=["N", "S2", "L2", "J2"], num_ev=0,
        Rpr_list=None, rho_name='/RHO', n_bot = 0, n_top = None):
    '''
    Get the local histogram with labels.
    '''
    from local_operator_factory import get_local_operators
    print(" Get local operators (N, J2 ...)")
    op_dict = get_local_operators(imp, f, op_list)

    N = op_dict["N"]
    if "S2" in op_list:
        S = op_dict["S2"]
    else:
        S = None
    if "L2" in op_list:
        L = op_dict["L2"]
    else:
        L = None
    if "J2" in op_list:
        J = op_dict["J2"]
    else:
        J = None

    rho = get_csr_matrix(f, "/Impurity_" + str(imp) + rho_name)

    n1 = N.shape[0]
    n2 = rho.shape[0]
    N = N[n1 - n2: n1, n1 - n2: n1]
    if S is not None:
        S = S[n1 - n2: n1, n1 - n2: n1]
    if L is not None:
        L = L[n1 - n2: n1, n1 - n2: n1]
    if J is not None:
        J = J[n1 - n2: n1, n1 - n2: n1]
    # one-based to zero-based
    ID_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1
    if n_bot == 0:
        n_bot = int(N[0, 0].real + 0.5) - \
              int(f["/Impurity_" + str(imp) + "/Sec_VAL"][0] + 0.5)
    if n_top == None:
        n_top = len(ID_block) - 1
    print(' n_bot = {}, n_top = {}'.format(n_bot, n_top))

    rho = rho[ID_block[n_bot] : ID_block[n_top],
            ID_block[n_bot] : ID_block[n_top]]
    N = N[ID_block[n_bot] : ID_block[n_top],
            ID_block[n_bot] : ID_block[n_top]]
    if S is not None:
        S = S[ID_block[n_bot] : ID_block[n_top],
                ID_block[n_bot] : ID_block[n_top]]
    if L is not None:
        L = L[ID_block[n_bot] : ID_block[n_top],
                ID_block[n_bot] : ID_block[n_top]]
    if J is not None:
        J = J[ID_block[n_bot] : ID_block[n_top],
                ID_block[n_bot] : ID_block[n_top]]
    if Rpr_list is not None:
        Rpr_list = Rpr_list[n_bot : n_top + 1]
    ID_block = ID_block[n_bot : n_top + 1] - ID_block[n_bot]

    return get_rho_histogram(rho, ID_block, N=N, S=S, L=L, J=J, num_ev=num_ev,
            Rpr_list=Rpr_list)


if __name__ == "__main__":
    generate_symmetrized_rho_by_hmexpand_f(imp=1, rpath="/RHO")
