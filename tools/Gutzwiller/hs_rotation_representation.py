import h5py
import numpy as np
from scipy.sparse import csr_matrix
from dataproc import h5write_csr_matrix


def h5get_hs_rotation_representation_N(imp):
    '''
    Calculate the Hilbert space rotation representation of each valence block
    of impurity imp in the hilbert_space_rotations.h5 file.
    '''

    f = h5py.File("init_ga_info.h5", 'r')
    # 0-based in python
    Lie_Jeven_params = f["/impurity_" + str(imp - 1) + "/Lie_Jeven"][...]
    Lie_Jodd_params = f["/impurity_" + str(imp - 1) + "/Lie_Jodd"][...]
    f.close()

    from local_operator_factory import get_local_operators

    f = h5py.File("glog.h5", 'r')
    op_dict = get_local_operators(imp, f, ["N", "Jx", "Jy", "Jz"])
    from dataproc import get_csr_matrix
    rho = get_csr_matrix(f, "/Impurity_" + str(imp) + "/RHO")
    # one-based to zero-based
    ID_N_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1
    VAL_N_block = f["/Impurity_" + str(imp) + "/Sec_VAL"][...]
    f.close()

    f_p = h5py.File("hilbert_space_rotations.h5", 'w')
    n1 = op_dict["Jx"].shape[0]
    n2 = rho.shape[0]
    base = n1 - n2

    from AtomSymmetry import get_representation
    num_N_block = len(ID_N_block) - 1
    dim_group_N = []
    for i in range(num_N_block):
        if ID_N_block[i + 1] - ID_N_block[i] <= 0:
            dim_group_N.append(0)
            continue
        Jx = op_dict["Jx"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        Jy = op_dict["Jy"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        Jz = op_dict["Jz"][base + ID_N_block[i] : base + ID_N_block[i + 1],
                           base + ID_N_block[i] : base + ID_N_block[i + 1]]
        if abs(np.mod(VAL_N_block[i], 2)) > 1.e-4:
            Rpr = get_representation([Jx, Jy, Jz], Lie_Jodd_params)
        else:
            Rpr = get_representation([Jx, Jy, Jz], Lie_Jeven_params)

        dim_group_N.append(len(Rpr))
        for j in range(dim_group_N[i]):
            Rpr[j] = csr_matrix(Rpr[j])
            h5write_csr_matrix(f_p, "Impurity_" + str(imp) + "/NBlock=" +
                    str(i) + "/rotation_" + str(j), Rpr[j])
    f_p["Impurity_" + str(imp) + "/dim_group_N"] = dim_group_N
    f_p.close()


def get_U_chi_sigma(imp):
    '''
    Get the U_list, chi_list and sigma_list in the many-body space.
    '''
    from dataproc import get_csr_matrix
    from AtomSymmetry import get_atom_U_sigma_chi
    f_p = h5py.File("hilbert_space_rotations.h5", 'r')
    list_U = []
    list_chi = []
    list_sigma = []
    for i, dim_g in enumerate(f_p["Impurity_" + str(imp) + "/dim_group_N"]):
        Rpr = []
        for j in range(dim_g):
            Rpr.append(get_csr_matrix(f_p, "Impurity_" + str(imp) +
                "/NBlock=" + str(i) + "/rotation_" + str(j)))
        U, sigma, chi = get_atom_U_sigma_chi(Rpr)
        list_U.append(U)
        list_chi.append(chi)
        list_sigma.append(sigma)
    f_p.close()
    return list_U, list_chi, list_sigma


def get_hs_rotations(f, imp):
    '''
    Get rotation representations in many-body space.
    '''
    Rpr_list = []
    dim_group_N = f["Impurity_" + str(imp) + "/dim_group_N"][...]
    from dataproc import get_csr_matrix
    for i, dim_group in enumerate(dim_group_N):
        Rpr_list.append([])
        for j in range(dim_group):
            Rpr_list[i].append(get_csr_matrix(f, "/Impurity_" + str(imp) + \
                    "/NBlock=" + str(i) + "/rotation_" + str(j)))
    return Rpr_list


if __name__ == "__main__":
    h5get_hs_rotation_representation_N(1)
