import h5py
import numpy as np
from matrix_util import matrix_basis_expand


def get_BS_matrix(B_vector, S_vector):
    '''
    Get B_vec . S_vec.
    '''
    V = np.zeros_like(S_vector[0])
    for i, B in enumerate(B_vector):
        V += B * S_vector[i]
    return V


def get_B_potential_matrix(B_list):
    '''
    Generate the local potential matrix due to the local magnetic field.
    '''
    f = h5py.File("init_ga_info.h5", 'r')
    matrix_basis_list = f["/matrix_basis_list"]
    S_vector_list = f["/S_vector_list"]
    v_list = []
    for i, S_vec in enumerate(S_vector_list):
        v = get_BS_matrix(B_list[i], S_vec)
        vp = matrix_basis_expand(v, matrix_basis_list[i])
        if np.max(np.abs(vp - v)) > 1.e-6:
            np.set_printoptions(precision=4, suppress=True, linewidth=150)
            print 'v_original_real:'
            print v.real
            print 'v_original_imag:'
            print v.imag
            print 'v_expanded_real:'
            print vp.real
            print 'V_expanded_imag:'
            print vp.imag
            quit()
        v_list.append(vp)
    f.close()
    return v_list


def write_B_potential_matrix(v_list):
    '''
    Write B_potential_matrix file WH_Bvec_Svec.INP.
    '''
    from fileio import write_slice
    dim_list = [v.shape[0] for v in v_list]
    dim_max = np.max(dim_list)
    with open("WH_Bvec_Svec.INP", 'w') as f:
        for v in v_list:
            mat = np.zeros((dim_max, dim_max), dtype=np.complex)
            dim = v.shape[0]
            mat[:dim, :dim] = v
            write_slice(f, v)
    f.close()


if __name__ == "__main__":
    # Create a local B field list
    B_list = [[0, 0, 0.1]]
    v_list = get_B_potential_matrix(B_list)
    write_B_potential_matrix(v_list)
