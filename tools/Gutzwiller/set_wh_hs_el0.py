'''
Manually set WH_HS_E.INP to have additional symmetrization for the one-body
part of the local oniste Hamiltonian.
'''
import numpy as np

def generate_wh_hs_e(sigma_list, file_name="WH_HS_E.INP"):
    '''
    Given the sigma_list in relativistic Harmonics basis,
    transform it to the current basis.
    '''
    import matrix_basis as mb
    matrix_basis_list = mb.ListSigmaToMatrixBasis(sigma_list)
    import fileio as fio
    u_list = fio.read_TRANS('WH_N2N.INP')

    # basis trnasformation
    for i, u in enumerate(u_list):
        udagger = np.conj(u.T)
        for j, matrix in enumerate(matrix_basis_list[i]):
            matrix_basis_list[i][j] = np.dot(udagger, np.dot(matrix, u))

    fio.write_Hs(file_name, matrix_basis_list)


if __name__=="__main__":
    sigma_list = [np.diag([1,1,1,1,1,1,2,2,2,2,2,2,2,2])]
    generate_wh_hs_e(sigma_list)
