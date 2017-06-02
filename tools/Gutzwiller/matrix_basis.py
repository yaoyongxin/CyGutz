"""
Given the matrix structure, generate the (Hermitian) matrix basis set.
"""
import numpy as np


def SigmaToMatrixBasis(sigma):
    matrix_basis = []
    sigma = np.asarray(sigma)
    for element in xrange(np.max(sigma), 0, -1):
        spots = np.argwhere(sigma == element)
        num_spots = len(spots)
        if num_spots == 0:
            continue
        spots = spots.T
        # Skip if located at lower trigonal block
        if spots[0][0] > spots[1][0]:
            continue
        if spots[0][0] == spots[1][0]:
            value = 1 / np.sqrt(float(num_spots))
            matrix = np.zeros_like(sigma, dtype=np.complex128)
            matrix[spots[0], spots[1]] = value
            matrix_basis.append(matrix)
        else:
            # non-zero element
            value = 1 / np.sqrt(float(num_spots * 2))
            matrix = np.zeros_like(sigma, dtype=np.complex128)
            matrix[spots[0], spots[1]] = value
            matrix[spots[1], spots[0]] = value
            matrix_basis.append(matrix)
            value = value * 1.j
            matrix = np.zeros_like(sigma, dtype=np.complex128)
            matrix[spots[0], spots[1]] = value
            matrix[spots[1], spots[0]] = -value
            matrix_basis.append(matrix)
    return matrix_basis


def ListSigmaToMatrixBasis(sigma_list):
    matrix_basis_list = []
    for sigma in sigma_list:
        matrix_basis_list.append(SigmaToMatrixBasis(sigma))
    return matrix_basis_list


def MatrixStructToBasis(M_struct):
    matrix_basis = []
    M_struct = np.asarray(M_struct)
    for element in xrange(np.max(M_struct), 0, -1):
        spots = np.argwhere(M_struct == element)
        num_spots = len(spots)
        if num_spots == 0:
            continue
        spots = spots.T
        value = 1 / np.sqrt(float(num_spots))
        matrix = np.zeros_like(M_struct, dtype=np.complex128)
        matrix[spots[0], spots[1]] = value
        matrix_basis.append(matrix)
    return matrix_basis


def ListMatrixStructToBasis(M_struct_list):
    matrix_basis_list = []
    for M_struct in M_struct_list:
        matrix_basis_list.append(MatrixStructToBasis(M_struct))
    return matrix_basis_list


def hermitian_matrix_basis(n, istart=0):
    '''
    Generate hermitian matrix basis set of dimention n.
    '''
    matrix_basis = []
    val = 1./np.sqrt(2.)
    for i in range(istart, n):
        for j in range(i, n):
            if i == j:
                matrix = np.zeros((n, n), dtype=np.complex128)
                matrix[i, j] = 1.0
                matrix_basis.append(matrix)
            else:
                matrix = np.zeros((n, n), dtype=np.complex128)
                matrix[i, j] = matrix[j, i] = val
                matrix_basis.append(matrix)
                matrix = np.zeros((n, n), dtype=np.complex128)
                matrix[i, j] = 1.j*val
                matrix[j, i] = -1.j*val
                matrix_basis.append(matrix)
    return matrix_basis
