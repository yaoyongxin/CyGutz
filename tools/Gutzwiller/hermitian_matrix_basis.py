"""
Given the self-energy structure sigma,
generate the Hermitian matrix basis set.
"""
import numpy as np


def SigmaToMatrixBasis(sigma):
    matrix_basis = []
    sigma = np.asarray(sigma)
    for element in xrange(np.max(sigma), 0, -1):
        spots = np.argwhere(sigma == element)
        # Terminate if reaching end
        num_spots = len(spots)
        if num_spots == 0:
            break
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
