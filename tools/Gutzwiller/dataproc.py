import h5py
import numpy as np
from scipy.sparse import csr_matrix

def get_csr_matrix(f, path):
    '''
    Read the csr_matrix located at path in the hdf5 file f.
    '''
    nrow = f[path + ".nrow"][0]
    ncol = f[path + ".ncol"][0]
    data = f[path + ".data"][...]
    base = f[path + ".base"][0]
    indices = f[path + ".indices"][...] - base # possible one-based to zero-based
    indptr = f[path + ".indptr"][...] - base # possible one-based to zero-based
    return csr_matrix((data, indices, indptr), shape=(nrow, ncol))

def h5write_csr_matrix(f, path, a):
    '''
    Read the csr_matrix located at path in the hdf5 file f.
    '''
    f[path + ".nrow"] = [a.shape[0]]
    f[path + ".ncol"] = [a.shape[1]]
    f[path + ".data"] = a.data
    f[path + ".base"] = [0]
    f[path + ".indices"] = a.indices
    f[path + ".indptr"] = a.indptr

def get_ed_lowest_vec(h_mat):
    '''
    Compute the lowest eigen-vector of the local many-body Hamiltonian.
    '''
    from scipy.sparse.linalg import eigsh
    vals, vecs = eigsh(h_mat, 3, which='SA')
    return vecs[:,0:1]

def get_ed_eigvec(h_mat):
    '''
    Compute the lowest eigen-vector of the local many-body Hamiltonian.
    '''
    from scipy.linalg import eigh
    vals, vecs = eigh(h_mat.todense())
    return vecs[:,0:1]

def get_ed_dm(evec):
    '''
    Compute the density matrix based on the given state vector.
    Dimension: evec[:,0:1]
    '''
    return np.dot(evec, np.conj(evec.T))

def trace_single_state(evec, op):
    '''
    Calculate Tr(|state><state| op).
    '''
    v_op = np.conj(evec.T) * op
    return np.sum(evec[:,0]*v_op[0,:])

def block_diagonalize(A, ID_block):
    '''
    Diagonalize the block diagonal matrix A characterized by chi.
    '''
    val_list = []; vec_list = []
    from scipy.linalg import eigh
    for i, istart in enumerate(ID_block[:-1]):
      iend = ID_block[i + 1]
      if iend <= istart:
        continue
      tmp = A[:istart, istart:iend]
      if (tmp.nnz > 0):
        assert np.amax(np.abs(tmp.data)) < 1.e-8, " Error in block_diagonalize!"
      tmp = A[istart:iend, :istart]
      if (tmp.nnz > 0):
        assert np.amax(np.abs(tmp.data)) < 1.e-8, " Error in block_diagonalize!"
      tmp = A[iend:, istart:iend]
      if (tmp.nnz > 0):
        assert np.amax(np.abs(tmp.data)) < 1.e-8, " Error in block_diagonalize!"
      tmp = A[istart:iend, iend:]
      if (tmp.nnz > 0):
        assert np.amax(np.abs(tmp.data)) < 1.e-8, " Error in block_diagonalize!"
      vals, vecs = eigh(A[istart:iend, istart:iend].todense())
      val_list.extend(vals.tolist())
      vec_list.append(vecs.T)
    from scipy.linalg import block_diag
    vec_list = block_diag(*vec_list)
    val_list = np.array(val_list)
    idx = val_list.argsort()[::-1]
    val_list = val_list[idx]
    vec_list = vec_list[idx, :]
    return val_list, vec_list

if __name__ == "__main__":
    '''
    Test.
    '''
