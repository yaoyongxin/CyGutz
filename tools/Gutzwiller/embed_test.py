import h5py
import numpy as np
from scipy.sparse import csr_matrix
from dataproc import get_csr_matrix, get_ed_lowest_vec, trace_single_state, get_ed_eigvec

def get_one_body_term(t, c):
    '''
    Get local one-body terms.
    '''
    res = csr_matrix(c[0].shape, dtype = np.complex64)
    indices = np.array(np.where(np.abs(t) > 1.e-10)).T
    for index in indices:
      i = index[0]; j = index[1]
      res += t[i, j]*c[i].getH()*c[j]
    return res

def get_two_body_term(U, c):
    '''
    Get local two-body terms.
    '''
    res = csr_matrix(c[0].shape, dtype = np.complex64)
    indices = np.array(np.where(np.abs(U) > 1.e-10)).T
    for index in indices:
      i = index[0]; j = index[1]; k = index[2]; l = index[3]
      if i ==j or k == l:
        continue
      res += U[i, j, k, l]/2.0*c[i].getH()*c[j].getH()*c[l]*c[k]
    return res

def get_D_term(D, c):
    '''
    Get D (cross) terms.
    '''
    res = csr_matrix(c[0].shape, dtype = np.complex64)
    indices = np.array(np.where(np.abs(D) > 1.e-10)).T
    N1 = D.shape[0]
    for index in indices:
      i = index[0]; j = index[1]
      I = i + N1
      res += D[i, j]*c[j].getH()*c[I] + np.conj(D[i, j])*c[I].getH()*c[j]
    return res

def get_Lc_term(Lc, c):
    '''
    Get Lc (auxillary) terms.
    '''
    res = csr_matrix(c[0].shape, dtype = np.complex64)
    indices = np.array(np.where(np.abs(Lc) > 1.e-10)).T
    N1 = Lc.shape[0]
    for index in indices:
      i = index[0]; j = index[1]
      I = i + N1; J = j + N1
      res += Lc[i, j]*c[J]*c[I].getH()
    return res

def get_var_nssp(c, gs):
    '''
    Get < c^{\dagger} c>_G.
    '''
    nssp = []
    N1 = len(c)/2
    for i in range(N1):
      I = i + N1
      nssp_row = []
      for j in range(N1):
        J = j + N1
        nssp_row.append(trace_single_state(gs, c[J]*c[I].getH()))
      nssp.append(nssp_row)
    return np.array(nssp)

def get_phy_nssp(c, gs):
    '''
    Get < c^{\dagger} c>_G.
    '''
    nssp = []
    N1 = len(c)/2
    for i in range(N1):
      nssp_row = []
      for j in range(N1):
        nssp_row.append(trace_single_state(gs, c[i].getH()*c[j]))
      nssp.append(nssp_row)
    return np.array(nssp)

def get_R(c, gs):
    '''
    Get the R-matrix.
    '''
    R = []
    N1 = len(c)/2
    for i in range(N1):
      R_row = []
      for j in range(N1):
        J = j + N1
        R_row.append(trace_single_state(gs, c[i].getH()*c[J]))
      R.append(R_row)
    return np.array(R)

def test_embed(param_file, op_file):
    '''
    Test the embedding formalism.
    '''
    from dataread import get_embed_parameters
    np.set_printoptions(precision=4)
    f = h5py.File(param_file, 'r')
    t, U, D, Lc, R = get_embed_parameters(1, f_p)
    f.close()
    t = t.T; D = D.T; Lc = Lc.T; R = R.T  # fortran to c convention
    f = h5py.File(op_file, 'r')
    c = get_c_operators(1, f)
    f.close()
    assert (len(t[0])*2 == len(c)), "len(t[0]) = " + str(len(t[0])) + " != " + str(len(c)) + "!"
    h_emb = get_one_body_term(t, c)
    h_emb += get_two_body_term(U, c)
    h_emb += get_D_term(D, c)
    h_emb += get_Lc_term(Lc, c)

    h_emb = h_emb[93:163, 93:163]
#    h_emb = h_emb[5:11, 5:11]
    state = get_ed_lowest_vec(h_emb)
    gs = np.zeros((c[0].shape[0], 1), dtype=np.complex64)
    gs[93:163,:] = state
    nssp = get_var_nssp(c, gs)
    print 'nssp_var'
    print nssp
    nssp = get_phy_nssp(c, gs)
    print 'nssp_phy'
    print nssp
    R = get_R(c, gs)
    print 'R'
    print R


if __name__ == "__main__":
    '''
    Test.
    '''
    test_embed("glog_2.h5", "glog_4.h5")
