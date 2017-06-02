'''
Tools to analyse the local many-body density matrix (multiplets).
'''

import numpy as np
from scipy.sparse import csr_matrix
import os
import h5py

def get_rho_histogram(rho, ID_block, N = None, J = None,
    num_ev = 0, Rpr_list = None):
  '''
  Get the histogram of the local density matrix with labels.
  '''
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
  n_label = get_label_list(N, vecs)
  j_label = get_label_list(J, vecs)
  if j_label is not None:
    j_label = np.sqrt(j_label + 0.25) -0.5

  from matrix_util import set_eigen_space, shrink_label
  idx = set_eigen_space(vals)
  vals    = shrink_label(vals   ,idx,method="sum")
  n_label = shrink_label(n_label,idx,method="average")
  j_label = shrink_label(j_label,idx,method="average")
  multiplet_degeneracies = np.array([idx[i+1]-idx[i] for i in range(len(idx)-1)])

  if Rpr_list is not None:
    # check commutation
    from AtomSymmetry import check_commute_G
    for ni in range(3):
      print 'Check commutation for ni = ',ni
      check_commute_G(rho[ID_block[ni]:ID_block[ni+1],ID_block[ni]:ID_block[ni+1]].todense(),Rpr_list[ni])

    from AtomSymmetry import get_characters_espace, check_sum_chi2_1
    chi_label = []
    for i,n in enumerate(n_label):
      ni = int(n+0.5)
      chi, repr = get_characters_espace(Rpr_list[ni],vecs[idx[i]:idx[i+1],ID_block[ni]:ID_block[ni+1]])
      chi_label.append(chi)
      check = check_sum_chi2_1(chi)
      if check != "OK":
        print "Warning: ",check, idx[i],idx[i+1]
    chi_label = np.array(chi_label)
  else:
    chi_label = None

  idx = vals.argsort()[::-1]
  vals = vals[idx]
  if j_label is not None:
    j_label = j_label[idx]
  if n_label is not None:
    n_label = n_label[idx]
  if chi_label is not None:
    chi_label = chi_label[idx]
  if multiplet_degeneracies is not None:
    multiplet_degeneracies = multiplet_degeneracies[idx]
  return vals, n_label, j_label, chi_label, multiplet_degeneracies

def generate_symmetrized_rho(rho, ID_block, VAL_N_block, Rpr_list):
  '''
  Explicitly symmetrize \rho.
  '''
  from scipy.sparse import block_diag
  from matrix_util import sym_csr_matrix
  rho_list = []
  num_block = len(ID_block) - 1
  for i in range(num_block):
    if ID_block[i+1] <= ID_block[i]: continue
    rho_list.append(sym_csr_matrix(rho[ID_block[i]:ID_block[i+1], ID_block[i]:ID_block[i+1]], Rpr_list[i]))
  return block_diag(rho_list).tocsr()

def generate_symmetrized_rho_f(imp):
  '''
  Generate and store the symmetrized \rho.
  '''
  f = h5py.File("hilbert_space_rotations.h5", 'r')
  from hs_rotation_representation import get_hs_rotations
  Rpr_list = get_hs_rotations(f, imp)
  f.close()
  f = h5py.File("glog.h5", 'r+')
  ID_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1 # one-based to zero-based
  VAL_N_block = f["/Impurity_" + str(imp) + "/Sec_VAL"][...]
  from dataproc import get_csr_matrix, h5write_csr_matrix
  rho = get_csr_matrix(f, "/Impurity_" + str(imp) + "/RHO")
  rho_sym = generate_symmetrized_rho(rho, ID_block, VAL_N_block, Rpr_list)
  try:
    del f["/Impurity_" + str(imp) + "/RHO_SYM.data"], \
        f["/Impurity_" + str(imp) + "/RHO_SYM.indices"], \
        f["/Impurity_" + str(imp) + "/RHO_SYM.indptr"], \
        f["/Impurity_" + str(imp) + "/RHO_SYM.nrow"], \
        f["/Impurity_" + str(imp) + "/RHO_SYM.ncol"]
  except:
    pass
  h5write_csr_matrix(f, "/Impurity_" + str(imp) + "/RHO_SYM", rho_sym)
  f.close()

def get_local_histogram(imp, f, num_ev = 0, Rpr_list = None, rho_name = '/RHO'):
  '''
  Get the local histogram with labels.
  '''
  from local_operator_factory import get_local_operators
  op_dict = get_local_operators(imp, f, ["N", "J2"])

  N = op_dict["N"]
  J = op_dict["J2"]

  from dataproc import get_csr_matrix
  rho = get_csr_matrix(f, "/Impurity_" + str(imp) + rho_name)

  n1 = J.shape[0]; n2 = rho.shape[0]
  N = N[n1 - n2 : n1, n1 - n2 : n1]
  J = J[n1 - n2 : n1, n1 - n2 : n1]
  ID_block = f["/Impurity_" + str(imp) + "/Sec_ID"][...] - 1 # one-based to zero-based
  return get_rho_histogram(rho, ID_block, N=N, J=J, num_ev=num_ev, Rpr_list=Rpr_list)


if __name__=="__main__":
  '''
  Test.
  '''
  np.set_printoptions(precision=8)
  f = h5py.File("glog.h5", 'r+')

  if "/Impurity_1/RHO_SYM.data" not in f:
    generate_symmetrized_rho_f(1)
    quit()

  Rpr_list = None
  try:
    f_p = h5py.File("hilbert_space_rotations.h5", 'r')
    from hs_rotation_representation import get_hs_rotations
    Rpr_list = get_hs_rotations(f_p, 1)
    f_p.close()
  except:
    pass
  vals, n_label, j_label, chi_label, multiplet_degeneracies \
      = get_local_histogram(1, f, num_ev=400, Rpr_list=Rpr_list, rho_name='/RHO_SYM')
  f.close()
  print " weights:"
  print vals
  print " sum-weight:"
  print np.sum(vals)
  print " n_label:"
  print n_label
  print " j_label:"
  print j_label
  print " chi_label:"
  print chi_label

  with open("multiplets.py", 'w') as f:
    print >> f, "weights = ", vals.tolist()
    print >> f, "n_label = ", n_label.tolist()
    print >> f, "j_label = ", j_label.tolist()
    print >> f, "deg_label = ", multiplet_degeneracies.tolist()
    if chi_label is not None:
      print >> f, "chi_label = ", chi_label.tolist()
  import matplotlib.pyplot as plt
  x = range(len(j_label))
  f, axarr = plt.subplots(4, sharex=True)
  axarr[0].plot(x, vals, 'o-')
  axarr[0].set_ylabel('weight')
  axarr[1].plot(x, n_label, 'o-')
  axarr[1].set_ylabel('N')
  axarr[2].plot(x, multiplet_degeneracies, 'o-')
  axarr[2].set_ylabel('mult. deg.')
  axarr[3].plot(x, j_label, 'o-')
  axarr[3].set_ylabel('J')
  axarr[3].set_xlabel('multiplets')
  plt.show()
