'''
Tools to calculate expectation values of the local observables.
'''

import numpy as np
from scipy.sparse import csr_matrix
import os
import h5py

def get_exp_val_local_observables(imp, f, ob_list):
  '''
  Calculate expectation values of the local observables.
  '''
  from local_operator_factory import get_local_operators
  op_dict = get_local_operators(imp, f, ob_list)

  from dataproc import get_csr_matrix
  rho = get_csr_matrix(f, "/Impurity_" + str(imp) + "/RHO")
  res = {}
  for key, op in op_dict.items():
    if key == "Nab":
      idx_start = op[0][0].shape[0]-rho.shape[0]
      num_op = len(op)
      Nab = np.zeros((num_op, num_op), dtype = np.complex128)
      for i in range(num_op):
        for j in range(i+1):
          Nab[i,j] = np.sum((rho*op[i][j][idx_start:,idx_start:]).diagonal())
          if i != j:
            Nab[j,i] = np.conj(Nab[i,j])
      res[key] = Nab
    else:
      idx_start = op.shape[0]-rho.shape[0]
      res[key] = np.sum((rho*op[idx_start:,idx_start:]).diagonal())
  return res


if __name__=="__main__":
  '''
  Test.
  '''
  np.set_printoptions(precision=2)
  f = h5py.File("glog.h5", 'r')
  ans = get_exp_val_local_observables(1, f, ["Sx", "Sy", "Sz", "N"])
  print ans
