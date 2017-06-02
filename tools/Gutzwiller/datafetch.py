'''
Functions to fetch various data from the Wien2K + Gutzwiller calculations.
'''
import os
import h5py
import numpy as np

def get_data_scf(path, label):
  '''
  Fetch data from case.scf.
  '''
  f = os.popen('grep ' + label + " " + path)
  line = f.readlines()
  assert len(line) > 0, "Error: No " + label + " found in " + path
  f.close()
  line = line[-1]
  line = line.split()
  return float(line[-1])

def get_data_list_scf(path_list, label):
  '''
  Fetch data list from a list of case.scf(m).
  '''
  res_list = []
  for path in path_list:
    res_list.append(get_data_scf(path, label))
  return res_list

def get_R_list(path_list, id):
  '''
  Get list of R of impurity id.
  '''
  R_list = []
  for path in path_list:
    f = h5py.File(path + "/glog.h5", 'r')
    R = f["/Impurity_" + str(id) + "/GA_R"][...]
    R_list.append(R.tolist())
    f.close()
  return R_list

def get_Nab_list(path_list, id):
  '''
  Get list of one-body density matrix Nab of impurity id.
  '''
  from one_body_density_matrix import get_one_body_density_matrix
  Nab_list = []
  for path in path_list:
    Nab_list.append(get_one_body_density_matrix(id, path))
  return Nab_list

def get_Z_evalues(R_list, n1, n2):
  '''
  Get eigen-values of Z = R^{\dagger} R.
  '''
  w_list = []
  from numpy import linalg as la
  for R in R_list:
    Z = np.dot(np.conj(np.array(R)[n1:n2, n1:n2].T), np.array(R)[n1:n2, n1:n2])
    w = la.eigvalsh(Z)
    w_list.append(w.tolist())
  return w_list


