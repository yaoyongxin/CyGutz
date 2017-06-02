'''
Common functions to read information from glog.h5.
'''
import h5py

def get_embed_parameters(imp, f):
  '''
  Get embedding parameters: D, Lc, t, U.
  '''
  D = f["/Impurity_" +str(imp) + "/GA_D"][...]
  Lc = f["/Impurity_" +str(imp) + "/GA_Lc"][...]
  t = f["/Impurity_" +str(imp) + "/One_body_local"][...]
  U = f["/Impurity_" +str(imp) + "/Two_body_local"][...]
  R = f["/Impurity_" +str(imp) + "/GA_R"][...]
  return t, U, D, Lc, R

def get_c_operators(imp, f):
  '''
  Get all the annihilation operators c.
  '''
  from dataproc import get_csr_matrix
  c = []
  i = 0
  while True:
    try:
      c.append(get_csr_matrix(f, "/Impurity_" + str(imp) + "/annihi.op._" + str(i + 1))) # one-based.
      i += 1
    except Exception:
      break
  return c


