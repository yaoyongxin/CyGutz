'''
Tools to analyse one-body density matrix.
'''

import numpy as np
from local_observables import get_exp_val_local_observables
import os
import h5py

def get_one_body_density_matrix(imp, path=''):
  f = h5py.File(path + "glog.h5", 'r+')
  try:
    Nab = f["/Impurity_" + str(imp) + "/GA_NC_PHY"][...]
  except:
    Nab = get_exp_val_local_observables(imp, f, ["Nab"])["Nab"]
    f["/Impurity_" + str(imp) + "/GA_NC_PHY"] = Nab
  f.close()
  return Nab


if __name__=="__main__":
  np.set_printoptions(precision=4, linewidth=150, suppress=True)
  Nab = get_one_body_density_matrix(1)
  print Nab.real
  print
  print Nab.imag
