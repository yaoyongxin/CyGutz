'''
Multiplets analysis for LDA+G+SOC calculations.
'''

import numpy as np
from scipy.sparse import csr_matrix
import os
import h5py
from multiplets_analysis_lib import get_local_histogram

def calc_save_atomic_states(imp=1, num_ev=100):
  '''
  Calculate and save the eigen-values of the local many-body density matrix
  and the labels.
  '''
  # Default hdf5-format log file of the LDA+G job.
  f = h5py.File("glog.h5", 'r+')
  # Get eigen-values and labels of the local atomic states.
  vals, n_labels, j_labels, chi_labels, multiplet_degeneracies \
      = get_local_histogram(imp, f, num_ev)
  f.close()
  # Store results to metadata.
  f = h5py.File("multiplets.h5", 'a')
  base = "/impurity_"+str(imp)
  try:
    del f[base]
  except:
    pass
  f[base+"/sum_weights"] = np.sum(vals)
  f[base+"/weights"] = vals
  f[base+"/n_labels"] = n_labels
  f[base+"/j_labels"] = j_labels
  f[base+"/deg_labels"] = multiplet_degeneracies
  f.close()

def plot_atomic_states(imp=1, num_label = 5):
  '''
  Plot the atomic state probabilities using the formula \rho = EXp(-F)
  '''
  # Read in metadata
  f = h5py.File("multiplets.h5", 'r')
  base = "/impurity_"+str(imp)
  vals = f[base+"/weights"][...]
  multiplet_degeneracies = f[base+"/deg_labels"][...]
  n_labels = f[base+"/n_labels"][...]
  j_labels = f[base+"/j_labels"][...]

  # Generate data corresponding to \rho = exp(-F)
  f_wt = -np.log(vals/multiplet_degeneracies)
  f_wt = f_wt - f_wt[0]

  # Labels of the atomic states
  labels = ['N=%1d,J=%3.1f'%(int(n+0.1),j) for n,j in zip(n_labels[:num_label],j_labels[:num_label])]

  # Scatter plot with annotation.
  import matplotlib.pyplot as plt
  plt.scatter(f_wt,vals)
  for label,x,y in zip(labels, f_wt[:num_label],vals[:num_label]):
    plt.annotate(
        label,
        xy = (x,y), xytext = (0,20),
        textcoords = 'offset points', ha = 'center', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
  plt.xlim(-0.1,12)
  plt.ylim(0,1)
  plt.xlabel("$f_{n}$")
  plt.ylabel("$e^{-f_{n}}d_{n}/\sum_{n}{e^{-f_{n}}d_{n}}$")
  plt.title("Eigen-values of the local many-body density matrix")
  plt.show()

if __name__=="__main__":
  #calc_save_atomic_states(imp=1, num_ev=400)
  plot_atomic_states(imp=1, num_label = 5)
