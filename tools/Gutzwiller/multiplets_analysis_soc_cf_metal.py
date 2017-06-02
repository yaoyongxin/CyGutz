'''
Multiplets analysis for the metallic phase of LDA+G calculations
with spin-orbit coupling (SOC) and crystal field effect (CF).
'''

import numpy as np
from scipy.sparse import csr_matrix
import os
import h5py
from hs_rotation_representation import h5get_hs_rotation_representation_N
from multiplets_analysis_lib import get_local_histogram


# Impurity site, one-based
imp = 1

if not os.path.isfile("hilbert_space_rotations.h5"):
  h5get_hs_rotation_representation_N(imp)

np.set_printoptions(precision=2)
f = h5py.File("glog.h5", 'r+')
f_p = h5py.File("hilbert_space_rotations.h5", 'r')
from hs_rotation_representation import get_hs_rotations
Rpr_list = get_hs_rotations(f_p, 1)
f_p.close()

# main analysis routine
vals, n_label, j_label, chi_label, multiplet_degeneracies \
    = get_local_histogram(1, f, num_ev=400, Rpr_list=Rpr_list)
f.close()

# brief info
print " weights: \n", vals
print " sum-weight: \n", np.sum(vals)
print " n_label: \n", n_label
print " j_label: \n", j_label
print " chi_label: \n", chi_label
print " deg_label: \n", multiplet_degeneracies

# file to store the results.
f = h5py.File("multiplets.h5", 'w')
f["/weights"] = vals
f["/n_label"] = n_label
f["/j_label"] = j_label
for i,chi in enumerate(chi_label):
  f["/chi_label"+str(i)] = chi
f["/deg_label"] = multiplet_degeneracies
f.close()

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
