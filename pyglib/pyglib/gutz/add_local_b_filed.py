from __future__ import print_function
from builtins import range
import h5py
import numpy as np
import itertools as it

with h5py.File("GPARAM.h5", 'r') as f:
    svec_list = []
    lvec_list = []
    db2sab_list = []
    for i in range(f["/num_imp"][0]):
        svec = []
        lvec = []
        svec.append(f["/IMPURITY_1/SX"][...].T)
        svec.append(f["/IMPURITY_1/SY"][...].T)
        svec.append(f["/IMPURITY_1/SZ"][...].T)
        lvec.append(f["/IMPURITY_1/LX"][...].T)
        lvec.append(f["/IMPURITY_1/LY"][...].T)
        lvec.append(f["/IMPURITY_1/LZ"][...].T)
        svec_list.append(svec)
        lvec_list.append(lvec)

# apply magnetic field
bvec = [0, 0, 0.2]
vb_list = []
for svec, lvec in it.izip(svec_list, lvec_list):
    vb = bvec[0]*(svec[0]) # + lvec[0])
    for b, s, l in it.izip(bvec[1:], svec[1:], lvec[1:]):
        vb += b*(s) # + l)
    vb_list.append(vb)

with h5py.File("WH_RL_INIT.h5", 'w') as f:
    f["/vshift"] = np.asarray([1], dtype=np.int32)
    for i, vb in enumerate(vb_list):
        f["/IMPURITY_" + str(i+1) + "/LAMBDA"] = vb.T
