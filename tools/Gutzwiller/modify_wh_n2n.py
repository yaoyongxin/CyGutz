import h5py
import numpy as np
import sys
sys.path.append("/home/ykent/BitBucket/pyscript/Gutzwiller")

f = h5py.File("glog.h5", 'r')
e_ab = f["/Impurity_1/One_body_local"][...]
f.close()

np.set_printoptions(precision = 3, suppress=True, linewidth=400)
print 'Original e_ab = \n', e_ab

# rotate locally
U = np.identity(14, dtype=complex)

w_1, v_1 = np.linalg.eigh(e_ab[4:6, 4:6])
print 'e_ab_1 = \n', w_1

U[4:6, 4:6] = U[6:8, 6:8] = U[8:10, 8:10] = U[10:12, 10:12] = v_1

e_ab = np.dot(np.conj(U.T), np.dot(e_ab, U))

print 'Rotated e_ab = \n', e_ab

from fileio import *
U_All = read_TRANS("WH_N2N.INP")
U_All[0][:,:] = np.dot(U_All[0][:,:], U)

write_TRANS('WH_N2N_NEW.INP', U_All)
