import h5py
import numpy as np
from scipy.linalg import eigh


f = h5py.File("glog.h5", 'r')
N = f["/Impurity_1/GA_NC_PHY"][...].T
R = f["/Impurity_1/GA_R"][...].T
Z = np.dot(np.conj(R.T),R)

w, v = eigh(Z[:2, 0:2])
print '\Gamma-7-Z:\n',w
w, v = eigh(Z[4:6, 4:6])
print '\Gamma-8-Z:\n',w
print '\Gamma-6-Z:\n',Z[12,12]

w, v = eigh(N[:2, 0:2])
print '\Gamma-7-N:\n',w
w, v = eigh(N[4:6, 4:6])
print '\Gamma-8-N:\n',w
print '\Gamma-6-N:\n',N[12,12]

