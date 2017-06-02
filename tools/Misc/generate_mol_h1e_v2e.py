import numpy as np
import h5py


with h5py.File('glog.h5', 'r') as f:
    _D = f['/Impurity_1/GA_D'][...].T
    _Lc = f['/Impurity_1/GA_Lc'][...].T
    _t = f['/Impurity_1/One_body_local'][...].T
    _u = f['/Impurity_1/Two_body_local'][...].T

n = _D.shape[0]
h1e = np.zeros((2*n, 2*n), dtype=_D.dtype)
h1e[:n, :n] = _t; h1e[n:, n:] = -_Lc
h1e[:n, n:] = _D.T; h1e[n:, :n] = _D.conj()
v2e = np.zeros((2*n, 2*n, 2*n, 2*n), dtype=_D.dtype)
v2e[:n, :n, :n, :n] = _u

print np.trace(_Lc)

with h5py.File('impurity_H.h5', 'w') as f:
    f['/h1e'] = h1e
    f['/v2e'] = v2e
