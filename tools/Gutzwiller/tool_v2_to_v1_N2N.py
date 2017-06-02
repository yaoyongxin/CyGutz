import h5py
import numpy as np

U_list = []
with h5py.File('GPARAM.h5', 'r') as f:
    for imp in range(f['/num_imp'][0]):
        _U = f['/IMPURITY_{}/DB_TO_SAB'.format(imp+1)][...].T
        U = np.zeros_like(_U)

        # for the case without soc
        n1 = U.shape[0]/2
        U[:n1,:n1] = _U[:n1,::2]
        U[n1:,n1:] = _U[n1:,1::2]
        U_list.append(U)

import fileio as fio
fio.write_TRANS("WH_N2N.INP", U_list)
