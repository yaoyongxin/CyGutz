import os, h5py
import numpy as np
from pyglib.model import semicir as bethlatt
import subprocess
from builtins import zip


# set range of Hubbard U.
u_list = np.arange(0.0, 5.1, 0.2)

# get *CyGutz* command. Choose option '-r -1' to skip the claculation of
# Gutzwiller renormalized electron density.
root = os.environ['WIEN_GUTZ_ROOT2']
cmd = [root+'/CyGutz', '-r', '-1']

# set chemical potential to 0.
mu=0

# set Hubbard U=0
u = 0.

# generate input files with updated with u=0, mu=0
bethlatt.gutz_model_setup(u=u, nmesh=5000, mu=mu)

# total energy list
e_list = []

# quasi-particle weight z list
z_list = []

# double occupancy d list
d_list = []

# loop over the list of U.
for u in u_list:

    print(' working on u = {}'.format(u))

    # modify the local Coulomb matrix
    with h5py.File('GPARAM.h5', 'a') as f:

        # note the transposition, which is transformation
        # from Fortran convention to c-convention.
        v2e = f['/IMPURITY_1/V2E'][()].T
        vdc2_list = f['VDC2_LIST'][()].T

        # now update the Coulom matrix
        v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
        f['/IMPURITY_1/V2E'][()] = v2e.T

        # update vdc2, which keeps the particle-hole symmetry of the model
        # for finite U.
        vdc2_list[0, 0:2] = -u/2 + mu
        f['VDC2_LIST'][()] = vdc2_list.T

    # perform the *CyGutz* calculation.
    subprocess.call(cmd)

    # get total energy
    with h5py.File('GLOG.h5', 'r') as f:
        nks = f['/IMPURITY_1/NKS_SYM'][()]
        n = nks[0, 0] + nks[1, 1]
        e = f['/etot_model'][0] + u/2*n
        e_list.append(e.real)

    # get Z = R^\dagger R
    with h5py.File('WH_RL_OUT.h5', 'r') as f:
        r = f['/IMPURITY_1/R'][0,0]

        # simple here since it is just a scalar
        z = r*r.conj()
        z_list.append(z.real)

    # to get double occupancy (of impurity 1), <n_up n_dn>_G,
    # we run analysis code *exe_spci_analysis*
    subprocess.call([root+'/exe_spci_analysis', '1'])

    # double occupancy is simply the local many-body density matrix element
    # in the valence=2 block.
    with h5py.File('EMBED_HAMIL_ANALYSIS_1.h5', 'r') as f:
        if '/valence_block_2/RHO' in f:
            d = f['/valence_block_2/RHO'][0, 0]
        else:
            d = 0.
        d_list.append(d.real)

with open('res_u.dat', 'w') as f:
    for u, e, z, d in zip(u_list, e_list, z_list, d_list):
        f.write('{} {} {} {}\n'.format(u, e, z, d))
