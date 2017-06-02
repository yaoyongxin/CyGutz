from __future__ import print_function
'''
Coherent density of states of UO2 in Mott insulating phase
with crystal field and spin-orbit interaction.
'''

import h5py
import numpy as np
from pyglib.basic import units
from pyglib.estructure import dos

f = h5py.File("glog.h5", 'r')
e_fermi = f["/E_FERMI"][0]
e_skn = (f["/BND_EK"][...] - e_fermi) * units.Ryd_eV
f_skn = f["/BND_FERWE"][...]
w_k = f["/KPT_WT"][...]
bnd_ne = f["/BND_NE"][...]
R = f["/Impurity_1/GA_R"][...].T
# ecdfdffdelta = f["/Impurity_1/ecdfdffdelta"][...]
num_orb = R.shape[0]
psi_sksn = f["/BND_VK"][:, :, :, :, : num_orb]
f.close()

# shift the flat bands to correct position.
e_skn[abs(e_skn-(30-e_fermi)*units.Ryd_eV)<0.01] = 0

psi_skn_w_ab = dos.get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne)

window = (-5., 15.)
_dos = dos.DOS(w_k, e_skn,  width=0.05, window=window, npts=101)
energies = _dos.get_energies()

# total quasi-particle
dos_t = _dos.get_dos_t()

# total quasi-particle f-component
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., :, :])
dos_f_qp = _dos.get_dos_component(psi_sksn_f)

# total f coherent part
psi_skn_w_ab = np.einsum('...ij,ik,jl', psi_skn_w_ab, np.conj(R), R)
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., :, :])
dos_f_ch = _dos.get_dos_component(psi_sksn_f)
#
#dos_t_nch = np.zeros((psi_skn_w_ab.shape[0], _dos.npts), dtype=np.complex)
#
#for i in range(1):
#
#    print(i)
#
#    dos_t_nch += _dos.get_dos_noncoherent(psi_sksn, ecdfdffdelta, f_skn,
#            bnd_ne, i)


#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

f, axarr = plt.subplots(1, 4, sharey=True)

axarr[0].plot(energies, dos_t[0].real, '-')
axarr[1].plot(energies, dos_f_qp[0].real, '-')
axarr[2].plot(energies, dos_f_ch[0].real, '-')
#axarr[3].plot(energies, dos_t_nch[0].real, '-')
axarr[0].set_ylabel('DOS')
axarr[1].set_xlabel('E (eV)')

#plt.savefig('dos_1.pdf')

plt.show()
