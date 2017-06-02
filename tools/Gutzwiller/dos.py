'''
Coherent density of states of UO2 in Mott insulating phase
with crystal field and spin-orbit interaction.
'''

import h5py
import numpy as np
from pyglib.basic import units
from pyglib.estructure import dos
from pyglib.symm import unitary

f = h5py.File("glog.h5", 'r')
e_fermi = f["/E_FERMI"][0]
e_skn = (f["/BND_EK"][...] - e_fermi) * units.Ryd_eV
w_k = f["/KPT_WT"][...]
bnd_ne = f["/BND_NE"][...]
R = f["/Impurity_1/GA_R"][...].T
jj_to_calc = f["/Impurity_1/PreferedBasisToCurrentBasis"][...].T
num_orb = R.shape[0]
psi_sksn = f["/BND_VK"][:, :, :, :, : num_orb]
f.close()

jj_to_cubic = unitary.jj_to_cubic_relativistic_harmonics('f')
calc_to_cubic = np.einsum('ij,ik', np.conj(jj_to_calc), jj_to_cubic)

psi_skn_w_ab = dos.get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne)

print(' Unitary transforming dos_ab_1')
psi_skn_w_ab = np.einsum('...ij,ik,jl', psi_skn_w_ab, np.conj(R), R)

print(' Unitary transforming dos_ab_2')
psi_skn_w_ab = np.einsum('...ij,ik,jl', psi_skn_w_ab,
        np.conj(calc_to_cubic), calc_to_cubic)

window = (-5., 15.)
_dos = dos.DOS(w_k, e_skn,  width=0.05, window=window, npts=1001)
energies = _dos.get_energies()
dos_t = _dos.get_dos()
# dos_f_G7_5/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 0:2, 0:2])
dos_f_G7_5 = _dos.get_dos(psi_sksn_f)
# dos_f_G81_5/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 2:4, 2:4])
dos_f_G81_5 = _dos.get_dos(psi_sksn_f)
# dos_f_G82_5/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 4:6, 4:6])
dos_f_G82_5 = _dos.get_dos(psi_sksn_f)
# dos_f_G6_7/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 6:8, 6:8])
dos_f_G6_7 = _dos.get_dos(psi_sksn_f)
# dos_f_G7_7/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 8:10, 8:10])
dos_f_G7_7 = _dos.get_dos(psi_sksn_f)
# dos_f_G81_7/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 10:12, 10:12])
dos_f_G81_7 = _dos.get_dos(psi_sksn_f)
# dos_f_G82_7/2
psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[..., 12:14, 12:14])
dos_f_G82_7 = _dos.get_dos(psi_sksn_f)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f, axarr = plt.subplots(3, 3, sharex=True, sharey=True)
axarr[0, 0].plot(energies, dos_t[0]/14, '-')
axarr[0, 1].plot(energies, dos_f_G7_5[0]/2, '-')
axarr[0, 2].plot(energies, dos_f_G81_5[0]/2, '-')
axarr[1, 0].plot(energies, dos_f_G82_5[0]/2, '-')
axarr[1, 1].plot(energies, dos_f_G6_7[0]/2, '-')
axarr[1, 2].plot(energies, dos_f_G7_7[0]/2, '-')
axarr[2, 0].plot(energies, dos_f_G81_7[0]/2, '-')
axarr[2, 1].plot(energies, dos_f_G82_7[0]/2, '-')
axarr[1, 0].set_ylabel('DOS per f-orbital')
axarr[2, 1].set_xlabel('E (eV)')

plt.savefig('dos_1.pdf')

f, axarr = plt.subplots(2, 2, sharex=True, sharey=True)
axarr[0, 0].plot(energies, dos_t[0]/14, '-')
axarr[0, 1].plot(energies, (dos_f_G7_5[0] + dos_f_G81_5[0] + \
        dos_f_G82_5[0])/6,'-')
axarr[1, 0].plot(energies, (dos_f_G6_7[0] + dos_f_G7_7[0] + \
        dos_f_G81_7[0] + dos_f_G82_7[0])/8, '-')
axarr[1, 0].set_ylabel('DOS per f-orbital')
axarr[1, 0].set_xlabel('E (eV)')

plt.savefig('dos_2.pdf')



