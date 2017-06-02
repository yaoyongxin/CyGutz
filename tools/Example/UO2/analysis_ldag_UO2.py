import numpy as np
import sys
sys.path.append("/home/ykent/BitBucket/pyscript/Gutzwiller/")
from unitconv import Bohr_A, Ryd_eV, eVA_GPa
from scipy.linalg import eigvalsh
import h5py

f = h5py.File('/home/ykent/WIEN_GUTZ/EXAMPLE/UO2/UO2_ldag_data.h5', 'a')

path = "/ldag_soc_cf/Ecut_12/U8.0/J0.60"

vi, vf = 0, 8
v_list = np.arange(vi, vf, 1)

path_list_scf = []
path_list = []
for i,vol in enumerate(v_list):
  path_list_scf.append('./' + '/UO2_' + str(vol) + '/UO2_' + str(vol) + '.scf')
  path_list.append('./' + '/UO2_' + str(vol) + '/')

from datafetch import get_data_list_scf
etot_list = get_data_list_scf(path_list_scf, ":ENE")
etot_list = np.array(etot_list)*Ryd_eV
try:
  f[path + "/etot_list"] = etot_list
except:
  f[path + "/etot_list"][...] = etot_list

v_list = get_data_list_scf(path_list_scf, ":VOL")
v_list = np.array(v_list)*Bohr_A**3
try:
  f[path + "/v_list"] = v_list
except:
  f[path + "/v_list"][...] = v_list

from datafetch import get_R_list, get_Z_evalues
R_list = get_R_list(path_list, 1)
Z_e_values = []
for i,z_list in enumerate(get_Z_evalues(R_list, 0, 2)):
  Z_e_values.append([])
  for z in z_list:
    Z_e_values[i].append(z)
for i,z_list in enumerate(get_Z_evalues(R_list, 4, 6)):
  for z in z_list:
    Z_e_values[i].append(z)
for i,z_list in enumerate(get_Z_evalues(R_list, 12, 13)):
  for z in z_list:
    Z_e_values[i].append(z)

Z_e_values = np.array(Z_e_values)
try:
  f[path + "/Z_eigen_values"] = Z_e_values
except:
  f[path + "/Z_eigen_values"][...] = Z_e_values

from datafetch import get_Nab_list
Nab_list = get_Nab_list(path_list, 1)
N_e_values = []
for i,Nab in enumerate(Nab_list):
  N_e_values.append([])
  vals = eigvalsh(Nab[0:2, 0:2])
  for val in vals:
    N_e_values[i].append(val.real)
  vals = eigvalsh(Nab[4:6, 4:6])
  for val in vals:
    N_e_values[i].append(val.real)
  N_e_values[i].append(Nab[12,12].real)
N_e_values = np.array(N_e_values)
try:
  f[path + "N_eigen_values"] = N_e_values
except:
  f[path + "N_eigen_values"][...] = N_e_values

f.close()

import matplotlib.pyplot as plt

plt.figure(0)
plt.plot(v_list, etot_list)

plt.figure(1)
plt.plot(v_list, Z_e_values[:,0])
plt.plot(v_list, Z_e_values[:,1])
plt.plot(v_list, Z_e_values[:,2])
plt.plot(v_list, Z_e_values[:,3])
plt.plot(v_list, Z_e_values[:,4])

plt.figure(2)
plt.plot(v_list, N_e_values[:,0])
plt.plot(v_list, N_e_values[:,1])
plt.plot(v_list, N_e_values[:,2])
plt.plot(v_list, N_e_values[:,3])
plt.plot(v_list, N_e_values[:,4])


plt.show()
