#!/usr/bin/env python
import numpy as np
import sys
import os
import h5py
import glob

help_msg ='''
No metadata will be saved! In order to save results, please rerun by specifying the metadata file to store the results:
  ${WIEN_GUTZ_ROOT}/tools/Gutzwiller/analysis_total_energy.py -p save_path -g datagroup
  e.g.,
  ${WIEN_GUTZ_ROOT}/tools/Gutzwiller/analysis_total_energy.py -p ../Actinides.h5 -g /Am/fcc/LDA/EV'''

if '-p' in sys.argv:
  save_path = sys.argv[sys.argv.index('-p')+1]
  save_data = True
else:
  print help_msg
  save_data = False
if save_data:
  if '-g' in sys.argv:
    datagroup = sys.argv[sys.argv.index('-g')+1]
  else:
    save_data = False
    print help_msg

sys.path.append(os.environ['WIEN_GUTZ_ROOT']+"/tools/Gutzwiller/")
from unitconv import Bohr_A, Ryd_eV, eVA_GPa
from datafetch import get_data_scf

scf_file_list = []

dir_list = [dir for dir in os.listdir('./') if os.path.isdir(dir)]
energies = []
volumes = []

for dir in dir_list:
  if len(glob.glob(dir+'/*scf'))==0: continue
  energy = get_data_scf(dir+'/*scf', ':ENE')
  volume = get_data_scf(dir+'/*scf', ':VOL')
  energies.append(energy)
  volumes.append(volume)

energies = np.array(energies)*Ryd_eV
volumes  = np.array(volumes)*Bohr_A**3

# Sort according to volume
idx = volumes.argsort()
volumes = volumes[idx]
energies = energies[idx]

# Store the results
if save_data:
  f = h5py.File(save_path, 'a')
  try:
    del f[datagroup]
  except:
    pass
  f[datagroup+'/E'] = energies
  f[datagroup+'/V'] = volumes
  f.close()

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

plt.figure(0)
plt.plot(volumes, energies, '-o')
plt.xlabel("V ($\AA^{3}/f.u.$)")
plt.ylabel("E ($eV/f.u.$)")

plt.savefig('EV.pdf')
