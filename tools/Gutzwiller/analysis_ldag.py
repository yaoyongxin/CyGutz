#!/usr/bin/env python
import numpy as np
import sys
import os
import h5py
import glob
sys.path.append(os.environ['WIEN_GUTZ_ROOT']+"/tools/Gutzwiller/")
from unitconv import Bohr_A, Ryd_eV, eVA_GPa
from datafetch import get_data_scf

def analysis_ldag():
  help_msg ='''
  No metadata will be saved! In order to save results, please rerun by specifying the metadata file to store the results:
    ${WIEN_GUTZ_ROOT}/tools/Gutzwiller/analysis_ldag.py -p save_path -g datagroup [total-energy-V|interaction-V|doble-counting-V]
    e.g.,
    ${WIEN_GUTZ_ROOT}/tools/Gutzwiller/analysis_ldag.py -p ../Actinides.h5 -g /Am/fcc/LDA/ total-energy-V '''

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

  dir_list = [dir for dir in os.listdir('./') if os.path.isdir(dir)]

  if 'total-energy-V' in sys.argv:
    analysis_total_energy_volume(dir_list, save_path, save_data, datagroup)
  if 'interaction-V' in sys.argv:
    analysis_energy_component_volume(dir_list, save_path, save_data, datagroup, 'E_POT2_U')
  if 'doble-counting-V' in sys.argv:
    analysis_energy_component_volume(dir_list, save_path, save_data, datagroup, 'E_DC2_U')

def analysis_total_energy_volume(dir_list, save_path, save_data, datagroup):
  '''
  Total energy vs volume data.
  '''
  energies = []; volumes = []

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
      del f[datagroup+'/EV/V'], f[datagroup+'/EV/E']
    except:
      pass
    f[datagroup+'/EV/E'] = energies
    f[datagroup+'/EV/V'] = volumes
    f.close()

  import matplotlib
  matplotlib.use('agg')
  import matplotlib.pyplot as plt

  plt.figure(0)
  plt.plot(volumes, energies, '-o')
  plt.xlabel("V ($\AA^{3}/f.u.$)")
  plt.ylabel(name+" ($eV/f.u.$)")

  plt.savefig('EV.pdf')

def analysis_energy_component_volume(dir_list, save_path, save_data, datagroup, name):
  '''
  Interaction energy vs volume data.
  '''
  energies = []; volumes = []

  for dir in dir_list:
    try:
      f = h5py.File(dir+'/glog.h5','r')
      energy = f["/"+name][0]
      volume = get_data_scf(dir+'/*scf', ':VOL')
      energies.append(energy)
      volumes.append(volume)
    except:
      pass

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
      del f[datagroup+'/EV/V'], f[datagroup+'/EV/'+name]
    except:
      pass
    f[datagroup+'/EV/'+name] = energies
    f[datagroup+'/EV/V'] = volumes
    f.close()

  import matplotlib
  matplotlib.use('agg')
  import matplotlib.pyplot as plt

  plt.figure(0)
  plt.plot(volumes, energies, '-o')
  plt.xlabel("V ($\AA^{3}/f.u.$)")
  plt.ylabel("E ($eV/f.u.$)")

  plt.savefig(name+'_V.pdf')


if __name__=="__main__":
  analysis_ldag()
