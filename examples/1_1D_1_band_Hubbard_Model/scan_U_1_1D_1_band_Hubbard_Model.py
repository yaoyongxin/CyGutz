import sys
import os
import glob
import numpy as np
import subprocess as sp
import shutil
import h5py
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools/Gutzwiller"))
from dataproc import get_csr_matrix

def scan_U_1D_1_band_Hubbard_Model(U_lower_bound, U_upper_bound, U_step):
  '''
  In this script, you will learn:

  * Perform Gutzwiller calculations with U in range(U_lower_bound, U_upper_bound, U_step).
  * Store the key quantities, including double_occupany d, Z(=R^\dagger R), interaction energy e_int, total energy e_tot as a function of U, and generate figures for them.

  Please click [source] on upper right to see detailed notes about the script. Note that here we have only one impurity.

  Type the following line to run the script::

    $ python ${WIEN_GUTZ_ROOT}/examples/1_1D_1_band_Hubbard_Model/scan_1D_1_band_Hubbard_Model.py

  '''
  # Prepare some lists to store the quantities of interest.
  d_list = []; Z_list = []; e_int_list = []; e_tot_list = []

  # Make a copy of GL.INP to GL.INP_TEMPLATE, since we will modify the U value in that file.
  shutil.copyfile("GL.INP", "GL.INP_TEMPLATE")

  # Get the content (line with u,j) to be modified.
  with open("GL.INP_TEMPLATE", 'r') as f:
    for line in f.readlines():
      if 'u,j' in line:
        target_str = line[:line.index('#')]
        break

  # begin loop over U
  U_list = np.arange(U_lower_bound, U_upper_bound, U_step)
  for U in U_list:
    # Modify parameter U in GL.INP
    with open("GL.INP_TEMPLATE", 'r') as source:
      with open("GL.INP", 'w') as target:
        target.write(source.read().replace(target_str, str(U) + '  0.0 '))
    # Do some cleaning
    for f in glob.glob("GL*TMP"):
      os.remove(f)
    # Execute CyGutz
    sp.call('${WIEN_GUTZ_ROOT}/CyGutz', shell=True)

    # Check convergence
    f = h5py.File("glog.h5", 'r')
    max_err = f["/GA_MAX_ERR"][0]
    if max_err > .001:
      # Try the default initial input.
      os.remove("WH_RLNEF.INP")
      # Execute CyGutz
      sp.call('${WIEN_GUTZ_ROOT}/CyGutz', shell=True)
    f.close()

    # Save the Gutzwiller solution for the next as better initial input.
    shutil.copyfile("WH_RLNEF.OUT", "WH_RLNEF.INP")

    # Store double occupancy, Z, e_int and e_tot
    f = h5py.File("glog.h5", 'r')

    # Interaction energy
    e_int_list.append(f["/E_POT2_U"][0])

    # Total energy
    e_tot_list.append(f["/E_TB_TOT"][0])

    # Z=R^\dagger R, proportional to identity in this case
    R = f["/Impurity_1/GA_R"][...].T; Z = np.dot(np.conj(R).T,R)
    Z_list.append(Z[0,0].real)

    # Double occupancy, the fourth (one-based) diagonal element in the local reduced density matrix \rho ( = \phi \phi^(\dagger)) in this case.
    # Note the order of the Fock basis is |0>, |1down>, |1up>, |1down, 1up>.
    # Read \phi in sparse matrix format
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    d_list.append(rho[3,3].real) # python uses zero-based convention.
    f.close()

  # Use matplotlib to plot figures.
  import matplotlib
  # Good if no display is installed
  matplotlib.use('agg')
  import matplotlib.pyplot as plt
  f, axarr = plt.subplots(2, 2, sharex = True)
  axarr[0, 0].plot(U_list, d_list, '-o')
  axarr[0, 0].set_ylabel('double occupancy')
  axarr[0, 1].yaxis.tick_right()
  axarr[0, 1].yaxis.set_label_position("right")
  axarr[0, 1].plot(U_list, Z_list, '-o')
  axarr[0, 1].set_ylabel('Z')
  axarr[1, 0].plot(U_list, e_int_list, '-o')
  axarr[1, 0].set_ylabel('interaction_energy/t')
  axarr[1, 0].set_xlabel('U/t')
  axarr[1, 1].yaxis.tick_right()
  axarr[1, 1].yaxis.set_label_position("right")
  axarr[1, 1].plot(U_list, e_tot_list, '-o')
  axarr[1, 1].set_ylabel('total_energy/t')
  axarr[1, 1].set_xlabel('U/t')

  plt.savefig('1DZEU.png')


if __name__=="__main__":
  scan_U_1D_1_band_Hubbard_Model(0.0, 12.0, 0.2)
