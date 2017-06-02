#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from scipy.sparse import csr_matrix
import os
import h5py
import sys
from hs_rotation_representation import h5get_hs_rotation_representation_N
from multiplets_analysis_lib import generate_symmetrized_rho_f, get_local_histogram


def multiplets_analysis_soc_cf_mott():
    '''
    Multiplets analysis for the Mott phase of LDA+G calculations
    with spin-orbit coupling (SOC) and crystal field effect (CF).
    '''
    # Impurity site, one-based
    if '-i' in sys.argv:
        imp = int(sys.argv[sys.argv.index('-i') + 1])
    else:
        imp = 1
    # Valence range
    if '-nbot' in sys.argv:
        n_bot = int(sys.argv[sys.argv.index('-nbot') + 1])
    else:
        n_bot = 0
    if '-ntop' in sys.argv:
        n_top = int(sys.argv[sys.argv.index('-ntop') + 1])
    else:
        n_top = None
    if '-nev' in sys.argv:
        num_ev = int(sys.argv[sys.argv.index('-nev') + 1])
    else:
        num_ev = 0

    print(" Multiplets analysis for impurity %2d" % (imp))

    if not os.path.isfile("hilbert_space_rotations.h5"):
        print(" Generate hilbert_space_rotations.")
        h5get_hs_rotation_representation_N(imp)

    np.set_printoptions(precision=2, suppress=True)
    f = h5py.File("glog.h5", 'r+')
    if "/Impurity_" + str(imp) + "/RHO_SYM.data" not in f:
        print(" Generate symmetrized \RHO.")
        generate_symmetrized_rho_f(imp, '/RHO')

    f_p = h5py.File("hilbert_space_rotations.h5", 'r')
    from hs_rotation_representation import get_hs_rotations
    Rpr_list = get_hs_rotations(f_p, imp)
    f_p.close()

    # main analysis routine
    print(" Get local histogram.")
    vals, n_label, j_label, chi_label, multiplet_degeneracies \
        = get_local_histogram(1, f, num_ev=num_ev, Rpr_list=Rpr_list,
                rho_name='/RHO_SYM', n_bot=n_bot, n_top=n_top)
    f.close()

    # brief info
    print(" weights: \n", vals)
    print(" sum-weight: \n", np.sum(vals))
    print(" n_label: \n", n_label)
    print(" j_label: \n", j_label)
    print(" chi_label: \n", chi_label[:min(20, len(vals))])
    print(" deg_label: \n", multiplet_degeneracies)

    # file to store the results.
    f = h5py.File("multiplets.h5", 'w')
    f["/weights"] = vals
    f["/n_label"] = n_label
    f["/j_label"] = j_label
    for i, chi in enumerate(chi_label):
        f["/chi_label" + str(i)] = chi
    f["/deg_label"] = multiplet_degeneracies
    f.close()

    import matplotlib.pyplot as plt
    x = range(len(j_label))
    f, axarr = plt.subplots(4, sharex=True)
    axarr[0].plot(x, vals, 'o-')
    axarr[0].set_ylabel('weight')
    axarr[1].plot(x, n_label, 'o-')
    axarr[1].set_ylabel('N')
    axarr[2].plot(x, multiplet_degeneracies, 'o-')
    axarr[2].set_ylabel('mult. deg.')
    axarr[3].plot(x, j_label, 'o-')
    axarr[3].set_ylabel('J')
    axarr[3].set_xlabel('multiplets')
    plt.show()


if __name__ == "__main__":
    multiplets_analysis_soc_cf_mott()
