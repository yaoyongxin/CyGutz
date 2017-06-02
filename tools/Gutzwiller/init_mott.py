#!/usr/bin/env python

import h5py
import numpy as np


def get_shrink_sigma(sigma):
    '''
    Shrink the structure of sigma if numbers skipped.
    '''
    elem_max = np.max(sigma)
    for elem in range(elem_max):
        if elem == np.max(sigma):
            break
        while elem + 1 not in sigma:
            spots = np.where(sigma > elem)
            sigma[spots] -= 1
    return sigma


def print_sigma(sigma, title):
    '''
    Print the structure of sigma with basis indices.
    '''
    print(title)
    print("index  " + ''.join("%4d " % (i) for i in range(len(sigma))) + '\n')
    for i, row in enumerate(sigma):
        print("%4d   " % (i) + ''.join("%4d " % (j) for j in row))


def init_mott(fname):
    '''
    Generate the input files for Mott phase calculation.
    '''
    sigma_L_list = []
    sigma_R_list = []
    orb_mott_list = []
    ne_mott_list = []
    f = h5py.File(fname, 'r')
    num_imp = f["/num_impurity"][...]
    for imp in range(num_imp):
        sigma = f["/impurity_" + str(imp) + "/sigma"][...]
        print "**********  Impurity " + str(imp) + "  **********"
        print_sigma(sigma, ' Sigma structure:')
        while True:
            orb_mott = raw_input(
                ' Please provide the indices of orbitals to be Mott localized \n (e.g., 0 2 ): ')
            yn = raw_input(
                ' You selected [' +
                orb_mott +
                '] to be Mott localized, right? (y/n):')
            if 'y' in yn or 'Y' in yn:
                orb_mott = orb_mott.split()
                ind_orb_mott = [int(s) for s in orb_mott]
                while True:
                    ne_mott = raw_input(
                        ' Please provide the total number of Mott localized electrons (per unit cell): ')
                    yn = raw_input(
                        ' Total ' +
                        ne_mott +
                        ' electrons will be Mott localized, right? (y/n):')
                    if 'y' in yn or 'Y' in yn:
                        ne_mott_list.append(int(ne_mott))
                        break
                break
        orb_mott_list.append(ind_orb_mott)
        for ind in ind_orb_mott:
            elem = sigma[ind, ind]
            sigma[ind, :] = 0
            sigma[ind, ind] = elem
        sigma = get_shrink_sigma(sigma)
        print_sigma(sigma, ' R structure:')
        sigma_R_list.append(np.copy(sigma))
        elem_mott = 0
        for ind in ind_orb_mott:
            elem = sigma[ind, ind]
            if elem_mott == 0:
                elem_mott = elem
            sigma[:, ind] = 0
            sigma[ind, ind] = elem_mott
        sigma = get_shrink_sigma(sigma)
        print_sigma(sigma, ' Lambda structure:')
        sigma_L_list.append(np.copy(sigma))
    from gl_inp import set_wh_hs
    set_wh_hs(sigma_R_list, file_name="WH_HS_R.INP", lsym=False)
    set_wh_hs(sigma_L_list, file_name="WH_HS_L.INP")
    from fileio import write_sigma_struct
    write_sigma_struct(sigma_R_list, file_name="WH_SIGMA_STRUCT_R.INP")
    write_sigma_struct(sigma_L_list, file_name="WH_SIGMA_STRUCT_L.INP")
    from fileio import write_frozen
    write_frozen(ne_mott_list, orb_mott_list)
    # Change LGPRJ=11 for the calculation on Mott phase.
    from shutil import copyfile
    copyfile('GL.INP', 'GL.INP_TEMPLATE')
    # Get the content (line with LGPRJ) to be modified.
    with open("GL.INP_TEMPLATE", 'r') as f:
        for line in f.readlines():
            if 'LGPRJ' in line:
                if '#' in line:
                    target_str = line[:line.index('#')]
                else:
                    target_str = line
                break
    with open("GL.INP_TEMPLATE", 'r') as source:
        with open("GL.INP", 'w') as target:
            target.write(source.read().replace(target_str, 'LGPRJ = 11   '))

if __name__ == "__main__":
    init_mott("init_ga_info.h5")
