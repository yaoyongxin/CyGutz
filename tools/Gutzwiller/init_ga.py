#!/usr/bin/env python

import os
import sys
import numpy as np
import shutil

def help():
    '''
    Print instructions.
    '''
    message = '''
    Usage:

    init_ga.py [OPTION]

    The script sets up the necessary input files for the Gutzwiller/Slave-Boson solver.
    By default, it will choose the rotated local orbital basis {hs} to work with,
    since it could be computationally more efficient than using the original basis.
    However, the script also prints WH_HS.INP_ORIG which contains {hs} in the original basis
    for consistence check.

    The script will generate file "WH_RLNELF.INP" if "WH_RLNEF.OUT_ORIG" exists.
    And it assumes that "WH_RLNELF.OUT_ORIG" is some previous converged solution
    with the original local basis.
    Tips: 1) you can use "init_ga.py [OPTION] < init_ga.input" to repeat the initialization,
             where init_ga.input was generated at the first run.
          2) init_ga.slog records the important information in the session.

    -c rcut_symmetry  Assign the cut-off distance for extracting the atom-centered"
                      cluster for symmetry analysis."
    -vasp             Set struct file in VASP format. By default it is in Wien2K format."
    '''

    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print message
        sys.exit(0)


def initialize():
    '''
    Initialization for the Gutzwiller/Slave-Boson solver.
    Store symmetry informations to init_ga_info.h5.
    '''
    log_file = open("init_ga.slog", 'w')
    print >> log_file, "The executed commant line is:"
    print >> log_file, " ".join(sys.argv[:])
    usr_input = open("input.slog", 'w')
    from structure import get_myatoms, check_material

    # get the material
    material = get_myatoms()

    # Q & A
    from usrqa import usr_qa_setup
    usr_qa_setup(material, log_file, usr_input)

    # check some info in log
    check_material(material, log_file)

    # write input files for CyGutz
    from gl_inp import set_gl_inp_material
    set_gl_inp_material(material, log_file, usr_input)

    # ga_init_dmft.py inputs.
    from usrqa import inp_ga_init_dmft
    inp_ga_init_dmft(material.get_mySOC(), log_file)

    log_file.close()
    usr_input.close()
    os.rename("input.slog", "init_ga.input")

    # save data for further process.
    from h5save_init_ga_info import h5save_init_ga_info
    h5save_init_ga_info(material)

if __name__ == '__main__':

    help()
    initialize()
