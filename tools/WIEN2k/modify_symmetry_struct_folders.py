#!/usr/bin/env python

import glob
import os
import sys
import numpy as np


def modify_symmetry_struct_foldes(file_ref):
    '''
    Change the symmetry of the structure file.
    '''
    dir_list = [_dir for _dir in os.listdir('./') if os.path.isdir(_dir)]
    for i, _dir in enumerate(dir_list):
        struct = glob.glob(_dir + '/*struct')
        if len(struct) != 1:
            print ' Skip ' + _dir + ': no or more than one struct files!'
            continue
        struct = struct[0]
        with open(struct, 'r') as f:
            line_3 = f.readlines()[3]
        with open(file_ref, 'r') as f_source:
            with open(struct, 'w') as f_target:
                for i, line in enumerate(f_source.readlines()):
                    if i == 3:
                        line = line_3
                    f_target.write(line)


if __name__ == "__main__":
    '''
    Call modify_symmetry_struct_foldes with inline argiments.
    '''
    if 'h' in sys.argv[1]:
        print "Usage: \n      modify_symmetry_struct_folders.py file_ref \n"
        quit()
    f_ref = sys.argv[1]
    modify_symmetry_struct_foldes(f_ref)
