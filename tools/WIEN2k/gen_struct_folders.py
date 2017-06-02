#!/usr/bin/env python

import os
import sys
import numpy as np


def gen_struct_foldes(v_list, case, v_ref, file_ref):
    '''
    Generate a list of folders for vol in v_list.
    '''
    for i, vol in enumerate(v_list):
        foldername = case + "_" + str(i)
        ratio = (vol / v_ref)**(1. / 3)
        print 'vol/atom:', i, vol
        os.mkdir(foldername)
        with open(file_ref, 'r') as f_source:
            with open(foldername + "/" + foldername + ".struct", 'w') \
                    as f_target:
                for i, line in enumerate(f_source.readlines()):
                    if i == 3:
                        a_ref, b_ref, c_ref = float(line[:10]), \
                                float(line[10:20]), float(line[20:30])
                        a, b, c = a_ref * ratio, b_ref * ratio, c_ref * ratio
                        line = "{:10.6f}{:10.6f}{:10.6f}".format(a, b, c) + \
                                line[30:]
                    f_target.write(line)


if __name__ == "__main__":
    '''
    Call gen_struct_foldes with inline argiments.
    '''
    if 'h' in sys.argv[1]:
        print "Usage: \n" + \
                "      gen_struct_folders.py " + \
                "case v_i v_f delta_v v_ref file_ref \n"
        quit()
    case, vi, vf, dv, vref, f_ref = [
        sys.argv[1]] + [float(x) for x in sys.argv[2:6]] + [sys.argv[6]]
    v_list = np.arange(vi, vf, dv)
    gen_struct_foldes(v_list, case, vref, f_ref)
