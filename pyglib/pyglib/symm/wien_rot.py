from __future__ import print_function
from builtins import range
import glob
import sys
from scipy.linalg import det


def get_rotations(case_fname):
    with open(str(case_fname), 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'SYMMETRY OPERATIONS' in line:
                break
        nops = int(lines[i].split()[0])
        lines = lines[i+1:]
        rot_list = []
        for i in range(0,nops*4,4):
            rot = []
            lok = True
            for j in range(3):
                line = lines[i+j]
                rot.append([int(line[:2]), int(line[2:4]), int(line[4:6])])
                if abs(float(line[6:17])) > 1.e-6:
                    lok = False
                    break
            if lok:
                if abs(det(rot)-1)<1.e-6:
                    rot_list.append(rot)
    return rot_list


if __name__ == "__main__":
    files = glob.glob('*.struct')
    if len(files) < 1:
        print(' no wien2k struct file exists.')
        sys.exit()
    for f in files:
        rot_list = get_rotations(f)
        print(' number of rotations = {}'.format(len(rot_list)))
        for rot in rot_list:
            print(rot)
