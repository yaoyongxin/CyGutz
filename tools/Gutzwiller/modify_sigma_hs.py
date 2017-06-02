import h5py
import numpy as np
import sys
sys.path.append("/home/ykent/BitBucket/pyscript/Gutzwiller")

sigma_list = []
sigma_list.append(np.zeros((14,14),dtype=int))
sigma_list = np.array(sigma_list)

sigma_list[0,0,0] = sigma_list[0,2,2] = 1
sigma_list[0,1,1] = sigma_list[0,3,3] = 2
sigma_list[0,0,1] = sigma_list[0,2,3] = 3
sigma_list[0,1,0] = sigma_list[0,3,2] = 4
sigma_list[0,4,4] = sigma_list[0,6,6] = sigma_list[0,8,8] = sigma_list[0,10,10] = 5
sigma_list[0,5,5] = sigma_list[0,7,7] = sigma_list[0,9,9] = sigma_list[0,11,11] = 6
sigma_list[0,12,12] = sigma_list[0,13,13] = 7

print 'Modified sigma = \n', sigma_list

from gl_inp import set_wh_hs
set_wh_hs(sigma_list)
from fileio import write_sigma_struct
write_sigma_struct(sigma_list)
