#!/usr/bin/env python
from __future__ import print_function
import h5py, numpy
from pyglib.math.matrix_util import sym_array
from pyglib.gutz.usrqa import get_usr_input
from pyglib.basic.units import Ryd_eV

def get_b_field_list():
    '''get the list of magnetic field by q&a.
    '''
    with h5py.File('ginit.h5', 'r') as f:
        unit = f['/usrqa/unit'][()]
    if 'ryd' in unit:
        unit = Ryd_eV
    else:
        unit = 1.0
    f = h5py.File('GPARAM.h5', 'r')
    num_imp = f["/num_imp"][0]
    imap_imp = f['/IMAP_IMP'][()]
    iso = f['/iso'][0]
    print(' total {} impurities with equivalence indices \n {}'.format( \
            num_imp, imap_imp))
    bvec_list = []
    for imp in range(num_imp):
        if imap_imp[imp] == imp + 1:
            print('\n IMPURITY {}'.format(imp+1))
            if iso == 2:
                while True:
                    vec = raw_input(
                            '\n please enter direction (to be normalized) \n'+\
                            ' of local magnetic field by components \n'+\
                            ' in global coordinate system: x y z\n ')
                    vec = vec.split()
                    vec = numpy.array(map(float, vec))
                    vec = vec/numpy.linalg.norm(vec)
                    if len(vec) == 3:
                        break
                    else:
                        print(' enter 3 float numbers with finger space only.')
            else:
                updn = get_usr_input(' enter spin up or down:', ['up', 'dn'])
                if updn == 'up':
                    vec = numpy.array([0., 0., -1.])
                else:
                    vec = numpy.array([0., 0.,  1.])
            bm = raw_input('\n please enter the magnitude of the field: \n'+\
                    ' in unit of eV per Bohr magneton: \n ')
            bm = float(bm)/unit
            bvec_list.append(bm*vec)
        else:
            bvec_list.append(bvec_list[imap_imp[imp] - 1])
    f.close()
    givext = get_usr_input( \
            '\n Is the external field applied only at initial step (0) \n'+ \
            ' or fixed through all the iterations (1).', ['0', '1'])
    givext = int(givext)
    return bvec_list, givext


def get_vext_list(bvec_list):
    '''get the list of magnetic potential matrix in the single particle space,
    given the list of local magnetic field.
    '''
    vext_list = []
    max_sym_err = 0.
    with h5py.File('GPARAM.h5', 'r') as f:
        for imp, bvec in enumerate(bvec_list):
            prepath = "/IMPURITY_" + str(imp + 1)
            sx = f[prepath+'/SX'][()].T
            sy = f[prepath+'/SY'][()].T
            sz = f[prepath+'/SZ'][()].T
            vext = (bvec[0]*sx + bvec[1]*sy + bvec[2]*sz)*2
            r_list = numpy.swapaxes(f[prepath+'/SP_ROTATIONS'][()], 1, 2)
            vext_sym = sym_array(vext, r_list)
            max_sym_err = max(max_sym_err, numpy.max(numpy.abs( \
                    vext - vext_sym)))
            vext_list.append(vext_sym)
    print(' maximal symmetrization error of vext = {}'.format(max_sym_err))
    return vext_list


def h5wrt_gmagnet(vext_list, g_ivext, fn='GVEXT.h5'):
    with h5py.File(fn, 'w') as f:
        for imp, vext in enumerate(vext_list):
            f['/IMPURITY_{}/VEXT'.format(imp+1)] = vext.T
        f['/givext'] = [g_ivext]


def init_magnet_conf():
    '''
    initialize the the magnetic configuration for magnetic calculation.
    '''
    bvec_list, g_ivext = get_b_field_list()
    vext_list = get_vext_list(bvec_list)
    h5wrt_gmagnet(vext_list, g_ivext)



if __name__ == "__main__":
    init_magnet_conf()
