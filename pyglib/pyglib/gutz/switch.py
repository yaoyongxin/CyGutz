from __future__ import print_function
import h5py, sys, os
from pyglib.gutz.init import initialize as ginit

def modify_gparam():
    '''convenient scipt to modify the settings in GPARAM.h5 file.
    '''
    if '-h' in sys.argv or len(sys.argv) == 1:
        print('\n inline argiment examples: \n' +
                ' -imix n -- change gimix to n \n' +
                ' -iembeddiag n -- change giembeddiag to n \n' +
                ' -dc_mode n -- change dc_mode to n \n' +
                ' -dc_nelf_list n1-n2-... -- change n_dc to n list \n' +
                ' -unique_j_ev j1-j2-... -- change unique_j_list_ev \n' +
                ' -unique_u_ev u1-u2-... -- change unique_u_list_ev \n' +
                ' -maxiter n -- change gmaxiter to n \n')
        return
    with h5py.File('GPARAM.h5', 'a') as f:
        if '-imix' in sys.argv:
            f['/gimix'][()] = [int(sys.argv[sys.argv.index( \
                    '-imix') + 1])]
        if '-iembeddiag' in sys.argv:
            f['/giembeddiag'][()] = [int(sys.argv[sys.argv.index( \
                    '-iembeddiag') + 1])]
        if '-dc_mode' in sys.argv:
            f['/dc_mode'][()] = [int(sys.argv[sys.argv.index( \
                    '-dc_mode') + 1])]
        if '-maxiter' in sys.argv:
            f['/gmaxiter'][()] = [int(sys.argv[sys.argv.index( \
                    '-maxiter') + 1])]
        if '-dc_nelf_list' in sys.argv:
            stmp = sys.argv[sys.argv.index('-dc_nelf_list')+1]
            f['/dc_nelf_list'][()] = [float(s) for s in
                    stmp.split('-')]

    # change GMOTT.h5 if necessary
    if os.path.isfile('GMOTT.h5'):
        with h5py.File('GMOTT.h5', 'a') as f:
            if '-iembeddiag' in sys.argv:
                f['/giembeddiag'][()] = [int(sys.argv[sys.argv.index( \
                        '-iembeddiag') + 1])]
    # change settings in 'ginit.h5', re-initialize if becessary
    re_init = False
    with h5py.File('ginit.h5', 'a') as f:
        if '-imix' in sys.argv:
            f['/usrqa/lnewton'][()] = int(sys.argv[sys.argv.index( \
                    '-imix') + 1])
        if '-iembeddiag' in sys.argv:
            f['/usrqa/iembeddiag'][()] = int(sys.argv[sys.argv.index( \
                    '-iembeddiag') + 1])
        if '-dc_mode' in sys.argv:
            f['/usrqa/ldc'][()] = int(sys.argv[sys.argv.index( \
                    '-dc_mode') + 1])
        if '-unique_j_ev' in sys.argv:
            stmp = sys.argv[sys.argv.index('-unique_j_ev')+1]
            f['/usrqa/unique_j_list_ev'][()] = [float(s) for s in
                    stmp.split('-')]
            re_init = True
        if '-unique_u_ev' in sys.argv:
            stmp = sys.argv[sys.argv.index('-unique_u_ev')+1]
            f['/usrqa/unique_u_list_ev'][()] = [float(s) for s in
                    stmp.split('-')]
            re_init = True

    if re_init:
        ginit()
