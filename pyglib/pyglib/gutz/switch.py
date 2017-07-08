from __future__ import print_function
import h5py,sys


def modify_gparam():
    '''convenient scipt to modify the settings in GPARAM.h5 file.
    '''
    if '-h' in sys.argv or len(sys.argv) == 1:
        print('\n inline argiment examples: \n' +
                ' -imix n -- change gimix to n \n' +
                ' -iembeddiag n -- change giembeddiag to n \n' +
                ' -dc_mode n -- change dc_mode to n \n' +
                ' -maxiter n -- chnage gmaxiter to n \n')
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
