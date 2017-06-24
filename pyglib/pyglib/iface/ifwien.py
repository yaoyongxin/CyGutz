import numpy as np
import scipy
import h5py
from pyglib.io.fio import file_exists


def trafoso(L):
    '''
    trafoso provides transformation matrices from
    |L,1/2,mL,mS> (L=0,1,2,3, mS=-1/2,1/2) basis to
    basis |J,L,S,mJ>, J=L-1/2, L+1/2
    H. Watanabe 'Operator Methods in Ligand Field Theory'
    Prentice Hall, 1966, Table 1.8-1.
    ordering because of the convention used in WIEN is:
                       mS=1/2        mS=-1/2
                     -L .......L  -L ...... L     (2*(2L+1) columns)
            -(L-1/2)
               .
    J=L-1/2    .
               .
             (L-1/2)
             -L-1/2
               .
    J=L+1/2    .
               .
              L+1/2
    '''
    cf = np.zeros((2*(2*L + 1), 2*(2*L + 1)))
    if L == 0:
        cf[0, 1] = 1.0
        cf[1, 0] = 1.0
    else:
        k1 = -1
        for ms in range(-1, 2, 2):
            ams = -ms / 2.
            for ml in range(-L, L + 1):
                k1 = k1 + 1
                k2 = -1
                for mj in range(-2*L + 1, 2*L, 2):  # L-1/2 states
                    amj = mj / 2.
                    k2 = k2 + 1
                    d = amj - ml - ams
                    if abs(d) < 0.0001:
                        if ms == 1:
                            cf[k2, k1] = - \
                                np.sqrt((L + 0.5 + amj) / (2*L + 1))
                        else:
                            cf[k2, k1] = np.sqrt((L + 0.5 - amj) / (2*L + 1))
                for mj in range(-2*L - 1, 2*L + 2, 2):  # L+1/2 states
                    amj = mj/2.
                    k2 = k2 + 1
                    d = amj - ml - ams
                    if abs(d) < 0.0001:
                        if ms == 1:
                            cf[k2, k1] = np.sqrt((L + 0.5 - amj) / (2*L + 1))
                        else:
                            cf[k2, k1] = np.sqrt((L + 0.5 + amj) / (2*L + 1))
    return cf


def get_orbital_transformation(l, qsplit):
    '''get the unitary transformation according to qsplit, which can simply
    be 3 for identity, and 4 for complex spherical Harmonics to relativisitc
    Harmonics.
    '''

    if qsplit == 3:
        u_trans = np.identity((2*l+1), dtype=complex)
    elif qsplit == 4:
        u_trans = trafoso(l)
    return u_trans


def h4set_indmfl(emin=-10., emax=10., projtype=-2):
    '''create case.indmfl for the interface given file ginit.h5.
    '''
    finit = h5py.File('ginit.h5', 'r')
    if '/struct/case' not in finit:
        raise ValueError('path /struct/case does not exist in ginit.h5!')

    case = finit['/struct/case'][()].split('.')[0]

    if file_exists(case+'.indmfl'):
        print(case+'.indmfl exists! \n' +
                'Please rm it if you intend to generate a new one.')
        return

    findm = open(case+'.indmfl', 'w')
    findm.write(('{:6.2f} {:6.2f} {:3d}  # hybrid. emin/max w.r.t. FS,' + \
            ' projector type\n').format(emin, emax, projtype))

    atom_symbols = finit['/struct/symbols'][()]
    unique_corr_atom_symbols = finit['/usrqa/unique_corr_symbol_list'][()]. \
            tolist()
    n_corr_atoms = 0
    for symbol in atom_symbols:
        if symbol in unique_corr_atom_symbols:
            n_corr_atoms += 1
    findm.write('{:2d}                 # number of correlated atoms\n'. \
            format(n_corr_atoms))

    if 'y' in finit['/usrqa/spin_orbit_coup'][()]:
        qsplit = 4; nspin = 2
    else:
        qsplit = 3; nspin = 1
    unique_df_list = finit['/usrqa/unique_df_list'][()]

    from pyglib.symm.angular_momentum_1p import get_l_list_from_string
    cix = 0
    norbsmax = 0
    l_super_list = []
    for i, symbol in enumerate(atom_symbols):
        if symbol not in unique_corr_atom_symbols:
            continue
        df = unique_df_list[unique_corr_atom_symbols.index(symbol)]
        l_list = get_l_list_from_string(df)
        findm.write('{:2d} {:2d}              # iatom, num_L\n'. \
                format(i+1, len(l_list)))
        l_super_list.extend(l_list)
        for l in l_list:
            cix += 1
            findm.write(' {:2d} {:2d} {:2d}          # L, qsplit, cix\n'.\
                    format(l, qsplit, cix))
            norbsmax = max(norbsmax, (2*l+1)*nspin)

    findm.write('----- unitary transformation for correlated orbitals -----\n')
    findm.write(('{:2d} {:2d}              # num indep. kcix blocks,' + \
            ' max dimension\n').format(cix, norbsmax))
    for i, l in enumerate(l_super_list):
        norbs = (2*l+1)*nspin
        findm.write('{:2d} {:2d}              # cix, dimension\n'. \
                format(i+1, norbs))
        findm.write('----- unitary transformation matrix -----\n')
        u_trans = get_orbital_transformation(l, qsplit)
        for u1 in u_trans:
            for u12 in u1:
                findm.write('{:20.16f}{:20.16f} '.format(u12.real, u12.imag))
            findm.write('\n')


if __name__=='__main__':
    h4set_indmfl()
