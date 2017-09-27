import h5py, numpy
from builtins import range
from pyglib.estructure.dos import get_bands

def get_greek_label(kname):
    '''return the possible greek label.
    '''
    if kname.upper() in ['LAMBDA', 'GAMMA', 'DELTA', 'SIGMA', 'THETA', \
            'XI', 'PI', 'UPSILON', 'PHI', 'PSI', 'OMEGA']:
        return '$\{}{}$'.format(kname[0].upper(), kname[1:].lower())
    else:
      return kname


def get_k_info():
    '''get the k distance list, special k-point label and position for
    the band structure plot.
    '''
    with h5py.File('GPARAMBANDS.h5', 'r') as f:
        kx = f['/kptx'][()]
        ky = f['/kpty'][()]
        kz = f['/kptz'][()]
        kn = f['/kptname'][()]
        # reciprocal lattice vectors ordered as colomn vectors
        br = f['/br2_car_dir'][()]

    ktick_pos = []
    ktick_label = []
    for i in range(len(kn)):
        if i == 0:
            klist = [0.]
        else:
            dk = numpy.array([kx[i]-kx[i-1], ky[i]-ky[i-1], kz[i]-kz[i-1]])
            dl = dk.dot(br.T).dot(br).dot(dk)
            klist.append(klist[i-1]+numpy.sqrt(dl))
        klabel = kn[i].strip()
        if klabel != '':
            ktick_pos.append(klist[i])
            ktick_label.append(get_greek_label(klabel))
    return klist, ktick_pos, ktick_label


def get_estimated_gap(nocc, fband='GBANDS_0.h5'):
    '''Given the number of occupied bands and the band energy file,
    estimate the band gap.
    '''
    # valence band maximum
    evmax = -1.e10
    # conduction band minimum
    ecmin = 1.e10
    with h5py.File(fband, 'r') as f:
        for ik in range(f['/IKP_START'][()], f['/IKP_END'][()]+1):
            ek = f['/ISPIN_1/IKP_{}/ek'.format(ik)][()]
            evmax = max(ek[nocc-1], evmax)
            ecmin = min(ek[nocc], ecmin)
    return max(ecmin - evmax, 0.)


def plot_band_sturture():
    '''plot band structure with overall correlated orbital character.
    '''
    klist, ktick_pos, ktick_label = get_k_info()
    e_skn, psi_sksna = get_bands()

