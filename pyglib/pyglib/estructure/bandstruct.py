import h5py


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
