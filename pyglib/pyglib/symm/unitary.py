import numpy as np


def jj_to_cubic_relativistic_harmonics(orbital='f'):
    if 'f' == orbital:
        jj_to_cubic = np.zeros((14,14))
        jj_to_cubic[0,8] = -np.sqrt(1./6.) # |5/2, -5/2>
        jj_to_cubic[4,8] =  np.sqrt(5./6.) # |5/2, +3/2> G7, 5/2, +
        jj_to_cubic[5,10] = -np.sqrt(1./6.) # |5/2, +5/2>
        jj_to_cubic[1,10] =  np.sqrt(5./6.) # |5/2, -3/2> G7, 5/2, -

        jj_to_cubic[4,0] =  np.sqrt(1./6.) # |5/2, +3/2>
        jj_to_cubic[0,0] =  np.sqrt(5./6.) # |5/2, -5/2> G81, 5/2, +
        jj_to_cubic[1,2] =  np.sqrt(1./6.) # |5/2, -3/2>
        jj_to_cubic[5,2] =  np.sqrt(5./6.) # |5/2, +5/2> G81, 5/2, -

        jj_to_cubic[3,4] =  1. # |5/2, +1/2> G82, 5/2, +
        jj_to_cubic[2,6] =  1. # |5/2, -1/2> G82, 5/2, -

        jj_to_cubic[13,12] =  np.sqrt(5./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,12] =  np.sqrt(7./12.) # |7/2, -1/2> G6, 7/2, +
        jj_to_cubic[ 6,13] =  np.sqrt(5./12.) # |7/2, -7/2>
        jj_to_cubic[10,13] =  np.sqrt(7./12.) # |7/2, +1/2> G6, 7/2, -

        jj_to_cubic[12,11] = -np.sqrt(3./4.) # |7/2, +5/2>
        jj_to_cubic[ 8,11] =  np.sqrt(1./4.) # |7/2, -3/2> G7, 7/2, +
        jj_to_cubic[ 7,9] =  np.sqrt(3./4.) # |7/2, -5/2>
        jj_to_cubic[11,9] = -np.sqrt(1./4.) # |7/2, +3/2> G7, 7/2, -

        jj_to_cubic[13,7] =  np.sqrt(7./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,7] = -np.sqrt(5./12.) # |7/2, -1/2> G81, 7/2, +
        jj_to_cubic[ 6,5] = -np.sqrt(7./12.) # |7/2, -7/2>
        jj_to_cubic[10,5] =  np.sqrt(5./12.) # |7/2, +1/2> G81, 7/2, -

        jj_to_cubic[12,3] = -np.sqrt(1./4.)  # |7/2, +5/2>
        jj_to_cubic[ 8,3] = -np.sqrt(3./4.)  # |7/2, -3/2> G82, 7/2, +
        jj_to_cubic[ 7,1] =  np.sqrt(1./4.)  # |7/2, -5/2>
        jj_to_cubic[11,1] =  np.sqrt(3./4.)  # |7/2, +3/2> G82, 7/2, -
    else:
        raise ValueError('UndefinedFunction')
    return jj_to_cubic


def comp_sph_harm_to_real_harm(orbital='d'):
    if 'd' == orbital:
        # Unitary transformation from complex spherical harmonics to
        # real harmonics {xy, yz, z2, zx, x2-y2}
        # i.e., matrix < cmp_sph_harm | real_harm >
        csh2rh = np.zeros((5, 5), dtype=complex)
        csh2rh[0, 0] =  1.j/np.sqrt(2.); csh2rh[4, 0] = -1.j/np.sqrt(2.)
        csh2rh[1, 1] =  1.j/np.sqrt(2.); csh2rh[3, 1] =  1.j/np.sqrt(2.)
        csh2rh[2, 2] = 1.0
        csh2rh[1, 3] = 1./np.sqrt(2); csh2rh[3, 3] = -1./np.sqrt(2.)
        csh2rh[0, 4] = 1./np.sqrt(2.); csh2rh[4, 4] = 1./np.sqrt(2.)
    else:
        raise ValueError('UndefinedFunction')
    return csh2rh
