import numpy as np


def jj_to_cubic_relativistic_harmonics(orbital='f'):
    if 'f' in orbital:
        jj_to_cubic = np.zeros((14,14))
        jj_to_cubic[0,0] = -np.sqrt(1./6.) # |5/2, -5/2>
        jj_to_cubic[4,0] =  np.sqrt(5./6.) # |5/2, +3/2> G7, 5/2, +
        jj_to_cubic[5,1] = -np.sqrt(1./6.) # |5/2, +5/2>
        jj_to_cubic[1,1] =  np.sqrt(5./6.) # |5/2, -3/2> G7, 5/2, -

        jj_to_cubic[4,2] =  np.sqrt(1./6.) # |5/2, +3/2>
        jj_to_cubic[0,2] =  np.sqrt(5./6.) # |5/2, -5/2> G81, 5/2, +
        jj_to_cubic[1,3] =  np.sqrt(1./6.) # |5/2, -3/2>
        jj_to_cubic[5,3] =  np.sqrt(5./6.) # |5/2, +5/2> G81, 5/2, -

        jj_to_cubic[3,4] =  1. # |5/2, +1/2> G82, 5/2, +
        jj_to_cubic[2,5] =  1. # |5/2, -1/2> G82, 5/2, -

        jj_to_cubic[13,6] =  np.sqrt(5./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,6] =  np.sqrt(7./12.) # |7/2, -1/2> G6, 7/2, +
        jj_to_cubic[ 6,7] =  np.sqrt(5./12.) # |7/2, -7/2>
        jj_to_cubic[10,7] =  np.sqrt(7./12.) # |7/2, +1/2> G6, 7/2, -

        jj_to_cubic[12,8] = -np.sqrt(3./4.) # |7/2, +5/2>
        jj_to_cubic[ 8,8] =  np.sqrt(1./4.) # |7/2, -3/2> G7, 7/2, +
        jj_to_cubic[ 7,9] =  np.sqrt(3./4.) # |7/2, -5/2>
        jj_to_cubic[11,9] = -np.sqrt(1./4.) # |7/2, +3/2> G7, 7/2, -

        jj_to_cubic[13,10] =  np.sqrt(7./12.) # |7/2, +7/2>
        jj_to_cubic[ 9,10] = -np.sqrt(5./12.) # |7/2, -1/2> G81, 7/2, +
        jj_to_cubic[ 6,11] = -np.sqrt(7./12.) # |7/2, -7/2>
        jj_to_cubic[10,11] =  np.sqrt(5./12.) # |7/2, +1/2> G81, 7/2, -

        jj_to_cubic[12,12] = -np.sqrt(1./4.)  # |7/2, +5/2>
        jj_to_cubic[ 8,12] = -np.sqrt(3./4.)  # |7/2, -3/2> G82, 7/2, +
        jj_to_cubic[ 7,13] =  np.sqrt(1./4.)  # |7/2, -5/2>
        jj_to_cubic[11,13] =  np.sqrt(3./4.)  # |7/2, +3/2> G82, 7/2, -
    else:
        raise Exception('UndefinedFunction')
    return jj_to_cubic


