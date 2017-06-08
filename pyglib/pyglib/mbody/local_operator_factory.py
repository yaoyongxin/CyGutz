'''
Get lcoal operators given the {n_ij} operators.
'''

import h5py
import itertools as it
import numpy as np
from scipy.sparse import csc_matrix

def get_linear_operators(nij, uij):
    '''
    Calculate \sum_{i,j} nij_{i,j}*uij(i,j).
    '''
    res = None
    for i, n1, u1 in it.izip(it.count(), nij, uij):
        for j, n12, u12 in it.izip(it.count(), n1, u1):
            if n12 is None:
                continue
            if abs(u12)>1.e-12:
                if res is None:
                    res = u12*n12
                else:
                    res += u12*n12
            if i==j:
                continue
            if abs(uij[j,i])>1.e-12:
                if res is None:
                    res = uij[j,i]*n12.getH()
                else:
                    res += uij[j,i]*n12.getH()
    if res is None:
        res = csc_matrix((1,1))
    return res


def get_am_op_from_nij(svec, lvec, nij, op_list):
    '''
    Get angular momentum operators given the {n_ij} operators and
    coiefficient matrice of spin and orbital momentum operators.
    '''
    if op_list is None:
        return None

    s_z = get_linear_operators(nij, svec[2])
    s_p = get_linear_operators(nij, svec[0]+1.j*svec[1])
    l_z = get_linear_operators(nij, lvec[2])
    l_p = get_linear_operators(nij, lvec[0]+1.j*lvec[1])
    j_z, j_p = s_z + l_z, s_p + l_p
    res = {}

    if "Sx" in op_list:
        res["Sx"] = (s_p.getH() + s_p) / 2.0
    if "Sy" in op_list:
        res["Sy"] = (s_p - s_p.getH()) / 2.0j
    if "Sz" in op_list:
        res["Sz"] = s_z
    if "Lx" in op_list:
        res["Lx"] = (l_p.getH() + l_p) / 2.0
    if "Ly" in op_list:
        res["Ly"] = (l_p - l_p.getH()) / 2.0j
    if "Lz" in op_list:
        res["Lz"] = l_z
    if "Jx" in op_list:
        res["Jx"] = (j_p.getH() + j_p) / 2.0
    if "Jy" in op_list:
        res["Jy"] = (j_p - j_p.getH()) / 2.0j
    if "Jz" in op_list:
        res["Jz"] = j_z
    if "S2" in op_list:
        res["S2"] = s_p.getH() * s_p + s_z * s_z + s_z
    if "L2" in op_list:
        res["L2"] = l_p.getH() * l_p + l_z * l_z + l_z
    if "J2" in op_list:
        res["J2"] = j_p.getH() * j_p + j_z * j_z + j_z
    for key in res.keys():
        res[key] = res[key].tocsc()
    return res


def get_local_operators(imp, ival, op_list):
    '''
    Get the local operators, like Sx, Sy, Sz, Lx, Ly, Lz, S^2, L^2
    and J^2 operators for impurity imp.
    '''
    from pyglib.io.h5io import get_coo_matrix

    # Read in coefficient matrices for S and L.
    with h5py.File('GPARAM.h5', 'r') as f:
        svec = []
        svec.append(f['/IMPURITY_'+str(imp)+'/SX'][...].T)
        svec.append(f['/IMPURITY_'+str(imp)+'/SY'][...].T)
        svec.append(f['/IMPURITY_'+str(imp)+'/SZ'][...].T)
        lvec = []
        lvec.append(f['/IMPURITY_'+str(imp)+'/LX'][...].T)
        lvec.append(f['/IMPURITY_'+str(imp)+'/LY'][...].T)
        lvec.append(f['/IMPURITY_'+str(imp)+'/LZ'][...].T)

    # Read n_ij operators
    norb = svec[0].shape[0]
    nij = []
    with h5py.File('EMBED_HAMIL_ANALYSIS_'+str(imp)+'.h5', 'r') as f:
        for i in range(norb):
            n_i = []
            for j in range(i+1):
                path = '/valence_block_{}/NP_{}_{}'.format(ival,i+1,j+1)
                if path in f:
                    n_i.append(get_coo_matrix(f, path).tocsc())
                else:
                    n_i.append(None)
            nij.append(n_i)

    return get_am_op_from_nij(svec, lvec, nij, op_list)


def get_label_list(J, vecs):
    '''
    Get list of expection values of J for vecs.
    '''
    if J is None:
        label = None
    else:
        label = []
        for vec in vecs:
            res = np.vdot(vec, J.dot(vec))
            assert np.abs(res.imag) < 1.e-8, " error in get_label_list!"
            label.append(res.real)
        label = np.array(label)
    return label
