#!/usr/bin/env python

from __future__ import print_function
import sys
import subprocess
import h5py
import itertools as it
import numpy as np
import numpy.linalg as LA
from scipy.linalg import expm
import timeit
import pyglib.math.matrix_basis as mb
from pyglib.io.h5io import get_csr_matrix
import Model.write_bndu as bndu


def get_hamilt_1d_nn_chain(n, e0, tp):
    v = -1.0 * np.ones(n - 1)
    v[0] = tp
    hamilt = np.zeros((n, n))
    hamilt[0, 0] = e0
    hamilt[np.arange(n - 1), np.arange(1, n)] = v
    hamilt[np.arange(1, n), np.arange(n - 1)] = v
    return hamilt


def get_1st_site_double_ocupancy(f):
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    c1up = get_csr_matrix(f, "/Impurity_1/annihi.op._1")
    c1dn = get_csr_matrix(f, "/Impurity_1/annihi.op._2")
    d1op = c1up.getH()*c1up*c1dn.getH()*c1dn
    nbase = d1op.shape[0]-rho.shape[0]
    return np.sum((rho * d1op[nbase:,nbase:]).diagonal())


def gutzwiller_energy(xs, hmat_basis, evals, wf, t_matrix):
    norb = len(evals)
    u = np.zeros([norb, norb], dtype=complex)
    for x, h in it.izip(xs, hmat_basis):
        if np.abs(x) < 1.e-16:
            continue
        else:
            u += x*h
#    u = expm(-1.j*u)
#    uwf = np.asarray([[np.dot(u, wf)]])

    uwf = np.asarray([[wf.T]])

    # Missing two additional dimensions: spin and k-points
    print(" Writing bndu...")
    bndu.write_bndu([[evals]], uwf.real, uwf.imag, ni=0)
    print(" Running CyGutz...")
    out, err = subprocess.Popen(
            'CyGutz', shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate()
    print(" Analysing results...")
    f = h5py.File('glog.h5', 'r')
    res = f['/E_TB_TOT'][0]

    # R
    R = f['/Impurity_1/GA_R'][...].T
    # spin degeneracy
    norb_c = R.shape[0]/2
    R = np.array([R[2*i, 2*j] for i in range(norb_c) \
            for j in range(norb_c)]).reshape((norb_c, norb_c))

    # impurity density matrix
    dm_imp = f['/Impurity_1/GA_NC_PHY'][...].T
    dm_imp = np.array([dm_imp[2*i, 2*j] for i in range(norb_c) \
            for j in range(norb_c)]).reshape((norb_c, norb_c))

    psi = f['/BND_VK'][0, 0, 0, :, :].T
    psi_wt = f['/BND_FERWE'][0, 0, :]/2

#    psi[:norb_c, :] = np.dot(np.conj(R.T), psi[:norb_c, :])

    dm = np.einsum('ij, j, kj->ki', psi, psi_wt, np.conj(psi))


    print(dm[:3, :3])

    dm[:norb_c, :norb_c] = dm_imp
    d_occ = get_1st_site_double_ocupancy(f)
    e_U = d_occ*0
    print(" Interaction energy = {}".format(e_U))
    e_t = np.einsum('ij, ij', t_matrix, dm*2)

    print(dm[:3, :3])
    print(t_matrix[:7, :7])

    print(" Total energy = {} vs {}".format(e_U + e_t, res))


def calc_1d_chain(nsites, e0, tp):
    print(" Get bare Hamiltonian...")
    t_matrix = get_hamilt_1d_nn_chain(nsites, e0, tp)
    print(" Diagonalize bare Hamiltonian...")
    evals, wf = LA.eigh(t_matrix)
    print(" Generating Hermitian matrix basis...")
    t0 = timeit.default_timer()
    hmat_basis = mb.hermitian_csc_matrix_basis(nsites, istart=1)
    t1 = timeit.default_timer()
    print(" Time used {} seconds.".format(t1 - t0))
    # initial guess
    xs = np.zeros((len(hmat_basis)))
    gutzwiller_energy(xs, hmat_basis, evals, wf, t_matrix)



if __name__=='__main__':
    if '-h' in sys.argv:
        print(' Usage: finite_1d_chain.py [-n num_sites -t tp] \n' +
                '        finite_1d_chain.py -n 500 -t -1.0')
        sys.exit(0)

    if '-n' in sys.argv:
        nsites = int(sys.argv[sys.argv.index('-n') + 1])
    else:
        nsites = 500
    if '-t' in sys.argv:
        tp = float(sys.argv[sys.argv.index('-t') + 1])
    else:
        tp = -1.0
    if '-e' in sys.argv:
        e0 = float(sys.argv[sys.argv.index('-e') + 1])
    else:
        e0 = 0.0
    calc_1d_chain(nsites, e0, tp)

