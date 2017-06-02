'''
Transformation of the GA/Slave-boson solutions.
'''
import numpy as np
import itertools as it
from scipy.linalg import block_diag
from fileio import read_RLNEF, write_RLNEF


def get_gauge_solution(U_All, fname='WH_RLNEF.OUT'):
    '''
    Generate another equivalent solution due to gauge invariance,
    given the old solution in fname. Note Hs might be a concern
    if symmetry is included.
    '''
    R_All, LA1_All, NKS_All, E_Fermi = read_RLNEF(fname)
    for i,U,R,LA1,NKS in it.izip(it.count(), U_All, R_All, LA1_All, NKS_All):
        Udagger = np.conj(U.T)
        R_All[i] = np.dot(Udagger, R)
        LA1_All[i] = np.dot(Udagger, np.dot(LA1, U))
        NKS_All[i] = np.dot(U.T, np.dot(NKS, Udagger.T))
    write_RLNEF("WH_RLNEF.INP", R_All, LA1_All, NKS_All, E_Fermi)


def get_gauge_solution_nks(fname='WH_RLNEF.OUT'):
    '''
    Generate another equivalent solution due to gauge invariance,
    given the old solution in fname.
    '''
    from numpy.linalg import eigh
    R_All, LA1_All, NKS_All, E_Fermi = read_RLNEF(fname)
    for i,R,LA1,NKS in it.izip(it.count(), R_All, LA1_All, NKS_All):
        _, U = eigh(NKS[:2,:2])
        U = np.conj(U)
        U = block_diag(U, U)

        Udagger = np.conj(U.T)
        R_All[i] = np.dot(Udagger, R)
        LA1_All[i] = np.dot(Udagger, np.dot(LA1, U))

        NKS_All[i] = np.dot(U.T, np.dot(NKS, Udagger.T))
    write_RLNEF("WH_RLNEF.INP", R_All, LA1_All, NKS_All, E_Fermi)


def trans_solution_n2n():
    '''
    Apply unitary trnasformation to the solution according to WH_N2N.INP.
    '''
    U_All = read_TRANS("WH_N2N.INP")
    get_gauge_solution(U_All, "WH_RLNEF.OUT_ORIG")


def trans_WH_Hs_n2n():
    '''
    Apply unitary trnasformation to the WH_HS.INP_ORIG according to WH_N2N.INP.
    '''
    U_All = read_TRANS("WH_N2N.INP")
    HsAll = read_Hs("WH_HS.INP_ORIG")
    for i, U in enumerate(U_All):
        Udagger = np.conj(U.T)
        for iHs in xrange(len(HsAll[i])):
            HsAll[i][iHs] = np.dot(Udagger, np.dot(HsAll[i][iHs], U))
    write_Hs("WH_HS_NEW.INP", HsAll)

if __name__ == "__main__":
    '''
    Test run.
    '''
#    from matrix_util import get_random_unitary_matrix_spin_deg
#    N1 = 2
#    U_All = []
#    U_All.append(get_random_unitary_matrix_spin_deg(N1, spin_fast=False))
#    get_gauge_solution(U_All)
    get_gauge_solution_nks(fname="WH_RLNEF.OUT")
