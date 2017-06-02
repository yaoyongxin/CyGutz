'''
Transformation of the GA/Slave-boson solutions.
'''
from fileio import *
import numpy as np

def get_gauge_solution(U_All, fname):
  '''
  Generate another equivalent solution due to gauge invariance,
  given the old solution in fname. Note Hs might be a concern
  if symmetry is included.
  '''
  R_All, LA1_All, NKS_All, E_Fermi = read_RLNEF(fname)
  for i in xrange(len(U_All)):
    Udagger = np.conj(U_All[i].T)
    R_All[i] = np.dot(Udagger, np.dot(R_All[i], U_All[i]))
    LA1_All[i] = np.dot(Udagger, np.dot(LA1_All[i], U_All[i]))
    NKS_All[i] = np.dot(U_All[i].T, np.dot(NKS_All[i], Udagger.T))
  write_RLNEF("WH_RLNEF_NEW.OUT", R_All, LA1_All, NKS_All, E_Fermi)

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
#  from matrix_util import get_random_unitary_matrix_spin_deg
#  N1 = 2
#  U_All = []
#  U_All.append(get_random_unitary_matrix_spin_deg(N1))
#  get_gauge_solution(U_All, "WH_RLNEF.INP")
  trans_solution_n2n()
#  trans_WH_Hs_n2n()
