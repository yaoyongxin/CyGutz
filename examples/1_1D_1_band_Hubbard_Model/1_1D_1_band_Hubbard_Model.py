import numpy
import pickle
from scipy.linalg import block_diag

import glob
import shutil
import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools/Model/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools/Gutzwiller"))

#import ase module
from ase.dft import kpoints
#import tbBase
from tbASE import *
from tbGutz import *

def tutorial_1_1D_1_band_Hubbard_Model():
    '''
    A 1-dim single-band Hubbard model interface for CyGutz.
    To run this example

    1) Create and cd to your working directory
    2) type::

         $ python ${WIEN_GUTZ_ROOT}/examples/1_1D_1_band_Hubbard_Model.py

       and answer a short list of questions::

         LHUB = 1: Slater-Condo parametrization.
         LHUB = 2: Kanamori parametrization.
         LHUB = 0: U_{i,j,k,l} (NO SPIN INDEX) = \int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} will be provided by file V2AO.INP
         LHUB =-1: U_{i,j,k,l} (INCLUDING SPIN INDEX) = \int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} will be provided by file V2H.INP
         Please select LHUB:  Pick one from [-1, 0, 1]...1

         LDC = 12:  Recommended. Standard double counting (updating Vdc at each charge iteration, initial n0 to be provided.)
         LDC =  2:  Fix double counting potential (keep same Vdc/n0 at each charge iteration, n0 to be provided.)
         LDC =  1:  Standard double counting potential (n0 self-consistently determined.)
         LDC =  0:  No double counting.
         Please select LDC:  Pick one from [0, 1, 2, 12]...0

         LCLUSTER = 0: Single-atom impurity.
         LCLUSTER = 1: Multi-atom (cluster) impurity.
         Please select LCLUSTER:  Pick one from [0, 1]...0

         Solution embedding system:
         LEIGV = 0:  Choose automatically solver depending on the size of the problem (DEFAULT)
                 1:  Exact diagonalization (ZHEEV) in LAPACK
                 2:  Lanczos (zhdrv1) in ARPACK
                 3:  Exact diagonalization (ZHEEVX, selective lowest two eigen-vectors) in LAPACK
                 5:  PRIMEE (Iterative MultiMethod Eigensolver)
         Please select LEIGV:  Pick one from [0, 1, 2, 3, 5]...1

         INFORMATION FOR g ELECTRONS OF X :
         Please provide interaction parameters U,J separated by a space: 7.0 0.0
         Please provide N1,N2 defining valence range [N1,N2] separated by a space ( 0 < N1 < N2 < 2 ): 0 2
    3) type::

       $ ${WIEN_GUTZ_ROOT}/CyGutz

       The text format result is in :download:`GUTZ.LOG <../../examples/1_1D_1_band_Hubbard_Model/reference/GUTZ.LOG>`.
       The key results in metadata format are also stored in ``glog.h5`` for easy post-processing using h5py.
       One can easily view the content by using the hdf5 tools::

         $ h5ls -r glog.h5

       which shows a list of metadata::

         /                        Group
         /BND_EK                  Dataset {1, 100, 2}   # band energy E[i_spin, i_kpoint, i_band]
         /BND_NE                  Dataset {100, 3}      # band_index[i_kpoint,0/1/2]: total number of band/ botom/top band index at each k-point to define the sets of bands to expand the local d/f projector.
         /BND_VK                  Dataset {1, 100, 1, 2, 2} # <basis orbitals | Gutzwiller band > [i_spin, i_kpoint, i_symmtry_op, i_orbital, i_band]
         /Distinct_impurity_indices Dataset {1}         # Indices of distinct impurities.
         /E_DC2_U                 Dataset {1}           # Total double counting energy of the Hubbard interactions, summing over all the impurity sites in the unit cell.
         /E_FERMI                 Dataset {1}           # Gutzwiller Fermi energy
         /E_GAMMA                 Dataset {1}           # Total onsite energy including one- and two-body interactions.
         /E_POT2_U                Dataset {1}           # Total onsite Hubbard interaction energy.
         /E_TB_TOT                Dataset {1}           # Total energy (Band energy - double counting + onsite energy.)
         /GA_MAX_ERR              Dataset {1}           # The final error of the Gutzwiller solution.
         /Impurity_1              Group                 # Information about the 1st impurity (one-based)
         /Impurity_1/ComplexHarmonicsToCurrentBasis Dataset {2, 2} # Transposition(Transformation matrix: <Complex Harmonics | Current local orbital basis in the calculation>.) Not relevant for models.
         /Impurity_1/GA_D         Dataset {2, 2}        # Transposition(D)[alpha_physical basis index, a_quasi_particle ]. D as defined in PHYS. REV. X 5, 011008 (2015).
         /Impurity_1/GA_La        Dataset {2, 2}        # Transposition(\lambda), \lambda as defined PHYS. REV. X 5, 011008 (2015).
         /Impurity_1/GA_Lc        Dataset {2, 2}        # Transposition(\lambda^{c}), \lambda^{c} as defined PHYS. REV. X 5, 011008 (2015).
         /Impurity_1/GA_NC_PHY    Dataset {2, 2}        # Transposition(physical one-body density matrix)
         /Impurity_1/GA_NC_VAR    Dataset {2, 2}        # Transposition(variational one-body density matrix)
         /Impurity_1/GA_NKS       Dataset {2, 2}        # Transposition(quasi-particle density matrix)
         /Impurity_1/GA_NKS_BARE  Dataset {2, 2}        # Transposition(bare one-body density matrix, for debug only)
         /Impurity_1/GA_R         Dataset {2, 2}        # Transposition(R)[alpha_physical basis index, a_quasi_particle ]. R as defined in PHYS. REV. X 5, 011008 (2015).
         /Impurity_1/H.base       Dataset {1}           # Local Hamiltonian including one- and two-body terms in spare matrix (csr) format.
         /Impurity_1/H.data       Dataset {1}
         /Impurity_1/H.indices    Dataset {1}
         /Impurity_1/H.indptr     Dataset {5}
         /Impurity_1/H.ncol       Dataset {1}
         /Impurity_1/H.nrow       Dataset {1}
         /Impurity_1/One_body_local Dataset {2, 2}      # The bare one-body part of the local Hamiltonian.
         /Impurity_1/PreferedBasisToCurrentBasis Dataset {2, 2} # Transposition(Transformation matrix: < prefered basis | Current local orbital basis in the calculation>). Prefered basis can be wither complex Harmonics or relativistic Harmonics, depending on the choice when executing "ga_init_dmft.py". Not relevant for models.
         /Impurity_1/RHO.base     Dataset {1}           # Local reduced many-body density matrix in the spare matrix (csr) format.
         /Impurity_1/RHO.data     Dataset {6}
         /Impurity_1/RHO.indices  Dataset {6}
         /Impurity_1/RHO.indptr   Dataset {5}
         /Impurity_1/RHO.ncol     Dataset {1}
         /Impurity_1/RHO.nrow     Dataset {1}
         /Impurity_1/Sec_ID       Dataset {4}           # N-block of the local Fock states, e.g., for calculations with {f0,f1,f2}, the states of f0 runs from Sec_ID[0] to Sec_ID[1]-1
         /Impurity_1/Sec_VAL      Dataset {3}           # with Sec_VAL[0]=0; ... fi from Sec_ID[i] to Sec_ID[i+1]-1 with Sec_VAL[i]=i.
         /Impurity_1/Two_body_local Dataset {2, 2, 2, 2}  # U_{i,j,k,l} (INCLUDING SPIN INDEX) = \int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}}
         /Impurity_1/annihi.op._1.base Dataset {1}      # c_{i_spin_orbital} in Fock state basis.
         /Impurity_1/annihi.op._1.data Dataset {2}
         /Impurity_1/annihi.op._1.indices Dataset {2}
         /Impurity_1/annihi.op._1.indptr Dataset {5}
         /Impurity_1/annihi.op._1.ncol Dataset {1}
         /Impurity_1/annihi.op._1.nrow Dataset {1}
         /Impurity_1/annihi.op._2.base Dataset {1}
         /Impurity_1/annihi.op._2.data Dataset {2}
         /Impurity_1/annihi.op._2.indices Dataset {2}
         /Impurity_1/annihi.op._2.indptr Dataset {5}
         /Impurity_1/annihi.op._2.ncol Dataset {1}
         /Impurity_1/annihi.op._2.nrow Dataset {1}
         /Impurity_1/phi.base     Dataset {1}           # \phi_{\Gamma, n} as defined in PHYS. REV. X 5, 011008 (2015) in the spare matrix (csr) format.
         /Impurity_1/phi.data     Dataset {6}
         /Impurity_1/phi.indices  Dataset {6}
         /Impurity_1/phi.indptr   Dataset {5}
         /Impurity_1/phi.ncol     Dataset {1}
         /Impurity_1/phi.nrow     Dataset {1}
         /Impurity_1/weight_nf_0  Dataset {1}           # weight of valence block 0 (f0)
         /Impurity_1/weight_nf_1  Dataset {1}
         /Impurity_1/weight_nf_2  Dataset {1}
         /KPT_WT                  Dataset {100}         # weight (i_kpoint)
    '''

    aTB=TB.gallery("Chain_nn").add_spindegeneracy()
    kps_size=(100,1,1)
    kps=kpoints.monkhorst_pack(kps_size)
    # a TB model on a Chain
    gTB=tbGutz(aTB.Atoms,aTB.Hr)
    # electron filling can be changed here
    gTB.output_CyGutz(kps, num_electrons = 1.0)

    # write WH_HS.INP
    sigma_list = []; U_list = []
    norbitals = gTB.Atoms.nspinorbitals/2
    for i in range(1):
        sig_half = (numpy.arange(norbitals*norbitals)+1).reshape(norbitals,norbitals)
        sigma_list.append(block_diag(sig_half, sig_half))  # assuming Sz conservation
        U_list.append(numpy.identity(norbitals*2, dtype = complex))

    from gl_inp import set_wh_hs, set_gl_inp
    set_wh_hs(sigma_list, U_list)
    # write GL.INP
    spin_pol = 'n'
    SOC = ['y']; CF = ['y']
    NTYPE = 1; NIONS = 1; ITYPE_list = [1]
    NI0_list = [1]; NIMAP_list = [1]; corr_atom_type = ["Cy"]
    type_1atom = [0]; df_list = ["s"]; dim_list = [norbitals*2]
    log = open("init_ga_a.slog", 'w')
    set_gl_inp(spin_pol, SOC, CF, NTYPE, NIONS, ITYPE_list, NI0_list, NIMAP_list, corr_atom_type, type_1atom, df_list, dim_list, log)
    # write ga_init_info.h5
    num_sigma = len(sigma_list)
    from h5save_init_ga_info import h5write_init_ga_info
    h5write_init_ga_info(num_sigma=num_sigma, sigma_list=sigma_list)
    log.close()

if __name__=="__main__":
    tutorial_1_1D_1_band_Hubbard_Model()
