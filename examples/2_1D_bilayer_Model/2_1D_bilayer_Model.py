import numpy
import pickle
from scipy.linalg import block_diag

import glob
import shutil
import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools/Model/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../../tools/Gutzwiller"))

# import ase module
from ase.dft import kpoints
# import tbBase
from tbASE import *
from tbGutz import *

def tutorial_2_1D_bilayer_Model():
    '''
    A 1-dim bilayer model interface for CyGutz.
    To run this example

    1) Create and cd your working directory.
    2) type::

         $ python ${WIEN_GUTZ_ROOT}/examples/2_1D_bilayer_Model.py

       and answer a short list of questions::

         LHUB = 1: Slater-Condo parametrization.
         LHUB = 0: U_{i,j,k,l} (NO SPIN INDEX) = \int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} will be provided by file V2AO.INP
         LHUB =-1: U_{i,j,k,l} (INCLUDING SPIN INDEX) = \int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} will be provided by file V2H.INP
         Please select LHUB:  Pick one from [-1, 0, 1]...0

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

         INFORMATION FOR s ELECTRONS OF Cy :
         Please provide N1,N2 defining valence range [N1,N2] separated by a space ( 0 < N1 < N2 < 4 ): 0 4

       Since we choose LHUB = 0, we also need create the file :download:`V2AO.INP <../../examples/2_1D_bilayer_Model/V2AO.INP>`.
    3) type::

         $ ${WIEN_GUTZ_ROOT}/CyGutz

       The text format result is again in ``GUTZ.LOG``,
       The key results in metadata format  are also stored in ``glog.h5`` for easy post-processing using h5py.

    4) You may further experiment it with different U, J, electron filling parameters.
    '''
    a = AtomsTB("N",[(0,0,0)],cell=(1,1,1))
    a.set_orbitals_spindeg(orbitals=[("s","p")])
    aTB = TB(a)
    aTB.set_hop([(( 1,0,0),0,0,-1),
                 ((-1,0,0),0,0,-1),
                 (( 1,0,0),1,1,-1),
                 ((-1,0,0),1,1,-1),
                 ((0,0,0),0,1,0.25),
                 ((0,0,0),1,0,0.25)])
    aTB = aTB.add_spindegeneracy()
    kps_size=(100,1,1)
    kps=kpoints.monkhorst_pack(kps_size)
    # unit cell
    gTB=tbGutz(aTB.Atoms,aTB.Hr)
    # electron filling can be changed here
    gTB.output_CyGutz(kps, num_electrons = 2.0)

    ### write WH_HS.INP
    sigma_list = []; U_list = []
    norbitals = gTB.Atoms.nspinorbitals/2
    for i in range(1):
        sig_half = (numpy.arange(norbitals*norbitals)+1).reshape(norbitals,norbitals)
        sigma_list.append(block_diag(sig_half, sig_half))  # assuming Sz conservation and orbital index as fast index
        U_list.append(numpy.identity(norbitals*2, dtype = complex))

    from gl_inp import set_wh_hs, set_gl_inp
    set_wh_hs(sigma_list, U_list)
    #### write GL.INP
    spin_pol = 'n'
    SOC = ['y'] # Necessary to keep the orbital index (instead of the spin) as the fast index.
    CF = ['y']
    NTYPE = 1; NIONS = 1; ITYPE_list = [1]
    NI0_list = [1]; NIMAP_list = [1]; corr_atom_type = ["X"]
    type_1atom = [0]; df_list = ["g"]; dim_list = [norbitals*2]
    log = open("init_ga_a.slog", 'w')
    set_gl_inp(spin_pol, SOC, CF, NTYPE, NIONS, ITYPE_list, NI0_list, NIMAP_list, corr_atom_type, type_1atom, df_list, dim_list, log)
    from h5save_init_ga_info import h5write_init_ga_info
    num_sigma = len(sigma_list)
    h5write_init_ga_info(num_sigma=num_sigma, sigma_list=sigma_list)
    log.close()

if __name__=="__main__":
    tutorial_2_1D_bilayer_Model()
