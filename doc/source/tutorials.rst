Tutorial
========

Single Band Hubbard Model
-------------------------
.. automodule:: 1_1D_1_band_Hubbard_Model
   :members:


Bilayer Model
---------------------
.. automodule:: 2_1D_bilayer_Model
   :members:


:math:`\alpha`-Ce with spin-orbit interaction
---------------------------------------------
1) Finish a self-consistent lapw + so run for :math:`\alpha`-Ce using WIEN2k. 
   Here is the structure file :download:`Ce.struct <../../examples/3_Ce_soc/Ce.struct>`.
2) type::

     $ ~/WIEN_GUTZ/bin/tools/Gutzwiller/init_ga.py 

   Read the list of questions and pick your options::

     User inputs to initialize the ga job.
     Do you want to BREAK SPIN-SYMMETRY? Pick one from [y, n]...n
     Do you want to COMPLETELY BREAK ORBITAL-SYMMETRY? Pick one from [y, n]...n
     atom  0   Ce
     Is this atom correlated? Pick one from [y, n]...y
     Enter correlated shells? Pick one or combinations separated by blank space from [s, p, d, f]...f
     Do you want to take into account the SPIN-ORBIT interaction? Pick one from [y, n]...y
     Do you want to take into account the CRYSTAL FIELD effect? Pick one from [y, n]...n
    
     LHUB = 1: Slater-Condo parametrization.
     LHUB = 0: U_{i,j,k,l} (NO SPIN INDEX) = int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} 
               will be provided by file V2AO.INP
     LHUB =-1: U_{i,j,k,l} (INCLUDING SPIN INDEX) = int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j) phi(r_k) phi(r'_l)}} 
               will be provided by file V2H.INP
     Please select LHUB:  Pick one from [-1, 0, 1]...1
    
     LDC = 12:  Recommended. Standard double counting (updating Vdc at each charge iteration, initial n0 to be provided.) 
     LDC =  2:  Fix double counting potential (keep same Vdc/n0 at each charge iteration, n0 to be provided.) 
     LDC =  1:  Standard double counting potential (n0 self-consistently determined.) 
     LDC =  0:  No double counting. 
     Please select LDC:  Pick one from [0, 1, 2, 12]...12
    
     LCLUSTER = 0: Single-atom impurity.
     LCLUSTER = 1: Multi-atom (cluster) impurity.
     Please select LCLUSTER:  Pick one from [0, 1]...0
    
     Solution embedding system:
     LEIGV = 0:  Choose automatically solver depending on the size of the problem (DEFAULT) 
             1:  Exact diagonalization (ZHEEV) in LAPACK
             2:  Lanczos (zhdrv1) in ARPACK 
             3:  Exact diagonalization (ZHEEVX, selective lowest two eigen-vectors) in LAPACK
             5:  PRIMEE (Iterative MultiMethod Eigensolver)
     Please select LEIGV:  Pick one from [0, 1, 2, 3, 5]...5
     INFORMATION FOR f ELECTRONS OF Ce :
     Please provide interaction parameters U,J separated by a space: 6.0 0.7
     Please provide N1,N2 defining valence range [N1,N2] separated by a space ( 0 < N1 < N2 < 14 ): 0 4
     Please provide guess n0 for valence (0 < n0 < 4): 1.0
     Please run ga_init_dmft.py with parameters given in  init_ga.slog

   Thus it finishes the setup of all the necessary files for CyGutz.
3) To use the w2k+dmft priofector, we will type::

     $ ~/WIEN_GUTZ/bin/ga_init_dmft.py

   Read the list of questions and use the answers provided in ``init_ga.slog``.
4) Now we are ready to run the w2k+Gutzwiler/Slave-Boson self-cponsistently by typing::

     $ ~/WIEN_GUTZ/bin/ga_run_dmft.py

   Like WIEN2k jobs, the total energy is stored in ``case.scf`` file like :download:`Ce.scf <../../examples/3_Ce_soc/Ce.scf>`. 
   The main CyGutz output (last electron density step) is located in :download:`GL_LOG.OUT <../../examples/3_Ce_soc/GL_LOG.OUT>`
   and :download:`glog.h5 <../../examples/3_Ce_soc/glog.h5>`.


Band structure calculation for UO2
---------------------------------------------
1) Finish a self-consistent LDA+G calculation. (Make a copy of the job by using e.g. `~/BitBucket/cygutz_release/tools/Gutzwiller/save_lapwg`)
2) Set `LSCF = 2` and add a new line `LEL0 = 1` in GL.INP
3) type::

     $ cp WH_EL0.OUT WH_EL0.INP && cp EFLDA.OUT EFLDA.INP

4) Prepare the `UO2_10.klist_band` file for the high-symmetry k-path of the bands. (Check Wien2k `user_guide`)
5) type::

     $ mv UO2_10.klist_band UO2_10.klist
     $ $(WIEN_GUTZ_ROOT)/ga_run_dmft.py -e CyGutz

6) Use :download:`plot_bands.py <../../examples/4_UO2_bands/plot_bands.py>` to dump the data to `BANDS.dat` and use xmgrace to plot it. 


DOS plot for UO2
---------------------------------------------
1) Finish a self-consistent LDA+G calculation.
2) Use :download:`dos.py <../../examples/5_UO2_dos/dos.py>` to plot the (partial) dos.


Crystal field for CeB6
---------------------------------------------
1) Finish a self-consistent LDA+G calculation.
2) Use :download:`dos.py <../../examples/6_CeB6_Crystal_Field/analysis_cf.py>` to analyse the crystal field.

