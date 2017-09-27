:math:`\alpha`-Ce with spin-orbit interaction
---------------------------------------------
In this example, you will learn how to perform a 
DFT+G calculation using Wien2k plus CyGutz packages,
including some typical post-analyses.

1) Finish a self-consistent lapw + spin-orbit calculation 
    for :math:`\alpha`-Ce using WIEN2k.
    Please follow the following steps:

    (a) create a directory with name 'Ce', change to that directory,
        and download the structure file 
        :download:`Ce.struct <_files/Ce.struct>`.

    (b) type the following command to initialize the wien2k job::

            $ ${WIENROOT}/init_lapw

        use the answers below for the questions:

        * Enter reduction in %: enter (choose default)
        * Use old or new scheme (o/N): enter
        * Do you want to accept these radii... (a/d/r): enter
        * nn  (14:13:54)  specify nn-bondlength factor: enter
        * continue with sgroup or edit the Ce.struct file (c/e): enter
        * continue with symmetry (old case.struct) ...: enter
        * continue with lstart or edit the Ce.struct_st file (c/e/x): enter
        * Eventually specify switches for instgen_lapw: enter
        * lstart ...  SELECT XCPOT: 5
        * SELECT ENERGY to separate core and valence states: -6.0
        * continue with kgen or edit the Ce.inst file ... (c/e): enter  
        * in 2nd line of Ce.in1_st file, change 7.00 (R-MT*K-MAX) to 9.00
        * kgen NUMBER OF K-POINTS IN WHOLE CELL: 5000
        * continue with dstart or execute kgen again or exit (c/e/x): enter
        * do you want to perform a spinpolarized calculation ? (n/y): n

    (c) type the following command to run lapw job::
        
            $ ${WIENROOT}/run_lapw

        check the last line of `Ce.scf` file, which should give a LDA total 
        energy close to -17717.65311 Ryd.

    (d) type the following command to add spin-orbit interaction::

            $ ${WIENROOT}/initso_lapw

        use the answers below for the questions:

        * Please select the direction of the moment ( h k l ): enter
        * for which you would NOT like to add SO interaction: enter
        * Please enter EMAX(default 5.0 Ryd): 7.5
        * Add RLO for NONE, ALL, CHOOSE elements? (N/a/c) : enter
        * in Ce.inso file, change the 3rd line by replacing 1.5 with 4.5 (Emax)
        * Do you have a spinpolarized case ...: enter

    (e) type the following command to run lda calculation with spin-orbit::

            $ ${WIENROOT}/run_lapw -so

        check the last line of `Ce.scf` file, which should give a LDA+so 
        total energy close to -17717.67370 Ryd.
        Use the following command to save the result::

            $ ${WIENROOT}/save_lapw -d lapwso

2) We are ready to initialize the DFT+G calculation. 
    Type::

        $ ${WIEN_GUTZ_ROOT2}/init_ga.py 

    Use the following answers or the questions:
        
    * Do you want to BREAK SPIN-SYMMETRY: n
    * Do you want to COMPLETELY break orbital-symmetry: n
    * Do you want to take into account the SPIN-ORBIT interaction: y
    * Do you want to take into account the CRYSTAL FIELD effect: n
    * Please select the method to parametrize Coulomb U-matrix: 1
    * Please select method for U-interaction double counting: 12
    * Symmetrically-equivalent atom indices ...: y
    * atom 0 Ce:

      * Is this atom correlated: y
      * Enter correlated shells: f
      * Please provide interaction parameters U,J ... (eV): 6.0 0.7
      * Please provide initial guess ... of localized f-electrons: 1.0

    * Please select the method to solve G-RISB equations: 0
    * Please select the method to solve embedding Hamiltonian: -1
        
    Thus it finishes the setup of all the necessary files for CyGutz.
    The main output file of the initialization is stored 
    in `init_ga.slog` file. 
    The user-provided entries is stored in HDF5 file `ginit.h5` 
    for the convenience of re-initialization.

3) To run the DFT+G calculation, 
    type::

        $ ${WIEN_GUTZ_ROOT2}/run_ga.py

    Check the last line of `Ce.scf` file, which should give a LDA+so+G
    total energy close to -17717.60065 Ryd.
    The main output message of CyGutz calculation is printed 
    in `GUTZ.LOG` file in text format. 
    The DFT+G calculation can be saved using the following commands::

        $ ${WIEN_GUTZ_ROOT2}/save_ldag -d lapwg

4) Let us learn how to plot the density of states of Ce 
    after the self-consistent DFT+G calculation.  

    (a) get the simple script 
        :download:`plot_dos_ce.py<./_files/plot_dos_ce.py>`

    (b) type::

            $ python ./plot_dos_ce.py

        you will get dos like

        .. image:: _images/ce_dos.png
           :alt: alternate text
           :scale: 100 %
           :align: center

    The above scirpt calls two predefined functions, 
    which serves a s a template to be adapted by users
    for specific purposes.

    .. autofunction:: pyglib.estructure.dos.h5get_dos
    .. autofunction:: pyglib.estructure.dos.plot_dos_tf

5) Another important analysis is the eigen-values of 
    the local reduced many-body density matrix 
    using the exponential form :math:`\rho=e^{-F}`, 
    as shown in the figure below:

    .. image:: _images/ce_hist_jj.png
       :alt: alternate text
       :scale: 100 %
       :align: center

    Follow the steps below to get the analysis done for the first impurity.

    (a) type::

            $ ${WIEN_GUTZ_ROOT2}/exe_spci_j2_mott_analysis 1

    (b) Use :download:`plot_hist_ce.py <_files/plot_hist_ce.py>` 
        to plot by typing::

            $ python ./plot_hist_ce.py

    The above scirpt calls two predefined functions, 
    which serves a s a template to be adapted by users
    for specific purposes.

    .. autofunction:: pyglib.mbody.multiplets_analysis_soc.calc_save_atomic_states
    .. autofunction:: pyglib.mbody.multiplets_analysis_soc.plot_atomic_states

6) To calculate the bands structure along selected k-path, 
    follow the steps below:

    (a) Prepare the Ce.klist_band file for the high-symmetry k-path 
        of the primitive Brillouin Zone. 
        The SRC_templates directory of Wien2k has some examples.
        For instance, we can use `fcc.klist` file.
        Type the command to get the file::

            $ cp ${WIENROOT}/SRC_templates/fcc.klist Ce.klist_band

    (b) Type the following command 
        to calculate the band structure::

            $ ${WIEN_GUTZ_ROOT2}/run_ga.py -band

    (c) Download the simple script 
        :download:`plot_band_ce.py<_files/plot_band_ce.py>`
        and type the following command to plot the band structure::

            $ python ./plot_band_ce.py

        You will see the band structure like the following

        .. image:: _images/ce_bands.png
          :alt: alternate text
          :scale: 100 %
          :align: center
      
        The above scirpt calls a predefined function,
        which serves a s a template to be adapted by users
        for specific purposes.

        .. autofunction:: pyglib.estructure.bandstruct.plot_band_sturture
