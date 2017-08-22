DFT+G batch calculations
------------------------

Following the previous tutorial on DFT batch calculations,
Here we present how to perform DFT+G batch calculations.


step a: initialize one calculation at the minimal volume point
==============================================================

Go to the folder V0 and initialize for DFT+G calculation. Type::

    $ cd V0/case
    $ ${WIEN_GUTZ_ROOT2}/init_ga.py

The response to the questions is listed below:

* Do you want to BREAK SPIN-SYMMETRY: n 
  (keep spin-symmetry for paramagnetic phase)
* Do you want to COMPLETELY break orbital-symmetry: n
  (there are still site symmetries which reduces complexility)
* Do you want to take into account the SPIN-ORBIT interaction: y
  (as we have already added SOI in the previous tutorial)
* Do you want to take into account the CRYSTAL FIELD effect: y
  (use site symmetry)
* Please select the method to parametrize Coulomb U-matrix: 1
  (Slater-Condo parametrization is the most widely used option
  and also a very good approximation for the real screened 
  Coulomb matrix.)
* Please select method for U-interaction double counting: 12
  (most stable option for Fully-localized-limit DC.)
* Symmetrically-equivalent atom indices: y
* atom 0 Ce

  * Is this atom correlated: y
  * Enter correlated shells: f
    (f-shell is selected as correlated orbitals)
  * Please provide interaction parameters U,J: 6.0 0.7
    (adjustable, but they are usually kept fixed.)
  * Please provide initial guess of the number of localized f-electrons: 0.5
    (estimate according to the valence)
* atom 1 O

  * Is this atom correlated: n
    (not correlated for usual s(p)-bands.)

* Please select the method to solve G-RISB equations: 0
  (recommended Modified Powell hybrid method)
* Please select the method to solve embedding Hamiltonian: -1
  (valence truncation-based exact diagonalization)

If there are mistakes poping up, please consult an expert.
Copy the main CyGutz input files to the template directory::

    $ cp init_ga.slog ginit.h5 GPARAM.h5 case.indmfl ../../template/


step b: batch initialize DFT+G jobs
===================================

Type the following command to automitically initialize 
the series of DFT+G jobs::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py batch_run_ga

The script simply copies the previously generated DFT+G initialization files
(ginit.h5, GPARAM.h5 and case.indmfl) to each job directory.


step c: run a series of DFT+G calculations
==========================================

Use the following command to directly run a series DFT+G calculations.
(In case that it is preferable to submit these jobs to the clusters,
one has to write own machine-dependent script.) Type::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py batch_run_ga -p 4

It will use 4 processors to run 4 DFT+G jobs each time
until all the jobs are done.
Depending on machines, the whole jobs may take several hours or one day. 


step d: save the DFT+G calculations
===================================

Save the DFT+G results by typing::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py batch_gsave_u6.0j0.7

It will loop over all the job directories and 
save the main results to a subfolder named ''u6.0j0.7''.
Note here we use the u/j values in the DFT+G calculations
to name the subfolder.


step d: energy-volume curve from DFT+G
======================================

We can easily check the energy vs volume curve by typing::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py ev_u6.0j0.7

The figure is plotted in a pdf file ''ev_u6.0j0.7.pdf''. 
The numerical data are also stored in metadata file ''results.h5''.
Get the pressure-volume curve by typing::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py eos_fit_u6.0j0.7

The energy-volume curve and pressure-volume curve from fitting 
to the Murnaghan equation of state are saved as ''u6.0j0.7_evfit.pdf'' 
and ''u6.0j0.7_pvfit.pdf'', 
the numerical results are also stored in the metadata file ''results.h5''.
By typing::

    $ h5ls -r results.h5

you will see the new data::

    /u6j0.7                  Group
    /u6j0.7/eosfit           Group
    /u6j0.7/eosfit/b0        Dataset {SCALAR}
    /u6j0.7/eosfit/bp        Dataset {SCALAR}
    /u6j0.7/eosfit/e0        Dataset {SCALAR}
    /u6j0.7/eosfit/e_list    Dataset {129}
    /u6j0.7/eosfit/p_list    Dataset {129}
    /u6j0.7/eosfit/v0        Dataset {SCALAR}
    /u6j0.7/eosfit/v_list    Dataset {129}
    /u6j0.7/etot_list        Dataset {13}
 
