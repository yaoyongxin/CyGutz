DFT energy-volume curve
-----------------------

Since DFT+Gutzwiller method is very efficient, 
it is a routine to perform the calculations for a series of volume points.
It can not only get the typical energy-volume (E-V) curve 
or equivalent equation of state (EOS), 
which can be directly compared with experiments,
but also the possible accompanying phase transition, 
e.g., metal-insulator transition, low spin-high spin transition, etc.
As a first step, here we provide a detailed instruction
on how to perform such DFT calculations using wien2k.
We choose CeO\ :sub:`2` as an example, 
which is known to be paramegnetic experimentally.


step 0: get the cif structure file
==================================

A convenient and free place to locate the cif structures is 
materialsproject.org. 
After login, on the webpage of https://materialsproject.org/#apps/phasediagram,
select ''Phase Diagram'' to explore, 
and choose the binary systems of ''Ce'' and ''O''. 
In the convex hull diagram, you can locate 
the stable structure of CeO\ :sub:`2`.
Click the id (mp-20194), and you will reach the page of 
https://materialsproject.org/materials/mp-20194/. 
You will see that it is a cubic structure with m-3m point group. 
The lattice parameters from calculation and ICSD database (experiments) 
are also listed at right. 
The structure can be cross-checked at other online database, 
for instance, http://mio.asminternational.org/apd/.

In your working directory, create a subdirectory named ''template'' by typing::

    $ mkdir template
    $ cd template

The symmetrized cif file can now be downloaded 
at the right ''Final Structure'' section 
and stored in your ''template'' directory.
Rename the cif file to be ''case.cif'' by typing::

    $ rename *cif case.cif


step 1: initialize and run the minimal volume point
===================================================

Go back to the working directory and create a minimal volume point,
with a default volume fraction 0.7 of the reference structure volume 
(V0) obtained previously. 
One can use the predifined script by typing::

    $ cd ..
    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py Vmin

Following the message printed on the screen, go to directory ''Ref_Vmin/case''
and manually initialize and run the DFT calculation. 
This step is to make sure that the parameters can be properly set 
and DFT calculation can be well converged. Typing::

    $ cd Ref_Vmin/case
    $ init_lapw

Here are the answers for the questions: 

* Enter reduction in %: 1
* Use old or new scheme (o/N): N
* Do you want to accept these radii: a
* specify nn-bondlength factor: 2
* continue with sgroup: c
* use/edit case.struct_sgroup: e
* Do you want to use the new struct file: y
* specify nn-bondlength factor: 2
* continue with sgroup: c
* continue with symmetry: c
* continue with lstart: c
* Eventually specify switches for instgen_lapw: -up
* SELECT XCPOT: 5
* SELECT ENERGY to separate core and valence states: -6.0
* continue with kgen: c
* NUMBER OF K-POINTS IN WHOLE CELL: 2000
* continue with dstart: y
* do you want to perform a spinpolarized calculation: n

Note that the some of the above parameters are for testing only.
There are default parameters setup by the ''stepin_wien_gutz.py'' script.
Test the convergence of the DFT calculations by typing::

    $ run_lapw

In a few minutes or so, you will see the convergence of the lapw calculation.
If not, please consult an expert.
Now we have a working DFT example of the minimal volume point 
stored at the directory of ''Ref_Vmin/case''.


step 2: update struct file in the template directory
====================================================

Go back to the top of the working directory and update the structure file 
in template directory, including muffin-tin radius and symmetry. Typing::

    $ cd ../../
    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py update_case_ref


step 2: create a list of volume samples
=======================================

Type the following command to create a list of jobs 
of uniform volume sampling::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py Vlist

By default, the volume fraction with respect to V0 goes from 0.7 to 1.3 
at step size of 0.05. 


step 3: batch initialize the jobs
=================================

Type the following command to automitically initialize 
the series of jobs just created::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py batch_init_lapw

It will take a few minutes. 
Log file ''binit_lapw.log'' is created in the working directory,
which records the main output, including possible errors or warnings.


step 4: run a series of lapw calculations
=========================================

Use the following command to directory run the series DFT calculations, 
since they are usually cheap. 
Otherwise, one should use job script file to submit these jobs 
to the clusters. Type::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py run_lapw

If there are several cores available, e.g., 4,
use the following command to save time::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py run_lapw -p 4

This can take up to a few hour.


step 5: save the lapw calculations
==================================

It is a good idea to save the main calculation results, 
which can be used for analysis or starting point for new calculations.
Type::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py batch_save_lapw

It saves the main results to s subfolder named ''lapw''.


step 6: energy-volume curve from lapw 
=====================================

We can easily check the energy vs volume curve by typing::

    $ ${WIEN_GUTZ_ROOT2}/stepin_wien_gutz.py ev_lapw

The figure is plotted in a pdf file ''ev_lapw.pdf''.
The numerical data are also stored in metadata file ''results.h5''.


step 7: adding spin-orbit interaction
=====================================


