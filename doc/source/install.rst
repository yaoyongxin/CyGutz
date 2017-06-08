Installation of CyGutz
======================

Prerequisites
-------------

CyGutz consists of programs, executables, and scripts, 
written in Fortran90, c (c++) and Python2.7. 
Before you start the installation, 
you must make sure that the following packages 
are installed in your system.

* Linux
* Fortran, C, CXX compiler and blas/lapack library. 
  The following have been tested.

  * ifort, icc, icpc,  and mkl (recommended)
  * gcc, g++, gfortran >= 4.8 and openblas. 

* MPI 
  
  * Intel MPI, (mpiifort, mpiicc, mpiicpc, mpirun, etc. 
    Check https://software.intel.com/en-us/qualify-for-free-software)
  * open MPI (mpif90, mpicc, mpicxx, mpirun, etc. 
    Check https://www.open-mpi.org/)

* HDF5 library (sequential version, 
  compiled using the **same Fortran90 compiler** 
  with flags **enable-fortran and enable-fortran2003**. 
  To configure, using line like
  ./configure --prefix=prefix CC=icc FC=ifort CXX=icpc --enable-fortran 
  --enable-fortran2003. 
  Check for details at https://www.hdfgroup.org/HDF5/release/obtain5.html)
* python = 2.7
* numpy >= 1.5.0
* scipy >= 0.9.0
* h5py  >= 2.5.0
* matplotlib >= 1.4.3
* ase >= 3.3.1 (https://wiki.fysik.dtu.dk/ase/install.html)
* pyspglib >= 1.7.2 
  (https://atztogo.github.io/spglib/python-spglib.html#python-spglib)
* pymatgen >= 3.3.5 
  (http://pymatgen.org/#guided-install)

A very convenient way to install python and the related packages is Anaconda (https://www.continuum.io/downloads).

Build and install
-----------------

The compilation ans installation consist of three simple steps: 

* Download CyGutz here at https://github.com/yaoyongxin/CyGutz, unzip it 
  and cd to that directory (``CyGutz-master``). 
  You may simply type::

  $ wget https://codeload.github.com/yaoyongxin/CyGutz/zip/master 
  $ unzip master
  $ cd CyGutz-master

* Choose a ``Makefile.in`` file in directory ``template_makefile.in``. 
  Check and make sure the path to compilers, their options, 
  and libraries are correctly. For example, you may type::

  $ cp template_makefile.in/Makefile.in.intelmpi_ifort_hal2011 Makefile.in
  $ vi ./Makefile.in

* Run the following command at the top directory to install it::

    $ # Optional, clean any previous installation.
    $ make clean && make clean_bin 
    $ ./install.sh 
    $ # Get the environment variable ${WIEN_GUTZ_ROOT2} 
    $ # for the first-time installation.
    $ source ~/.bashrc 

Explanations on the scripts and Makefile.in
-------------------------------------------

* By default, ``${WIEN_GUTZ_ROOT2}=${HOME}/WIEN_GUTZ/bin2`` will be set 
  in your .bashrc file.
* Key executables, such as ``CyGutz``, ``exe_spci``, and python scripts
  are located at ``${WIEN_GUTZ_ROOT2}``.
* In ``Makefile.in`` file:

  * Use sequential version of mkl, so set ``MKL_LIB = -mkl=sequential``
    in ``Makefile.in``.
  * The base path of the hdf5 library has to be set correctly.
    For instance, ``HDF5_BASE = /home/ykent/OPT/lib/hdf5-1.8.15-intel/``
