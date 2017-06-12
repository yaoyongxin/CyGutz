Installation of CyGutz with WIEN2k and model interfaces
=======================================================

Prerequisites
-------------

CyGutz consists of programs, executables, and scripts, written in Fortran90, c (c++) and Python2.7. Before you start the installation, you must make sure that the following packages are installed in your system.

* Linux
* WIEN2k (http://www.wien2k.at/)
* Fortran, C, CXX compiler and blas/lapack library. The following have been tested.

  * ifort, icc, icpc,  and mkl (recommended)
  * gcc, g++, gfortran >= 4.8 and openblas. 

* MPI 
  
  * Intel MPI, (mpiifort, mpiicc, mpiicpc, mpirun, etc. Check https://software.intel.com/en-us/qualify-for-free-software)
  * open MPI (mpif90, mpicc, mpicxx, mpirun, etc. Check https://www.open-mpi.org/)

* HDF5 library (sequential version, compiled using the **same Fortran90 compiler** with flags **enable-fortran and enable-fortran2003**. Check for details at https://www.hdfgroup.org/HDF5/release/obtain5.html)
* python = 2.7
* numpy >= 1.5.0
* scipy >= 0.9.0
* h5py  >= 2.5.0
* matplotlib >= 1.4.3
* ase >= 3.3.1 (https://wiki.fysik.dtu.dk/ase/download.html)
* pyspglib >= 1.7.2 (http://spglib.sourceforge.net/python-spglib.html#python-spglib)
* pymatgen >= 3.3.5 (http://pymatgen.org/#guided-install)

A convenient way to install python and the related packages is Anaconda (https://www.continuum.io/downloads).

Build and install
-----------------

The compilation ans installation consist of three simple steps: 

* Download CyGutz here at 
  https://codeload.github.com/yaoyongxin/CyGutz/zip/0-1-alpha.
  unzip it and cd to that directory (``CyGutz-0-1-alpha``). 
* Check and make sure the path to compilers, their options, and libraries are correctly set by typing::

  $ vi ./Makefile.in

or try to find the suitable ``Makefile.in`` in directory ``template_makefile.in``.

* Run the following command at the top directory to install it::

    $ make clean  # Optional, clean any previous installation.
    $ ./install.sh 
    $ source ~/.bashrc  # Get the environment variable ${WIEN_GUTZ_ROOT}

Explanations on the scripts and Makefile.in
-------------------------------------------

* By default, ``${WIEN_GUTZ_ROOT}=${HOME}/WIEN_GUTZ/bin`` will be set in your .bashrc file.
* Key executables ``CyGutz`` and ``dmft, dmft2`` for the WIEN2k interface are located at ``${WIEN_GUTZ_ROOT}``.
* Useful initialization and analysis scripts are located at ``${WIEN_GUTZ_ROOT}/tools``.
* We use sequential version of mkl, so set ``WLIBS = -mkl=sequential``
* It is a must to set correctly the base path of the hdf5 library. For instance, ``HDF5_BASE = /home/ykent/OPT/lib/hdf5-1.8.15-intel/``. 
