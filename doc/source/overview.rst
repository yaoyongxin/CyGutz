Overview of CyGutz with WIEN2k and model interfaces
===================================================

DFT-LDA
-------

The Density Functional Theory implementation within WIEN2k consists of several steps depicted in the figure (Check http://www.wien2k.at/) 

.. image:: /images/WIEN2kFC.png
   :alt: alternate text
   :scale: 100 %
   :align: center

The steps are:

* init_lapw:

  * NN: lists nearest neighbor distances and determine atomic sphere radii.
  * SGROUP: determines the spacegroup of the structure deﬁned in your case.struct ﬁle.
  * SYMMETRY: generates from a raw case.struct ﬁle the space group symmetry operations, determines the point group of the individual atomic sites, generates the LM expansion for the lattice harmonics and determines the local rotation matrices.
  * LSTART: generates free atomic densities and determines how the different orbitals are treated in the band structure calculations.
  * KGEN: generates a k-mesh in the irreducible part of the BZ.
  * DSTART: generates a starting density for the scf cycle by a superposition of atomic densities generated in LSTART.

* run_lapw:
  
  * LAPW0: (POTENTIAL) generates potential from density.
  * LAPW1: (BANDS) calculates valence bands (eigenvalues and eigenvectors).
  * LAPWSO: adds spin-orbit coupling.
  * LAPW2: (RHO) computes valence densities from eigenvectors.
  * LCORE: computes core states and densities.
  * MIXER: mixes input and output densities.

LDA + Gutzwiller-Slave-Boson
----------------------------

The LDA+G inserts two new steps (dmft1 and CyGutz) and replaces lapw2 with dmft2 step, following vlosely the DFT+DMFT code (Check http://hauleweb.rutgers.edu/). The process looks like that: 

.. image:: /images/WIENGAFC.png
   :alt: alternate text
   :scale: 100 %
   :align: center

The LDA+G steps (``ga_run_dmft.py``) are:

* LAPW0:  (POTENTIAL) generates potential from density.
* LAPW1:  (BANDS) calculates valence bands (eigenvalues and eigenvectors).
* LAPWSO: adds spin-orbit coupling.
* DMFT1:  calculates the local correlated orbitals (projectors) using the Kohn-Sham bands as the basis.
* CYGUTZ: Solve the generic Kohn-Sham-Hubbard Hamiltonian using the Gutzwiller Slave-Boson method.
* dmft2:  calculates LDA+G renormalized valence electron density.
* LCORE:  computes core states and densities.
* MIXER:   mixes input and output densities.

CyGutz
------

CyGutz is a stand-alone program to solve a generic Hubbard model (including local correlated orbitals and nonlocal orbitals) using Gutzwiller-Slave-Boson method (See http://journals.aps.org/prx/abstract/10.1103/PhysRevX.5.011008). The steps to reach self-consistent solutions are as follows:


.. image:: /images/CyGutzFC.png
   :alt: alternate text
   :scale: 100 %
   :align: center

CyGutz starts with the setup of a generic [Kohn-Sham]-Hubbard model based on the input bare band energies 
and the specified local correlated orbitals, 
which are expressed in terms of the bare band wave functions. 
Here the bare bands can be obtained from LDA calculations or tight-binding models. 
The solution of the generic Hubbard model is cast to a root-finding problem within the Gutzwiller-Slave-Boson method. 
Check http://journals.aps.org/prx/abstract/10.1103/PhysRevX.5.011008 for details.
