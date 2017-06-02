#!/usr/bin/evn python
#  A TB model interface for Gutz
import numpy
import pickle

import glob
import shutil
import os

## import ase module
from ase.dft import kpoints

# import tbBase
from tbASE import *

# import CyGutz interface
import gl_interface


class tbGutz(TB):
    """
    TB model with routines for Gutzwiller.

    Attributes:
      * *Hr*: bare hopping
      * *iscorrelated*: a list with same structure as orbitals. For example, for a two band system with orbital [("s","p")], iscorrelated= [(0,1)] means the second orbital is correlated. This is for further interface and not used currently.
    """
    def __init__(self,Atoms,Hr=None,iscorrelated=None):
        """
        Init object. As in TB object, but including also the correlated orbitals

        Args:
          * *Atoms*: AtomsTB objects with structural information, see AtomsTB in tbASE
          * *Hr*: dict of hopping matrix. the bare hopping part is stored in self.Hr.
          * *iscorrelated*: (optional) set flags of correlated orbitals.

        Note h_onsite is a block-diagonalized matrix, with non-zero elements only for the corresponding block for each atom.
        """
        super(tbGutz,self).__init__(Atoms,Hr)
        self.set_correlatedorbitals(iscorrelated)

    def set_correlatedorbitals(self,iscorrelated=None):
        """
        Set the attribute iscorrelated.
        """
        if iscorrelated is None:  # by default all orbitals are correlated.
            self.iscorrelated=[tuple([1 for j in i]) for i in self.Atoms.orbitals]
        else:
            self.iscorrelated=iscorrelated

    def get_correlatedorbitals(self):
        return self.iscorrelated

    def output_CyGutz(self,kps,num_electrons=None,mpiinfo=None):
        """
        This is output for the gl_interface to CyGutz codel according to Yao's format.

        Args:
          * *kps*: kpoints, numpy.ndarry (N,3), in unit of reciprocal lattice vectors.
          * *num_electrons*: number of electrons
          * *mpiinfo*: (nprocs, master), by default, nprocs=1,master=0

        Outputs: generating three files:
          * *GUTZ1.INP*       :
          * *GUTZ2.INP*       :
          * *GUTZ3.INP*       :
          * *GUTZ4.INP*       :
          * *GUTZ5.INP*       :
          * *BNDU_0.TXT*      : <bare band wave functions | local correlated orbitals>
          * *U2H.INP*         : Coulomb matrix elements including spin-orbit indices.
        """
        ## GMPI_?.INP
        if mpiinfo is None:
            nprocs=1
            master=0
        else:
            nprocs=mpiinfo[0]
            master=mpiinfo[1]
        for myrank in xrange(nprocs):
            gl_interface.write_gl_mpi(myrank,nprocs,master)

        ##
        num_atoms=len(self.Atoms)
        ## Gutz1.INP
        units=0
        gl_interface.write_gutz1(num_atoms,units)

        # Gutz2.INP
        num_kpts=len(kps)
        index_spin_orbit= 2 if self.Atoms.spinorbit else 1
        index_spin_bare = 2 if (self.Atoms.spindeg and not self.Atoms.spinorbit) else 1
        # tricky thing, if spin-up and spin-down  re not supposed to be treated sepaarted, we set index_spin_orbit=2
        if index_spin_bare == 2: index_spin_orbit= 2

        max_num_bands=self.Atoms.nspinorbitals
        gl_interface.write_gutz2(index_spin_orbit,index_spin_bare,max_num_bands,num_kpts)

        # Gutz3.INP
        gl_interface.write_gutz3(0,unitary_trans=None,translations=None)

        # Gutz4.INP
        ### a bit tricky here since in principle cluster has to be defined.
        #So here I force here that all the correlated atom in the unit cell form a cluster and correspond to different orbitals, that is, the num_corr_atoms is always 1.
        num_corr_atoms=1
        max_dim_sorbit=self.Atoms.nspinorbitals
        ## not max_dim_sorbit consider spin index only when spinorbit is not considered.
        #if (not self.Atoms.spinorbit) and self.Atoms.spindeg:
        #    max_dim_sorbit/=2
        U_CH_to_local_basis=None
        gl_interface.write_gutz4(max_dim_sorbit,num_corr_atoms,U_CH_to_local_basis)

        ## Gutz5.INP
        weight_kpts=1.0/num_kpts*numpy.ones((num_kpts))
        ##
        if not self.Atoms.spindeg:
            weight_kpts*=2.0
        index_smear=0
        delta=1e-2
        if num_electrons is None: # by default, half filled.
            num_electrons=self.Atoms.nspinorbitals*1.0/len(self.Atoms.spin)
        index_bands=numpy.zeros((num_kpts,3),dtype=numpy.int)
        index_bands[:,0]=self.Atoms.nspinorbitals
        index_bands[:,1]=1
        index_bands[:,2]=self.Atoms.nspinorbitals
        gl_interface.write_gutz5(num_kpts,weight_kpts,index_smear,delta,num_electrons,index_bands)

        ### BNDU_
        # define and diganolize Hk
        Norb=self.Atoms.nspinorbitals
        for myrank in xrange(nprocs):
            nkt=(len(kps)+nprocs-1) // nprocs   # get a ceiling division
            k_list=[i for i in xrange(nkt*myrank,min(nkt*myrank+nkt,len(kps)))]
            ek_list=numpy.zeros((len(k_list),Norb),dtype=numpy.float)
            Uk_list=numpy.zeros((len(k_list),Norb,Norb),dtype=numpy.complex)
            for ik in k_list:
                hk,ikpscart=self.Hk([kps[ik]])
                ## add h_on_site to hk
                ekt,Ukt=numpy.linalg.eigh(hk[0])
                ### sort ekt Ukt
                sorta=ekt.real.argsort()
                ekt_sorted=ekt[sorta].real
                Ukt_sorted=Ukt[:,sorta]
                ek_list[ik]=ekt_sorted
                Uk_list[ik]=Ukt_sorted.transpose().conjugate()   # get the complex conjugate

            gl_interface.write_gutz_bndu(myrank,ek_list,Uk_list)
            #import  write_bndu
            #write_bndu.write_bndu([ek_list], [Uk_list.real],
            #        [Uk_list.imag], myrank)


if __name__=="__main__":
    #Example . SquareLattice.
    ### a normal square lattice. default in gallery
    aTB=TB.gallery().add_spindegeneracy()

    # unit cell
    ### a Gutz TB model on a square lattice.
    gTB=tbGutz(aTB.Atoms,aTB.Hr)
    ###print gTB.get_Hloc()[0]
    kps_size=(10,10,1)
    kps=kpoints.monkhorst_pack(kps_size)
    gTB.output_CyGutz(kps)
