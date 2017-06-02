'''
A TB model interface, using ASE Atoms as bases.
Following the convention in Wannier90, H(R) is based on unit cell, not on orbitals. H(k)=sum_R H(R) exp(-ikR).

A tight-binding model is defined by a set of hoppoing matrix H(R). H(R) is a matrix with dimension of the number of orbitals in the unit cell. A usefull technique is to construct supercell, which needs to construct H(R) for the supercell
Dimension: all cells are three dimension in ASE, however one can mimick two dimensional model by setting hopping along the third axis zero.
'''
import numpy
import pickle

import glob
import shutil
import os

from ase import Atoms
from ase.dft import kpoints

class AtomsTB(Atoms):
    """
    ASE ATOMS object with orbitals.

    Attributes:
      * For attributes of ASE ATOMS object, please check ASE documentation.
      * *orbitals*:  orbital names/index, list with each item is a list of orbital names on each atom, for example, [["s"],["s","p"]]
      * *spin*:      name of spins, by default, spin=["up"], only one spin component is considered. Can be also ["up","down"].
      * nspinorbitals*: number of orbitals with possible spin degeneracy
    """
    def set_orbitals_spindeg(self,orbitals=None,spindeg=False,spinorbit=False):
        """
        Set the orbitals in a ASE ATOMS object.

        Args:
          * *orbitals*: orbitals names,list with each item is a list of orbital names on each atom, for example, [["s"],["s","p"]]
          * *spindeg*: spin degeneracy, Bool, if TRUE spin=["up"]; if False, spin=["up","down"]
          * *spinorbit*: True, orbitals are in fact spin-orbitals, if this is true, spin is set to ["so"], spindeg=False. False, orbitals are without spinorbit coupling.

        Attributes created:
          * *spindeg*:
          * *spin*:
          * *nspinorbitals*:
          * *idx_sao_spinorbital*: index(ispin,iatom,iorbital) is reduced to a number (ispinorbital):
          * *idx_spinorbital_sao*: spinorbital is transform to index(ispin,iatom,iorbital).

        Example::

        >>> # 1. set a atom with two bands "s" and spin degenaracy
        >>> a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
        >>> a.set_orbitals_spindeg(orbitals=[("s","p")],spindeg=True)
        """

        if orbitals is None: ## one s  orbital per atoms by default
            self.orbitals=[("s",) for i in xrange(len(self.positions))]
        else:
            self.orbitals=orbitals[:]   ## make a copy
        assert type(self.orbitals) == type([]), "orbitals should be a list of tuples"

        self.spindeg=spindeg
        self.spinorbit=spinorbit
        if spindeg:
            self.spin=["up","down"]
        elif spinorbit:
            self.spin=["so"]
        else:
            self.spin=["up"]

        self.nspinorbitals=0
        for i in self.orbitals:
            self.nspinorbitals+=len(i)
        self.nspinorbitals*=len(self.spin)

        ### build a index
        idx_spinorbital_sao={}
        idx_sao_spinorbital={}
        idx=0
        for ispin in self.spin:
            for iatom in xrange(len(self.positions)):
                for iorb in xrange(len(self.orbitals[iatom])):
                    idx_spinorbital_sao[idx]=(ispin,iatom,iorb)
                    idx_sao_spinorbital[(ispin,iatom,iorb)]=idx
                    idx+=1
        self.idx_spinorbital_sao=idx_spinorbital_sao
        self.idx_sao_spinorbital=idx_sao_spinorbital

class TB(object):
    """
    A tight-binding mode is defined with an ASE object Atoms and Hopping.
    """
    def __init__(self,AtomsTB,Hr=None):
        """
        Init TB object

        Args:
          * *AtomsTB*: ASE object define the unit cell, the atoms, with orbitals set.
          * *Hr*: Hr is a dict of (R, hop), with R is tranlational vector and hop is hopping matrix between unit cells.
            len(Hr[R])=self.AtomsTB.nspinorbitals.
        """
        self.Atoms=AtomsTB
        if Hr is None:
            self.Hr={}
        else:
            self.Hr=Hr


    @staticmethod
    def gallery(name="SimpleSquareLattice"):
        """
        A predefined gallery of predefined lattice that frequentely used.

        Args:
          * *name*: model name; currently "SimpleSquareLattice" implemented, a singe band mode on square lattice.

        Returns:
          * TB object

        """
        if name == "SimpleSquareLattice":
            a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
            a.set_orbitals_spindeg()
            aTB=TB(a)
            aTB.set_hop([((0,1,0),0,0,-1),
                 ((1,0,0),0,0,-1),
                 ((0,-1,0),0,0,-1),
                 ((-1,0,0),0,0,-1)])
            return aTB
        if name =="Chain_nn":
            a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
            a.set_orbitals_spindeg()
            aTB=TB(a)
            aTB.set_hop([((1,0,0),0,0,1),
                 ((-1,0,0),0,0,1)])
            return aTB

    def set_hop(self,hoppings=None):
        """
        Set hopping mannually. In addition to a Hr matrix. Hopping is a tupe with the form (R, iorb,jorb,t_hop).
        hoppings could be a list of hoppings.

        Args:
          * *hoppings*: a list of hopping defined as tuple(R, iorb,jorb,t_hop). For example, [((0,0),0,0,1)], or a tuple(R, hr), where hr is a matrix of hoppings

        Example::

        >>>a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
        >>>a.set_orbitals_spindeg()
        >>>a.set_hop([((0,1),0,0,1])   #note other directions has to be set too
        """
        if type(hoppings) == type((1,)): # if only set only one term, change it to list for consistency.
            hoplist=[hoppings,]
        else:
            hoplist=hoppings
        for ihop in hoplist:             #
            if len(ihop)>2:
                R,iorb,jorb,t_hop=ihop
                if R not in self.Hr:  # the matrix for R is not set yet
                    self.Hr[R]=numpy.zeros((self.Atoms.nspinorbitals,self.Atoms.nspinorbitals),dtype=numpy.complex)
                self.Hr[R][iorb,jorb]=t_hop
            else:
                R,hr=ihop
                self.Hr[R]=hr

    def supercell(self,extent=None,pbc=None):
        """
        Generate a supercell with respect to original unit cell.
        The main problem is to generate hopping matrix between super cells.
        Following ASE atoms, the supercell can only be constructed by repeating the unit cell in three directions.

        Args:
          * *extent*: a tuple(m,n,l) specifying the multiplicity of the unit cell in each direction along it lattice vectors, which consistute a supercell. If it is a number then extent=(m,m,m), thus, same multiplicity along each direction.
          * *pbc*: peoriodic boundary conditions. reserved flag for further usage in slabs, not used now.

        Returns:
          * *sTB*: a TB object on the supercell.
        """
        atoms=self.Atoms*extent
        m=extent
        if type(extent) is int:
            m=(extent,extent,extent)
        trans=[numpy.array(i) for i  in numpy.ndindex(m)]  # way to trans unit cell to supercell.
        dup=len(trans)          # num of duplicate of unit cell
        atoms.set_orbitals_spindeg(self.Atoms.orbitals*dup,self.Atoms.spindeg,self.Atoms.spinorbit)

        ex=numpy.diag(m)
        def reduceVector(R,trans):
            # for vector R, rewrite in the unit cell of supercell
            # relative position respect to supercell
            tpos=numpy.dot(numpy.array(R),numpy.linalg.inv(ex))
            # set tpos to nearest int if it is really close to the integer
            tpostoint=numpy.rint(tpos)
            select=numpy.where(abs(tpos-tpostoint)<1e-6)
            tpos[select]=tpostoint[select]
            # vector shift will shit tpos in the unit cell of supercell
            shift=numpy.floor(tpos).astype(int)
            # tau is the relative position in the unit cell of supercell
            tau=numpy.rint((numpy.dot(tpos-shift,ex))).astype(int)
            td=-1
            #print tpos, shift, tau,
            for i in xrange(len(trans)):
                if numpy.linalg.norm(trans[i]-tau)<1e-6:
                    td=i
            #print shift,tau,td
            assert td !=-1,"shift vector is not in the unit cell"
            return tuple(shift),tau,td

        # define hopping matrix. This is done by iterating all possible hopping matrix from cells inside the supercell.
        Norb=self.Atoms.nspinorbitals
        Hr={}
        R=list(self.Hr)
        for i in xrange(dup):
            icord=trans[i]
            for iR in R:
                shiftR=numpy.array(iR)+icord
                shift,tau,td=reduceVector(shiftR,trans)
                if shift not in Hr:
                    Hr[shift]=numpy.zeros((Norb*dup,Norb*dup),dtype=type(self.Hr[R[0]][0,0]))
                if not self.Atoms.spindeg:   # spin degree of freedom is true
                    Hr[shift][i*Norb:(i+1)*Norb,td*Norb:(td+1)*Norb]=self.Hr[iR]
                else:  # make sure the spin index is the slowest one.
                    Hr[shift][i*Norb/2:(i+1)*Norb/2,td*Norb/2:(td+1)*Norb/2]=self.Hr[iR][:Norb/2,:Norb/2]
                    Hr[shift][i*Norb/2+dup*Norb/2:(i+1)*Norb/2+dup*Norb/2,td*Norb/2:(td+1)*Norb/2]=self.Hr[iR][Norb/2:,:Norb/2]
                    Hr[shift][i*Norb/2:(i+1)*Norb/2,td*Norb/2+dup*Norb/2:(td+1)*Norb/2+dup*Norb/2]=self.Hr[iR][:Norb/2,Norb/2:]
                    Hr[shift][i*Norb/2+dup*Norb/2:(i+1)*Norb/2+dup*Norb/2,td*Norb/2+dup*Norb/2:(td+1)*Norb/2+dup*Norb/2]=self.Hr[iR][Norb/2:,Norb/2:]

        sTB=TB(AtomsTB=atoms,Hr=Hr)
        return sTB

    def Hk(self,kps=numpy.zeros((1,3))):
        """
        Construct Hamiltonian for given k points.Fourier transform of Hr.

        Args:
          * *kps*: array of size (n,3), n is number of k points. k point is required to be in unit of reciprocal periodic lattice vector.

        Return:
          * *Hk*: ndarray(n,norbitals,norbitals),with each element the Hamiltonian for the corresponding kpoints.
          * *kcart*: k points in cartesian coordinates.
        """
        reciprocalvector=self.Atoms.get_reciprocal_cell()
        hk=numpy.zeros((len(kps),self.Atoms.nspinorbitals,self.Atoms.nspinorbitals),dtype=numpy.complex)
        ikcart=[]
        for ik in xrange(len(kps)):
            ikcartesian=numpy.dot(kps[ik],reciprocalvector)
            ikcart.append(ikcartesian)
            #Sites_cart=numpy.dot(self.Atoms.positions,numpy.array(self.Lattice.LatticeVector))
            for ir in self.Hr:
                R=numpy.array(ir)
                Rcartesian=numpy.dot(R,numpy.array(self.Atoms.cell))
                expk=numpy.exp(-1j*numpy.dot(ikcartesian,Rcartesian)*2.0*numpy.pi)
                #print hk.shape,self.Hr[ir].shape,ikcartesian
                hk[ik]+=expk*self.Hr[ir]
        return hk,numpy.array(ikcart)

    def eigens(self,kps=numpy.zeros((1,3)),keepwf=False):
        """
        Eigenvalue of a given list of kpoints. It calls Hk and diagonalize the hamiltonian.

        Args:
          * *kps*: kpoints, numpy.ndarry (N,3), in unit of reciprocal lattice vectors.
          * *keepwf*: keep wave functions (True) or not (False). default: False

        Returns:
          * *(eks,Uk) or (eks,)*: eks, eigenvalues, taken to be real Uk, wavefunctions. Both are sorted according to the order of eks. If keepwf is True: return (eigenvalues, wavefunctions) else return (eigenvalues,)

        """
        Norb=self.Atoms.nspinorbitals
        ek=numpy.zeros((len(kps),Norb),dtype=numpy.float)
        if keepwf:
            Uk=numpy.zeros((len(kps),Norb,Norb),dtype=numpy.complex)
        kpscart=[]
        for ik in xrange(len(kps)):
            hk,ikpscart=self.Hk([kps[ik]])
            #if ik%50 ==0:
            #    print ik
            ekt,Ukt=numpy.linalg.eigh(hk[0])
                ### sort ekt Ukt
            sorta=ekt.real.argsort()
            ekt_sorted=ekt[sorta].real
            Ukt_sorted=Ukt[:,sorta]
            ek[ik]=ekt_sorted
            if keepwf:
                Uk[ik]=Ukt_sorted
            kpscart.append(ikpscart[0])
        if keepwf:
            #return ek,Uk,numpy.array(kpscart)
            return ek,Uk
        else:
            return ek,

    def get_bandstructure(self,kps,saveto=None,with_weights=False):
        """
        Get the band structure for given kpath.

        Args:
          * *kps*: kpoints, kpath object given by ase.dft.kpoints.
          * *saveto*: save data to specific file with_weights:

        """
        if with_weights:
            eks,Uks,kpscart=self.eigens(kps[0],keepwf=with_weights)
        else:
            eks,=self.eigens(kps[0])
        if saveto is None:
            filename="band.dat"
        else:
            filename=saveto
        with open(filename,"w") as f:
            for iorb in xrange(self.Atoms.nspinorbitals):
                for ik in xrange(len(kps[0])):
                    f.write("%f  %f "%(kps[1][ik],eks[ik,iorb]))
                    if with_weights:
                        for ilay in xrange(self.Atoms.nspinorbitals):
                            weight=Uks[ik,ilay,iorb]*Uks[ik,ilay,iorb].conjugate()
                            f.write("%f  "%(weight.real))
                    f.write("\n")
                f.write("\n")


    def get_dos(self,  kps_size=None,saveto=None,dos_mesh=None,eta=1e-2):
        """
        Get the density of states. Note: symmetry is not used

        Args:
          * *dos_mesh*: (optional),the mesh of dos, numpy.array. default: numpy.linspace(Emax,Emin,500), (Emax, Emin)= (max,min) of all the eigenvalues
          * *kps_size*: (optional),the size of Monkhorst_pack mesh, tuple. default: (8,8,8)
          * *saveto*: (optional), data is saved to a file.
          * *eta*: broadening, current only Gaussian braodening available, default: 1e-2

        Returns:
          * *dos*: density of states

        Examples::

        >>># dos of 2D square lattice
        >>>aTB=TB.gallery()
        >>>aTB.get_dos(kps_size=(400,400,1))

        """
        if kps_size is None:
            kps_size=(8,8,8)
        kps=kpoints.monkhorst_pack(kps_size)
        eks,=self.eigens(kps)

        if dos_mesh is None:
            dos_mesh=numpy.linspace(eks.min(),eks.max(),500)
        if saveto is None:
            filename="dos_ek.dat"
        else:
            filename=saveto
        dos=numpy.zeros(len(dos_mesh),dtype=numpy.float)
        for ek in numpy.nditer(eks):
            dos+=numpy.exp(-(ek-dos_mesh)**2/2.0/eta**2)
        dos*=1.0/numpy.sqrt(2*numpy.pi)/eta
        dos/=len(kps)

        with open(filename,"w") as f:
            for idx in xrange(len(dos_mesh)):
                f.write("%f  %f \n"%(dos_mesh[idx],dos[idx]))


    def add_spindegeneracy(self):
        """
        Add spin degeneracy to a non-spin TB hamiltonian.
        Note simply increase the number of orbitals and increase the size of Hr to 2x2 block diagonal form. up and down spin in different block

        Returns:
          * A new TB object with spin degeneracy
        """
        atoms=self.Atoms.copy()
        assert not self.Atoms.spindeg,"Spin degeneray is already considered.! ERROR!"
        assert not self.Atoms.spinorbit,"Cann't add spin degeneracy to spin-orbit orbitals! ERROR!"
        atoms.set_orbitals_spindeg(self.Atoms.orbitals,spindeg=True)
        norb=atoms.nspinorbitals
        Hr={}
        for iR in self.Hr:
            Hr[iR]=numpy.zeros((norb,norb),dtype=type(self.Hr[iR][0,0]))
            Hr[iR][0:norb/2,0:norb/2]=self.Hr[iR][:,:]
            Hr[iR][norb/2:,norb/2:]=self.Hr[iR][:,:]

        return TB(atoms,Hr)


    def trans_Nambubasis(self):
        """
        Transform the TB to nambu basis. Only hopping is changed.
        This transpose is correct only if the one-site hopping is set to zero.

        Returns:
          * A new TB project with modified Hr.
        """
        atoms=self.Atoms.copy()
        assert self.Atoms.spindeg, "Orbital has no spin degeneracy! Add spin degenaracy fisrt! ERROR!"
        assert not self.Atoms.spinorbit,"Cann't transform to Nambubasis with spin-orbit orbitals! ERROR!"
        atoms.set_orbitals_spindeg(orbitals=self.Atoms.orbitals,spindeg=self.Atoms.spindeg)
        norb=atoms.nspinorbitals
        Hr={}
        for iR in self.Hr:
            Hr[iR]=numpy.zeros((norb,norb),dtype=type(self.Hr[iR][0,0]))
            Hr[iR][0:norb/2,0:norb/2]=self.Hr[iR][0:norb/2,0:norb/2]
            ### for down spin, Hr[R]=-(Hr[-R]).transpose
            minusR=tuple([-i for i in iR])
            assert minusR in self.Hr, "Error!, inverse R is not in the hopping matrix"  # create Hr matrix
            Hr[iR][norb/2:,norb/2:]=-self.Hr[minusR].transpose()[norb/2:,norb/2:]

        return TB(atoms,Hr)


if __name__=="__main__":
    # The following is a simple test for the above codes.

    # AtomsTB object
    ### set up an AtomsTB object. one band without spin degeneracy.
    a=AtomsTB("N",[(0,0,0)],cell=(1,1,1))
    a.set_orbitals_spindeg()

    # sqare lattice
    ### set up a TB object and the hoppings. This corresponds to a 2D squared lattice with nearest-neighbor hopping. This is the default one in TB.gallery()
    aTB=TB(a)
    aTB.set_hop([((0,1,0),0,0,-1),
                 ((1,0,0),0,0,-1),
                 ((0,-1,0),0,0,-1),
                 ((-1,0,0),0,0,-1)])

    # bands and dos
    ### set special k points for bands
    kG=[0,0,0]
    kX=[0.5,0,0]
    kM=[0.5,0.5,0]
    ### set up a ase.dft.kpoints kpath object
    kps=kpoints.get_bandpath([kG,kX,kM,kG],a.cell)
    ### get band structure of a square lattice
    aTB.get_bandstructure(kps,saveto="pcell_band.dat")
