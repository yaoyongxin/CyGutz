import itertools as it
import numpy as np
import h5py, os
from pyglib.gutz.atoms import Atoms
from pyglib.basic.units import Ryd_eV
try:
    import spglib
except:
    import pyspglib.spglib as spglib


class gAtoms(Atoms):

    '''
    A class to describe the calculated system with all the relevant
    informations, such as coordinates, correlated atom indices,
    symmetry information, etc.
    '''
    def __init__(self, symbols=None, positions=None, numbers=None,
                 masses=None, magmoms=None, scaled_positions=None,
                 cell=None, pbc=None, wiencase=None, locrot_list=None):
        Atoms.__init__(self, symbols=symbols, positions=positions,
                numbers=numbers, masses=masses, magmoms=magmoms,
                scaled_positions=scaled_positions,
                cell=cell, pbc=pbc)
        self.wiencase = wiencase
        self.locrot_list = locrot_list

    def h5set(self):
        from pyglib.io.h5io import h5auto_read
        with h5py.File('ginit.h5', 'r') as f:
            dist_cut = h5auto_read(f, \
                    '/usrqa/dist_cut', default=-1.0)
            self.set_sym_dist_cut(dist_cut)

            self.unit = h5auto_read(f, \
                    '/usrqa/unit', default='rydberg')

            spin_polarization = h5auto_read(f, \
                    '/usrqa/spin_polarization', default='n')
            self.set_ispin(spin_polarization)

            orbital_polarization = h5auto_read(f, \
                    '/usrqa/full_orbital_polarization', default='n')
            self.set_orbital_polarization(orbital_polarization)

            spin_orbit_coup = h5auto_read(f, \
                    '/usrqa/spin_orbit_coup', default='n')
            self.set_iso(spin_orbit_coup)

            crystal_field = h5auto_read(f, \
                    '/usrqa/crystal_field', default='y')
            self.set_crystal_field(crystal_field)

            lhub = h5auto_read(f, \
                    '/usrqa/u_matrix_type', default=1)
            self.set_lhub(lhub)

            ldc = h5auto_read(f, \
                    '/usrqa/ldc', default=1)
            self.set_ldc(ldc)

            idx_equivalent_atoms = h5auto_read(f, \
                    '/usrqa/idx_equivalent_atoms')
            idx_equivalent_atoms = idx_equivalent_atoms.tolist()
            self.set_idx_equivalent_atoms(idx_equivalent_atoms)

            if self.iso != self.ispin == 2:
                self.updn_full_list = h5auto_read(f, '/usrqa/updn_full_list')
                self.b_field = h5auto_read(f,
                        '/usrqa/bfield_ev_per_bohr_magneton', 0.5)
                if 'ryd' in self.unit:
                    self.b_field /= Ryd_eV
                self.updn_full_list = self.updn_full_list*self.b_field

            unique_corr_symbol_list = h5auto_read(f, \
                    '/usrqa/unique_corr_symbol_list')
            self.unique_corr_symbol_list = unique_corr_symbol_list.tolist()
            self.unique_df_list = h5auto_read(f, \
                    '/usrqa/unique_df_list')
            if self.ldc > 1:
                self.unique_nf_list = np.asfarray(h5auto_read(f, \
                        '/usrqa/unique_nf_list'))
            if self.lhub > 0:
                self.unique_u_list = np.asfarray(h5auto_read(f, \
                        '/usrqa/unique_u_list_ev'))
                self.unique_j_list = np.asfarray(h5auto_read(f, \
                        '/usrqa/unique_j_list_ev'))
                if 'ryd' in self.unit:
                    self.unique_u_list /= Ryd_eV
                    self.unique_j_list /= Ryd_eV

            self.set_nval_range_list()
            self.set_CorrAtoms()
            self.set_na2_list()

            lnewton = h5auto_read(f, \
                    '/usrqa/lnewton', default=0)
            self.set_gimix(lnewton)

            iembeddiag = h5auto_read(f, \
                    '/usrqa/iembeddiag', default=-1)
            self.set_giembeddiag(iembeddiag)


    # corr_list = list of indeces of correlated atoms listed in self.symbols
    def set_CorrAtoms(self):
        corr_list = []
        ityp_list = []
        imap_list = []
        df_list = []
        nf_list = []
        u_list = []
        j_list = []
        updn_list = []
        for i, s in enumerate(self.symbols):
            if s in self.unique_corr_symbol_list:
                corr_list.append(i)
                ityp_list.append(self.unique_corr_symbol_list.index(s))
                df_list.append(self.unique_df_list[ityp_list[-1]])
                if self.ldc > 1:
                    nf_list.append(self.unique_nf_list[ityp_list[-1]])
                if self.lhub > 0:
                    u_list.append(self.unique_u_list[ityp_list[-1]])
                    j_list.append(self.unique_j_list[ityp_list[-1]])
                idx_equ = self.idx_equivalent_atoms[i]
                imap = self.idx_equivalent_atoms.index(idx_equ)
                imap_list.append(imap)
                if self.iso != self.ispin == 2:
                    updn_list.append(self.updn_full_list[i])
        self.corr_list = corr_list
        self.ityp_list = ityp_list
        self.imap_list = imap_list
        self.df_list = df_list
        self.updn_list = updn_list
        self.u_list = u_list
        self.j_list = j_list

        if self.modify_mode and self.ldc > 1:
            with h5py.File('GPARAM.h5', 'r') as f:
                if '/dc_nelf_list' in f:
                    self.nelf_list = f['/dc_nelf_list'][()]
                else:
                    self.nelf_list = nf_list
        else:
            self.nelf_list = nf_list


    def set_lhub(self, lhub):
        self.lhub = lhub


    def set_ldc(self, ldc):
        self.ldc= ldc


    def set_modify_mode(self):
        self.modify_mode = os.path.exists('GPARAM.h5')


    def set_orbital_polarization(self, orbital_polarization):
        self.orbital_polarization = orbital_polarization


    def set_crystal_field(self, crystal_field):
        self.crystal_field = crystal_field


    def set_ispin(self, spin_polarization):
        self.ispin = 2 if 'y' in spin_polarization else 1


    def set_iso(self, spin_orbit_coup):
        self.iso = 2 if 'y' in spin_orbit_coup else 1


    def set_imap_list(self, imap_list):
        self.imap_list = imap_list


    def set_gimix(self, gimix):
        self.gimix = gimix


    def set_giembeddiag(self, giembeddiag):
        self.giembeddiag = int(giembeddiag)


    def set_sym_dist_cut(self, dist_cut):
        self.sym_dist_cut = dist_cut


    def get_Rotations(self, iat, Nmax=10):
        symbols = self.symbols
        cell = self.cell
        scaled_positions = self.scaled_positions
        if self.locrot_list is not None:
            locrot = self.locrot_list[iat]
        else:
            locrot = None

        equivalent_indices = self.idx_equivalent_atoms
        from pyglib.gutz.molecule import xtal_get_local_rot
        return xtal_get_local_rot(symbols, scaled_positions, cell, iat,
                self.sym_dist_cut, equivalent_indices=equivalent_indices,
                locrot=locrot, Nmax=Nmax)


    def get_EquivalentAtoms(self):
        # Getting general information about symmetry of the lattice
        if self.wiencase is not None:
            from pyglib.dft.wien import get_equivalent_atom_indices
            idx_equivalent_atoms = get_equivalent_atom_indices(
                    self.wiencase)
        else:
            dataset = spglib.get_symmetry_dataset(self, symprec=1e-5)
            idx_equivalent_atoms = dataset['equivalent_atoms'].tolist()
        return idx_equivalent_atoms


    def set_idx_equivalent_atoms(self, idx_equivalent_atoms):
        self.idx_equivalent_atoms = idx_equivalent_atoms


    def get_llist(self):
        if not hasattr(self, 'l_list'):
            from pyglib.symm.angular_momentum_1p import get_l_list_from_string
            l_list = []
            AM_list = self.df_list
            for i, s in enumerate(self.corr_list):
                l_sublist = get_l_list_from_string(AM_list[i])
                l_list.append(l_sublist)
        return l_list


    def set_SL_vector_list(self):
        '''
        Set S(L)_vector_list in the working local basis.
        '''
        l_list = self.get_llist()
        sx_list = []; sy_list = []; sz_list = []
        lx_list = []; ly_list = []; lz_list = []
        from pyglib.symm.angular_momentum_1p import get_S_vector,get_L_vector
        for i, _ls in enumerate(l_list):
            S_vector = get_S_vector(_ls, self.iso)
            L_vector = get_L_vector(_ls, self.iso)
            utrans = self.utrans_list[i]
            for j, _S, _L in it.izip(it.count(), S_vector, L_vector):
                S_vector[j] = utrans.T.conj().dot(_S).dot(utrans)
                L_vector[j] = utrans.T.conj().dot(_L).dot(utrans)
            sx_list.append(S_vector[0])
            sy_list.append(S_vector[1])
            sz_list.append(S_vector[2])
            lx_list.append(L_vector[0])
            ly_list.append(L_vector[1])
            lz_list.append(L_vector[2])
        self.sx_list = sx_list
        self.sy_list = sy_list
        self.sz_list = sz_list
        self.lx_list = lx_list
        self.ly_list = ly_list
        self.lz_list = lz_list


    def get_J(self, basis):
        '''
        Get J vector in basis of 'CH' (complex harmonics) if spin-orbit
        interaction is neglected or 'JJ' (J^2-J_z basis) if spin-orbit
        interaction is present.
        '''
        from pyglib.symm.angular_momentum_1p import get_J_vector
        J_list = []
        l_list = self.get_llist()
        for i,l in enumerate(l_list):
            J_list.append(get_J_vector(l, basis))
        return J_list


    def set_SelfEnergy(self):
        '''
        Get the self-energy structure and the corresponding unitary
        transformation of the local basis.
        '''
        from pyglib.gutz.self_energy import get_self_energy
        utrans_list = []
        sigma_list = []
        l_list = self.get_llist()
        rotations_list = []
        jgenerator_list = []
        idx_equivalent_atoms = self.idx_equivalent_atoms
        if self.modify_mode:
            f = h5py.File('GPARAM.h5', 'r')
        for i, l in enumerate(self.corr_list):
            if self.modify_mode:
                utrans_list.append(f['/IMPURITY_'+str(i+1)+ \
                        '/DB_TO_SAB'][()].T)
                sigma_list.append(f['/IMPURITY_'+str(i+1)+ \
                        '/SIGMA_STRUCT'][()].T)
                if 'y' in self.crystal_field and 'n' in \
                        self.orbital_polarization:
                    rotations_list.append(f['/IMPURITY_'+str(i+1)+ \
                            '/rotations'][()])
                    jgenerator_list.append(f['/IMPURITY_'+str(i+1)+ \
                            '/JGENERATOR'][()])
                    jgenerator_list[-1] = np.swapaxes(jgenerator_list[-1],1,2)
            else:
                if i > 0:
                    idx_equ = idx_equivalent_atoms[i]
                    if idx_equ in idx_equivalent_atoms[:i]:
                        imap = idx_equivalent_atoms.index(idx_equ)
                        utrans_list.append(utrans_list[imap])
                        sigma_list.append(sigma_list[imap])
                        if len(rotations_list) > 0:
                            rotations_list.append(rotations_list[imap])
                            jgenerator_list.append(jgenerator_list[imap])
                        continue
                if 'y' in self.crystal_field and 'n' in \
                        self.orbital_polarization:
                    rotations = self.get_Rotations(l)
                    rotations_list.append(rotations)
                else:
                    rotations = None
                jgen, utrans, sigma = get_self_energy(l_list[i], self.ispin,
                        self.orbital_polarization,
                        self.crystal_field,self.iso,
                        rotations=rotations)
                utrans_list.append(utrans)
                sigma_list.append(sigma)
                if jgen is not None:
                    jgenerator_list.append(jgen)
        self.utrans_list = utrans_list
        self.sigma_list = sigma_list
        if len(rotations_list) > 0:
            self.rotations_list = rotations_list
            self.jgenerator_list = jgenerator_list
        else:
            self.rotations_list = None
            self.jgenerator_list = None


    def set_LieParameters(self):
        '''
        Set the Lie parameters of rotation operators for odd and even J.
        '''
        import pyglib.symm.atom_symm as atsym
        Lie_Jeven_list = []
        Lie_Jodd_list = []
        if 'y' in self.crystal_field and 'n' in self.orbital_polarization:
            for i, l in enumerate(self.corr_list):
                rotations = self.get_Rotations(l)
                Lie = atsym.get_Lie_parameters(rotations, plus2pi=False)
                Lie_Jeven_list.append(Lie)
                Lie = atsym.get_Lie_parameters(rotations, plus2pi=True)
                Lie_Jodd_list.append(Lie)
            self.Lie_Jeven_list = Lie_Jeven_list
            self.Lie_Jodd_list = Lie_Jodd_list
        else:
            self.Lie_Jeven_list = None
            self.Lie_Jodd_list = None


    def set_one_particle_rotation_list(self):
        '''Setup the rotation matrix in the single-particle basis
        '''
        import pyglib.symm.atom_symm as atsym
        if self.Lie_Jodd_list is not None:
            sp_rotations_list = []
            for J, lies in it.izip(self.jgenerator_list, self.Lie_Jodd_list):
                sp_rotations_list.append(atsym.get_representation(J, lies))
            self.sp_rotations_list = sp_rotations_list
        else:
            self.sp_rotations_list = None


    def set_na2_list(self):
        l_list = self.get_llist()
        na2_list = [np.sum(2 * (2 * np.array(ls) + 1)) for ls in l_list]
        self.na2_list = na2_list


    def set_nval_range_list(self):
        self.nval_bot_list=[]
        self.nval_top_list=[]
        for s in self.unique_corr_symbol_list:
            self.nval_bot_list.append(nval_range_list[s][0])
            self.nval_top_list.append(nval_range_list[s][1])


    def set_v2e_list(self):
        from pyglib.mbody.coulomb_matrix import U_matrix
        from pyglib.math.matrix_util import trans_JJ_to_CH_sup_sdn
        mode_list = ['manual', 'slater-condon', 'kanamori']
        if self.lhub > 0:
            self.v2e_list = []
            self.u_avg_list = []
            self.j_avg_list = []
            l_list = self.get_llist()

            for i,imap in enumerate(self.imap_list):
                if i > imap:
                    self.v2e_list.append(self.v2e_list[imap])
                    self.u_avg_list.append(self.u_avg_list[imap])
                    self.j_avg_list.append(self.j_avg_list[imap])
                    continue
                assert len(l_list[i]) == 1, " more than one l with lhub>0!"
                l_imp = l_list[i][0]
                utrans = self.utrans_list[i]

                if self.iso == 2:
                    u_cmplx_harm_to_rel_harm =  \
                            trans_JJ_to_CH_sup_sdn(l_imp).T # convention
                    utrans = u_cmplx_harm_to_rel_harm.dot(utrans)
                v2e, u_avg, j_avg = U_matrix(mode_list[self.lhub], l_imp,
                        U_int=self.u_list[i],
                        J_hund=self.j_list[i], T=utrans)
                self.v2e_list.append(v2e)
                self.u_avg_list.append(u_avg)
                self.j_avg_list.append(j_avg)
        else:
            self.v2e_list = None
            self.u_avg_list = None
            self.j_avg_list = None

nval_range_list = {
    "H":[0,2],
    "Li":[0,2],
    "Na":[0,2],
    "Sc":[0,10],
    "Ti":[0,10],
    "V":[0,10],
    "Cr":[0,10],
    "Mn":[0,10],
    "Fe":[0,10],
    "Co":[0,10],
    "Ni":[0,10],
    "Cu":[0,10],
    "Y":[0,10],
    "Zr":[0,10],
    "Nb":[0,10],
    "Mo":[0,10],
    "Tc":[0,10],
    "Ru":[0,10],
    "Rh":[0,10],
    "Pd":[0,10],
    "Ag":[0,10],
    "Ce":[0,3],
    "Pr":[0,4],
    "Nd":[0,5],
    "Pm":[0,8],
    "Sm":[0,9],
    "Eu":[0,10],
    "Gd":[0,11],
    "Tb":[0,12],
    "Dy":[0,13],
    "Ho":[0,14],
    "Er":[1,14],
    "Tm":[2,14],
    "Yb":[3,14],
    "Hf":[0,10],
    "Ta":[0,10],
    "W":[0,10],
    "Re":[0,10],
    "Os":[0,10],
    "Ir":[0,10],
    "Pt":[0,10],
    "U":[0,4],
    "Np":[0,14],
    "Pu":[0,14],
    "Am":[0,14],
    "Cm":[0,14],
    "Bk":[0,14],
    "Cf":[0,14],
    "Es":[0,14],
    "Fm":[0,14],
    "Md":[0,14]
        }
