import itertools as it
import numpy as np
from atoms import Atoms
try:
    import spglib
except:
    import pyspglib.spglib as spglib


class myAtoms(Atoms):

    '''
    A class to describe the calculated system with all the relevant
    informations, such as coordinates, correlated atom indices,
    symmetry information, etc.
    '''

    # corr_list = list of indeces of correlated atoms listed in self.symbols
    def set_myCorrAtoms(self, corr_list):
        N = len(self.symbols)
        assert(len(corr_list) <= N)
        self.corr_list = corr_list

    def get_myCorrAtoms(self):
        return self.corr_list

    def set_1p_U_list(self, U_list):
        '''
        Save the one-particle unitary transformation from eithter CH
        (no spin-orbit) or JJ (with spin-orbit)
        to the basis used in the real calculation.
        '''
        self.one_body_U_list = U_list

    def get_1p_U_list(self):
        return self.one_body_U_list

    def set_sigma_list(self, sigma_list):
        self.sigma_list = sigma_list

    def get_sigma_list(self):
        return self.sigma_list

    # df_list = list of angular momentum correlated atoms (s, p, d or f)
    def set_myAM(self, df_list):
        N = len(self.symbols)
        assert(len(df_list) <= N)
        self.df_list = df_list

    def get_myAM(self):
        return self.df_list

    # soc_list = list of yes/no for correlated atoms (yes = soc, no = no soc)
    def set_mySOC(self, soc_list):
        N = len(self.symbols)
        assert(len(soc_list) <= N)
        self.soc_list = soc_list

    def get_mySOC(self):
        return self.soc_list

    def set_spin_polarization(self, spin_polarization):
        self.spin_polarization = spin_polarization

    def get_spin_polarization(self):
        return self.spin_polarization

    def set_orbital_polarization(self, orbital_polarization):
        self.orbital_polarization = orbital_polarization

    # soc_list = list of yes/no for correlated atoms (yes = soc, no = no soc)
    def set_myCF(self, cf_list):
        N = len(self.symbols)
        assert(len(cf_list) <= N)
        self.cf_list = cf_list

    def set_Mvec(self, Mvec_list):
        '''
        Set magnetization direction list.
        '''
        self.Mvec_list = Mvec_list

    def get_Mvec_list(self):
        return self.Mvec_list

    def set_sym_dist_cut(self, dist_cut):
        self.sym_dist_cut = dist_cut

    def get_sym_dist_cut(self):
        return self.sym_dist_cut

    def get_myCF(self):
        return self.cf_list

    def get_myRotations(self, iat, Nmax=10):
        symbols = self.symbols
        cell = self.cell
        scaled_positions = self.scaled_positions
        Mvec_list = self.get_Mvec_list()
        from molecule import xtal_get_local_rot
        return xtal_get_local_rot(symbols, scaled_positions, cell, iat,
                self.sym_dist_cut, Mvec_list=Mvec_list, Nmax=Nmax)

    def get_myEquivalentAtoms(self):
        # Getting general information about symmetry of the lattice
        dataset = spglib.get_symmetry_dataset(self, symprec=1e-5)
        return dataset

    def set_idx_equivalent_atoms(self, idx_equivalent_atoms):
        self.idx_equivalent_atoms = idx_equivalent_atoms

    def get_idx_equivalent_atoms(self):
        return self.idx_equivalent_atoms

    def get_myllist(self):
        from angular_momentum_1p import get_l_list_from_string
        l_list = []
        myAM_list = self.get_myAM()
        for i, s in enumerate(self.get_myCorrAtoms()):
            l_sublist = get_l_list_from_string(myAM_list[i])
            l_list.append(l_sublist)
        return l_list

    def get_SL_vector_list(self):
        '''
        Get S(L)_vector_list in the working local basis.
        '''
        from matrix_util import trans_orbital_fast_to_spin_fast
        l_list = self.get_myllist()
        soc_list = self.get_mySOC()
        S_vector_list = []
        L_vector_list = []
        from angular_momentum_1p import get_S_vector, get_L_vector
        for i, _ls in enumerate(l_list):
            S_vector = get_S_vector(_ls, soc_list[i])
            L_vector = get_L_vector(_ls, soc_list[i])
            if 'y' in soc_list[i]:
                U = self.one_body_U_list[i]
            else:
                U = trans_orbital_fast_to_spin_fast(self.one_body_U_list[i])
            for j, _S, _L in it.izip(it.count(), S_vector, L_vector):
                S_vector[j] = np.dot(np.conj(U.T), np.dot(_S, U))
                L_vector[j] = np.dot(np.conj(U.T), np.dot(_L, U))
            S_vector_list.append(S_vector)
            L_vector_list.append(L_vector)
        return S_vector_list, L_vector_list

    def get_myJ(self, basis):
        '''
        Get J vector in basis of 'CH' (complex harmonics) if spin-orbit
        interaction is neglected or 'JJ' (J^2-J_z basis) if spin-orbit
        interaction is present.
        '''
        from angular_momentum_1p import get_J_vector
        J_list = []
        l_list = self.get_myllist()
        for i in range(len(self.get_myCorrAtoms())):
            J_list.append(get_J_vector(l_list[i], basis))
        return J_list

    def get_mySelfEnergy(self):
        '''
        Get the self-energy structure and the corresponding unitary
        transformation of the local basis.
        '''
        from self_energy import get_self_energy
        U_list = []
        sigma_list = []
        soc_list = self.get_mySOC()
        cf_list = self.get_myCF()
        l_list = self.get_myllist()
        idx_equivalent_atoms = self.get_idx_equivalent_atoms()
        for i, l in enumerate(self.get_myCorrAtoms()):
            if l > 0 and idx_equivalent_atoms[l - 1] == idx_equivalent_atoms[l]:
                U_list.append(U_list[-1])
                sigma_list.append(sigma_list[-1])
                continue
            if 'y' in cf_list[i] and 'n' in self.orbital_polarization:
                rotations = self.get_myRotations(l)
            else:
                rotations = None
            U, sigma = get_self_energy(l_list[i], self.spin_polarization,
                    self.orbital_polarization, cf_list[i], soc_list[i],
                    rotations=rotations)
            U_list.append(U)
            sigma_list.append(sigma)
        self.set_1p_U_list(U_list)
        self.set_sigma_list(sigma_list)
        return U_list, sigma_list

    def get_myLieParameters(self):
        '''
        Get the Lie parameters of rotation operators for odd and even J.
        '''
        from angular_momentum_1p import get_J_vector
        import AtomSymmetry as atsym
        cf_list = self.get_myCF()
        Lie_Jeven_list = []
        Lie_Jodd_list = []
        for i, l in enumerate(self.get_myCorrAtoms()):
            if 'y' in cf_list[i] and 'n' in self.orbital_polarization:
                rotations = self.get_myRotations(l)
                Lie = atsym.get_Lie_parameters(rotations, plus2pi=False)
                Lie_Jeven_list.append(Lie)
                Lie = atsym.get_Lie_parameters(rotations, plus2pi=True)
                Lie_Jodd_list.append(Lie)
            else:
                Lie_Jeven_list.append([])
                Lie_Jodd_list.append([])
        return Lie_Jeven_list, Lie_Jodd_list

    def get_myLabels(self):
        corr_atoms = self.get_myCorrAtoms()
        # list symbols correlated atoms (not repeated, e.g., Ce1 and Ce2 have
        # the same symbol, Ce)
        corr_atom_type = []
        # correlated-atom label such that: corr_atom_type[k] =
        # self.symbols[type_1atom[k]]
        type_1atom = []
        for i, l in enumerate(corr_atoms):
            a = self.symbols[l]
            if not a in corr_atom_type:
                corr_atom_type.append(a)
                type_1atom.append(i)
        NI0_list = []  # NI0 of GL.INP
        NIMAP_list = []  # NIMAP of GL.INP
        X = self.get_idx_equivalent_atoms()
        NI0 = 0
        for Y in X:
            if Y in corr_atoms:
                if Y + 1 not in NIMAP_list:
                    NI0 = NI0 + 1
                NIMAP_list.append(Y + 1)
                NI0_list.append(NI0)
        ITYPE_list = []  # ITYPE of GL.INP
        for l in corr_atoms:
            a = self.symbols[l]
            ITYPE_list.append(corr_atom_type.index(a) + 1)
        NTYPE = max(ITYPE_list)
        NIONS = len(ITYPE_list)
        return NTYPE, NIONS, ITYPE_list, NI0_list, NIMAP_list, \
               corr_atom_type, type_1atom
