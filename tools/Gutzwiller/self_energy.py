import sys
import numpy as np
from scipy.linalg import block_diag
from matrix_util import trans_orbital_fast_to_spin_fast


def get_self_energy(l_list, spin_pol, orbital_pol, cf, soc, rotations=None):
    '''
    Dispatcher.
    '''
    perturbation = None
    if '-p' in sys.argv:
        perturbation = sys.argv[sys.argv.index('-p') + 1]

    if 'y' in orbital_pol and 'y' in soc:
        return get_self_energy_op_soc(l_list)
    elif 'y' in cf and 'y' in soc and perturbation is None:
        return get_self_energy_cf_soc(l_list, rotations)
    elif 'y' in cf and 'y' in soc and perturbation == 'cf':
        return get_self_energy_soc_perturb_cf(l_list, rotations)
    elif 'y' in soc:
        return get_self_energy_soc(l_list, spin_pol)
    elif 'y' in orbital_pol:
        return get_self_energy_op_nosoc(l_list, spin_pol)
    elif 'y' in cf:
        return get_self_energy_cf_nosoc(l_list, spin_pol, rotations)
    else:
        return get_self_energy_average(l_list, spin_pol)


def get_self_energy_op_soc(l_list):
    '''
    Get the self energy structure in the case of fully orbital symmetry
    breaking and with spin-orbit interaction.
    '''
    dim_tot = np.sum(2 * (2 * np.array(l_list) + 1))
    self_energy = np.arange(dim_tot * dim_tot, dtype=int) + 1
    self_energy = self_energy.reshape((dim_tot, dim_tot))
    U = np.identity(dim_tot, dtype=complex)
    return U, self_energy


def get_self_energy_op_nosoc(l_list, spin_pol):
    '''
    Get the self energy structure in the case of fully orbital symmetry
    breaking but with negligible spin-orbit interaction.
    '''
    dim_t = int(np.sum(2 * np.array(l_list) + 1) + 0.5)
    m_half = np.arange(dim_t * dim_t, dtype=int) + 1
    m_half = m_half.reshape((dim_t, dim_t))
    if 'y' in spin_pol:
        shift = np.max(m_half)
    else:
        shift = 0
    from scipy.linalg import block_diag
    self_energy = block_diag(m_half, m_half + shift)
    U = np.identity(dim_t*2, dtype=complex)
    return U, trans_orbital_fast_to_spin_fast(self_energy)


def get_self_energy_cf_soc(l_list, rotations):
    '''
    Get the self energy structure in the case of crystal field splitting
    and with spin-orbit interaction.
    '''
    from angular_momentum_1p import get_J_vector
    Jorig = get_J_vector(l_list, 'JJ')
    import AtomSymmetry as atsym
    J, U, self_energy = atsym.get_atom_Jnew(rotations, Jorig)
    return U, self_energy


def get_self_energy_soc_perturb_cf(l_list, rotations):
    '''
    Get the self energy structure in the case of adding crystal field
    perturbatively with spin-orbit interaction.
    '''
    from angular_momentum_1p import get_J_vector
    assert len(l_list) == 1, " Not tested for multiple-l!"
    Jorig = get_J_vector(l_list, 'JJ')
    l = (len(Jorig[0]) - 2) / 4
    import AtomSymmetry as atsym
    J1_orig = [Jc[:2 * l, :2 * l] for Jc in Jorig]
    J1, U1, self_energy1 = atsym.get_atom_Jnew(rotations, J1_orig)
    J2_orig = [Jc[2 * l:, 2 * l:] for Jc in Jorig]
    J2, U2, self_energy2 = atsym.get_atom_Jnew(rotations, J2_orig)
    self_energy2[np.where(self_energy2 > 0)] += np.max(self_energy1)
    return block_diag(U1, U2),  block_diag(self_energy1, self_energy2)


def get_self_energy_cf_nosoc(l_list, spin_pol, rotations):
    '''
    Get the self energy structure in the case of crystal field splitting
    and without spin-orbit interaction.
    '''
    from angular_momentum_1p import get_J_vector
    Jorig = get_J_vector(l_list, 'CH')
    import AtomSymmetry as atsym
    Jhalf, Uhalf, self_energy_half = atsym.get_atom_Jnew(rotations, Jorig)
    from scipy.linalg import block_diag
    U = block_diag(Uhalf, Uhalf)
    if 'y' in spin_pol:
        shift = np.max(self_energy_half)
    else:
        shift = 0

    self_energy = block_diag(
        self_energy_half, shift_self_energy(self_energy_half, shift))
    self_energy = trans_orbital_fast_to_spin_fast(self_energy)
    return U, self_energy


def get_self_energy_average(l_list, spin_pol):
    '''
    Get the self energy structure in the case of negligible crystal field
    splitting and without spin-orbit interaction.
    '''
    diag_elem = []
    elem_base = 1
    for l in l_list:
        diag_elem += [elem_base for i in range(2 * l + 1)]
        elem_base += 1
    elem_base -= 1
    if 'y' in spin_pol:
        diag_elem += [i + elem_base for i in diag_elem]
    else:
        diag_elem += diag_elem
    self_energy = np.diag(diag_elem)
    self_energy = trans_orbital_fast_to_spin_fast(self_energy)
    U = np.identity(len(self_energy), dtype=complex)
    return U, self_energy


def get_self_energy_soc(l_list, spin_pol):
    '''
    Get the self energy structure in the case that spin-orbit interaction
    is dominant.
    '''
    diag_elem = []
    elem_base = 1
    for l in l_list:
        if l == 0:
            if 'n' in spin_pol:
                elem = [elem_base, elem_base]
            else:
                elem = [elem_base, elem_base + 1]
        else:
            elem = []
            if 'n' in spin_pol:
                elem += [elem_base for i in range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base for i in range(int(2 * (l + 0.51)) + 1)]
            else:
                elem += [elem_base + i for i in range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base + i for i in range(int(2 * (l + 0.51)) + 1)]
        diag_elem += elem
        elem_base = max(elem) + 1
    self_energy = np.diag(diag_elem)
    U = np.identity(len(self_energy), dtype=complex)
    return U, self_energy


def shift_self_energy(self_energy, shift):
    '''
    Shift non-zero self-energy elements by shift.
    '''
    res = np.zeros_like(self_energy)
    spots = np.where(self_energy > 0.1)
    res[spots] = self_energy[spots] + shift
    return res


if __name__ == "__main__":
    '''
    A test.
    '''
    l_list = [1, 1]
    spin_pol = 'n'
    orbital_pol = 'n'
    cf = 'y'
    soc = 'n'
    rotations = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
    U, self_energy = get_self_energy(
        l_list, spin_pol, orbital_pol, cf, soc, rotations)
    print " self_energy = "
    print self_energy

    # H cluster
    symbols = ['H', 'H', 'H', 'H', 'H', 'H']
    positions = [[0., 0., 0.], [2.0, 0., 0.], [0., 2.0, 0.0],
                 [0., -2.0, 0.], [0., 0., 2.0], [0., 0., -2.0]]
    from molecule import xyz_get_rot_list
    rotations = xyz_get_rot_list(symbols, positions, log='screen')
    print "H-dimer rotations:"
    print rotations
    l_list = [0, 0, 1]
    spin_pol = 'n'
    orbital_pol = 'n'
    cf = 'y'
    soc = 'n'
    U, self_energy = get_self_energy(
        l_list, spin_pol, orbital_pol, cf, soc, rotations)
    print " self_energy = "
    print self_energy
