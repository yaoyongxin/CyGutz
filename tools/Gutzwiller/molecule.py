'''
Help funtions for extracting molecule from crystal and determine local rotational operations.
'''

import numpy as np

def xtal_get_local_rot(
        symbols, scaled_positions, cell, iat, dist_cut_,
        Mvec_list=None, Nmax=4, tol=1.e-7, log=None):
    '''
    Get rotation operations of center atom iat in crystal.
    '''
    molecule, mol_Mvec_list = xtal_extract_mol(
            symbols, scaled_positions, cell, iat,
            dist_cut_, Nmax=Nmax, Mvec_list=Mvec_list)
    return mol_get_rot_list(
            molecule, Mvec_list=mol_Mvec_list, tol=tol, log=log)


def xtal_extract_mol(
        symbols, scaled_positions, cell, iat, dist_cut_,
        Nmax=10, Mvec_list=None):
    '''
    Extracting molecule of center atom iat from crystal.
    atom iat is the molecule center
    '''
    center = np.copy(scaled_positions[iat])
    new_scaled_positions = np.copy(scaled_positions)
    new_scaled_positions = new_scaled_positions - center
    molecule_positions = []
    molecule_symbols = []
    pair_dist = []
    if Mvec_list is not None:
        molecule_Mvec_list = []
    else:
        molecule_Mvec_list = None
    # Get the (2*Nmax+1)x(2*Nmax+1)x(2*Nmax+1) block
    Nmax_list = np.arange(-Nmax, Nmax)
    for i in Nmax_list:
        for j in Nmax_list:
            for k in Nmax_list:
                for jat, sp in enumerate(new_scaled_positions):
                    n_sp = sp + np.array([i, j, k])
                    n_p =   n_sp[0] * np.array(cell[0]) + \
                            n_sp[1] * np.array(cell[1]) + \
                            n_sp[2] * np.array(cell[2])
                    pair_dist.append(np.linalg.norm(n_p))
                    molecule_positions.append(n_p)
                    molecule_symbols.append(symbols[jat])
                    if Mvec_list is not None:
                        molecule_Mvec_list.append(Mvec_list[jat])
    # Get reasnoable dist_cut
    dist = pair_dist[:]
    dist.sort()
    if dist_cut_ > dist[1]:
        dist_cut = dist_cut_
    else:
        dist_cut = dist[12] + 0.1
        print " The default dist_cut for extracting a centered cluster" +\
                "\n for symmetry evaluation = ", dist_cut
    molecule_positions = [molecule_positions[i]
            for i in xrange(len(pair_dist)) if pair_dist[i] < dist_cut]
    molecule_symbols = [molecule_symbols[i]
            for i in xrange(len(pair_dist)) if pair_dist[i] < dist_cut]
    if Mvec_list is not None:
        molecule_Mvec_list = [molecule_Mvec_list[i]
                for i in xrange(len(pair_dist)) if pair_dist[i] < dist_cut]

    from itertools import groupby
    mol_name = ''
    mol_symbols = sorted(molecule_symbols)
    for key, group in groupby(mol_symbols):
        mol_name += key + '_' + str(len(list(group))) + ' '

    print(" Molecule extracted " + mol_name + ':')
    print(" atom   x      y       z    distance")
    for symbo, position in zip(molecule_symbols, molecule_positions):
        print(" %3s %6.2f %6.2f %6.2f  %8.4f" %
                ((symbo,) + tuple(position)+
                (np.sqrt(np.dot(position, position)),)))

    from pymatgen import Molecule
    return Molecule(molecule_symbols, molecule_positions), molecule_Mvec_list


def xyz_get_rot_list(symbols, positions, Mvec_list=None, tol=1.e-7, log=None):
    '''
    Get the rotation element list of a molecule,
    given symboal and position list.
    '''
    from pymatgen import Molecule
    mol = Molecule(symbols, positions)
    return mol_get_rot_list(mol, Mvec_list=Mvec_list, tol=tol, log=log)


def check_rot_magnetism(rot, mol, Mvec_list, tol=1.e-7):
    coords = mol.cart_coords
    for i, coord in enumerate(coords):
        coordp = np.dot(rot, coord)
        diff = coords - coordp
        ind = np.where(np.all(np.abs(diff) < tol, axis=1))[0]
        assert len(ind) == 1, " Error in check_rot_magnetism!"
        if mol.species[i] != mol.species[ind[0]]:
            return False
        Mvec = np.dot(rot, Mvec_list[i])
        if np.max(np.abs(Mvec - Mvec_list[ind[0]])) > tol:
            return False
    return True


def mol_get_rot_list(mol, Mvec_list=None, tol=1.e-7, log=None):
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer,\
            generate_full_symmops
    from scipy.linalg import det
    analyzer = PointGroupAnalyzer(mol)
    if log == 'screen':
        print " sch_symbol = ", analyzer.sch_symbol
    elif log != None:
        print >> log, " sch_symbol = ", analyzer.sch_symbol
    # Pick rotations only.
    symmops = [o for o in analyzer.symmops
            if abs(det(o.rotation_matrix)-1)<1.e-6]
    symmops = generate_full_symmops(symmops, tol)
    rot_list=[o.rotation_matrix for o in symmops]
    # Check for magnetism.
    if Mvec_list is not None:
        rot_list_tmp = rot_list
        rot_list = []
        for rot in rot_list_tmp:
            if check_rot_magnetism(rot, mol, Mvec_list, tol):
                rot_list.append(rot)
        assert np.mod(len(rot_list_tmp), len(rot_list)) < tol, \
                " Error in mol_get_rot_list!"
    return rot_list


if __name__ == "__main__":
    '''
    A test.
    '''
    # H cluster
    symbols = ['H', 'H', 'H', 'H', 'H', 'H']
    positions = [[0., 0., 0.], [2.0, 0., 0.], [0., 2.0, 0.0],
                 [0., -2.0, 0.], [0., 0., 2.0], [0., 0., -2.0]]
    rotations = xyz_get_rot_list(symbols, positions, log='screen')
    print "H-dimer rotations:"
    print rotations
    # extract molecule
    symbols = ['Ce']
    scaled_positions = [[0, 0, 0]]
    cell = [[0.5, 0.5, 0.], [0.5, 0., 0.5], [0., 0.5, 0.5]]
    iat = 0
    dist_cut = -1.0
    molecule = xtal_extract_mol(symbols, scaled_positions, cell, iat, dist_cut)
    print molecule.species
    print molecule.sites
