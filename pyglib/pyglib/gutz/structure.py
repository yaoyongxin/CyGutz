from __future__ import print_function
import os
import h5py
from pyglib.gutz.gatom import gAtoms


def read_wien_structure():
    from getwiencase import get_wiencase
    structure_file = get_wiencase() + '.struct'
    from ase.io.wien2k import read_struct
    material = read_struct(structure_file)
    from pyglib.dft.wien import get_local_rotations
    locrot_list = get_local_rotations(structure_file)
    return material, structure_file, locrot_list


def read_vasp_contcar():
    from ase.io.vasp import read_vasp
    material = read_vasp()
    return material


def h5save_structure_params():
    import sys
    if '-vasp' in sys.argv:
        material = read_vasp_contcar()
        case=None
    else:
        material, case, locrot_list = read_wien_structure()
    with h5py.File('ginit.h5', 'a') as f:
        f['/struct/symbols'] = material.get_chemical_symbols()
        f['/struct/cell'] = material.get_cell()
        f['/struct/scaled_positions'] = material.get_scaled_positions()
        if case is not None:
            f['/struct/case'] = case
            f['/struct/locrot_list'] = locrot_list

def get_gatoms():
    if not os.path.isfile('ginit.h5'):
        h5save_structure_params()

    with h5py.File('ginit.h5', 'r') as f:
        symbols = f['/struct/symbols'][()]
        cell = f['/struct/cell'][()]
        scaled_positions = f['/struct/scaled_positions'][()]
        if '/struct/case' in f:
            case = f['/struct/case'][()]
        else:
            case = None
        if '/struct/locrot_list' in f:
            locrot_list = f['/struct/locrot_list'][()]
        else:
            locrot_list = None
    return gAtoms(symbols=symbols, cell=cell,
            scaled_positions=scaled_positions, pbc=True,
            wiencase=case, locrot_list=locrot_list)


def check_material(material, f):
    if f is None:
        return
    print(" Cell: \n", material.get_cell(), file=f)
    print(" Symbols: ", material.get_chemical_symbols(), file=f)
    print(" Scaled_positions: \n", material.get_scaled_positions(), file=f)
    print(" Mapping-table equivalent atoms: ",
            material.idx_equivalent_atoms, file=f)
    print(" Labels of correlated atoms: ", material.corr_list, file=f)
    print(" Orbital angular momentums of correlated electrons: ",
            material.df_list, file=f)
    print(" Take into accout SOC or not: ",
            material.iso, file=f)
    print(" Take into accout crystal-fields or not: ",
            material.crystal_field, file=f)
    print(" Cut-off distance for symmetry analysis: ",
            material.sym_dist_cut, file=f)
    for i, sigma in enumerate(material.sigma_list):
        print(" self energy structure ", i, file=f)
        print(sigma, file=f)

