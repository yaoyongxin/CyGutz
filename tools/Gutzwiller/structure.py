def read_wien_structure():
  from getwiencase import get_wiencase
  structure_file = get_wiencase() + '.struct'
  from ase.io.wien2k import read_struct
  material = read_struct(structure_file)
  return material

def read_vasp_contcar():
  from ase.io.vasp import read_vasp
  material = read_vasp()
  return material

def get_myatoms():
  import sys
  if '-vasp' in sys.argv:
    material = read_vasp_contcar()
  else:
    material = read_wien_structure()
  from myatom import myAtoms
  symbols = material.get_chemical_symbols()
  cell = material.get_cell()
  scaled_positions = material.get_scaled_positions()
  return myAtoms(symbols = symbols, cell = cell, scaled_positions = scaled_positions, pbc = True )

def check_material(material, log):
  print >> log, " Cell: ", '\n', material.get_cell()
  print >> log, " Symbols: ", material.get_chemical_symbols()
  print >> log, " Scaled_positions: ", '\n', material.get_scaled_positions()
  print >> log, " Mapping-table equivalent atoms: ", material.get_idx_equivalent_atoms()
  print >> log, " Labels correlated atoms: ", material.get_myCorrAtoms()
  print >> log, " Orbital angular momentums of correlated electrons: ", material.get_myAM()
  print >> log, " Take into accout SOC or not: ", material.get_mySOC()
  print >> log, " Take into accout crystal-fields or not: ", material.get_myCF()
  print >> log, " Cut-off distance for symmetry analysis: ", material.get_sym_dist_cut()

