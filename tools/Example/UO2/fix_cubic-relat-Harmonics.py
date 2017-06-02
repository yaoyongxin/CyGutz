'''
Using the standard cubic relativistic Harmonics for WH_N2N and
WH_SIGMA.INP.
'''


import h5py
import pyglib.symm.unitary as un
import pyglib.math.matrix_basis as mb
import pyglib.io.fio as fio

# Typical cubic relativistic Harmonics in the literature.
jj_to_cubic_basis = un.jj_to_cubic_relativistic_harmonics('f')
U_list = [jj_to_cubic_basis]

# The corresponding self-energy structure.
Sigma_list = [[
        [ 1,  2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 3,  4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 3, 4, 0, 0, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 5, 6, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 7, 8, 0, 0, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 6, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 8, 0, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0],
        [ 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9]]]

matrix_basis_list = mb.ListSigmaToMatrixBasis(Sigma_list)
fio.write_Hs("WH_HS.INP", matrix_basis_list)
fio.write_TRANS("WH_N2N.INP", U_list)
fio.write_sigma_struct(Sigma_list)

f = h5py.File("init_ga_info.h5", 'a')
del f["/impurity_0/matrix_basis"], f["/impurity_0/sigma"]
f["/impurity_0/matrix_basis"] = matrix_basis_list[0]
f["/impurity_0/sigma"] = Sigma_list[0]
f.close()
