'''
jj to cubic relativistic Harmonics.
'''

import h5py
import itertools as it
import numpy as np
from numpy.linalg import eigh, inv
import pyglib.basic.prints as pt
import pyglib.symm.unitary as un

jj_to_cubic_basis = un.jj_to_cubic_relativistic_harmonics('f')

pt.print_2d_array(
        np.dot(jj_to_cubic_basis, np.conj(jj_to_cubic_basis.T)), 'UU\dagger')
pt.print_2d_array(
        np.dot(np.conj(jj_to_cubic_basis.T), jj_to_cubic_basis), 'U\daggerU')

f = h5py.File("glog.h5", 'r')
dm_calc_basis = f["/Impurity_1/GA_NC_PHY"][...].T
jj_to_calc_basis = f["/Impurity_1/PreferedBasisToCurrentBasis"][...].T

pt.print_2d_array(dm_calc_basis, 'Density matrix in calculation basis:')

dm_jj_basis = np.dot(
        np.dot(jj_to_calc_basis, dm_calc_basis), np.conj(jj_to_calc_basis.T))
pt.print_2d_array(dm_jj_basis, 'Density matrix in jj basis:')

dm_cubic_basis = np.dot(
        np.dot(np.conj(jj_to_cubic_basis.T), dm_jj_basis), jj_to_cubic_basis)
pt.print_2d_array(dm_cubic_basis, 'Density matrix in cubic basis:')

#la_calc_basis = f["/Impurity_1/GA_La"][...].T
#pt.print_2d_array(la_calc_basis, 'Lambda in calculation basis:')
#
#la_jj_basis = np.dot(np.dot(jj_to_calc_basis,la_calc_basis),np.conj(jj_to_calc_basis.T))
#pt.print_2d_array(la_jj_basis, 'Lambda in jj basis:')
#
#la_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), la_jj_basis), jj_to_cubic_basis)
#pt.print_2d_array(la_cubic_basis, 'Lambda in cubic basis:')

R_calc_basis = f["/Impurity_1/GA_R"][...].T

# R_inv_calc_basis = inv(R_calc_basis)
# Rdagger_inv_calc_basis = inv(np.conj(R_calc_basis.T))
# sigma0_calc_basis = np.dot(R_inv_calc_basis, np.dot(la_calc_basis, Rdagger_inv_calc_basis))
# pt.print_2d_array(sigma0_calc_basis, 'Sigma0 in calculation basis:')
#
# sigma0_jj_basis = np.dot(np.dot(jj_to_calc_basis,sigma0_calc_basis),np.conj(jj_to_calc_basis.T))
# sigma0_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), sigma0_jj_basis), jj_to_cubic_basis)
# pt.print_2d_array(Sigma0_cubic_basis, 'Sigma0 in cubic basis:')

Z_calc_basis = np.dot(np.conj(R_calc_basis.T), R_calc_basis)
pt.print_2d_array(Z_calc_basis, 'Z matrix in calculation basis:')

Z_jj_basis = np.dot(
    np.dot(jj_to_calc_basis, Z_calc_basis), np.conj(jj_to_calc_basis.T))
pt.print_2d_array(Z_jj_basis, 'Z matrix in jj basis:')

Z_cubic_basis = np.dot(
        np.dot(np.conj(jj_to_cubic_basis.T), Z_jj_basis), jj_to_cubic_basis)
pt.print_2d_array(Z_cubic_basis, 'Z matrix in cubic basis:')


Z = np.array(
        [Z_cubic_basis[[2, 2], [2, 13]], Z_cubic_basis[[13, 13], [2, 13]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array(
        [Z_cubic_basis[[3, 3], [3, 12]], Z_cubic_basis[[12, 12], [3, 12]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array(
        [Z_cubic_basis[[5, 4], [5, 10]], Z_cubic_basis[[10, 10], [5, 10]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array(
        [Z_cubic_basis[[4, 5], [4, 11]], Z_cubic_basis[[11, 11], [4, 11]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v
