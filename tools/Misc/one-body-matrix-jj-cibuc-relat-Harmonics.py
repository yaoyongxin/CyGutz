'''
jj to cubic relativistic Harmonics.
'''
import h5py
import itertools as it
import numpy as np
from numpy.linalg import eigh, inv

def print_array(dm, msg):
  print msg
  print 'real part'
  for row in dm:
    print ''.join('%9.4f'%(x.real) for x in row)
  print 'imag part'
  for row in dm:
    print ''.join('%9.4f'%(x.imag) for x in row)

jj_to_cubic_basis = np.zeros((14,14))
jj_to_cubic_basis[0,0] = -np.sqrt(1./6.) # |5/2, -5/2>
jj_to_cubic_basis[4,0] =  np.sqrt(5./6.) # |5/2, +3/2> G7, 5/2, +
jj_to_cubic_basis[5,1] = -np.sqrt(1./6.) # |5/2, +5/2>
jj_to_cubic_basis[1,1] =  np.sqrt(5./6.) # |5/2, -3/2> G7, 5/2, -

jj_to_cubic_basis[4,2] =  np.sqrt(1./6.) # |5/2, +3/2>
jj_to_cubic_basis[0,2] =  np.sqrt(5./6.) # |5/2, -5/2> G81, 5/2, +
jj_to_cubic_basis[1,3] =  np.sqrt(1./6.) # |5/2, -3/2>
jj_to_cubic_basis[5,3] =  np.sqrt(5./6.) # |5/2, +5/2> G81, 5/2, -

jj_to_cubic_basis[3,4] =  1. # |5/2, +1/2> G82, 5/2, +
jj_to_cubic_basis[2,5] =  1. # |5/2, -1/2> G82, 5/2, -

jj_to_cubic_basis[13,6] =  np.sqrt(5./12.) # |7/2, +7/2>
jj_to_cubic_basis[ 9,6] =  np.sqrt(7./12.) # |7/2, -1/2> G6, 7/2, +
jj_to_cubic_basis[ 6,7] =  np.sqrt(5./12.) # |7/2, -7/2>
jj_to_cubic_basis[10,7] =  np.sqrt(7./12.) # |7/2, +1/2> G6, 7/2, -

jj_to_cubic_basis[12,8] = -np.sqrt(3./4.) # |7/2, +5/2>
jj_to_cubic_basis[ 8,8] =  np.sqrt(1./4.) # |7/2, -3/2> G7, 7/2, +
jj_to_cubic_basis[ 7,9] =  np.sqrt(3./4.) # |7/2, -5/2>
jj_to_cubic_basis[11,9] = -np.sqrt(1./4.) # |7/2, +3/2> G7, 7/2, -

jj_to_cubic_basis[13,10] =  np.sqrt(7./12.) # |7/2, +7/2>
jj_to_cubic_basis[ 9,10] = -np.sqrt(5./12.) # |7/2, -1/2> G81, 7/2, +
jj_to_cubic_basis[ 6,11] = -np.sqrt(7./12.) # |7/2, -7/2>
jj_to_cubic_basis[10,11] =  np.sqrt(5./12.) # |7/2, +1/2> G81, 7/2, -

jj_to_cubic_basis[12,12] = -np.sqrt(1./4.)  # |7/2, +5/2>
jj_to_cubic_basis[ 8,12] = -np.sqrt(3./4.)  # |7/2, -3/2> G82, 7/2, +
jj_to_cubic_basis[ 7,13] =  np.sqrt(1./4.)  # |7/2, -5/2>
jj_to_cubic_basis[11,13] =  np.sqrt(3./4.)  # |7/2, +3/2> G82, 7/2, -

print_array(np.dot(jj_to_cubic_basis, np.conj(jj_to_cubic_basis.T)), 'UU\dagger')
print_array(np.dot(np.conj(jj_to_cubic_basis.T),jj_to_cubic_basis), 'U\daggerU')

f = h5py.File("glog.h5", 'r')
dm_calc_basis = f["/Impurity_1/GA_NC_PHY"][...].T
jj_to_calc_basis = f["/Impurity_1/PreferedBasisToCurrentBasis"][...].T

print_array(dm_calc_basis, 'Density matrix in calculation basis:')

dm_jj_basis = np.dot(np.dot(jj_to_calc_basis,dm_calc_basis),np.conj(jj_to_calc_basis.T))
print_array(dm_jj_basis, 'Density matrix in jj basis:')

dm_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), dm_jj_basis), jj_to_cubic_basis)
print_array(dm_cubic_basis, 'Density matrix in cubic basis:')

la_calc_basis = f["/Impurity_1/GA_La"][...].T
print_array(la_calc_basis, 'Lambda in calculation basis:')

la_jj_basis = np.dot(np.dot(jj_to_calc_basis,la_calc_basis),np.conj(jj_to_calc_basis.T))
print_array(la_jj_basis, 'Lambda in jj basis:')

la_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), la_jj_basis), jj_to_cubic_basis)
print_array(la_cubic_basis, 'Lambda in cubic basis:')

R_calc_basis = f["/Impurity_1/GA_R"][...].T

# R_inv_calc_basis = inv(R_calc_basis)
# Rdagger_inv_calc_basis = inv(np.conj(R_calc_basis.T))
# sigma0_calc_basis = np.dot(R_inv_calc_basis, np.dot(la_calc_basis, Rdagger_inv_calc_basis))
# print_array(sigma0_calc_basis, 'Sigma0 in calculation basis:')
#
# sigma0_jj_basis = np.dot(np.dot(jj_to_calc_basis,sigma0_calc_basis),np.conj(jj_to_calc_basis.T))
# sigma0_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), sigma0_jj_basis), jj_to_cubic_basis)
# print_array(Sigma0_cubic_basis, 'Sigma0 in cubic basis:')

Z_calc_basis = np.dot(np.conj(R_calc_basis.T), R_calc_basis)
print_array(Z_calc_basis, 'Z matrix in calculation basis:')

Z_jj_basis = np.dot(np.dot(jj_to_calc_basis,Z_calc_basis),np.conj(jj_to_calc_basis.T))
print_array(Z_jj_basis, 'Z matrix in jj basis:')

Z_cubic_basis = np.dot(np.dot(np.conj(jj_to_cubic_basis.T), Z_jj_basis), jj_to_cubic_basis)
print_array(Z_cubic_basis, 'Z matrix in cubic basis:')


Z = np.array([Z_cubic_basis[[2,2], [2,13]], Z_cubic_basis[[13,13], [2,13]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array([Z_cubic_basis[[3,3], [3,12]], Z_cubic_basis[[12,12], [3,12]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array([Z_cubic_basis[[5,4], [5,10]], Z_cubic_basis[[10,10], [5,10]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v

Z = np.array([Z_cubic_basis[[4,5], [4,11]], Z_cubic_basis[[11,11], [4,11]]])
print Z
w, v = eigh(Z)

print 'w = ', w
print 'v: \n', v


