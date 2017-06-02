from __future__ import print_function

'''
Interface to CyGutz.
'''

import numpy as np


def write_gl_mpi(myrank, nprocs, master):
    '''
    GMPI_${myrank}.INP file.

    Args:
      * *myrank*: the rank of the current processor.
      * *nprocs*: total number of processors.
      * *master*: master node.

    Outputs:
      * Create GMPI_${myrank}.INP file.
    '''
    with open("GMPI_" + str(myrank) + ".INP", 'w') as f:
        f.write('{} {} {} {} {}'.format(myrank, nprocs, master, ".FALSE.", 1))


def write_gutz1(num_atoms, units):
    '''
    GUTZ1.INP file.

    Args:
      * *num_atoms*: total number of atoms (correlated and uncorrelated)
      * *units*

        *  *1*: Wien2k convention, Rydberg/Bohr.
        *  *0*: VASP convention, eV/A.

    Outputs:
      * Create GUTZ1.INP file.
    '''
    with open("GUTZ1.INP", 'w') as f:
        f.write('{} {}'.format(num_atoms, units))


def write_gutz2(index_spin_orbit, index_spin_bare, max_num_bands, num_kpts):
    '''
    GUTZ2.INP file.

    Args:
      * *index_spin_orbit*:

        * *1*: without spin-orbit interaction.
        * *2*: with spin-orbit interaction.
      * *index_spin_bare*:

        * *1*: spin-degenerate bare band dispersion. (e.g., LDA + G)
        * *2*: spin-polarized  bare band dispersion. (e.g.,
        if one wants LSDA + G)
      * *max_num_bands*: maximal number of bands over all the k-points.
      * *num_kpts*: total number of k-points.

    Outputs:
      * Create GUTZ2.INP file.
    '''
    with open("GUTZ2.INP", 'w') as f:
        f.write('{} {} {} {}'.format(index_spin_orbit, index_spin_bare,
                                     max_num_bands, num_kpts))


def write_gutz3(num_symop=0, unitary_trans=None, translations=None):
    '''
    GUTZ3.INP file. It contains the symmetry operations.

    Args:
      * *num_symop*:

        * *<0*: no symmetry operations.
        * *>0*: number of symmetry operation.
      * *unitary_trans*: unitary transformation part.
      * *translations*: transplational part.

    Outputs:
      * Create GUTZ3.INP file.
    '''
    with open("GUTZ3.INP", 'w') as f:
        f.write('{}\n'.format(num_symop))
        if num_symop > 0:
            for i, trans in enumerate(len(translations)):
                for row in np.array(unitary_trans).T:
                    f.write(' '.join(map(str, row)) + '\n')
                f.write(' '.join(map(str, trans)) + '\n')


def write_gutz4(max_dim_sorbit, num_corr_atoms, U_CH_to_local_basis=None):
    '''
    GUTZ4.INP file.

    Args:
      * *max_dim_sorbit*: maximal dimension of local correlated orbitals
      over all the atoms. It includes the spin-factor only if spin-orbit
      interaction is present.
      * *num_corr_atoms*: total number of correlated atoms.
      * *U_CH_TO_local_basis*: (complex) transformations from complex
      spherical Harmonics to the desried local basis, e.g.,
      real Harmonics or relativistic Harmonics, with dimension =
      num_corr_atoms * max_dim_sorbit * max_dim_sorbit.

    Outputs:
      * Create GUTZ4.INP file.
    '''
    with open("GUTZ4.INP", 'w') as f:
        f.write('{} {}\n'.format(max_dim_sorbit, num_corr_atoms))
        if U_CH_to_local_basis == None:
            return
        for U in U_CH_to_local_basis:
            U_real = np.real(U.T)
            for row in U_real:
                f.write(' '.join(["%20.16f" % (elem) for elem in row]) + '\n')
            U_imag = np.imag(U.T)
            for row in U_imag:
                f.write(' '.join(["%20.16f" % (elem) for elem in row]) + '\n')


def write_gutz5(num_kpts, weight_kpts, index_smear, delta, num_electrons, index_bands):
    '''
    GUTZ5.INP file.

    Args:
      * *num_kpts*: Total number of k-points.
      * *weight_kpts*: k-points weight.
      * *index_smear*

        * *-5*: linear tetrahedron method.
        * *-1*: Fermi smearing.
        * * 0*: Gaussian smearing.
      * *delta*: smearing factor (take the energy unit specified
      by units in GUTZ1.INP).
      * *num_electrons*: total number of electrons.
      * *index_bands*: dimention num_kpts * 3.

        * *index_bands[:,0]*: total number of bands for each k-point.
        * *index_bands[:,1]*: starting correlated band index for each k-point.
        * *index_bands[:,2]*: ending correlated band index for each k-point.

    Outputs:
      * Create GUTZ5.INP file.
    '''
    with open("GUTZ5.INP", 'w') as f:
        f.write('{}\n'.format(num_kpts))
        if num_kpts < 0:  # impurity model
            f.write('{}\n'.format(num_electrons))
        else:
            for i in range(0, num_kpts, 10):
                f.write(' '.join( \
                        map(str, weight_kpts[i: min(i + 10, num_kpts)])) + '\n')
            f.write('{} {}\n{}\n'.format(index_smear, delta, num_electrons))
            for indice in index_bands:
                f.write(' '.join(map(str, indice)) + '\n')


def write_gutz_bndu(myrank, ek_list, Uk_list):
    '''
    BNDU_${myrank}.TXT file.
    It contains the bare eigen-values and the correlated orbitals
    expanded by bare eigen-vectors.

    Args:
      * *myrank*: the rank of the current processor.
      * *ek_list*: dimension = num_local_k_points * num_bands.
      list of bare band eigen-values for each local k-point
      (considering k-point parallelization in round-robin fashion).
      * *Uk_list*: dimension = num_local_k_points *
      num_local_correlated_orbitals * num_correlated_bands.
      list of expansion coefficients ($<\psi_nk | \phi_i\alpha>$)
      of local orbitals with the basis of the bare bare eigen-vectors
      for each local k-point (considering k-point parallelization).
      (It could just be the complex conjugate of the coefficient
      matrix of the bare eigen-vector if the correlated orbitals
      were used as tight-binding basis.)

    Outputs:
      * Create BNDU_${myrank}.TXT file.
    '''
    with open("BNDU_" + str(myrank) + ".TXT", 'w') as f:
        for ik, ek in enumerate(ek_list):
            num_bands = len(ek)
            # dump eigen-values
            for ib in range(0, num_bands, 10):
                f.write(' '.join(
                        ["%20.12f" % (elem) for elem in
                         ek[ib: min(ib + 10, num_bands)]]) + '\n')
            # dump < \psi_nk | \phi_i\alpha >
            num_corr_bands = len(Uk_list[ik][0])
            # dump column first
            for phi in Uk_list[ik].transpose():
                for ib in range(0, num_corr_bands, 10):
                    f.write(' '.join(
                            ["%20.16f" % (np.real(elem)) for elem in
                             phi[ib: min(ib + 10, num_corr_bands)]]) + '\n')
                for ib in range(0, num_corr_bands, 10):
                    f.write(' '.join(
                            ["%20.16f" % (np.imag(elem)) for elem in
                             phi[ib: min(ib + 10, num_corr_bands)]]) + '\n')


def get_symmetrized_Umatrix(U_m):
    '''
    Get the fully symmetrized U matrix.
    '''
    num_orbs = len(U_m)
    for j1 in range(num_orbs):
        for j2 in range(num_orbs):
            for j3 in range(num_orbs):
                for j4 in range(num_orbs):
                    U1 = U_m[j1][j2][j3][j4]
                    if (np.abs(U1) < 1.e-10):
                        continue
                    if (j1 == j2 and j1 == j3 and j1 == j4):
                        continue
                    # 1
                    if (np.abs(U_m[j2][j1][j4][j3]) > 1.e-10):
                        assert np.abs(U_m[j2][j1][j4][j3] - U1) < 1.e-10
                    else:
                        U_m[j2][j1][j4][j3] = U1
                    # 2
                    if (np.abs(U_m[j3][j4][j1][j2]) > 1.e-10):
                        assert np.abs(
                            U_m[j3][j4][j1][j2] - np.conj(U1)) < 1.e-10
                    else:
                        U_m[j3][j4][j1][j2] = np.conj(U1)
                    # 3
                    if (np.abs(U_m[j4][j3][j2][j1]) > 1.e-10):
                        assert np.abs(
                            U_m[j4][j3][j2][j1] - np.conj(U1)) < 1.e-10
                    else:
                        U_m[j4][j3][j2][j1] = np.conj(U1)
    return U_m


def write_coulomb_full(U_matrix_full_list):
    '''
    V2H.INP file.

    Args:
      * *U_matrix_full_list*: Full two-body U matrix including orbital
        and spin indices. U_matrix_full[1,2,3,4] =
        \int{dr \int {dr' phi_1^{*}(r) phi_2^{*}(r') U(r - r')
        phi_4(r') phi_3(r)}}

    Outputs:
      * Create V2H.INP file.
    '''
    with open("V2H.INP", 'w') as f:
        for i, U_matrix in enumerate(U_matrix_full_list):
            U_matrix = get_symmetrized_Umatrix(U_matrix)
            print >> f, "NT=", i + 1  # one-based
            num_orbs = len(U_matrix)
            for j1 in range(num_orbs):
                for j2 in range(num_orbs):
                    for j3 in range(num_orbs):
                        for j4 in range(num_orbs):
                            U1 = U_matrix[j1][j2][j3][j4]
                            if np.abs(U1) < 1.e-10:
                                continue
                            f.write(" %3d %3d %3d %3d %20.12f %20.12f \n" % (
                                    j1 + 1, j2 + 1, j3 + 1, j4 + 1, U1.real,
                                    U1.imag))

if __name__ == "__main__":
    '''
    Simple test.
    '''
    # serial jobs.
    myrank = 0
    nprocs = 1
    master = 0
    write_gl_mpi(myrank, nprocs, master)

    # for super-cell cluster version, units eV/A.
    num_atoms = 1
    units = 0
    write_gutz1(num_atoms, units)

    # for superconductivity?
    index_spin_orbit = 2
    index_spin_bare = 2
    max_num_bands = 8
    num_kpts = 10
    write_gutz2(index_spin_orbit, index_spin_bare, max_num_bands, num_kpts)

    # for models
    write_gutz3()

    # Unitary transformations
    max_dim_sorbit = max_num_bands
    num_corr_atoms = 1
    write_gutz4(max_dim_sorbit, num_corr_atoms)

    #
    weight_kpts = np.empty(num_kpts)
    weight_kpts.fill(1. / num_kpts)
    index_smear = 0  # Gaussian smearing
    delta = 0.01  # smearing factor
    num_electrons = 2
    index_bands = []
    for k in range(num_kpts):
        index_bands.append([max_num_bands, 1, max_num_bands])
    write_gutz5(
        num_kpts, weight_kpts, index_smear, delta,
        num_electrons, index_bands)
