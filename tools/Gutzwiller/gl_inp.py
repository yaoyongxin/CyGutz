import numpy as np


def set_gl_inp_material(material, log, usr_input):
    '''
    Write GL.INP file given material
    '''
    NTYPE, NIONS, ITYPE_list, NI0_list, NIMAP_list, corr_atom_type, \
            type_1atom = material.get_myLabels()
    SOC = material.get_mySOC()
    CF = material.get_myCF()
    l_list = material.get_myllist()
    dim_list = [np.sum(2 * (2 * np.array(ls) + 1)) for ls in l_list]
    df_list = material.get_myAM()
    spin_pol = material.get_spin_polarization()
    # GL.INP
    set_gl_inp(
            spin_pol, SOC, CF, NTYPE, NIONS, ITYPE_list, NI0_list,
            NIMAP_list, corr_atom_type, type_1atom, df_list, dim_list,
            log, usr_input)
    # U_list and self energy
    U_list, Sigma_list = material.get_mySelfEnergy()
    for i, sigma in enumerate(Sigma_list):
        print >> log, "self_energy ", i
        print >> log, sigma
    # WH_RLNEF.INP
    set_wh_rlnef(U_list, log)
    # WH_HS.INP
    set_wh_hs(Sigma_list, U_list)
    # parameters.dat
    set_scf_parameters()
    # WH_N2N.INP
    set_wh_n2n(U_list)
    # WH_SIGMA_STRUCT.INP
    from fileio import write_sigma_struct
    write_sigma_struct(Sigma_list)


def set_gl_inp(spin_polarization, SOC, CF, NTYPE, NIONS, ITYPE_list,
        NI0_list, NIMAP_list, corr_atom_type, type_1atom, df_list,
        dim_list, log, usr_input=None):
    '''
    Write GL.INP file
    '''
    from usrqa import get_usr_input
    # Coulomb U
    print
    print " LHUB = 1: Slater-Condo parametrization."
    print "        2: Kanamori parametrization (useful for models)."
    print "        0: U_{i,j,k,l} (NO SPIN INDEX) = "
    print "           int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j)"
    print "           V(|r-r'|) phi(r_k) phi(r'_l)}} "
    print "           will be provided by file V2AO.INP"
    print "       -1: U_{i,j,k,l} (INCLUDING SPIN INDEX) = "
    print "           int_{dr int_{dr' phi^{*}(r_i) phi^{*}(r'_j)"
    print "           phi(r_k) phi(r'_l)}} "
    print "           will be provided by file V2H.INP"
    LHUB = get_usr_input(" Please select LHUB: ", ['1', '-1', '0', '2'])
    if usr_input is not None:
        usr_input.write(LHUB + '\n')
    # Here I assume that we want to generate and reuse the variational setup
    # files (GLXXX.VSP) and run self-consistently
    LSCF = '1'

    # Choose double counting functional
    print
    print " LDC = 12: Recommended. Standard double counting "
    print "           (updating Vdc at each charge iteration, "
    print "           initial n0 to be provided.) "
    print "        2: Fix double counting potential "
    print "           (keep same Vdc/n0 at each charge iteration,"
    print "           n0 to be provided.) "
    print "        1: Standard double counting potential "
    print "           (n0 self-consistently determined.) "
    print "        0:  No double counting. "
    LDC = get_usr_input(" Please select LDC: ", ['12', '0', '1', '2'])
    if usr_input is not None:
        usr_input.write(LDC + '\n')

    # Choose tolerance of {R, \lambda} for the solver
    print
    while True:
        try:
            ans = raw_input(
                    " Please enter the tolerance of solving the equations of" +
                    " {R, \lambda}\n (recommend 1.e-5 or smaller)...")
            r_tol = float(ans)
            if r_tol > 0 and r_tol < 1.e-3:
                break
        except:
            pass
    if usr_input is not None:
        usr_input.write(ans + '\n')

    # Choose spin factor
    print
    print " LNEWTON = 0: Recommended. Modified Powell hybrid method (HYDRD1)."
    print "          -1: Broyden method. Faster for solutions with Z "
    print "              much larger than 0, e.g., magnetic calculations."
    LNEWTON = get_usr_input(" Please select LNEWTON: ", ['-1', '0'])
    if usr_input is not None:
        usr_input.write(LNEWTON + '\n')

    if 'y' in spin_polarization:
        ISPIN = '2'
    else:
        ISPIN = '1'

    NMAX_ITER = '100000'

    # Cluster-GA?
    print
    print " LCLUSTER = 0: Single-atom impurity."
    print "            1: Multi-atom (cluster) impurity."
    LCLUSTER = get_usr_input(" Please select LCLUSTER: ", ['0', '1'])
    if usr_input is not None:
        usr_input.write(LCLUSTER + '\n')

    # Here it is assumed that the same type of variational setup
    # is used for all of the correlated atoms in the material
    # (it might be good to generalize it in the future)
    if '1' in ISPIN and 'n' in CF:
        print >> log, " LGPRJ = 1:  Assuming  $[\phi,N] = [\phi,J_k] = 0$" +\
                " for all $k=1,2,3$ (averaging over the crystal field," +\
                " but SOC fully taken into account)"
        LGPRJ = '1'
    elif 'n' in SOC:
        print >> log, " LGPRJ = 14:  Assuming  $[\phi,N] = [\phi,S_z] = 0$"+\
                "  (SOC neglected, magnetism and crystal field allowed)"
        LGPRJ = '14'
    else:
        print >> log, " LGPRJ = 11:  Assuming only that  $[\phi,N] = 0$" +\
                " (allowing simultaneously SOC, magnetism and crystal field)"
        LGPRJ = '11'

    # Solver embedding system
    print "\n Solution embedding system:"
    print " LEIGV = 0: Choose automatically solver depending "
    print "            on the size of the problem (DEFAULT) "
    print "         1: Exact diagonalization (ZHEEV) in LAPACK"
    print "         2: Lanczos (zhdrv1) in ARPACK "
    print "         3: Exact diagonalization (ZHEEVX, "
    print "            selective lowest two eigen-vectors) in LAPACK"
    print "         5: PRIMME (Recommended for large dimension.)"
    LEIGV = get_usr_input(" Please select LEIGV: ", ['5', '0', '1', '2', '3'])
    if usr_input is not None:
        usr_input.write(LEIGV + '\n')

    if LEIGV == '2':
        NEV_F = '20'
        nev_f = int(NEV_F)
        NCV_F_est = str(36 * nev_f)
        print >> log, ' LEIGV = ' + LEIGV + '   (Lanczos)'
        print >> log, ' NEV_F = ' + NEV_F + \
            ' (default number of eigenvectors calculated by Lanczos,' +\
            ' change in GL.INP if necessary)'
        print >> log, ' NCV_F = NEV_F x 36 = ' + NCV_F_est + \
            ' (default Number of Lanczos basis vectors, change in GL.INP' +\
            ' if necessary)'

    # Solver (R,\lambda) fix-point problem
    LSOLVER = '1'
    UJ_list = []
    VR_list = []
    n0_list = []
    for k in range(NTYPE):
        print '\n INFORMATION FOR ' + df_list[type_1atom[k]] + \
                ' ELECTRONS OF ' + corr_atom_type[k] + ' :'
        if LHUB == '1' or LHUB == '2':
            # UJ
            while True:
                UJ = raw_input(
                        ' Please provide interaction parameters U,J ' +\
                        '\n separated by a space: ')
                if len(UJ.split(" ")) != 2:
                    print " Wrong format, please try again!"
                else:
                    break
            if usr_input is not None:
                usr_input.write(UJ + '\n')
        else:
            UJ = "0 0"
        UJ_list.append(UJ.split())
        # Valence range
        while True:
            VR = raw_input(
                    ' Please provide N1,N2 defining valence range [N1,N2]' +
                    '\n separated by a space ( 0 < N1 < N2 < ' +
                    str(dim_list[type_1atom[k]]) + ' ): ')
            if len(VR.split(" ")) != 2:
                print " Wrong format, please try again!"
            else:
                break
        VR_list.append(VR)
        if usr_input is not None:
            usr_input.write(VR + '\n')
        N1, N2 = VR.split(" ")

        if not (LDC == '0' or LDC == '1'):
            # n0
            while True:
                ans = raw_input(
                        ' Please provide guess n0 for valence\n (' +
                        N1 + ' < n0 < ' + N2 +
                        ' overwritten by GL_NELF1.INP): ')
                n0 = ans.split()[0]
                if float(n0) > float(N1) and float(n0) < float(N2):
                    break

            n0_list.append(n0)
            if usr_input is not None:
                usr_input.write(ans + '\n')
            print >> log, '( the value of n0 will be overwritten ' +\
                    'by GL_NELF1.INP & GL_NELF1.OUT during the calculation )'

    def GL_print(key, value, comments):
        GL_INP.write("{:<12}{:<3}{:<10}{:<50}\n".format(key, '=', value,
                '# ' + comments))

    def GL_matprint(Mat):
        for i in range(Mat.shape[0]):
            GL_print(' '.join(map(str, Mat[i])))

    GL_INP = open("GL.INP", "w")
    GL_INP.write("# One key word per line.\n")
    GL_print('LSCF', LSCF, 'Always generate variational setup for a new run.')
    GL_print('RTOL', r_tol, 'Tolerance for solving {R, \lambda} equations.')
    GL_print('LGPRJ', LGPRJ, 'Variational setup.')

    if LEIGV == '1':
        GL_print('LEIGV', LEIGV, 'Exact diagonalization, ZHEEV.')
    elif LEIGV == '2':
        GL_print('LEIGV', LEIGV, 'Lanczos, ARPACK.')
        GL_print('NEV_F', NEV_F, 'Number of eigenvectors by lanczos.')
        GL_print('NCV_F', NCV_F_est, 'Number of Lanczos basis vectors.')
    elif LEIGV == '3':
        GL_print('LEIGV',LEIGV, 'Exact diagonalization, ZHEEVX.')
    elif LEIGV == '5':
        GL_print('LEIGV', LEIGV, 'PRIMEE (Iterative MultiMethod Eigensolver).')

    if LDC == '12':
        GL_print('LDC', LDC, 'Standard D.C., n0 updated in charge loop.')
    elif LDC == '2':
        GL_print('LDC', LDC, 'Fixed V_{dc}.')
    elif LDC == '1':
        GL_print('LDC', LDC, 'Standard D.C., n0 updated in CyGutz loop.')
    elif LDC == '0':
        GL_print('LDC', LDC, 'No double counting.')

    if ISPIN == '2':
        GL_print('ISPIN', ISPIN, 'Magnetism allowed in z-direction.')

    if LNEWTON == '-1':
        GL_print('LNEWTON', LNEWTON, 'Broyden method.')
    elif LNEWTON == '0':
        GL_print('LNEWTON', LNEWTON, 'Modified Powell hybrid method (HYDRD1).')

    if LHUB == '1':
        GL_print('LHUB', LHUB, 'Slater-Condo parametrization.')
    elif LHUB == '2':
        GL_print('LHUB', LHUB, 'Kanamori parametrization.')
    elif LHUB == '0':
        GL_print('LHUB', LHUB, 'V2AO.INP read in (no spin index).')
    elif LHUB == '-1':
        GL_print('LHUB', LHUB, 'V2H.INP read in (including spin index).')

    if LCLUSTER == '0':
        GL_print('LCLUSTER', LCLUSTER, 'Single-atom impurity.')
    elif LCLUSTER == '1':
        GL_print('LCLUSTER', LCLUSTER, 'Multi-atom (cluster) impurity.')

    if LSOLVER == '1':
        GL_print('LSOLVER', LSOLVER, 'Solve {R, \lambda} equations.')
        GL_print('NMAX_ITER', NMAX_ITER, 'Maximum CyGutz iterations.')

    GL_INP.write('\nATOM TYPE INFO\n')
    GL_INP.write("{:<25}{:<50}\n".format(NTYPE,
            '# Number of types of correlated atoms.'))

    for k in range(NTYPE):
        GL_print('NT', k + 1, 'The type of correlated atom.')
        _type = type_1atom[k]
        GL_INP.write("{:<8}{:<17}{:<50}\n".format(
                UJ_list[k][0], UJ_list[k][1],
                '# U,J of ' + df_list[_type] +
                ' electrons of ' + corr_atom_type[k]))
        GL_INP.write("{:<25}{:<50}\n".format(
                dim_list[type_1atom[k]] / 2,
                '# Number of ' + df_list[_type] +
                ' orbitals (no spin).'))

        valences = VR_list[k].split()
        if len(n0_list) > 0:
            GL_INP.write("{:<8}{:<8}{:<9}{:<50}\n".format(
                    valences[0], valences[1], n0_list[k],
                    '# Valence range chosen & n0.'))
        else:
            GL_INP.write("{:<8}{:<17}{:<50}".format(
                    valences[0], valences[1],
                    '# Valence range chosen & n0.'))

    GL_INP.write('\nATOM INFO\n')
    GL_INP.write("{:<25}{:<50}\n".format(
            NIONS, '# NIONS (total correlated atoms, equivalent or not).'))
    for i in range(NIONS):
        GL_INP.write("{:<8}{:<8}{:<9}{:<50}\n".format(
                ITYPE_list[i], NI0_list[i], NIMAP_list[i],
                '# itype, ni_global, nimap'))
    GL_INP.close()


def set_wh_rlnef(U_list, log):
    '''
    Write file WH_RLNEF.INP in the rotated basis given WH_RLNEF.INP_ORIG
    in the original basis.
    '''
    import os.path
    import fileio as fio
    if os.path.isfile("WH_RLNEF.OUT_ORIG"):
        log.write(" Found WH_RLNEF.INP.\n")
        try:
            R_All, LA1_All, NKS_All, E_Fermi = fio.read_RLNEF("WH_RLNEF.OUT_ORIG")
            for i, U in enumerate(U_list):
                Udagger = np.conj(U.T)
                R_All[i] = np.dot(Udagger, np.dot(R_All[i], U))
                LA1_All[i] = np.dot(Udagger, np.dot(LA1_All[i], U))
                NKS_All[i] = np.dot(Udagger.T, np.dot(NKS_All[i], U.T))
            import fileio as fio
            fio.write_RLNEF("WH_RLNEF.INP", R_All, LA1_All, NKS_All, E_Fermi)
            log.write(" WH_RLNEF.INP generated.\n")
        except:
            log.write(" WH_RLNEF.INP generation failed.\n")


def set_wh_hs(Sigma_list, U_list=None, file_name="WH_HS.INP", lsym=True):
    '''
    Write file WH_HS.INP in the rotated basis given WH_HS.INP_ORIG
    in the original basis. It stores the Hermitian matrix basis set.
    '''
    import matrix_basis as mb
    if lsym:
        matrix_basis_list = mb.ListSigmaToMatrixBasis(Sigma_list)
    else:
        matrix_basis_list = mb.ListMatrixStructToBasis(Sigma_list)
    import fileio as fio
    fio.write_Hs(file_name, matrix_basis_list)

    if U_list is None:
        return
    for i, U in enumerate(U_list):
        Udagger = np.conj(U.T)
        for iHs in xrange(len(matrix_basis_list[i])):
            matrix_basis_list[i][iHs] = np.dot(
                U, np.dot(matrix_basis_list[i][iHs], Udagger))
    fio.write_Hs(file_name + "_ORIG", matrix_basis_list)


def set_scf_parameters():
    with open("params.dat", 'w') as f:
        print >> f, "finish = 20"
        print >> f, "cc = 1e-3"
        print >> f, "ec = 1e-5"


def set_wh_n2n(U_list):
    '''
    The original representation of J is obtained by the basis,
    which by default, is the complex spherical Harmonics basis (CH)
    if the spin-orbit interaction is considered,
    or the relativistic Harmonics basis (JJ)
    if spin-orbit interaction is present.
    This means that  $J_{ij} \equiv <Y_{i}|\hat{J}|Y_{j}>$.
    The outputs of the symmetry code contain the matrices $J'_{ij}$
    and $U_{ij}$ satisfying the relation $J = U J' U^\dagger$.
    Using the original definition of $J_{ij}$ in terms of
    spherical harmonics, we have:
    $J'_{ij} = [U^\dagger J U]_{ij} =
    U^\dagger_{ii'} <Y_{i'}|\hat{J}|Y_{j'}> U_{j'j} \equiv
    <Y'_{i}|\hat{J}|Y'_{j}>$,
    where we have defined $|Y'_{j}> \equiv U_{j'j} |Y_{j'}>$
    ( or, equivalently, $|Y'_{j}> \equiv \hat{U} |Y_{j}>$ ).
    Thus, the first-quantization interpretation of $U$ is
    that $U_{ij} = <Y_{i}|Y'_{j}>$.
    The LDA+GA interface requires that we write in "WH_N2N.INP"
    the unitary transformations for the correlated atoms.
    '''
    import fileio as fio
    fio.write_TRANS("WH_N2N.INP", U_list)


if __name__ == "__main__":
    '''
    A simple test of gl_inp.py
    '''
    # GL.INP
    spin_pol = 'n'
    SOC = ['y', 'y', 'y']
    CF = ['n', 'n', 'n']
    NTYPE = 2
    NIONS = 3
    ITYPE_list = [1, 1, 2]
    NI0_list = [1, 1, 2]
    NIMAP_list = [1, 1, 3]
    corr_atom_type = ["C", "C", "Si"]
    type_1atom = [0, 0, 1]
    df_list = ["p", "p", "p"]
    dim_list = [6, 6, 6]
    log = open("init_ga_.log", 'w')
    set_gl_inp(
            spin_pol, SOC, CF, NTYPE, NIONS, ITYPE_list, NI0_list,
            NIMAP_list, corr_atom_type, type_1atom, df_list, dim_list, log)

    # WH_HS.INP
    sigma_list = []
    U_list = []
    for i in range(NIONS):
        sigma_list.append((np.arange(6 * 6) + 1).reshape(6, 6))
        U_list.append(np.identity(6, dtype=complex))
    set_wh_hs(sigma_list, U_list)

    # WH_N2N.INP
    set_wh_n2n(U_list)

    log.close()
