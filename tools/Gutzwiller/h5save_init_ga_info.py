import h5py


def h5write_init_ga_info_help(items, f, label):
    if items is not None:
        for i, item in enumerate(items):
            f["/impurity_" + str(i) + '/' + label] = item


def h5write_init_ga_info(
        num_sigma=0, sigma_list=None, Lie_Jeven_list=None,
        Lie_Jodd_list=None, matrix_basis_list=None,
        S_vector_list=None, Mvec_list=None):
    f = h5py.File("init_ga_info.h5", 'w')
    f["/num_impurity"] = num_sigma
    h5write_init_ga_info_help(sigma_list, f, 'sigma')
    h5write_init_ga_info_help(Lie_Jeven_list, f, 'Lie_Jeven')
    h5write_init_ga_info_help(Lie_Jodd_list, f, 'Lie_Jodd')
    h5write_init_ga_info_help(matrix_basis_list, f, 'matrix_basis')
    h5write_init_ga_info_help(S_vector_list, f, 'S_vector')
    h5write_init_ga_info_help(Mvec_list, f, 'M_vector')
    f.close()


def h5save_init_ga_info(material):
    '''
    Save data for post-processing.
    '''
    Lie_Jeven_list, Lie_Jodd_list = material.get_myLieParameters()
    sigma_list = material.get_sigma_list()
    import matrix_basis as mb
    matrix_basis_list = mb.ListSigmaToMatrixBasis(sigma_list)
    S_vector_list, L_vector_list = material.get_SL_vector_list()
    Mvec_list = material.get_Mvec_list()

    num_sigma = len(sigma_list)
    h5write_init_ga_info(
            num_sigma=num_sigma, sigma_list=sigma_list,
            Lie_Jeven_list=Lie_Jeven_list,
            Lie_Jodd_list=Lie_Jodd_list,
            matrix_basis_list=matrix_basis_list,
            S_vector_list=S_vector_list,
            Mvec_list=Mvec_list)

    from fileio import write_WH_SL_VEC
    write_WH_SL_VEC(S_vector_list, L_vector_list)
