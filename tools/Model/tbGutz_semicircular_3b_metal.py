#!/usr/bin/evn python
#  A TB model interface for Gutz
import numpy
import matplotlib.pyplot as plt

# import CyGutz interface
import gl_interface

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), \
        "/home/ykent/bitbucket/cygutz_release/tools/Gutzwiller"))


def get_e_list():
    dos = lambda e: 2./numpy.pi * numpy.sqrt(1-e**2)
    cdos = lambda e: ( e*numpy.sqrt(1-e**2) + numpy.arcsin(e) ) / numpy.pi + 0.5

    '''
    e_list = numpy.linspace(-1,1,100)
    plt.plot(e_list,dos(e_list))
    plt.plot(e_list,cdos(e_list))
    plt.show()
    quit()
    '''

    from scipy.optimize import bisect

    nmesh = 2000
    cdos_list = numpy.linspace(0,1,nmesh+1)
    e_list = [bisect(lambda x: cdos(x)-a, -1 ,1) for a in cdos_list]
    e_list = numpy.asarray(e_list)
    e_list = (e_list[1:] + e_list[0:-1])/2

    '''
    plt.plot(e_list,cdos(e_list),'o')
    plt.plot(e_list,cdos_list,'-')
    plt.show()
    quit()
    '''
    return e_list


def setup_gutzwiller(e_list):
    nprocs=1
    master=0
    myrank=0
    gl_interface.write_gl_mpi(myrank,nprocs,master)

    # Gutz1.INP
    num_atoms=1
    units=0
    gl_interface.write_gutz1(num_atoms,units)

    # Gutz2.INP
    num_kpts = len(e_list)
    index_spin_orbit= 2
    index_spin_bare = 1
    max_num_bands = 6
    gl_interface.write_gutz2(index_spin_orbit,index_spin_bare,max_num_bands,
            num_kpts)

    # Gutz3.INP
    gl_interface.write_gutz3(0, unitary_trans=None,translations=None)

    # Gutz4.INP
    num_corr_atoms = 1
    max_dim_sorbit = max_num_bands
    U_CH_to_local_basis = None
    gl_interface.write_gutz4(max_dim_sorbit,num_corr_atoms,U_CH_to_local_basis)

    # Gutz5.INP
    weight_kpts = 1.0/num_kpts*numpy.ones((num_kpts))
    # spin degenerate
    weight_kpts *= 2.0
    index_smear = 0
    delta = 1e-2
    num_electrons = 3
    index_bands = numpy.zeros((num_kpts,3),dtype=numpy.int)
    index_bands[:,0] = max_dim_sorbit
    index_bands[:,1] = 1
    index_bands[:,2] = max_dim_sorbit
    gl_interface.write_gutz5(num_kpts,weight_kpts,index_smear,delta,
            num_electrons,index_bands)

    # sigma matrix structure
    from scipy.linalg import block_diag
    sig_half = numpy.arange(1,(max_dim_sorbit/2)**2+1).reshape(
            (max_dim_sorbit/2,max_dim_sorbit/2))
    sigma = numpy.zeros((max_dim_sorbit, max_dim_sorbit),dtype=numpy.int)
    sigma[:max_dim_sorbit/2,:max_dim_sorbit/2] = sig_half
    sigma[max_dim_sorbit/2:,max_dim_sorbit/2:] = sig_half
    sigma_list = [sigma]
    U_list = [numpy.identity(max_dim_sorbit, dtype = numpy.complex)]
    from gl_inp import set_wh_hs, set_gl_inp
    # write WH_HS.INP
    set_wh_hs(sigma_list, U_list)

    # write GL.INP
    spin_pol = 'n'
    SOC = ['y']; CF = ['y']
    NTYPE = 1; NIONS = 1; ITYPE_list = [1]
    NI0_list = [1]; NIMAP_list = [1]; corr_atom_type = ["X"]
    type_1atom = [0]; df_list = ["g"]; dim_list = [max_dim_sorbit]

    log = open("init_ga_a.slog", 'w')
    set_gl_inp(spin_pol, SOC, CF, NTYPE, NIONS, ITYPE_list, NI0_list,
            NIMAP_list, corr_atom_type, type_1atom, df_list, dim_list, log)
    log.close()


def setup_bare_hamiltonian(e_list, U):
    # BNDU_
    ek_list = [[0., e - U/2., 0., 0., e - U/2., 0.] for e in e_list]
    Uk_list = [numpy.diag([1., 1+0.j, 1., 1., 1., 1.]) for e in e_list]
    gl_interface.write_gutz_bndu(0,ek_list,Uk_list)
    with open('V2AO.INP', 'w') as f:
        f.write("NT=  1 \n")
        f.write("2 2 2 2 {}".format(U))

if __name__=="__main__":

    # Get energy mesh.
    e_list = get_e_list()

    # setup_gutzwiller done once and for all.
    setup_gutzwiller(e_list)

    # Modify U
    U = 3.0
    setup_bare_hamiltonian(e_list, U)
