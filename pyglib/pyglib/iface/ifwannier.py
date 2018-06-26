from __future__ import print_function

from mpi4py import MPI
from pythtb import w90
import numpy as np
import h5py
from pyglib.symm.unitary import get_u_csh2rh_all


def if_gwannier(corbs_list, delta_charge=0., wpath="./",
        prefix="wannier", k_grid=None, lrot_list=None,
        ncpus=1, iso=1, ispin=1, ismear=0, delta=0.005,
        icycle=0):
    # read output from wannier90.
    gwannier = w90(wpath, prefix)
    gmodel = gwannier.model(zero_energy=0.0, min_hopping_norm=1.e-8,
            ignorable_imaginary_part=1.e-8)
    # uniform grid of k-points always contains the origin.
    if k_grid is None:
        k_grid = gwannier.kgrid
    elif isinstance(k_grid, np.int):
        k_grid = np.asarray(gwannier.kgrid)*k_grid
    else:
        assert(len(k_grid)==3), "Dim(k_grid) is not 3!"

    k_mesh = gmodel.k_uniform_mesh(k_grid)
    numk = k_mesh.shape[0]
    wk = 1./numk
    nbmax = gmodel.get_num_orbitals()

    # GMPI_x.h5 file
    comm = MPI.COMM_WORLD
    ncpu = comm.Get_size()
    nk_per_cpu = numk//ncpu
    if nk_per_cpu*ncpu < numk:
        nk_per_cpu += 1
    myrank = comm.Get_rank()
    kvec = [[]]
    kvec[0].append(myrank)
    k_start = nk_per_cpu*myrank
    kvec[0].append(min(nk_per_cpu, numk-k_start))
    kvec[0].append(k_start)
    with h5py.File("GMPI_{}.h5".format(myrank), "w") as f:
        f["/nvec"] = [1]
        f["/KVEC"] = np.asarray(kvec).T

    # List of one-body parts of Hamiltonian.
    h1e_list = []
    # basis transformation matrix which merges the corelated orbitals
    # to the first places, followed by the uncorrelated orbitals.
    ubasis = np.zeros((nbmax,nbmax), dtype=np.complex)
    # the spin-orbital index remapping accordingly.
    orbs_map = []
    for i,corbs in enumerate(corbs_list):
        norbs = len(corbs)
        if ispin == 1:
            h1e_list.append([np.zeros((norbs,norbs), dtype=np.complex)])
        else:
            raise ValueError("Not implemented for ispin = {}!".format(ispin))

        for j1 in range(norbs):
            _j1 = corbs[j1]
            orbs_map.append(_j1)
            for j2 in range(norbs):
                _j2 = corbs[j2]
                h1e_list[-1][0][j1,j2] = gwannier.ham_r[(0,0,0)]["h"][_j1,_j2]\
                        /float(gwannier.ham_r[(0,0,0)]["deg"])
    # appending the uncorrelated orbitals
    for i in range(nbmax):
        if i not in orbs_map:
            orbs_map.append(i)
    # ubasis: <original basis | correlated orbs first basis>
    for i in range(nbmax):
        ubasis[orbs_map[i],i] = 1.

    # Unitary transformation from complex Harmonics to real Harmonics.
    u_csh2rh = get_u_csh2rh_all([len(corbs) for corbs in corbs_list])

    # get the transformation from wannier basis to correlated
    # orbital-ordered complex spherical Harmonics basis.
    ncorbs = u_csh2rh.shape[0]
    u_wan2csh = ubasis.copy()
    u_wan2csh[:,:ncorbs] = ubasis[:,:ncorbs].dot(u_csh2rh.T.conj())

    nelectron = 0.
    with h5py.File('BAREHAM_{}.h5'.format(myrank), 'w') as f:
        for ik in range(kvec[0][2],kvec[0][2]+kvec[0][1]):
            kpt = k_mesh[ik]
            hmat = gmodel._gen_ham(kpt)
            # from wannier basis to correlated orbital-ordered csh basis.
            hmat = u_wan2csh.T.conj().dot(hmat).dot(u_wan2csh)
            f['/IKP_{}/ISYM_1/HK0'.format(ik+1)] = hmat.T
            evals, evecs = gmodel._sol_ham(hmat, eig_vectors=True)
            nelectron += np.count_nonzero(evals < 0.)*wk
            f['/IKP_{}/ek0'.format(ik+1)] = evals
            # evec is actually evec.T
            f['/IKP_{}/T_PSIK0_TO_HK0_BASIS'.format(ik+1)] = \
                    evecs.T.conj()

    nelectron = comm.reduce(nelectron, op=MPI.SUM)
    if myrank == 0:
        if icycle <= 1:
            with h5py.File("ginit.h5", "w") as f:
                f['/struct/symbols'] = gwannier.symbols
                f['/struct/cell'] = gwannier.lat
                f['/struct/scaled_positions'] = gwannier.atomic_positions
                if lrot_list is not None:
                    f['/struct/locrot_list'] = lrot_list
        with h5py.File('GPARAMBANDS.h5', 'w') as f:
            f['/iso/'] = [iso]
            f['/ispin'] = [ispin]
            f['/kptdim'] = [numk]
            f['/nbmax'] = [nbmax]
            f['/NE_LIST'] = [[nbmax, 1, nbmax] for k in range(numk)]
            f['/kptwt'] = [wk for k in range(numk)]
            f['/ismear'] = [ismear]
            f['/delta'] = [delta]
            with open("gwannier.log", "w") as flog:
                flog.write("Est. vs used num. of elecns per unit cell: {} {}".\
                        format(nelectron, int(nelectron+0.5)+delta_charge))
            f['/nelectron'] = [int(nelectron+0.5)+delta_charge]
            f['/symnop'] = [1]
            f['/symie'] = [1]
            for i, h1es in enumerate(h1e_list):
                for isp, h1e in enumerate(h1es):
                    f['/IMPURITY_{}/H1E_SPIN{}'.format(i+1, isp+1)] = h1e.T


