from __future__ import print_function
from mpi4py import MPI
from pythtb import w90
import numpy as np
import h5py, pickle, glob
from pyglib.symm.unitary import get_u_csh2rh_all
from scipy.linalg import block_diag
from scipy.io import FortranFile
from pyglib.estructure.gwannier import get_gmodel, mpiget_bndev, \
        get_wannier_den_matrix_risb
from pyglib.estructure.fermi import get_fermi_weight, get_fermi_level


def if_gwannier(corbs_list, delta_charge=0., wpath="./",
        prefix="wannier", k_grid=None, lrot_list=None,
        iso=1, ispin=1, ismear=0, delta=0.0258,
        icycle=0):
    # read output from wannier90.
    gwannier = w90(wpath, prefix)
    # save wannier90 model.
    with open("wannier90.pkl", "wb") as f:
        pickle.dump(gwannier, f, protocol=pickle.HIGHEST_PROTOCOL)
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
    # spin degeneracy
    spin_deg = 3-max(ispin, iso)

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # GMPI_x.h5 file
    kvec = set_kvec_para(numk)
    # wannier basis to complex spherical basis tranformation.
    u_wan2csh = get_wann2csh(nbmax, corbs_list)
    h1e_list = get_h1e_list_wannier(gwannier, corbs_list)
    nelectron = 0.
    with h5py.File('BAREHAM_{}.h5'.format(myrank), 'w') as f:
        for ik in range(kvec[0][2],kvec[0][2]+kvec[0][1]):
            kpt = k_mesh[ik]
            hmat = gmodel._gen_ham(kpt)
            # from wannier basis to correlated orbital-ordered csh basis.
            hmat = u_wan2csh.T.conj().dot(hmat).dot(u_wan2csh)
            f['/IKP_{}/ISYM_1/HK0_SPIN1'.format(ik+1)] = hmat.T
            evals, evecs = gmodel._sol_ham(hmat, eig_vectors=True)
            nelectron += np.count_nonzero(evals < 0.)*wk*spin_deg
            f['/IKP_{}/ek0_spin1'.format(ik+1)] = evals
            # evec is actually evec.T
            f['/IKP_{}/T_PSIK0_TO_HK0_BASIS_SPIN1'.format(ik+1)] = \
                    evecs.T.conj()
    nelectron = comm.reduce(nelectron, op=MPI.SUM)
    if myrank == 0:
        if icycle <= 1:
            wrt_ginit(gwannier.symbols, gwannier.lat,
                    gwannier.atomic_positions, u_wan2csh, lrot_list=lrot_list)

        ne_list = [[nbmax, 1, nbmax] for k in range(numk)]
        wk_list = [wk for k in range(numk)]
        nelectron = int(nelectron+0.5)+delta_charge
        wrt_gparambands(numk, nbmax, ne_list, wk_list, k_mesh, nelectron,
                h1e_list, iso=iso, ispin=ispin, ismear=ismear, delta=delta)


def if_gwannier90(corbs_list, delta_charge=0., wpath="./",
        prefix="wannier", k_grid=None, lrot_list=None,
        iso=1, ispin=1, ismear=0, delta=0.0258,
        icycle=0):
    # read output from wannier90.
    gwannier = w90(wpath, prefix)
    # save wannier90 model.
    with open("wannier90.pkl", "wb") as f:
        pickle.dump(gwannier, f, protocol=pickle.HIGHEST_PROTOCOL)
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
    # spin degeneracy
    spin_deg = 3-max(ispin, iso)

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()

    # GMPI_x.h5 file
    kvec = set_kvec_para(numk)
    # wannier basis to complex spherical basis tranformation.
    u_wan2csh = get_wann2csh(nbmax, corbs_list)
    h1e_list = get_h1e_list_wannier(gwannier, corbs_list)
    nelectron = 0.
    with h5py.File('BAREHAM_{}.h5'.format(myrank), 'w') as f:
        for ik in range(kvec[0][2],kvec[0][2]+kvec[0][1]):
            kpt = k_mesh[ik]
            hmat = gmodel._gen_ham(kpt)
            # from wannier basis to correlated orbital-ordered csh basis.
            hmat = u_wan2csh.T.conj().dot(hmat).dot(u_wan2csh)
            f['/IKP_{}/ISYM_1/HK0_SPIN1'.format(ik+1)] = hmat.T
            evals, evecs = gmodel._sol_ham(hmat, eig_vectors=True)
            nelectron += np.count_nonzero(evals < 0.)*wk*spin_deg
            f['/IKP_{}/ek0_spin1'.format(ik+1)] = evals
            # evec is actually evec.T
            f['/IKP_{}/T_PSIK0_TO_HK0_BASIS_SPIN1'.format(ik+1)] = \
                    evecs.T.conj()
    nelectron = comm.reduce(nelectron, op=MPI.SUM)
    if myrank == 0:
        if icycle <= 1:
            wrt_ginit(gwannier.symbols, gwannier.lat,
                    gwannier.atomic_positions, u_wan2csh, lrot_list=lrot_list)

        ne_list = [[nbmax, 1, nbmax] for k in range(numk)]
        wk_list = [wk for k in range(numk)]
        nelectron = int(nelectron+0.5)+delta_charge
        wrt_gparambands(numk, nbmax, ne_list, wk_list, k_mesh, nelectron,
                h1e_list, iso=iso, ispin=ispin, ismear=ismear, delta=delta)


def set_kvec_para(numk):
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
    return kvec


def get_h1e_list_wannier(gwannier, corbs_list):
    # List of one-body parts of Hamiltonian.
    h1e_list = []
    for i, corbs in enumerate(corbs_list):
        norbs = len(corbs)
        h1e_list.append([np.zeros((norbs,norbs), dtype=np.complex)])
        for j1 in range(norbs):
            _j1 = corbs[j1]
            for j2 in range(norbs):
                _j2 = corbs[j2]
                h1e_list[-1][0][j1,j2] = gwannier.ham_r[(0,0,0)]["h"][_j1,_j2]\
                        /float(gwannier.ham_r[(0,0,0)]["deg"])
    # Unitary transformation from complex Harmonics to real Harmonics.
    u_csh2rh_list = get_u_csh2rh_all([len(corbs) for corbs in corbs_list])

    for i, u_csh2rh in enumerate(u_csh2rh_list):
        h1e_list[i][0] = u_csh2rh.dot(h1e_list[i][0]).dot(u_csh2rh.T.conj())
    return h1e_list


def get_wann2csh(nbmax, corbs_list):
    '''
    get transformation from wannier basis to complex spherical harmonics basis.
    '''
    # basis transformation matrix which merges the corelated orbitals
    # to the first places, followed by the uncorrelated orbitals.
    ubasis = np.zeros((nbmax,nbmax), dtype=np.complex)
    # the spin-orbital index remapping accordingly.
    orbs_map = []
    for i, corbs in enumerate(corbs_list):
        norbs = len(corbs)
        for j1 in range(norbs):
            orbs_map.append(corbs[j1])
    # appending the uncorrelated orbitals
    for i in range(nbmax):
        if i not in orbs_map:
            orbs_map.append(i)
    # ubasis: <original basis | correlated orbs-first basis>
    for i in range(nbmax):
        ubasis[orbs_map[i],i] = 1.

    # Unitary transformation from complex Harmonics to real Harmonics.
    u_csh2rh_list = get_u_csh2rh_all([len(corbs) for corbs in corbs_list])
    # get the transformation from wannier basis to correlated
    # orbital-ordered complex spherical Harmonics basis.
    u_csh2rh = block_diag(*u_csh2rh_list)
    ncorbs = u_csh2rh.shape[0]
    u_wan2csh = ubasis.copy()
    u_wan2csh[:,:ncorbs] = ubasis[:,:ncorbs].dot(u_csh2rh.T.conj())
    return u_wan2csh


def wrt_ginit(symbols, cell, scaled_positions, u_wan2csh, lrot_list=None):
    with h5py.File("ginit.h5", "w") as f:
        f['/struct/symbols'] = symbols
        f['/struct/cell'] = cell
        f['/struct/scaled_positions'] = scaled_positions
        if lrot_list is not None:
            f['/struct/locrot_list'] = lrot_list
        f["/u_wan2csh"] = u_wan2csh


def wrt_gparambands(numk, nbmax, ne_list, wk_list, kpoints, nelectron,
        h1e_list, iso=1, ispin=1, ismear=0, delta=0.0258):
    with h5py.File('GPARAMBANDS.h5', 'w') as f:
        f['/iso/'] = [iso]
        f['/ispin'] = [ispin]
        f['/kptdim'] = [numk]
        f['/nbmax'] = [nbmax]
        f['/NE_LIST'] = ne_list
        f['/kptwt'] = wk_list
        f["/kpoints"] = kpoints
        f['/ismear'] = [ismear]
        f['/delta'] = [delta]
        f['/nelectron'] = [nelectron]
        f['/symnop'] = [1]
        f['/symie'] = [1]
        for i, h1es in enumerate(h1e_list):
            for isp, h1e in enumerate(h1es):
                f['/IMPURITY_{}/H1E_SPIN{}'.format(i+1, isp+1)] = h1e.T



def get_wannier_dat(path="./"):
    '''get the contents in wannier.dat Fortran binary file.
    '''
    with FortranFile("{}/wannier.dat".format(path), "r") as f:
        reals_lat = f.read_reals()
        recip_lat = f.read_reals()
        num_bands = f.read_ints()[0]
        num_wann = f.read_ints()[0]
        ndiv = f.read_ints()
        nqdiv = ndiv[0]*ndiv[1]*ndiv[2]
        kpts = f.read_reals().reshape((nqdiv, 3))
        include_bands = f.read_ints()
        # no swapaxes(1,2) because it will be written in fortran file
        wfwannier_list = f.read_reals().view(np.complex).reshape(\
                (nqdiv, num_wann, num_bands))
        bnd_es = f.read_reals().reshape((1, nqdiv, num_bands))
    return kpts, include_bands, wfwannier_list, bnd_es


def wannier_den_matrix(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    kpts, include_bands, wfwannier_list, bnd_es_in = get_wannier_dat(
            path=wannier_path)
    wfwannier_list = wfwannier_list.swapaxes(1, 2)
    gmodel = get_gmodel()
    bnd_es, bnd_vs = mpiget_bndev(kpts, gmodel, mode="risb")
    nktot = len(kpts)
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        delta = f["/delta"][0]
        ismear = f["/ismear"][0]
        iso = f["/iso"][0]
        num_elec = f["/nelectron"][0]
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    if rank == 0:
        efermi = get_fermi_level(bnd_es, wklist, num_elec, delta=delta, \
                ismear=ismear, iso=iso)
    else:
        efermi = None
    efermi = comm.bcast(efermi, root=0)
    ncpu = comm.Get_size()
    nk_cpu = nktot//ncpu
    if nk_cpu*ncpu < nktot:
        nk_cpu += 1
    # reduce bnd_es to local part only
    if rank == 0:
        bnd_es = bnd_es[:nk_cpu]
        wklist = wklist[:nk_cpu]
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, delta=delta,
            ismear=ismear, iso=iso)
    # calculate wannier_den_matrix
    wan_den = get_wannier_den_matrix_risb(bnd_vs, ferwes, wklist, nktot)
    if rank == 0:
        fwrt_wan_den(wan_den, wfwannier_list, include_bands)


def wannier_den_matrix_lda_chk(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    kpts, include_bands, _, bnd_es = get_wannier_dat(
            path=wannier_path)
    nktot = len(kpts)
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        num_elec = f["/nelectron"][0]

    # chop bnd_es
    nbnd = bnd_es.shape[2]
    nwan = int(num_elec/2+3)
    bnd_es = bnd_es[:,:,:nwan]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    efermi = get_fermi_level(bnd_es, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    wfwannier_list = []
    vmat = np.zeros((nbnd, nwan), dtype=np.complex)
    np.fill_diagonal(vmat, 1.0)
    for ik in range(nktot):
        wfwannier_list.append(vmat)
        wan_den[0].append(np.diag(ferwes[0][ik]*nktot/(2+0.j)))
    wan_den = np.asarray(wan_den)
    wfwannier_list = np.asarray(wfwannier_list)
    fwrt_wan_den(wan_den, wfwannier_list, include_bands)


def wannier_den_matrix_lda_chk2(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    kpts, include_bands, wfwannier_list, bnd_es = get_wannier_dat(
            path=wannier_path)
    wfwannier_list = wfwannier_list.swapaxes(1, 2)
    nktot = len(kpts)
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        num_elec = f["/nelectron"][0]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    efermi = get_fermi_level(bnd_es, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    for ik in range(nktot):
        dm = np.diag(ferwes[0][ik]*nktot/(2+0.j))
        wan_den[0].append(wfwannier_list[ik].T.conj().dot(dm).dot(
            wfwannier_list[ik]))
    wan_den = np.asarray(wan_den)
    fwrt_wan_den(wan_den, wfwannier_list, include_bands)


def wannier_den_matrix_lda_chk3(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    kpts, include_bands, wfwannier_list, bnd_es = get_wannier_dat(
            path=wannier_path)
    wfwannier_list = wfwannier_list.swapaxes(1, 2)
    nktot = len(kpts)
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        num_elec = f["/nelectron"][0]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    # get eigen-vector from interpolation
    bnd_es2 = [[]]
    bnd_ev2 = [[]]
    for ik in range(nktot):
        hk = wfwannier_list[ik].T.conj().dot(np.diag(bnd_es[0][ik])).dot(
                wfwannier_list[ik])
        w, v = np.linalg.eigh(hk)
        bnd_es2[0].append(w)
        bnd_ev2[0].append(v)
    bnd_ev2 = np.asarray(bnd_ev2)
    with h5py.File("ev_lda_ref.h5", "w") as f:
            f["e"] = bnd_es2
            f["v"] = bnd_ev2
    efermi = get_fermi_level(bnd_es2, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es2, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    for ik in range(nktot):
        dm = np.diag(ferwes[0][ik]*nktot/(2+0.j))
        wan_den[0].append(bnd_ev2[0][ik].dot(dm).dot(
            bnd_ev2[0][ik].T.conj()))
    wan_den = np.asarray(wan_den)
    fwrt_wan_den(wan_den, wfwannier_list, include_bands)


def fwrt_wan_den(wan_den, wfwannier_list, include_bands):
    """write wannier_den_matrix.dat file in fortran binary format.
    for spin-restricted case only as of now.
    note the index order wfwannier_list[isp, ibnd, iwann].
    """
    nktot = wan_den.shape[1]
    wfwannier_list = wfwannier_list.swapaxes(1, 2)
    nband = wfwannier_list.shape[2]
    nwann = wfwannier_list.shape[1]
    with FortranFile('wannier_den_matrix.dat', 'w') as f:
        f.write_record(300.0)
        f.write_record(nktot)
        f.write_record([nband for k in range(nktot)])
        f.write_record([nwann for k in range(nktot)])
        f.write_record([include_bands for k in range(nktot)])
        f.write_record(wfwannier_list)
        f.write_record(wan_den[0].swapaxes(1, 2))


def get_risb_bndes(path="./"):
    num_list = [int(x.split("_")[1].split(".")[0]) \
            for x in glob.glob(path+"/GBANDS_*.h5")]
    num_list.sort()
    bnd_es = []
    for isp in range(2):
        for i in num_list:
            fband = "GBANDS_{}.h5".format(i)
            with h5py.File(fband, "r") as f:
                if "/ISPIN_{}".format(isp+1) in f:
                    if len(bnd_es) == isp:
                        bnd_es.append([])
                    ik_start = f["/IKP_START"][0]
                    ik_end = f["/IKP_END"][0]
                    for ik in range(ik_start, ik_end+1):
                        bnd_es[isp].append(f["/ISPIN_{}/IKP_{}/ek".format(\
                                isp+1, ik)][()])
    bnd_es = np.asarray(bnd_es)
    return bnd_es
