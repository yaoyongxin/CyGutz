import h5py, pickle, numpy, warnings, sys
from scipy.linalg import block_diag
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin
from builtins import zip,range


def get_csh2sab():
    '''get transformation from complex spherical harmonics basis to
    g-risb symmetry adapted basis.
    '''
    csh2sab_list = []
    with h5py.File("GPARAM.h5", "r") as f:
        imap_list = f["/IMAP_IMP"][()]
        for i,imap in enumerate(imap_list):
            if i == imap-1:
                csh2sab_list.append(\
                        f["/IMPURITY_{}/DB_TO_SAB".format(i+1)][()].T)
            else:
                csh2sab_list.append(csh2sab_list[imap-1])
    return csh2sab_list


def get_wan2sab():
    '''get transformation from wannier basis to cygutz symmetry adapted basis.
    '''
    csh2sab_list = get_csh2sab()
    with h5py.File("ginit.h5", "r") as f:
        wan2csh = f["/u_wan2csh"][()]
        spin_orb = f["/usrqa/spin_orbit_coup"][()]
    if 'n' in spin_orb.lower():
        iso = 1
    else:
        iso = 2
    csh2sab = []
    for u1 in csh2sab_list:
        n1 = u1.shape[0]/(3-iso)
        csh2sab.append(u1[:n1, ::3-iso])
    csh2sab = block_diag(*csh2sab)
    wan2sab = wan2csh.copy()
    n1 = csh2sab.shape[0]
    wan2sab[:,:n1] = wan2csh[:,:n1].dot(csh2sab)
    return wan2sab


def get_rlambda_in_wannier_basis():
    '''get the r-matrix and lambda-matrix from grisb calculation.
    '''
    # read from cygutz output
    with h5py.File("GLOG.h5", "r") as f:
        lambda_list = f["/BND_LAMBDA"][()]
        r_list = f["/BND_R"][()]

    # get transformation from wannier basis to cygutz symmetry adapted basis.
    wan2sab = get_wan2sab()
    # convert to wannier basis in each spin block.
    lam2 = []
    r2 = []
    n2 = wan2sab.shape[0]
    for rmat,lam in zip(r_list, lambda_list):
        rmat = rmat.T
        lam = lam.T
        n1 = rmat.shape[0]
        if n2 > n1:
            rmat = block_diag(rmat, numpy.eye(n2-n1, dtype=numpy.complex))
            lam = block_diag(lam, numpy.zeros((n2-n1,n2-n1), \
                    dtype=numpy.complex))
        r2.append(wan2sab.dot(rmat).dot(wan2sab.T.conj()))
        lam2.append(wan2sab.dot(lam).dot(wan2sab.T.conj()))
    return numpy.asarray(r2), numpy.asarray(lam2)


def get_h1e_in_wannier_basis():
    '''get the h1e-matrix.
    '''
    with h5py.File("GPARAM.h5", "r") as f:
        num_imp = f["/num_imp"][0]
        ispin = f["/ispin"][0]
    # read wannier to csh-basis transformation.
    with h5py.File("ginit.h5", "r") as f:
        wan2csh = f["/u_wan2csh"][()]
    n2 = wan2csh.shape[0]
    h1e_list = []
    with h5py.File("GPARAMBANDS.h5", "r") as f:
        for isp in range(ispin):
            h1e = []
            for i in range(num_imp):
                h1e.append(f["/IMPURITY_{}/H1E_SPIN{}".format(\
                        i+1, isp+1)][()].T)
            h1e = block_diag(*h1e)
            n1 = h1e.shape[0]
            if n2 > n1:
                h1e = block_diag(h1e, numpy.zeros((n2-n1,n2-n1), \
                        dtype=numpy.complex))
            h1e_list.append(wan2csh.dot(h1e).dot(wan2csh.T.conj()))
    return h1e_list


def get_structure():
    with h5py.File("ginit.h5", "r") as f:
        lattice = f["/struct/cell"][()]
        species = f["/struct/symbols"][()]
        coords = f["/struct/scaled_positions"][()]
    return Structure(lattice=lattice, species=species, coords=coords)


def get_bands(kpoints, gmodel, mode="tb"):
    if mode == "risb":
        r_mat, lam_mat = get_rlambda_in_wannier_basis()
        h1_mat = get_h1e_in_wannier_basis()
        ispin = r_mat.shape[0]
    else:
        ispin = 1
    bnd_es = []
    bnd_vs = []
    for isp in range(ispin):
        if mode == "risb":
            ispp = min(isp, len(h1_mat))
        bnd_es.append([])
        bnd_vs.append([])
        for kpt in kpoints:
            hmat = gmodel._gen_ham(kpt)
            if mode == "risb":
                hmat -= h1_mat[ispp]
                hmat = r_mat[isp].dot(hmat).dot(r_mat[isp].T.conj())
                hmat += lam_mat[isp]
            evals, evecs = gmodel._sol_ham(hmat, eig_vectors=True)
            bnd_es[isp].append(evals)
            bnd_vs[isp].append(evecs)
    return numpy.asarray(bnd_es), numpy.asarray(bnd_vs)


def get_symkpath(atol=1.e-6):
    struct = get_structure()
    kpath = HighSymmKpath(struct)
    # check warning and perform transformation if needed.
    if not numpy.allclose(kpath._structure.lattice.matrix,
            kpath._prim.lattice.matrix, atol=atol):
        warnings.warn("Input structure does not match expected standard "
                "primitive! Try k-path transformation.")
        ktrans = kpath._prim.lattice.reciprocal_lattice.matrix.dot(\
                numpy.linalg.inv(kpath._structure.lattice.\
                reciprocal_lattice.matrix))
        for kname in kpath.kpath["kpoints"]:
            kpath.kpath["kpoints"][kname] = \
                    kpath.kpath["kpoints"][kname].dot(ktrans)
    return kpath


def get_gmodel():
    with open("wannier90.pkl", "rb") as f:
        wannier90 = pickle.load(f)
    gmodel = wannier90.model(zero_energy=0.0, min_hopping_norm=1.e-8,
            ignorable_imaginary_part=1.e-8)
    return gmodel


def get_bands_symkpath(efermi=0., mode="tb"):
    gmodel = get_gmodel()
    kpath = get_symkpath()
    nk = (len(kpath.kpath["kpoints"])-1)*7
    k_vec, k_dist, k_node = gmodel.k_path(kpath.kpath["kpoints"].values(), nk)
    bnd_es, bnd_vs = get_bands(k_vec, gmodel, mode=mode)
    # prepare the args for pymatgen bs class.
    eigenvals = {}
    eigenvals[Spin.up] = bnd_es[0].T
    if len(bnd_es) == 2:
        eigenvals[Spin.down] = bnd_es[1].T
    bs = BandStructureSymmLine(k_vec, eigenvals, \
            kpath._structure.lattice.reciprocal_lattice, \
            efermi, kpath.kpath["kpoints"])
    return bs


def plot_bandstructure():
    if "-h" in sys.argv:
        print("usage: complot_bands.py [-g] [-f fname] [-el emin] [-eu emax]")
        sys.exit()

    if "-g" in sys.argv:
        mode = "risb"
    else:
        mode = "tb"
    bs = get_bands_symkpath(mode=mode)
    bsplot = BSPlotter(bs)

    if "-f" in sys.argv:
        fname = sys.argv[sys.argv.index("-f")+1]
    else:
        fname = "bndstr.pdf"
    if "-el" in sys.argv:
        emin = float(sys.argv[sys.argv.index("-el")+1])
    else:
        emin = numpy.min(bs.bands.values())
    if "-eu" in sys.argv:
        emax = float(sys.argv[sys.argv.index("-eu")+1])
    else:
        emax = numpy.max(bs.bands.values())
    bsplot.save_plot(fname, img_format="pdf", ylim=(emin, emax))



if __name__ == "__main__":
    plot_bandstructure()
