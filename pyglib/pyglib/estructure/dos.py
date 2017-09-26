import os, h5py, glob
from math import pi, sqrt
import numpy as np
from future_builtins import zip
from pyglib.basic import units
import matplotlib.pyplot as plt


class DOS:

    def __init__(self, w_k, e_skn,  width=0.1, window=None, npts=201):
        """Electronic Density Of States object.

        width: float
          Width of guassian smearing.
        window: tuple of two float
          Use ``window=(emin, emax)``.  If not specified, a window
          big enough to hold all the eigenvalues will be used.
        npts: int
          Number of points.

        """

        self.npts = npts
        self.width = width
        self.w_k = w_k
        self.e_skn = e_skn

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = np.linspace(emin, emax, npts)


    def delta(self, energy):
        """
        Return a delta-function centered at energy.
        """
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)


    def get_energies(self):
        """"
        return energy mesh.
        """
        return self.energies


    def get_dos_component(self, psi_skn_w):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn, psi_kn_w in zip(self.e_skn,psi_skn_w):
            dos = np.zeros(self.npts, dtype=np.complex)
            # k-point
            for w, e_n, psi_n_w in zip(self.w_k, e_kn, psi_kn_w):
                # band index
                for e, psi_w in zip(e_n, psi_n_w):
                    dos += w * self.delta(e) * psi_w
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_t(self):
        """
        Get array of DOS values.
        """
        dos_list = []

        # spin
        for e_kn in self.e_skn:
            dos = np.zeros(self.npts)
            # k-point
            for w, e_n in zip(self.w_k, e_kn):
                # band index
                for e in e_n:
                    dos += w * self.delta(e)
            dos_list.append(dos)
        return np.array(dos_list)


def get_all_psi_skn_w(e_skn, psi_sksn, bnd_ne):
    '''
    Reduce psi_sksn to psi_skn_w.
    '''
    psi_skn_w = np.zeros(list(e_skn.shape) + [psi_sksn.shape[-1]], dtype=float)
    for isp in range(e_skn.shape[0]):
        for k in range(e_skn.shape[1]):
            for isy in range(psi_sksn.shape[2]):
                for n in range(bnd_ne[k, 1], bnd_ne[k, 2]):
                    for a in range(psi_sksn.shape[-1]):
                        psi_skn_w[isp, k, n, a] += (
                                psi_sksn[isp, k, isy, n - bnd_ne[k, 1], a] *
                                np.conj(
                                psi_sksn[isp, k, isy, n - bnd_ne[k, 1], a])).\
                                real
    return psi_skn_w / psi_sksn.shape[2]


def get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne):
    '''
    Reduce psi_sksn to psi_skn_w_ab.
    '''
    n_orb = psi_sksn.shape[-1]
    psi_skn_w_ab = np.zeros(list(e_skn.shape) + [n_orb, n_orb], dtype=complex)
    for k in range(e_skn.shape[1]):
        for n in range(bnd_ne[k, 1], bnd_ne[k, 2]):
            psi_skn_w_ab[:, k, n, :, :] = np.einsum('ijk, ijl->ikl',
                    psi_sksn[:, k, :, n - bnd_ne[k, 1], :],
                    np.conj(psi_sksn[:, k, :, n - bnd_ne[k, 1], :]))
    return psi_skn_w_ab / psi_sksn.shape[2]


def h5get_dos(ewin=(-3., 5.), delta=0.05, npts=1001):
    '''
    get total dos and the total correlated orbital component.
    '''
    with h5py.File("GPARAMBANDS.h5", 'r') as f:
        # list of k-point weight.
        w_k = f["/kptwt"][...]
        # band index list specifying the range of bands used
        # for the construction of correlated orbitals.
        bnd_ne = f["/NE_LIST"][...]
        # number of symmetry operations.
        nsym = f["/symnop"][0]
        # number of k-points
        nkp = f["/kptdim"][0]
        # maximal number of bands calculated over the k-points.
        nbmax = f["/nbmax"][0]

    with h5py.File("GLOG.h5", 'r') as f:
        # Gutzwiller fermi level
        e_fermi = f["/e_fermi"][0]
        # numbe of spin components
        nspin = f["/nspin"][0]
        # total number of correlated orbitals with spin included
        # for cases with spin-orbit interaction.
        nasotot = f["/nasotot"][0]

    if os.path.isfile('ginit.h5'):
        with h5py.File('ginit.h5', 'r') as f:
            use_rydberg = 'ryd' in f['/usrqa/unit'][()]
    else:
        use_rydberg = True

    # band energy array.
    e_skn = np.zeros((nspin,nkp,nbmax), dtype=np.float)
    # maximal number of bands for the construction of correlated orbitals.
    nbmaxin = np.max(bnd_ne[:,2]-bnd_ne[:,1]+1)

    # expansion coefficients of the correlated orbitals in terms of the
    # band wavefunctions, i.e., <\psi_{sks, n}|\phi_{sks, \alpha}>
    # with sks := ispin, ikpt, isym.
    psi_sksna = np.zeros((nspin,nkp,nsym,nbmaxin,nasotot),
            dtype=np.complex)

    # including the case with MPI run.
    for fname in glob.glob('GBANDS_*h5'):
        with h5py.File(fname, 'r') as f:
            for isp in range(nspin):
                for k in range(f['/IKP_START'][0],f['/IKP_END'][0]+1):
                    e_n = f['/ISPIN_{}/IKP_{}/ek'.format(isp+1,k)][...]
                    e_n = e_n - e_fermi
                    # convert to eV
                    if use_rydberg:
                        e_n *= units.Ryd_eV
                    nbands = len(e_n)
                    e_skn[isp,k-1,:nbands] = e_n
                    # for bands not available, push it to high
                    # irrelevant value.
                    e_skn[isp,k-1,nbands:] = 100.
                    for isym in range(nsym):
                        v = f['/ISPIN_{}/IKP_{}/ISYM_{}/EVEC'.format( \
                                isp+1,k,isym+1)][:,:nasotot]
                        psi_sksna[isp,k-1,isym,:v.shape[0],:] = v

    # get total dos
    dos = DOS(w_k, e_skn,  width=delta, window=ewin, npts=npts)
    energies = dos.get_energies()
    dos_t = dos.get_dos_t()

    # get total correlated orbital component.
    psi_skn_w_ab = get_all_psi_skn_w_ab(e_skn, psi_sksna, bnd_ne)
    psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[...,:,:])
    dos_f = dos.get_dos_component(psi_sksn_f)

    return energies, dos_t, dos_f


def plot_dos_tf(energies, dos_t, dos_f):
    '''plot total dos and total correlated component.
    '''
    fig, ax = plt.subplots()
    ax.fill_between(
        energies, 0, dos_t, facecolor='grey', alpha=0.5)
    ax.plot(energies, dos_t, color='grey', label='tot')
    ax.plot(energies, dos_f, color='red', label='$f$')
    ax.axvline(x=0, ls='--')
    ax.set_ylim(0.)
    ax.set_ylabel("DOS (states/f.u.)")
    ax.set_xlabel("E (eV)")
    plt.title("Density of States with $f$-component")
    plt.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig('gap_j.pdf')



if __name__ == "__main__":
    '''
    Test run.
    '''
    energies, dos_t, dos_f = h5get_dos()
    # pick up-component.
    plot_dos_tf(energies, dos_t[0], dos_f[0])
