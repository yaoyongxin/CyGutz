from math import pi, sqrt
import numpy as np
from future_builtins import zip
from pyglib.basic import units


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



if __name__ == "__main__":
    '''
    Test run.
    '''
    import h5py
    import glob

    with h5py.File("GPARAMBANDS.h5", 'r') as f:
        w_k = f["/kptwt"][...]
        bnd_ne = f["/NE_LIST"][...]
        nsym = f["/symnop"][0]
        nkp = f["/kptdim"][0]
        nbmax = f["/nbmax"][0]

    with h5py.File("GLOG.h5", 'r') as f:
        e_fermi = f["/e_fermi"][0]
        nspin = f["/nspin"][0]
        nasotot = f["/nasotot"][0]

    e_skn = np.zeros((nspin,nkp,nbmax), dtype=np.float)
    nbmaxin = np.max(bnd_ne[:,2]-bnd_ne[:,1]+1)

    psi_sksna = np.zeros((nspin,nkp,nsym,nbmaxin,nasotot), dtype=np.complex)

    for fname in glob.glob('GBANDS_*h5'):
        with h5py.File(fname, 'r') as f:
            for isp in range(nspin):
                e_kn = []
                for k in range(f['/IKP_START'][0],f['/IKP_END'][0]+1):
                    e_n = f['/ISPIN_{}/IKP_{}/ek'.format(isp+1,k)][...]
                    e_n = e_n - e_fermi
                    e_n *= units.Ryd_eV
                    nbands = len(e_n)
                    e_skn[isp,k-1,:nbands] = e_n
                    e_skn[isp,k-1,nbands:] = 100.
                    for isym in range(nsym):
                        v = f['/ISPIN_{}/IKP_{}/ISYM_{}/EVEC'.format( \
                                isp+1,k,isym+1)][:,:nasotot]
                        psi_sksna[isp,k-1,isym,:v.shape[0],:] = v


    window = (-10., 15.)
    dos = DOS(w_k, e_skn,  width=0.4, window=window, npts=1001)
    energies = dos.get_energies()
    dos_t = dos.get_dos_t()

    psi_skn_w_ab = get_all_psi_skn_w_ab(e_skn, psi_sksna, bnd_ne)
    psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[...,:,:])
    dos_f = dos.get_dos_component(psi_sksn_f)

#    import matplotlib
#    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure(0)
    plt.plot(energies, dos_t[0].real, '-', label='t')
    plt.plot(energies, dos_f[0].real, 'o', label='f')
    plt.legend()
    plt.xlabel("E (eV)")
    plt.ylabel("DOS (states/f.u.)")
#    plt.savefig('dos.pdf')
    plt.show()
