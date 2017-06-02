from __future__ import print_function
from math import pi, sqrt
import numpy as np
from future_builtins import zip
from itertools import product,count


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


    def get_dos_noncoherent(self, psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
            i, j):
        """
        Non-coherent part of dos.
        """

        # symmetry operations
        s_range = range(psi_sksn.shape[2])
        norb = ecdfdffdelta.shape[0]

        dos_list = []
        # spin
        for e_kn, f_kn, psi_ksn in zip(self.e_skn, f_skn, psi_sksn):
            dos = np.zeros(self.npts, dtype=np.complex)

            icount = 0
            # k-point
            for w1, e_n1, f_n1, psi_sn1, nb1 in zip(self.w_k, e_kn, f_kn,\
                    psi_ksn, bnd_ne):
                # k-point

                icount += 1; print(icount)

                for w2, e_n2, f_n2, psi_sn2, nb2 in zip(self.w_k, e_kn, \
                        f_kn, psi_ksn, bnd_ne):
                    # k-point
                    for w3, e_n3, f_n3, psi_sn3, nb3 in zip(self.w_k, \
                            e_kn, f_kn, psi_ksn, bnd_ne):
                        for n1, e1, f1 in zip(count(),
                                e_n1[nb1[1]-1 : nb1[2]],
                                f_n1[nb1[1]-1 : nb1[2]]):
                            for n2, e2, f2 in zip(count(),
                                    e_n2[nb2[1]-1 : nb2[2]],
                                    f_n2[nb2[1]-1 : nb2[2]]):
                                    for n3, e3, f3 in zip(count(),
                                            e_n3[nb3[1]-1 : nb3[2]],
                                            f_n3[nb3[1]-1 : nb3[2]]):
                                        f123 = (1-f1)*(1-f2)*f3 + f1*f2*(1-f3)
                                        if np.abs(f123) < 1.e-12:
                                            continue
                                        res1 = w1*w2*w3*f123*self.delta(
                                                e1+e2-e3)
                                        for s1, s2, s3 in product(s_range,
                                                s_range, s_range):
                                            resi = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[i],
                                                    psi_sn1[s1,n1,:norb]\
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n2,:norb])
                                            resj = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[j],
                                                    psi_sn1[s1,n1,:norb] \
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n3,:norb])
                                            dos += res1*resi.conj()*resj
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_noncoherent_test1(self, f_skn):
        """
        Non-coherent part of dos.
        """
        dos_list = []
        # spin
        for e_kn, f_kn in zip(self.e_skn, f_skn):
            dos = np.zeros(self.npts, dtype=np.complex)

            icount = 0
            # k-point
            for w1, e_n1, f_n1 in zip(self.w_k, e_kn, f_kn):
                # k-point

                icount += 1; print(icount)

                for w2, e_n2, f_n2 in zip(self.w_k, e_kn, \
                        f_kn):
                    # k-point
                    for w3, e_n3, f_n3 in zip(self.w_k, \
                            e_kn, f_kn):
                        for n1, e1, f1 in zip(count(),
                                e_n1, f_n1):
                            for n2, e2, f2 in zip(count(),
                                    e_n2, f_n2):
                                    for n3, e3, f3 in zip(count(),
                                            e_n3, f_n3):
                                        f123 = (1-f1)*(1-f2)*f3 + f1*f2*(1-f3)
                                        if np.abs(f123) < 1.e-12:
                                            continue
                                        res1 = w1*w2*w3*f123*self.delta(
                                                e1+e2-e3)
                                        dos += res1
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_noncoherent_test2(self, psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
            i, j):
        """
        Non-coherent part of dos.
        """

        # symmetry operations
        s_range = range(psi_sksn.shape[2])
        norb = ecdfdffdelta.shape[0]

        dos_list = []
        # spin
        for e_kn, f_kn, psi_ksn in zip(self.e_skn, f_skn, psi_sksn):
            dos = np.zeros(self.npts, dtype=np.complex)

            icount = 0
            # k-point
            for w1, e_n1, f_n1, psi_sn1, nb1 in zip(self.w_k, e_kn, f_kn,\
                    psi_ksn, bnd_ne):
                # k-point

                icount += 1; print(icount)

                for w2, e_n2, f_n2, psi_sn2, nb2 in zip(self.w_k, e_kn, \
                        f_kn, psi_ksn, bnd_ne):
                    # k-point
                    for w3, e_n3, f_n3, psi_sn3, nb3 in zip(self.w_k, \
                            e_kn, f_kn, psi_ksn, bnd_ne):
                        for n1, e1, f1 in zip(count(),
                                e_n1[nb1[1]-1 : nb1[2]],
                                f_n1[nb1[1]-1 : nb1[2]]):
                            for n2, e2, f2 in zip(count(),
                                    e_n2[nb2[1]-1 : nb2[2]],
                                    f_n2[nb2[1]-1 : nb2[2]]):
                                    for n3, e3, f3 in zip(count(),
                                            e_n3[nb3[1]-1 : nb3[2]],
                                            f_n3[nb3[1]-1 : nb3[2]]):
                                        f123 = (1-f1)*(1-f2)*f3 + f1*f2*(1-f3)
                                        if np.abs(f123) < 1.e-12:
                                            continue
                                        res1 = w1*w2*w3*f123*self.delta(
                                                e1+e2-e3)
                                        for s1, s2, s3 in product(s_range,
                                                s_range, s_range):
                                            dos += res1
            dos_list.append(dos)
        return np.array(dos_list)


    def get_dos_noncoherent_test3(self, psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
            i, j):
        """
        Non-coherent part of dos.
        """

        # symmetry operations
        s_range = range(psi_sksn.shape[2])
        norb = ecdfdffdelta.shape[0]

        dos_list = []
        # spin
        for e_kn, f_kn, psi_ksn in zip(self.e_skn, f_skn, psi_sksn):
            dos = np.zeros(self.npts, dtype=np.complex)

            icount = 0
            # k-point
            for w1, e_n1, f_n1, psi_sn1, nb1 in zip(self.w_k, e_kn, f_kn,\
                    psi_ksn, bnd_ne):
                # k-point

                icount += 1; print(icount)

                for w2, e_n2, f_n2, psi_sn2, nb2 in zip(self.w_k, e_kn, \
                        f_kn, psi_ksn, bnd_ne):
                    # k-point
                    for w3, e_n3, f_n3, psi_sn3, nb3 in zip(self.w_k, \
                            e_kn, f_kn, psi_ksn, bnd_ne):
                        for n1, e1, f1 in zip(count(),
                                e_n1[nb1[1]-1 : nb1[2]],
                                f_n1[nb1[1]-1 : nb1[2]]):
                            for n2, e2, f2 in zip(count(),
                                    e_n2[nb2[1]-1 : nb2[2]],
                                    f_n2[nb2[1]-1 : nb2[2]]):
                                    for n3, e3, f3 in zip(count(),
                                            e_n3[nb3[1]-1 : nb3[2]],
                                            f_n3[nb3[1]-1 : nb3[2]]):
                                        f123 = (1-f1)*(1-f2)*f3 + f1*f2*(1-f3)
                                        if np.abs(f123) < 1.e-12:
                                            continue
                                        res1 = w1*w2*w3*f123*self.delta(
                                                e1+e2-e3)
                                        for s1, s2, s3 in product(s_range,
                                                s_range, s_range):
                                            resi = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[i],
                                                    psi_sn1[s1,n1,:norb]\
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n2,:norb])
                                            resj = np.einsum('abc,b,c,a',
                                                    ecdfdffdelta[j],
                                                    psi_sn1[s1,n1,:norb] \
                                                    .conj(),
                                                    psi_sn2[s2,n2,:norb] \
                                                    .conj(),
                                                    psi_sn3[s3,n3,:norb])
                                            dos += res1*resi.conj()*resj
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
        for n in range(bnd_ne[k, 1], bnd_ne[k, 2]+1):
            psi_skn_w_ab[:, k, n - bnd_ne[k, 1], :, :] = \
                    np.einsum('ijk, ijl->ikl',
                    psi_sksn[:, k, :, n - bnd_ne[k, 1], :],
                    np.conj(
                    psi_sksn[:, k, :, n - bnd_ne[k, 1], :]))
    return psi_skn_w_ab / psi_sksn.shape[2]



if __name__ == "__main__":
    '''
    Test run.
    '''
    import h5py

    with h5py.File('glog.h5', 'r') as f:
        e_fermi = f["/E_FERMI"][0]
        bnd_ne = f["/BND_NE"][...]
        e_skn = (f["/BND_EK"][...] - e_fermi)
        f_skn = f["/BND_FERWE"][...]
        psi_sksn = f["/BND_VK"][...]
        w_k = f["/KPT_WT"][...]
        R = f["/Impurity_1/GA_R"][...].T
        ecdfdffdelta = f["/Impurity_1/ecdfdffdelta"][...]

    print(' number of electron = {}'.format(np.sum(f_skn)))
    # different convention of f_skn
    for ik, w in enumerate(w_k):
        f_skn[:,ik,:] /= w

    print(' fermi weight: \n', f_skn)
    print(' number of electron = {}'.format(np.einsum('ijk, j->...', \
            f_skn, w_k)))

    psi_skn_w_ab = get_all_psi_skn_w_ab(e_skn, psi_sksn, bnd_ne)

    window = (-6., 6.)
    dos = DOS(w_k, e_skn,  width=0.1, window=window, npts=1001)
    energies = dos.get_energies()
    dos_t = dos.get_dos_t()

    # dos_f
    psi_sksn_f = np.einsum('...ii', psi_skn_w_ab[...,:,:])
    dos_f = dos.get_dos_component(psi_sksn_f)

    psi_skn_w_RabR = np.einsum('...ij,ik,jl', psi_skn_w_ab, np.conj(R), R)
    psi_sksn_f_qp = np.einsum('...ii', psi_skn_w_RabR[...,:,:])
    dos_f_qp = dos.get_dos_component(psi_sksn_f_qp)

    # test, symmetric
#    dos_test1 = dos.get_dos_noncoherent_test1(f_skn)

    dos_test1 =  dos.get_dos_noncoherent_test2(psi_sksn, ecdfdffdelta,
            f_skn, bnd_ne, 1, 1)

    # non-coherent
#    dos_nch = dos.get_dos_noncoherent(psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
#            0, 0)
#    dos_nch11 = dos.get_dos_noncoherent(psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
#            1, 1)
#    dos_nch22 = dos.get_dos_noncoherent(psi_sksn, ecdfdffdelta, f_skn, bnd_ne,
#            1, 1)

#    import matplotlib
#    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure(0)
    plt.plot(energies, dos_t[0].real, '-', label='t')
#    plt.plot(energies, dos_f_qp[0].real, 'o', label='qp')
#    plt.plot(energies, dos_nch[0].real, 'o', label='nch')
    plt.legend()
    plt.xlabel("E (eV)")
    plt.ylabel("DOS (states/f.u.)")

#    plt.figure(1)
#    plt.plot(energies, dos_nch[0].real, 'o', label='nch')
#    plt.plot(energies, dos_nch11[0].real, 'o', label='nch')
#    plt.plot(energies, dos_nch22[0].real, 'o', label='nch')
#    plt.legend()
#    plt.xlabel("E (eV)")
#    plt.ylabel("DOS (states/f.u.)")

    plt.figure(2)
    plt.plot(energies, dos_test1[0].real, 'o', label='test1')
    plt.legend()
    plt.xlabel("E (eV)")
    plt.ylabel("DOS (states/f.u.)")


    plt.show()

#    plt.savefig('dos.pdf')
