from math import pi, sqrt
import numpy as np
# Based on ase.dft.dos.

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
    """Return a delta-function centered at 'energy'."""
    x = -((self.energies - energy) / self.width)**2
    return np.exp(x) / (sqrt(pi) * self.width)

  def get_energies(self):
    """"return energy mesh."""
    return self.energies

  def get_dos(self, psi_skn_w=None):
    """Get array of DOS values."""

    dos_list = []
    for isp, e_skn in enumerate(self.e_skn):
      dos = np.zeros(self.npts)
      for k, w in enumerate(self.w_k):
        for n, e in enumerate(e_skn[k]):
          if n>0 and e<e_skn[k][n-1]:
            break
          res = w * self.delta(e)
          if psi_skn_w is None:
            dos += res
          else:
            dos += res*psi_skn_w[isp,k,n]
      dos_list.append(dos)
    return np.array(dos_list)

def get_all_psi_skn_w(e_skn, psi_sksn, bnd_ne):
  '''
  Reduce psi_sksn to psi_skn_w.
  '''
  psi_skn_w = np.zeros(list(e_skn.shape)+[psi_sksn.shape[-1]], dtype=float)
  for isp in range(e_skn.shape[0]):
    for k in range(e_skn.shape[1]):
      for isy in range(psi_sksn.shape[2]):
        for n in range(bnd_ne[k,1],bnd_ne[k,2]):
          for a in range(psi_sksn.shape[-1]):
            psi_skn_w[isp,k,n,a] += (psi_sksn[isp,k,isy,n-bnd_ne[k,1],a] \
                            *np.conj(psi_sksn[isp,k,isy,n-bnd_ne[k,1],a])).real
  return psi_skn_w/psi_sksn.shape[2]

if __name__=="__main__":
  '''
  Test run.
  '''
  import h5py
  from unitconv import Ryd_eV
  f = h5py.File("glog.h5", 'r')
  e_fermi = f["/E_FERMI"][0]
  e_skn = (f["/BND_EK"][...]-e_fermi)*Ryd_eV
  w_k = f["/KPT_WT"][...]
  bnd_ne = f["/BND_NE"][...]
  psi_sksn = f["/BND_VK"][...]
  f.close()

  psi_skn_w = get_all_psi_skn_w(e_skn, psi_sksn, bnd_ne)
  window = (-20., 20.)
  dos = DOS(w_k, e_skn,  width=0.05, window=window, npts=1001)
  energies = dos.get_energies()
  dos_t = dos.get_dos()
  psi_sksn_f = np.zeros(e_skn.shape, dtype=float)
  for a in range(14):
    psi_sksn_f += psi_skn_w[:,:,:,a]
  dos_f = dos.get_dos(psi_sksn_f)

  f = h5py.File("dos.h5", 'w')
  f["/energies"] = energies
  f["/dos_tot"] = dos_t[0]
  f["/dos_f"] = dos_f[0]
  f.close()

  import matplotlib.pyplot as plt
  plt.figure(0)
  plt.plot(energies, dos_t[0])
  plt.xlabel("E (eV)")
  plt.ylabel("DOS (states/f.u.)")
  plt.figure(1)
  plt.plot(energies, dos_f[0])
  plt.xlabel("E (eV)")
  plt.ylabel("DOS (states/f.u.)")
  plt.show()
