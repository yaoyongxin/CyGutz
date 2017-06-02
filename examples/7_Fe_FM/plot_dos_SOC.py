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
  Reduce psi_sksn to psi_skn_w by summing over symmetry operation indiex.
  '''
  psi_skn_w = np.zeros(list(e_skn.shape)+[psi_sksn.shape[-1]], dtype=float)
  for isp in range(e_skn.shape[0]):
    for k in range(e_skn.shape[1]):
      for isy in range(psi_sksn.shape[2]):
        for n in range(bnd_ne[k,1],bnd_ne[k,2]+1):
          for a in range(psi_sksn.shape[-1]):
            psi_skn_w[isp,k,n-1,a] += (psi_sksn[isp,k,isy,n-bnd_ne[k,1],a] \
                            *np.conj(psi_sksn[isp,k,isy,n-bnd_ne[k,1],a])).real
  return psi_skn_w/psi_sksn.shape[2]

if __name__=="__main__":
  '''
  Test run.
  '''
  import h5py
  Ryd_eV = 13.605698066
  f = h5py.File("glog.h5", 'r')
  # Fermi level
  e_fermi = f["/E_FERMI"][0]*Ryd_eV
  # Band energy to eV
  e_skn = f["/BND_EK"][...]*Ryd_eV
  # k-ponits weight
  w_k = f["/KPT_WT"][...]
  # Band indices for contructing the local projector.
  # bnd_ne[ik,0]: total number of bands
  # bnd_ne[ik,1:2]: bottom/top bands used to contruct the local projector
  # bnd_ne is one-based.
  bnd_ne = f["/BND_NE"][...]
  # expansion coefficients <\psi | f>, (spin, k-points, sym_ops, bands, orbitals), band index shifted by ne(:,1)
  psi_sksn = f["/BND_VK"][...]
  f.close()

  psi_skn_w = get_all_psi_skn_w(e_skn, psi_sksn, bnd_ne)
  nbmax_plot = np.min(bnd_ne[:,2])

  # Energy window
  emin = e_fermi - 8
  emax = e_fermi + 4
  window = (emin, emax)
  dos = DOS(w_k, e_skn,  width=0.05, window=window, npts=5001)
  energies = dos.get_energies()
  dos_t = dos.get_dos()
  psi_sksn_f = np.zeros(e_skn.shape, dtype=float)
  for a in range(10):
    psi_sksn_f += psi_skn_w[:,:,:,a]
  dos_f = dos.get_dos(psi_sksn_f)

  # store metadata for the plot.
  f = h5py.File("dos.h5", 'w')
  f["/energies"] = energies
  f["/dos_tot"] = dos_t[0]
  f["/dos_f"] = dos_f[0]
  f.close()

  energies -= e_fermi
  e_fermi = 0

  import matplotlib.pyplot as plt
  plt.figure(0)
  plt.fill_between(energies, 0, dos_t[0], facecolor = 'grey', alpha = 0.5, label = 'tot')
  plt.plot(energies, dos_t[0], color = 'grey', label = 'tot')
  plt.plot(energies, dos_f[0], color = 'red', label = '$Fe-3d$')
  plt.axvline(x = e_fermi, ls = '--')
  plt.title("Density of States with $Fe-3d$ component")
  plt.xlabel("E (eV)")
  plt.ylabel("DOS (states/f.u.)")
  plt.legend()
  plt.savefig("FeDOS_SOC.png")
#  plt.show()

