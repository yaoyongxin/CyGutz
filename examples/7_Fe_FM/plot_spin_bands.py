'''
Band structure plot for CyGutz by Yongxin Yao.
'''
import numpy as np
import h5py

def get_letter_special(kn):
  '''
  Return label of possible Greek letters.
  '''
  if len(kn) == 1:
    return kn
  elif kn == 'LAMBDA':
    return "$\Lambda$"
  elif kn == 'GAMMA':
    return "$\Gamma$"
  elif kn == 'DELTA':
    return "$\Delta$"
  elif kn == 'SIGMA':
    return "$\Sigma$"
  else:
    return "Unknown"

# Open the main output metadata file of CyGutz
f = h5py.File("glog.h5", 'r')

# Read in the array of x,y,z components of the k-points.
kx = f["/KPT_X"][...]
ky = f["/KPT_Y"][...]
kz = f["/KPT_Z"][...]
kname = f["/KPT_NAME"][...]

# Band indices for contructing the local projector.
# ne[ik,0]: total number of bands
# ne[ik,1:2]: bottom/top bands used to contruct the local projector
ne = f["/BND_NE"][...]

# Minimal index of the top bands included to expand the local correlated orbitals.
nbmax_plot = np.min(ne[:,0])

# expansion coefficients <\psi | f>, (spin, k-points, sym_ops, bands, orbitals), band index shifted by ne(:,1)
psia = f["/BND_VK"]

# k-path length
x_list = [0]
for i in range(1,len(kx)):
  x_list.append(x_list[-1]+np.sqrt((kx[i]-kx[i-1])**2+(ky[i]-ky[i-1])**2+(kz[i]-kz[i-1])**2))

# Convert to eV.
eks = f["/BND_EK"][...]*13.605698066
e_fermi = f["/E_FERMI"][0]*13.605698066

# Ticks for high symmetry point label
x_tick = []; label_tick = []
for i,kn in enumerate(kname):
  if kn == '':
    continue
  x_tick.append(x_list[i])
  label_tick.append(get_letter_special(kn))

# Local correlated orbital weight of each band.
band_f_wt_dn = []
band_f_wt_up = []
for i in range(nbmax_plot):
  f_wt_dn = []
  f_wt_up = []
  for k in range(eks.shape[1]):
    i_ = i - ne[k,1]+1
    if i_ >= 0 and i < ne[k,2]:
      # iSpin = 0, iSymOp = 0
      f_wt_dn.append(np.sum(psia[0,k,0,i_,:]*np.conj(psia[0,k,0,i_,:])).real)
      # iSpin = 1, iSymOp = 0
      f_wt_up.append(np.sum(psia[1,k,0,i_,:]*np.conj(psia[1,k,0,i_,:])).real)
    else:
      f_wt_dn.append(0)
      f_wt_up.append(0)
  band_f_wt_dn.append(f_wt_dn)
  band_f_wt_up.append(f_wt_up)

# Shift energy
eks -= e_fermi
e_fermi = 0

# Plot band structure
import matplotlib.pyplot as plt
plt.figure(0)

# Bands with character
band_f_wt_dn = np.array(band_f_wt_dn)*20
band_f_wt_up = np.array(band_f_wt_up)*20

for n in range(nbmax_plot):
  plt.plot(x_list, eks[0,:,n],'k-')
  plt.scatter(x_list, eks[0,:,n], s = band_f_wt_dn[n], c = 'r',
          edgecolors = 'r')
  plt.plot(x_list, eks[1,:,n],'k--')
  plt.scatter(x_list, eks[1,:,n], s = band_f_wt_up[n], c = 'b',
          edgecolors = 'b')

# Fermi level
plt.axhline(y = e_fermi, ls = ':', lw = 2)
# High-symmetry lines and labels
for x1 in x_tick:
  plt.axvline(x = x1, ls = '--')
plt.xticks(x_tick, label_tick)
plt.ylabel("E (eV)")
# Set range
plt.xlim(x_list[0], x_list[-1])

# Energy range for band structure plot.
emin = -9; emax = 10
plt.ylim(emin, emax)
plt.title("Band structure with $Fe-3d$ character")

plt.savefig("FeBands.png")

#plt.show()

