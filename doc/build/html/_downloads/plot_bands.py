import numpy as np
import h5py

f = h5py.File("glog.h5", 'r')
kx = f["/KPT_X"][...]
ky = f["/KPT_Y"][...]
kz = f["/KPT_Z"][...]

# band indices for contructing the local projector
# ne[ik,0]: total number of bands; ne[ik,1:2]: bottom/top bands used to contruct the local projector
ne = f["/BND_NE"][...]

nmax = np.min(ne[:,2])
psia = f["/BND_VK"] # (spin, k-points, sym_ops, bands, orbitals), band index shifted by ne(:,1)

# k-path length
x_list = [0]
for i in range(1,len(kx)):
  x_list.append(x_list[-1]+np.sqrt((kx[i]-kx[i-1])**2+(ky[i]-ky[i-1])**2+(kz[i]-kz[i-1])**2))

eks = (f["/BND_EK"][...] - f["/E_FERMI"])*13.605698066 # Rydberg to eV

#for k in range(eks.shape[1]):
#  idx = eks[0,k,:].argsort()[::-1]
#  eks[0,k,:] = eks[0,k,idx]
#

emin = np.min(eks); emax = np.max(eks)

with open('BANDS.dat', 'w') as f:
  for i in range(nmax):
    for k in range(eks.shape[1]):
      i_ = i - ne[k,1]+1
      if i_ >= 0 and i < ne[k,2]:
        wt = np.sum(psia[0,k,0,i_,:]*np.conj(psia[0,k,0,i_,:])).real
      else:
        wt = 0
      print >> f, x_list[k], '   ',eks[0,k,i].real,wt # k-length, band energy, weight
    print >> f, '  \n'
  for i in [0,20,40,60,70,80,90,100]:
    print >> f, x_list[i], emin, 0
    print >> f, x_list[i], emax, 0
    print >> f, '  \n'

