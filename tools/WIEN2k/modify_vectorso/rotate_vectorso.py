import itertools as it
import h5py


def rotate_vectorso(f_h5, psi_kn, bnd_uk, bnd_ne, bnd_ek):
    with h5py.File(f_h5, 'a') as f:
        for k, psi, uk, ne, ek in it.izip(it.count(), psi_kn, bnd_uk,
                bnd_ne, bnd_ek):
            a = f['/IKP_'+str(k+1)+'/A'][...].T
            n = ne[2]-ne[1]+1
            a[:,ne[1]-1:ne[2]] = a[:,ne[1]-1:ne[2]].dot(
                    uk[:n,:n].T).dot(psi[:n,:n].T)
            try:
                f['/IKP_'+str(k+1)+'/AP'] = a.T
            except:
                f['/IKP_'+str(k+1)+'/AP'][...] = a.T
            try:
                f['/IKP_'+str(k+1)+'/eep'] = ek[:ne[0]]
            except:
                f['/IKP_'+str(k+1)+'/eep'][...] = ek[:ne[0]]

from subprocess import call
call(["/home/ykent/WIEN_GUTZ/bin/tools/WIEN2k/modify_vectorso/vectorso2hdf5", \
        "Ce"])

with h5py.File('glog.h5', 'r') as f:
    sym_ie = f['/SYM_IE'][0]
    psi_kn = f['/BND_VK'][0,:,sym_ie-1,:,:]
    bnd_uk = f['/BND_UK0'][:,sym_ie-1,:,:]
    bnd_ne = f['/BND_NE'][...]
    bnd_ek = f['/BND_EK'][0,:,:]

rotate_vectorso('Ce_vectorso.h5', psi_kn, bnd_uk, bnd_ne, bnd_ek)
rotate_vectorso('Ce_vectorsodn.h5', psi_kn, bnd_uk, bnd_ne, bnd_ek)

call(["/home/ykent/WIEN_GUTZ/bin/tools/WIEN2k/modify_vectorso/hdf52vectorso", \
        "Ce"])

