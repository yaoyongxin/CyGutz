# tranform the GPARAM.h5 file to local basis of
# relativistic spherical harmonics.
import h5py,numpy,os
from pyglib.math.matrix_util import unitary_transform_coulomb_matrix
from pyglib.symm.angular_momentum_1p import get_JU_relat_sph_harm_cg



def trans_gparam(u_list):
    with h5py.File("GPARAM.h5", "a") as f:
        for i,u in enumerate(u_list):
            u_db2sab = f["/IMPURITY_{}/DB_TO_SAB".format(i+1)][()].T
            u_db2sab = u_db2sab.dot(u)
            f["/IMPURITY_{}/DB_TO_SAB".format(i+1)][...]=u_db2sab.T
            hs_list = f["/IMPURITY_{}/HS".format(i+1)][()].swapaxes(1,2)
            hs_list = [u.T.conj().dot(hs).dot(u) for hs in hs_list]
            f["/IMPURITY_{}/HS".format(i+1)][...]=numpy.asarray( \
                    hs_list).swapaxes(1,2)
            for entry in ["LX", "LY", "LZ", "SX", "SY", "SZ"]:
                op = f["/IMPURITY_{}/{}".format(i+1,entry)][()].T
                op = u.T.conj().dot(op).dot(u)
                f["/IMPURITY_{}/{}".format(i+1,entry)][...] = op.T
            # to be safe, set self-energy form to be random.
            f["/IMPURITY_{}/SIGMA_STRUCT".format(i+1)][...] = \
                    (numpy.arange(len(u)**2)+1).reshape(len(u),len(u))
            v2e = f["/IMPURITY_{}/V2E".format(i+1)][()].T
            v2e = unitary_transform_coulomb_matrix(v2e, u)
            f["/IMPURITY_{}/V2E".format(i+1)][...]=v2e.T
    if os.path.isfile("WH_RL_INP.h5"):
        with h5py.File("WH_RL_INP.h5", "a") as f:
            for i,u in enumerate(u_list):
                for entry in ["LAMBDA", "R"]:
                    a = f["/IMPURITY_{}/{}".format(i+1, entry)][()].T
                    a = u.T.conj().dot(a).dot(u)
                    f["/IMPURITY_{}/{}".format(i+1, entry)][...] = a.T



def trans_gparam_rsh():
    '''transform GPARAM.h5 to relativisitic spherical harmonics basis.
    '''
    l_from_dimms = {2:0, 6:1, 10:2, 14:3}
    with h5py.File("GPARAM.h5", "a") as f:
        nimp = f["/num_imp"][0]
        u_list = []
        for i in range(nimp):
            u = f["/IMPURITY_{}/DB_TO_SAB".format(i+1)][()].conj()
            l_list = [l_from_dimms[len(u)]]
            _, u_csh2rsh = get_JU_relat_sph_harm_cg(l_list)
            u = u.dot(u_csh2rsh)
            u_list.append(u)
    trans_gparam(u_list)



if __name__ == "__main__":
    trans_gparam_rsh()
