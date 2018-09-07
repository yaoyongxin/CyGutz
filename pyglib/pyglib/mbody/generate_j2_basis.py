# generate the j_square basis for each valence block f-shell.
from pyglib.mbody.local_operator_factory import get_nij_op, get_am_op_from_nij
from pyglib.symm.angular_momentum_1p import get_J_vector
from pyglib.symm.unitary import comp_sph_harm_to_relativistic_harm
import numpy, h5py


def h5generate_j2_basis(l="f"):
    '''generate j2 eigen-states with relativistic harmonics as the
    single-particle basis.
    '''
    norbs = {"d":10, "f":14}
    l_list = [(norbs[l]/2-1)/2]
    jvec = get_J_vector(l_list)
    ucsh2rh = comp_sph_harm_to_relativistic_harm(norbs[l])
    jvec = [ucsh2rh.T.conj().dot(jcomp).dot(ucsh2rh) for jcomp in jvec]
    with h5py.File("j2evec.h5", "a") as f:
        group = "/{}".format(l)
        if group in f:
            del f[group]
        for val in range(1,norbs[l]):
            nij = get_nij_op(l=l, ival=val)
            jop_list = get_am_op_from_nij(jvec=jvec, nij=nij,
                    op_list=['Jz', 'J2'])
            j2op = jop_list["J2"].todense()
            jzop = jop_list["Jz"].todense()
            w,v = numpy.linalg.eigh(j2op)
            j_list = map(lambda x: abs(round(numpy.sqrt(x+0.25)-0.5, 1)), w)
            j_unique = numpy.unique(j_list)
            for j in j_unique:
                istr = j_list.index(j)
                iend = len(j_list) - j_list[::-1].index(j)
                v_j = v[:,istr:iend]
                jzop_j = v_j.T.conj().dot(jzop).dot(v_j)
                _w,_v = numpy.linalg.eigh(jzop_j)
                v_j = v_j.dot(_v)
                f["{}/val_{}/j_{}/V".format(group,val,j)] = v_j



if __name__ == "__main__":
    h5generate_j2_basis(l="f")
