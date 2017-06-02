import h5py
import numpy as np
from scipy.sparse import csr_matrix
from dataproc import get_csr_matrix, get_ed_lowest_vec, trace_single_state

def get_ed_d_occ_1():
    '''
    Get double occupancy of the first (correlated) site based on the exact diagonalization (ED) of the many-body Hamiltonian of a cluster.
    '''
    f = h5py.File("glog.h5", "r")
    h = get_csr_matrix(f, "/Impurity_1/H")
    state = get_ed_lowest_vec(h)
    c_up = get_csr_matrix(f, "/Impurity_1/annihi.op._1")
    c_dn = get_csr_matrix(f, "/Impurity_1/annihi.op._2")
    n_updn = c_up.getH() * c_up * c_dn.getH() * c_dn
    n_updn = n_updn[n_updn.shape[0] - h.shape[0] : n_updn.shape[0],
        n_updn.shape[1] - h.shape[1] : n_updn.shape[1]]
    f.close()
    return trace_single_state(state, n_updn)

def get_n12_ed():
    '''
    Get <c^{\dagger}_{1, up} c_{2, up}>.
    '''
    f = h5py.File("glog.h5", "r")
    h = get_csr_matrix(f, "/Impurity_1/H")
    state = get_ed_lowest_vec(h)
    c1 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(1)) # 0 to 1-based
    c2 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(3))
    f.close()
    c1_dagger_c2 = c1.getH() * c2
    c1_dagger_c2 = c1_dagger_c2[c1_dagger_c2.shape[0] - h.shape[0] : c1_dagger_c2.shape[0],
                   c1_dagger_c2.shape[1] - h.shape[1] : c1_dagger_c2.shape[1]]
    return trace_single_state(state, c1_dagger_c2)

def get_w0_sb(num_spin_orbital):
    '''
    Get W^{0}_{R, \Psi}.
    '''
    f = h5py.File("glog.h5", "r")
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    w0 = []
    for i in range(0, num_spin_orbital, 2):
        c_up = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(i + 1)) # 0 to 1-based
        c_dn = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(i + 2))
        n_up = c_up.getH() * c_up
        n_up = n_up[n_up.shape[0] - rho.shape[0] : n_up.shape[0],
                    n_up.shape[1] - rho.shape[1] : n_up.shape[1]]
        n_dn = c_dn.getH() * c_dn
        n_dn = n_dn[n_dn.shape[0] - rho.shape[0] : n_dn.shape[0],
                    n_dn.shape[1] - rho.shape[1] : n_dn.shape[1]]
        exp_val_n_up = np.sum((rho * n_up).diagonal())
        exp_val_n_dn = np.sum((rho * n_dn).diagonal())
        exp_val_n_up_n_dn = np.sum((rho * n_up * n_dn).diagonal())
        w0.append(np.abs(exp_val_n_up_n_dn - exp_val_n_up*exp_val_n_dn))
    f.close()
    return w0

def get_w1_sb(num_spin_orbital):
    '''
    Get W^{1}_{R, \Psi}.
    '''
    f = h5py.File("glog.h5", "r")
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    w1 = []
    for i in range(0, num_spin_orbital - 2, 2):
        c1 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(i + 1)) # 0 to 1-based
        c2 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(i + 3))
        n1 = c1.getH() * c1
        n1 = n1[n1.shape[0] - rho.shape[0] : n1.shape[0],
                      n1.shape[1] - rho.shape[1] : n1.shape[1]]
        n2 = c2.getH() * c2
        n2 = n2[n2.shape[0] - rho.shape[0] : n2.shape[0],
                      n2.shape[1] - rho.shape[1] : n2.shape[1]]
        c1_dagger_c2 = c1.getH() * c2
        c1_dagger_c2 = c1_dagger_c2[c1_dagger_c2.shape[0] - rho.shape[0] : c1_dagger_c2.shape[0],
                                                c1_dagger_c2.shape[1] - rho.shape[1] : c1_dagger_c2.shape[1]]
        c1_c2_dagger = c1 * c2.getH()
        c1_c2_dagger = c1_c2_dagger[c1_c2_dagger.shape[0] - rho.shape[0] : c1_c2_dagger.shape[0],
                                                c1_c2_dagger.shape[1] - rho.shape[1] : c1_c2_dagger.shape[1]]
        exp_val_n1 = np.sum((rho * n1).diagonal())
        exp_val_n2 = np.sum((rho * n2).diagonal())
        exp_val_n1_n2 = np.sum((rho * n1 * n2).diagonal())
        exp_val_c1_dagger_c2 = np.sum((rho * c1_dagger_c2).diagonal())
        exp_val_c1_c2_dagger = np.sum((rho * c1_c2_dagger).diagonal())
        w1.append(np.abs(exp_val_n1_n2 - (exp_val_n1*exp_val_n2 + exp_val_c1_dagger_c2*exp_val_c1_c2_dagger)))
    f.close()
    return w1

def get_n12_sb():
    '''
    Get <c^{\dagger}_{1, up} c_{2, up}>_{G}.
    '''
    f = h5py.File("glog.h5", "r")
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    c1 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(1)) # 0 to 1-based
    c2 = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(3))
    f.close()
    c1_dagger_c2 = c1.getH() * c2
    c1_dagger_c2 = c1_dagger_c2[c1_dagger_c2.shape[0] - rho.shape[0] : c1_dagger_c2.shape[0],
                                            c1_dagger_c2.shape[1] - rho.shape[1] : c1_dagger_c2.shape[1]]
    return np.sum((rho * c1_dagger_c2).diagonal())

def get_d11_sb():
    '''
    Get <n_{1, up} n_{1, dn}>_{G}--the double occupancy of first impurity first orbital in Gutzwiller-Slave-boson calculations.
    '''
    f = h5py.File("glog.h5", "r")
    rho = get_csr_matrix(f, "/Impurity_1/RHO")
    c1up = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(1)) # 0 to 1-based
    c1dn = get_csr_matrix(f, "/Impurity_1/annihi.op._" + str(2))
    f.close()
    n1up = c1up.getH() * c1up
    n1dn = c1dn.getH() * c1dn
    nbase = n1dn.shape[0]-rho.shape[0]
    return np.sum((rho * n1up[nbase:,nbase:] * n1dn[nbase:,nbase:]).diagonal())


if __name__ == "__main__":
    '''
    Test.
    '''
    print 'd_occ_1_gsb = ', get_d11_sb()
    #print 'd_occ_1 = ', get_ed_d_occ_1()
    #print 'w0 = ', get_w0_sb(10)
    #print 'w1 = ', get_w1_sb(10)
    #print '<c^{\dagger}_{1, up} c_{2, up}>_{G} = ', get_n12_sb()
