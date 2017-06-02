#!/usr/bin/env python
import h5py
import numpy as np
from matrix_util import matrix_basis_expand


def get_BS_matrix(B_vector, S_vector):
    '''
    Get B_vec . S_vec.
    '''
    v = np.zeros_like(S_vector[0])
    for i, B in enumerate(B_vector):
        v += B * S_vector[i]
    return v


def get_B_potential_matrix(B_mag):
    '''
    Generate the local potential matrix due to the local magnetic field.
    '''
    f = h5py.File("init_ga_info.h5", 'r')
    v_list = []
    num_impurity = f["/num_impurity"][...]
    for imp in range(num_impurity):
        matrix_basis = f["/impurity_" + str(imp) + "/matrix_basis"]
        S_vector = f["/impurity_" + str(imp) + "/S_vector"][...]
        print "\n Impurity %2d"%(imp)
        # Get local B-field direction.
        while True:
            answer = raw_input(
                    " Please enter the local spin magnitization" +
                    " direction \n in global coordinate system" +
                    " (e.g., 0 0 1):  ")
            try:
                Mvec = [float(v) for v in answer.split()[:3]]
                break
            except:
                pass
        Mvec = np.array(Mvec)
        Mvec = Mvec / np.linalg.norm(Mvec)
        B_vector = Mvec * B_mag
        v = get_BS_matrix(B_vector, S_vector)
        vp = matrix_basis_expand(v, matrix_basis)
        if np.max(np.abs(vp - v)) > 1.e-6:
            np.set_printoptions(precision=4, suppress=True, linewidth=150)
            print 'v_original_real:'
            print v.real
            print 'v_original_imag:'
            print v.imag
            print 'v_expanded_real:'
            print vp.real
            print 'V_expanded_imag:'
            print vp.imag
            raise AssertionError("Matrix basis is not complete!")
        v_list.append(vp)
    f.close()
    return v_list


def write_B_potential_matrix(v_list):
    '''
    Write B_potential_matrix file WH_Bvec_Svec.INP.
    '''
    from fileio import write_slice
    dim_list = [v.shape[0] for v in v_list]
    dim_max = np.max(dim_list)
    with open("WH_Bvec_Svec.INP", 'w') as f:
        for v in v_list:
            mat = np.zeros((dim_max, dim_max), dtype=np.complex)
            dim = v.shape[0]
            mat[:dim, :dim] = v
            write_slice(f, v)
    f.close()


if __name__ == "__main__":
    # Enter B-field strength
    while True:
        answer = raw_input(
                " Please enter the magnetic field strength "+
                "\n (0.1 introduces splitting of ~1 eV)...")
        try:
            B_mag = float(answer)
            break
        except:
            pass

    v_list = get_B_potential_matrix(B_mag)
    write_B_potential_matrix(v_list)
    print ' Apply_initial_B_field completed.'
