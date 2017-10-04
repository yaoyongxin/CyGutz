######################################
# some double counting potential form.
######################################
import numpy


def get_vdc_hf(v2e, dm):
    '''get Hartree-Fock type double counting potential,
    given the full coulomb matrix (chemist convention, including spin index)
    and the one-particle density matrix.
    '''
    vdc = -numpy.einsum('ijkl,kl->ij', v2e, dm)
    vdc += numpy.einsum('ikjl,kl->ij', v2e, dm)
    return vdc
