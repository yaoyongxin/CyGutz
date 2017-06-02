'''
Multiplets analysis for LDA+G+SOC calculations.
'''

from multiplets_analysis_soc import calc_save_atomic_states, \
        plot_atomic_states
from multiplets_analysis_lib import generate_symmetrized_rho_by_hmexpand_f


def multiplets_analysis_soc_mott(imp=1, num_ev=400, rho_name='/RHO',
        num_label = 5):
    generate_symmetrized_rho_by_hmexpand_f(imp=imp, rpath=rho_name)
    calc_save_atomic_states(imp=imp, num_ev=num_ev, rho_name=rho_name+'_SYM')
    plot_atomic_states(imp=imp, num_label = num_label)


if __name__=='__main__':
    multiplets_analysis_soc_mott(imp=1)
