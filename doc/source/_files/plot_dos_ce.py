from pyglib.estructure.dos import h5get_dos, plot_dos_tf

# get dos
energies, dos_t, dos_f = h5get_dos()
# pick up-component to plot
plot_dos_tf(energies, dos_t[0], dos_f[0])
