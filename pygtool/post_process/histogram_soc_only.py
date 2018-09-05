#!/usr/bin/env python

'''Script works for both metallic phase and Mott phase with
spin-orbit interaction ONLY.'''

import os,sys

# impurity index
imp = 1

if not os.path.isfile('EMBED_HAMIL_ANALYSIS_{}.h5'.format(imp)):
    print("please run the proper multiplet analysis code first!")
    sys.exit()

if "-e" in sys.argv:
    num_ev = int(sys.argv[sys.argv.index("-e")+1])
else:
    num_ev = 50

if "-l" in sys.argv:
    num_label = int(sys.argv[sys.argv.index("-l")+1])
else:
    num_label = 5


from pyglib.mbody.multiplets_analysis_soc import calc_save_atomic_states, \
        plot_atomic_states

# Calculate local histograms for impurity 1
# Check the 100 dominant eigen-states.
calc_save_atomic_states(imp=imp, num_ev=num_ev)

# Plot the histogram
plot_atomic_states(imp=imp, num_label=num_label)
