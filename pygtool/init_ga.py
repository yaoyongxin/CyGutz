#!/usr/bin/env python

import sys
from pyglib.gutz.init import initialize
from pyglib.iface.ifwien import h4set_indmfl

# Energy window which defines the band wave functions used to expand the
# local projector.

emin = -10.
emax = 10.

# Projector type. The default one should be good.
projtype = -2

if '-emin' in sys.argv:
    emin = float(sys.argv[sys.argv.index('-emin') + 1])
if '-emax' in sys.argv:
    emax = float(sys.argv[sys.argv.index('-emax') + 1])
if '-projtype' in sys.argv:
    projtype = int(sys.argv[sys.argv.index('-projtype') + 1])

initialize()
h4set_indmfl()
