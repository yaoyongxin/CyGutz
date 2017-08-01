#!/usr/bin/env python

import sys
from pyglib.gutz.init import initialize
from pyglib.iface.ifwien import h4set_indmfl


initialize()
if '-no-indmfl' not in sys.argv or '-vasp' in sys.argv:
    h4set_indmfl()
