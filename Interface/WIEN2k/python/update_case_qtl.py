#!/usr/bin/env python

import os
import shutil
import ga_utils

w2k = ga_utils.W2kEnvironment()    # W2k filenames and paths
f = os.popen('grep "GA FERMI LEVEL" ./'+w2k.case+'.outputdmf1.0')
lines = f.readlines(); line=lines[-1]
ef_ga = float(line[16:33])
qtl_file = []

shutil.copy2('./'+w2k.case+'.qtl', 'W2KQTL.OUT')
with open('./W2KQTL.OUT','r') as f:
  for line in f.read().split('\n'):
    if 'FERMI ENERGY' in line:
      line=line[0:56]+'%11.6f'%(ef_ga)
    elif 'BAND' in line:
      break
    qtl_file.append(line)
with open('./GLQTL.OUT','r') as f:
  for line in f.read().split('\n'):
    qtl_file.append(line)
with open('./'+w2k.case+'.qtl','w') as f:
  for line in qtl_file:
    print >>f, line
