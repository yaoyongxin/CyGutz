#!/usr/bin/env python

import os

f = os.popen('grep "FERMI-LDA" ./GL_LOG.OUT')
lines = f.readlines(); line=lines[-1]
ef_lda = float(line[11:21])
with open('EFLDA.INP','w') as f:
  print >>f, '%12.6f ! EF_LDA'%(ef_lda)
