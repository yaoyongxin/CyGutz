#!/usr/bin/env python

log=open('ga_vdc_solver.log','w')
from sys import argv
script, nmin, nmax = argv
print>>log, nmin, nmax
log.close()

# Some script for "stable" solver based on fixing V_dc trick
# I will re-write it in a nice way, but currently let it work
import ast
import os
from scipy import optimize
import ga_utils
import sys
import getopt

nf=[]
f=open('GL.INP','r')
for line in f:
    if 'GL%LDC' in line:
        if '2' not in line.split()[0]:
            print>>log, ' Please set GL%LDC=2 for this stable Vdc solver!'
            quit()
f.close()

def fun(nf):
    log=open('ga_vdc_solver.log','a')
    print>>log, 'nf_in=',nf
    GLINP=''
    f=open('GL.INP','r')
    for line in f:
        if 'OCC' in line:
            line_=line.split()
            line_[2]=str(nf)
            line='  '.join(line_)+'\n'
        GLINP+=line
    f.close()
    f=open('GL.INP','w')
    print >> f, GLINP
    f.close()
    dmfe = ga_utils.DmftEnvironment()  # DMFT paths
    w2k = ga_utils.W2kEnvironment()    # W2k filenames and paths

    os.system(dmfe.ROOT+"/ga_run_dmft.py >& out.ga_run")
    fname=w2k.case+'.outputdmf1.0'
    f=open(fname,'r')
    lread=0
    for line in f:
        if lread == 1 :
            nf_out=float(line.split()[0])
            lread=0
        if 'NELE_LOC TOTAL' in line:
            lread=1
    f.close()
    diff=nf_out-nf
    print>>log, 'nf_out=',nf_out,'fun=',diff
    fname=w2k.case+'.scf'
    f=open(fname,'r')
    for line in reversed(f.readlines()):
        if ':ENE' in line:
            energy=line.split()[8]
            break
    f.close()
    print>>log, 'energy=',energy
    log.close()
    return diff

#sol = optimize.fsolve(fun, nf, xtol=1.e-4, epsfcn=1.e-5)
#print sol.x
sol=optimize.brentq(fun, float(nmax),float(nmin),xtol=1.e-6,rtol=1.e-6)
