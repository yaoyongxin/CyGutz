#!/usr/bin/env python

import ast
import os
from scipy import optimize
import ga_utils
import sys
import getopt

globiter=0

def fun(xin):
    global globiter
    globiter+=1
    print 'Iter=', globiter
    print 'xin=',xin
    # Set up GL.INP for Charge scf with fixed {R,LA1}
    GLINP=''
    fname='GL.INP'
    with open(fname, 'r') as f:
        for line in f:
            if 'GL%LSCF' in line:
                line='2           ! GL%LSCF : -1->HF; 1->GA'+'\n'
            GLINP+=line
    with open(fname, 'w') as f:
        print >> f, GLINP
    fname='GLX.INP'
    with open(fname, 'w') as f:
        f.writelines('  '.join(str(xin[i]) for i in range(len(xin)/2))+'\n')
        f.writelines('  '.join(str(xin[i]) for i in range(len(xin)/2,len(xin))))
    # rho scf with fixed x
    dmfe = ga_utils.DmftEnvironment()  # DMFT paths
    w2k = ga_utils.W2kEnvironment()    # W2k filenames and paths
    os.system(dmfe.ROOT+"/ga_run_dmft.py >> out.ga_rhoscf_run")
    # One-shot dmft1 to get xout
    os.system("sed -i 's/2           ! GL%LSCF/3           ! GL%LSCF/g' GL.INP")
    os.system(dmfe.ROOT+"/dmft dmft1.def >> dmft1_info.out")
    diff=[]
    fname=w2k.case+'.outputdmf1.0'
    with open(fname, 'r') as f:
        ldiff=0
        for line in f:
            if 'DIF_X' in line:
                ldiff=1
            elif ldiff==1:
                diff+=[float(i) for i in line.split()]; ldiff+=1
            elif ldiff==2:
                diff+=[float(i) for i in line.split()]; ldiff=0
    print 'diff=',diff
    return diff


# get initial guess xin
xin=[]
fname='GLX.INP'
with open(fname, 'r') as f:
    for line in f:
        line=[float(i) for i in line.split()]
        xin+=line
print 'xin_start=',xin
#sol = optimize.fsolve(fun, xin, xtol=1.e-4, epsfcn=1.e-6)
sol = optimize.fsolve(fun, xin, xtol=1.e-6,epsfcn=1.e-6)
print sol
