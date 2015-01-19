#!/usr/bin/env python
import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
utime = '/usr/bin/time'  # This time function is compatible with w2k timing
import ga_utils,ga_indmffile
import ga_convert_processes
import numpy
from scipy import *
import copy

class Params(object):
    """Class to store various control paramers which control the main flow of the program
       The class is constructed such that it can handle two types of parameters:
         - starting parameters which are read only at the instantiation
         - on fly parameters which can be refreshed (reread from the file) at any time
           by a member function refresh

       Input:
          f_sparams    -- filename of the parameters read only at the beginning
          sparams      -- a dictionary of pairs {name of the variable: default value} which are looked for in the f_sparams file
          f_params     -- filename of the parameters which can be refreshed several times by refresh member function
          params       -- a dictionary of {variable: default value} pairs which can be updated at runtime
       Output:
          None
       Variables which become members of the class are those listed in params and sparams dictionaries
    """
    def __init__(self, f_params, params, f_sparams=None):

        self.p ={}
        self.p.update(params) # default values from the __main__

        # remember params and its data file
        self.f_params = f_params

        # values from the 'params.dat' file
        # executes the parameters at startup
        if f_sparams is not None:
            if os.path.isfile(f_sparams):
                execfile(f_sparams)
                sp = locals()      # stores parameters
                del sp['self']     # deletes automatic locals
                self.p.update(sp)  # new or updates variables
        
    def refresh(self, fh_info):

        # executes the parameters file
        if os.path.isfile(self.f_params):
            execfile(self.f_params)
            sp = locals()   # stores locals
            del sp['self']  # deletes automatics

        # prints what was changed
            for c in sp.keys():
                if (c not in self.p.keys()) or (sp[c] != self.p[c]):
                    print >>fh_info, '********------ parameter change ------***********'
                    print >>fh_info, c, "=", sp[c]

        # updates the dictonary
            self.p.update(sp)
        
    def __getitem__(self, x):
        return self.p[x]

    def __setitem__(self, i, value):
        self.p[i]=value

    def keys(self):
        return self.p.keys()

class Struct:
    def __init__(self, filename, fh_info=sys.stdout):
        """ Reads Win2K structure file"""
        print >> fh_info, '***** Reading Structure file *****'
        fh = open(filename, 'r')

        self.title = fh.next()
        line = fh.next()
        lattic = line[:4]
        self.nat = int(line[27:30])
        print >> fh_info, 'Number of sorts, nat=', self.nat
        
        fh.next()
        line = fh.next()
        self.angle=zeros(3,dtype=float)

        (self.a,self.b,self.c,self.angle[0],self.angle[1],self.angle[2]) = map(float, [line[:10],line[10:20],line[20:30],line[30:40],line[40:50],line[50:60]])
        for i in range(3):
            if abs(self.angle[i])<1e-5: self.angle[i]=90.
        print >> fh_info, 'unit vectors (a,b,c) are ', self.a, self.b, self.c
        print >> fh_info, 'angles betwen unit vectors are', self.angle

        print >> fh_info, '---- over all types and atoms ----'
        self.iatnr=[]
        self.pos=[]
        self.mult=[]
        self.isplit=[]
        self.aname=[]
        self.jrj=[]
        self.r0=[]
        self.rmt=[]
        self.Znuc=[]
        self.RotLoc=[]
        for iat in range(self.nat):
            print >> fh_info, '---------- atom type', iat, '------------'
            line = fh.next()
            self.iatnr.append( int(line[4:8]) )
            print >> fh_info, 'iatnr=', self.iatnr[-1]
            tpos = map(float, [line[12:22], line[25:35],line[38:48]])
            self.pos.append(tpos)
            print >> fh_info, 'Position inside the unit cell, pos=', tpos
            line = fh.next()
            self.mult.append( int(line[15:17]) )
            print >> fh_info, 'Number of atoms of this type, mult=', self.mult[-1]
            self.isplit.append( int(line[34:37]) )
            print >> fh_info, 'Symmetry of the atom, isplit=', self.isplit[-1]

            if self.mult[-1]>1: print >> fh_info, '   *** Equivalent atoms ***'
            
            for mu in range(self.mult[-1]-1):
                line = fh.next()
                self.iatnr.append( int(line[4:8]) )
                print >> fh_info, '    iatnr=', self.iatnr[-1],
                self.pos.append( map(float, [line[12:22], line[25:35],line[38:48]]) )
                print >> fh_info, '    pos=', self.pos[-1]
                
            line = fh.next()
            self.aname.append( line[:10] )
            self.jrj.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( float(line[55:65]) )
            print >> fh_info, '--------- Information about the current type -------------'
            print >> fh_info, 'aname=', self.aname[-1], 'Znuc=', self.Znuc[-1], 'Rmt=', self.rmt[-1]
            rt=[]
            for i in range(3):
                line = fh.next()
                rt.append( map(float, [line[20:30], line[30:40],line[40:50]]) )
            self.RotLoc.append(rt)
            print >> fh_info, 'RotLoc=', self.RotLoc[-1]
            print >> fh_info, '-----------------------------------------------------------'

        self.aZnuc=[]
        for iat in range(self.nat):
            for mu in range(self.mult[iat]):
                self.aZnuc.append(self.Znuc[iat])
        
        print >> fh_info, 'Z=', self.aZnuc
        print >> fh_info

def dmft1(fday, case, fh_info, MPI, ROOT):
    name = 'dmft1'
    if not os.path.isfile(name+'.def'):
        print '  stop error: the required input file '+name+'.def for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+name+'.def for the next step could not be found!'
        
    tim = time.strftime("%H:%M:%S")

    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()
    
    print >> fh_info, 'Running ---- ksum -----'
    cmd = utime+' '+MPI+ ' ./dmft '+name+'.def >> dmft1_info.out '
    print >> fh_info, cmd
    fh_info.flush()
    print cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    print >> fh_info, out, err
    fh_info.flush()
    os.system('cat '+case+'.outputdmf1.0 >> dmft1_info.out')
    
    fe = open('dmft1_log','w'); print >> fe, err;  fe.close()
    print err.splitlines()[0]
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, err.splitlines()[-2])
    fday.flush()
    for fe in glob.glob('dmft1.error*'):
        if getsize(fe) !=0:
            print 'ERROR in dmft1 from file:', fe, open(fe,'r').read()
            sys.exit(1)

def CyGutz(fday, case, fh_info, MPI, ROOT):
    name = 'CyGutz'
    tim = time.strftime("%H:%M:%S")
    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()

    print >> fh_info, 'Running ---- CyGutz -----'
    cmd = utime+' '+MPI+ ' ./CyGutz '+'>> CyGutz_info.out '
    print >> fh_info, cmd
    fh_info.flush()
    print cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    print >> fh_info, out, err
    fh_info.flush()
    os.system('cat GUTZ.LOG >> '+case+'.outputdmf1.0')
    os.system('cat '+case+'.outputdmf1.0 >> dmft1_info.out')
    print err.splitlines()[0]
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, err.splitlines()[-2])
    fday.flush()

################################################
## The charge self-consistency part
################################################

def AStep(fday, case, name, inp_name, WIEN, para, fh_info):
    if not os.path.isfile(case+'.'+inp_name):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
    
    tim = time.strftime("%H:%M:%S")
    cmd = WIEN+'/x'+para+' -f '+case+' '+name
    print cmd
    print >> fh_info, (name+': '), cmd,'\n'
    fh_info.flush()
    stdin, stdout, stderr = os.popen3(cmd)
    err,out = stderr.read().strip(), stdout.read().strip()
    print err
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, out)
    fday.flush()

    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print 'ERROR in '+fe+' from file:', open(fe,'r').readlines()
            #sys.exit(1)
    
    #if os.path.exists(name+'.error.0') and os.path.getsize(name+'.error.0')>0:
    #    ferr = open(name+'.error.0', 'r')
    #    print ferr.readlines()

def lapw0(fday, case, WIEN, para, fh_info):
    AStep(fday, case, 'lapw0', 'in0', WIEN, para, fh_info)
        
def lapw1(fday, case, WIEN, para, fh_info):
    AStep(fday, case, 'lapw1', 'in1', WIEN, para, fh_info)

def lapwso(fday, case, WIEN, para, fh_info):
    AStep(fday, case, 'lapwso', 'inso', WIEN, para, fh_info)

def lcore(fday, case, WIEN, fh_info):
    AStep(fday, case, 'lcore', 'inc', WIEN,'', fh_info)

def mixer(fday, case, WIEN, fh_info):
    AStep(fday, case, 'mixer', 'inm', WIEN,'', fh_info)

def dmft2(fday, case, fh_info, WIEN, MPI ):
    name = 'dmft2'
    if not os.path.isfile(name+'.def'):
        print '  stop error: the required input file '+name+'.def for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+name+'.def for the next step could not be found!'
        
    tim = time.strftime("%H:%M:%S")
    cmd = utime+' '+MPI + ' ./dmft2 '+name+'.def >> dmft2_info.out '
    print >> fh_info, (name+': '), cmd
    fh_info.flush()
    
    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()

    print cmd    
    stdin, stdout, stderr = os.popen3(cmd)
    err,out = stderr.read().strip(), stdout.read().strip()
    print err.splitlines()[0]
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, err.splitlines()[-2])
    fday.flush()

    for fe in glob.glob(name+'.error*'):
        if getsize(fe) !=0:
            print 'ERROR in '+fe+' from file:', open(fe,'r').readlines()
            #sys.exit(1)

def scf(fday, case, WIEN):
    dat=[]
    for i in ['0','1','so','2']:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()
            
    if os.path.exists('Eorb.dat'):
        dat += '\n---------- Eorb.dat --------------\n'
        dat += open('Eorb.dat','r').readlines()
    
    for i in ['1s','2s','c']:
        if os.path.exists(case+'.scf'+i):
            dat += '\n------- '+case+'.scf'+i+' --------\n'
            dat += open(case+'.scf'+i,'r').readlines()
            
    fs = open(case+'.scf', 'a')
    fs.writelines(dat)
    fs.close()

def scfm(fday, case, WIEN):
    dat=[]
    for i in ['m']:
        dat += '\n------- '+case+'.scf'+i+' --------\n'
        if os.path.exists(case+'.scf'+i):
            dat += open(case+'.scf'+i,'r').readlines()
            
    fs = open(case+'.scf', 'a')
    fs.writelines(dat)
    fs.close()
    
def scf1(fday, case, WIEN):

    for i in ['clmsum','vsp','vns','vrespsum']:
        name = case+'.'+i
        if os.path.exists(name) and os.path.getsize(name)>0:
            shutil.copy2(name, name+'_old')		#save last cycle
            
        
def Diff(fday, case):
    #
    fs = open(case+'.scf', 'r')
    dat = fs.readlines()
    fs.close()

    Ene=[]
    DChr=[]
    for line in dat:
        if re.match(r':ENE', line):
            Ene.append(float(line.split()[-1]))
            #print line.split()[-1]
        if re.match(r':DIS', line):
            DChr.append(float(line.split()[-1]))
            #print line.split()[-1]
    drho = DChr[-1]
    if size(Ene)>1:
        dene = abs(Ene[-1]-Ene[-2])
    else:
        dene = 1.0
    print >> fday, ':ENERGY convergence: ', dene
    print >> fday, ':CHARGE convergence: ', drho
    return (drho, dene)

def SimplifySiginds(siginds):
    " Takes dictionary of Sigind's and creates list or non-zero columns"

    def union(data):
        " Takes a union of array or list"
        c = []
        for d in data:
            if d not in c:
                c.append(d)
        return c

    cols={}
    for icix in siginds.keys():
        Sigind = siginds[icix]
        col = sort(union(array(Sigind).flatten())).tolist()
        if 0 in col: col.remove(0)
        cols[icix] = col
    return cols

# Do not need here. Trafo has been modified instead.
#def YLM_convert_w2k_to_standard(CF, l, qsplit):
#    """WIEN2K uses a different convention for spherical harmonics than the ED code
#    in the case of d-orbitals.
#
#    Input: T(i,m) = 'Transformation matrix' printed inside case.indmfl
#     - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
#     - m runs over complex harmonics (-2, -1, 0, 1, 2)
#
#    Output: T(i,m')
#
#    The matrix converting between conventions is
#              1  0 0 0 0
#              0 -i 0 0 0
#    R(m,m') = 0  0 1 0 0
#              0  0 0 i 0
#              0  0 0 0 1
#    This function just does: T(i,m') = T(i,m) * R(m,m')
#    """
#    if l == 2:
#        # only apply conversion if this is NOT a spin-orbit run
#        if qsplit not in [1, -1, 4, 11]:
#            nspin = len(CF)/(2*l+1)  # double size of R if two spins
#            R = diag([1, -1j, 1, 1j, 1]*nspin)
#            CFnew = dot(CF, R)
#    else:
#        CFnew = CF
#
#    return CFnew

def combineud(case, ROOT, fh_info):
    cmd = ROOT+'/combineud '+case
    print >> fh_info, ('combineud: '), cmd
    fh_info.flush()
    stdin, stdout, stderr = os.popen3(cmd)
    err,out = stderr.read().strip(), stdout.read().strip()
    if err: print err
    if out: print out
    shutil.move(case+'.clmval', case+'.clmval_up')
    shutil.move(case+'.clmvaldn', case+'.clmval_dn')
    shutil.move(case+'.clmval_aver', case+'.clmval')

if __name__ == '__main__':
    # -------------- Parameters which are set through input files sparams.dat and params.dat -----------
    # Default values of some parameters
    params = {
              'DCs'           : 'fixn',    # the double counting scheme, which fixes Edc through n0
              'finish'        : 1000,       # number of iterations of the charge loop
              'mix_delta'     : 1.0,       # whether to mix delta, or not
              'riter'         : 100,       # How often to restart broyden for charge mixing
              'sleeptime'     : 2,         # If broyden file are present at the submit time, user has some time to kill the job
              'cc'            : 1e-3,      # the charge density precision to stop the LDA+DMFT run
              'ec'            : 1e-5,      # the energy precision to stop the LDA+DMFT run
              'so'            : False,     # spin-orbit coupling
              'rCF'           : None,      # Reduction of the crystal field splitting, if necessary
              'saver'         : 0.0,       # Average over the last few DMFT self-energies is used for charge. Should be between [0,1]
              'KeepPhiJac'    : False,     # Keep Jacobian of the Phi-solver 
              'StartProg'     : 'lapw0',   # The program to start with
              'EndProg'       : ' '        # Exit after program EndProg
              }

    p = Params('params.dat', params)
    if '-s' in sys.argv:
      p['StartProg']=sys.argv[sys.argv.index('-s')+1]
      print 'StartProg =',p['StartProg']
    if '-e' in sys.argv:
      p['EndProg']=sys.argv[sys.argv.index('-e')+1]
      print 'EndProg =',p['EndProg']
    if '-kj' in sys.argv:
      p['KeepPhiJac'] = True
      print 'KeepPhiJac =',p['KeepPhiJac']
    if '-cc' in sys.argv:
      p['cc']=sys.argv[sys.argv.index('-cc')+1]
      print 'cc =',p['cc']
    if '-ec' in sys.argv:
      p['ec']=sys.argv[sys.argv.index('-ec')+1]
      print 'ec =',p['ec']
    if '-iter' in sys.argv:
      p['finish']=sys.argv[sys.argv.index('-iter')+1]
      print 'finish =',p['finish']


    if len(sys.argv)>1 and sys.argv[1] in ['-h', '--help']:
        help="""
        The script runs LDA+DMFT. It is a wrapper, which calls other scripts and executables.
        Usuall execution goes through the following steps in a loop:
        
          x lapw0        -- computes LDA potential on current LDA+DMFT charge
          x lapw1        -- solves LDA eigenvalue problem
          [x lapwso]     -- adds spin-orbit
          [x dmft0]      -- computes dmft eigenvalues
          [x mu]         -- computes dmft chemical potential
          x dmft1        -- computes local green's function and hybridization
          run impurity   -- runs impurity problem
          x dmft2        -- computes LDA+DMFT valence charge, and the chemical potential
          x lcore        -- computes LDA core charge
          x mixer        -- mixes total charge with the previous charge

        The most common parameters, which should be given through 'params.dat' file, include:

        name       possible values   default     help
        --------------------------------------------------------------------------------------
        max_iterations [int]         1        # number of iteration of the dmft-loop only
        finish         [int]         100      # number of iterations of the charge+dmft loop
        cc             [float]       1e-5,    # the charge density precision to stop the LDA+DMFT run
        ec             [float]       1e-5,    # the energy precision to stop the LDA+DMFT run
        broyd          [True| False] True     # Are we using broyden for charge mixing
        riter          [int]         99       # How often to restart broyden for charge mixing
        sleeptime      [int]         2        # If broyden file are present at the submit time, user
                                                # has some time to kill the job
        so             [True| False] False    # spin-orbit coupling
        DCs            [fixn| default| S0| fixEimp]
                                     fixn     # double-counting scheme
        mix_delta      [float]       1.0      # linear mixing parameter for Delta.
        """
        print help
        sys.exit(0)

    para = ''
    if os.path.isfile('.machines') and os.path.getsize('.machines')>0 : para = ' -p'
    
    dmfe = ga_utils.DmftEnvironment()  # DMFT paths
    w2k = ga_utils.W2kEnvironment()    # W2k filenames and paths

    # Processing 'case.indmfl' file
    inl = ga_indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    toclean = glob.glob('*.scf')+glob.glob('Edc.dat')+glob.glob('dmft1.error*')+glob.glob('dmft2.error*')+glob.glob('*.cdos.*')+glob.glob('*.dlt.*')+glob.glob('*.gc1.*')+glob.glob('*.outputdmf1.*')+glob.glob('*.outputdmf2.*')
    for f in toclean: os.remove(f)
    
    # Info files
    fday = open(w2k.case+'.dayfile', 'w')
    print >> fday, 'Calculating '+w2k.case+' in '+os.getcwd()+'\n'+'on '+socket.gethostname()+' with PID '+str(os.getpid())
    print >> fday, '\n\n'
    
    fh_info = open('dmft_info.out','w')
    fh_pinfo = open('info.iterate', 'a')
    fh_pinfo.write('%3s.%3s %12s %12s %12s %12s %12s\n' % ('#','#','mu','Eimp','Eimp[-1]','Edc','nf'))
    fh_pinfo.flush()
    
    if len(glob.glob('dmft0_info.out'))!=0 : os.remove('dmft0_info.out')
    if len(glob.glob('dmft1_info.out'))!=0 : os.remove('dmft1_info.out')
    if len(glob.glob('dmft2_info.out'))!=0 : os.remove('dmft2_info.out')

    m_extn = 'dn' if os.path.exists(w2k.case+'.indmfl'+'dn') else ''
    if m_extn:
        print >> fh_info, 'INFO: case.indmfldn present => magnetic calculation with two dmft2 steps'
        inldn = ga_indmffile.Indmfl(w2k.case, 'indmfl'+m_extn)

    # Reading parameters from params.dat
    p.refresh(fh_info)

    # Reads structure file
    struct = Struct(w2k.case+'.struct', fh_info)

    # Check spin-orbit coupling
    if os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print 'Found '+w2k.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!'
        p['so'] = True

    # corelated indexes
    cixs = inl.siginds.keys()        # all columns for cix

    # Nuclear charges of all correlated blocks
    Znuc={} 
    for icix in cixs:
        atm = inl.cps[icix][0][0]                # correlated atoms
        Znuc[icix] = int(struct.aZnuc[atm-1]+0.5)  # nuclear charge of each correlated block
    print >> fh_info, 'Znucs=', Znuc
    
    shutil.copy2(dmfe.ROOT+'/dmft', '.')  # For parallel execution, executable has to be on the current directory
    shutil.copy2(dmfe.ROOT+'/dmft2', '.') # For parallel execution, executable has to be on the current directory
    shutil.copy2(dmfe.ROOT+'/CyGutz', '.') # For parallel execution, executable has to be on the current directory

    print >> fh_info, '--------- Preparing GF calculation ---------'
    cmd = dmfe.ROOT+'/ga_x_dmft.py -d'+para+' dmft1 >> dmft1_info.out 2>&1'
    print >> fh_info, cmd
    stdin, stdout, stderr = os.popen3(cmd)
    print >> fh_info, stdout.read(), stderr.read()

    print >> fh_info, '--------- Preparing GF calculation ---------'
    cmd = dmfe.ROOT+'/ga_x_dmft.py -d'+para+' dmft2 >> dmft1_info.out 2>&1'
    print >> fh_info, cmd
    stdin, stdout, stderr = os.popen3(cmd)
    print >> fh_info, stdout.read(), stderr.read()

    if p['finish']>0: # Charge self-consistent; 1: One shot
        
        if os.path.isfile(w2k.case+'.clmsum_old'): shutil.copy2(w2k.case+'.clmsum_old', w2k.case+'.clmsum')
        if not os.path.isfile(w2k.case+'.clmsum'):
            print 'no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            print >> fday, 'no '+w2k.case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            sys.exit(1)
        
        
        if os.path.isfile(w2k.case+'.broyd1'):
            print w2k.case+'.broyd* files present! You did not save_lapw a previous clculation.' 
            print 'You have '+str(p['sleeptime'])+' seconds to kill this job ( ^C   or   kill '+str(os.getpid())+' )'
            print 'or the script will rm *.broyd* and continue (use -NI to avoid automatic rm)'
            time.sleep(p['sleeptime'])
            cmd = 'rm -f *.broyd*'
            os.popen3(cmd)
            stdin, stdout, stderr = os.popen3(cmd)
            err,out = stderr.read().strip(), stdout.read().strip()
            print >> fday, w2k.case+'.broyd* files removed!'

        print >> fday, '\n   start'+(' '*8)+ time.asctime()+' with '+p['StartProg']+' ('+str(1)+'/'+str(p['riter'])+' to go)'
        print >> fday, '\n   cycle %s \t%s %s/%s to go\n' % (0,time.asctime(),p['finish'],(p['finish'])%p['riter'])

        if p['StartProg'] in 'lapw0':
          lapw0(fday,w2k.case,w2k.WIENROOT,para,fh_info);              fday.flush()
        if p['StartProg'] in 'lapw0 lapw1':
          lapw1(fday,w2k.case,w2k.WIENROOT,para,fh_info);              fday.flush()
        if p['StartProg'] in 'lapw0 lapw1 lapwso':
          if p['so']: lapwso(fday,w2k.case,w2k.WIENROOT,para,fh_info); fday.flush()

    if para : ga_convert_processes.convert(p['so'], fh_info)

    #####################
    # Major charge loop #
    #####################
    icycle = 0

    while icycle < abs(p['finish']):
    
        #####################
        # DMFT loop         #
        #####################
        if p['KeepPhiJac'] and os.path.isfile('GLP_FJAC.OUT'):
          shutil.copy2('GLP_FJAC.OUT', 'GLP_FJAC.INP')
        if icycle>0 or (icycle==0 and p['StartProg'] in 'lapw0 lapw1 lapwso dmft1'):
          dmft1(fday,w2k.case,fh_info,dmfe.MPI,w2k.WIENROOT);  fday.flush()
        if p['EndProg'] in 'dmft1': break
        if icycle>0 or (icycle==0 and p['StartProg'] in 'lapw0 lapw1 lapwso dmft1 CyGutz'):
          CyGutz(fday,w2k.case,fh_info,dmfe.MPI,w2k.WIENROOT);  fday.flush()
        if p['EndProg'] in 'CyGutz': break
        if os.path.isfile(w2k.case+'.outputdmf1.0'):
          shutil.copy2(w2k.case+'.outputdmf1.0', 'GL_LOG.OUT')
        dmft2(fday,w2k.case,fh_info,w2k.WIENROOT,dmfe.MPI);  fday.flush()
        if p['EndProg'] in 'dmft2': break

        lcore(fday,w2k.case,w2k.WIENROOT,fh_info);       fday.flush()
        scf  (fday,w2k.case,w2k.WIENROOT);               fday.flush()
        scf1 (fday,w2k.case,w2k.WIENROOT);               fday.flush()
        mixer(fday,w2k.case,w2k.WIENROOT,fh_info);       fday.flush()
        scfm (fday,w2k.case,w2k.WIENROOT);               fday.flush()
        drho, dene = Diff(fday, w2k.case)
        print >> fh_info, '------- charge step', icycle, 'done ---------'
        fday.flush()
        print 'cc=',drho,'->',p['cc'],'ec=',dene,'->',p['ec'],'icycle=',icycle
        if drho<p['cc'] and dene<p['ec'] and icycle>-1 : break

        #if (p['finish']-icycle)%p['riter'] == 0:
        #    cmd = 'rm -f *.broyd*'
        #    os.popen3(cmd)
        #    stdin, stdout, stderr = os.popen3(cmd)
        #    print stderr.read().strip(), stdout.read().strip()
        #    print >> fday, w2k.case+'.broyd* files removed!'
            
        lapw0(fday,w2k.case,w2k.WIENROOT,para,fh_info);       fday.flush()
        lapw1(fday,w2k.case,w2k.WIENROOT,para,fh_info);       fday.flush()

        if p['so']:  lapwso(fday,w2k.case,w2k.WIENROOT,para,fh_info)
        if os.path.exists('GLX.OUT'):
          shutil.copy2('GLX.OUT', 'GLX.INP')
        if os.path.exists('GL_NELF1.OUT'):
          shutil.copy2('GL_NELF1.OUT', 'GL_NELF1.INP')

        icycle += 1   # END CHARGE LOOP
        
#    toclean = glob.glob('dmft1.error*')+glob.glob('dmft2.error*')+glob.glob('*.outputdmf1.*')+glob.glob('*.outputdmf2.*')+glob.glob('*.broyd*')+glob.glob('dmft')+glob.glob('dmft2')
#    for f in toclean: os.remove(f)
