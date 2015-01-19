#!/usr/bin/env python
from sys import *
import os, socket, shutil, time, re
import glob, datetime

bin=''
utime='/usr/bin/time'
dmft2_exe = './dmft2'

def AStep(fday, case, name, inp_name):
    if not os.path.isfile(case+'.'+inp_name):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        
    tim = time.strftime("%H:%M:%S")
    cmd = bin+'x '+name
    print cmd
    stdin, stdout, stderr = os.popen3(cmd)
    err,out = stderr.read().strip(), stdout.read().strip()
    print err
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, out)
    if os.path.exists(name+'.error') and os.path.getsize(name+'.error')>0:
        ferr = open(name+'.error', 'r')
        print ferr.readlines()

def lapw0(fday, case):
    AStep(fday, case, 'lapw0', 'in0')
        
def lapw1(fday, case):
    AStep(fday, case, 'lapw1', 'in1')

def lapwso(fday, case):
    AStep(fday, case, 'lapwso', 'inso')

def lcore(fday, case):
    AStep(fday, case, 'lcore', 'inc')

def mixer(fday, case):
    AStep(fday, case, 'mixer', 'inm')

def dmft2(fday, case):
    inp_name = 'indmf2'
    name = 'dmft2'
    if not os.path.isfile(case+'.'+inp_name):
        print '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+case+'.'+inp_name+' for the next step could not be found!'
    if not os.path.isfile(name+'.def'):
        print '  stop error: the required input file '+name+'.def for the next step could not be found!'
        print >> fday, '  stop error: the required input file '+name+'.def for the next step could not be found!'
        
    tim = time.strftime("%H:%M:%S")
    cmd = utime+' '+dmft2_exe+' '+name+'.def'

    fl = open(':log', 'a')
    print >> fl, time.strftime("%a %b %d %H:%M:%S %Z %Y")+'>     '+name
    fl.close()
    
    print cmd    
    stdin, stdout, stderr = os.popen3(cmd)
    err,out = stderr.read().strip(), stdout.read().strip()
    print err.splitlines()[0]
    print >> fday, '>%-10s ( %s ) %s' % (name, tim, err.splitlines()[1])
    if os.path.exists(name+'.error') and os.path.getsize(name+'.error')>0:
        ferr = open(name+'.error', 'r')
        print ferr.readlines()

    for line in out.split('\n'):
        if re.match('Difference in the chemical potential', line):
            dEF = float(line.split()[5])
        
    fd = open(name+'.out', 'a');    fd.write(out);    fd.close()
    return dEF

def scf(fday, case):
    dat=[]
    for i in ['0','1','so','2','1s','2s','c', 'm']:
        dat += '\n------- '+case+'.scf'+i+' --------\n'
        if os.path.exists(case+'.scf'+i):
            dat += open(case+'.scf'+i,'r').readlines()
    fs = open(case+'.scf', 'a')
    fs.writelines(dat)
    fs.close()
             
def scf1(fday, case):

    for i in ['clmsum','vsp','vns','vrespsum']:
        name = case+'.'+i
        if os.path.exists(name) and os.path.getsize(name)>0:
            shutil.copy2(name, name+'_old')		#save last cycle
            
        
def Diff(fday, case, dEF):
    #
    fs = open(case+'.scf', 'r')
    dat = fs.readlines()
    fs.close()

    DEne=[]
    DChr=[]
    for line in dat:
        if re.match(r':ENE', line):
            DEne.append(float(line.split()[-1]))
            #print line.split()[-1]
        if re.match(r':DIS', line):
            DChr.append(float(line.split()[-1]))
            #print line.split()[-1]
    drho = DChr[-1]
    dene = abs(DEne[-1]-DEne[-2])
    print >> fday, ':ENERGY convergence: ', dene
    print >> fday, ':CHARGE convergence: ', drho
    print >> fday, ':EF     convergence: ', dEF
    return (drho, dene)
    
if __name__ == '__main__':

    broyd = True
    riter=99
    iter=100
    sleeptime=2
    
    pwd = os.getcwd()    
    case = os.path.splitext(glob.glob("*.struct")[0])[0] # finds file xxx.struct and sets case to xxx
 
    fday = open(case+'.dayfile', 'w')

    print >> fday, 'Calculating '+case+' in '+pwd+'\n'+'on '+socket.gethostname()+' with PID '+str(os.getpid())
    print >> fday, '\n\n'

    if os.path.isfile(case+'.clmsum'):
        if os.path.isfile(case+'.clmsum_old'):
            shutil.copy2(case+'.clmsum_old', case+'.clmsum')
        else:
            print 'no '+case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            print >> fday, 'no '+case+'.clmsum(_old) file found, which is necessary for lapw0 !'
            exit(1)

  
    if broyd:
        if os.path.isfile(case+'.broyd1'):
            print case+'.broyd* files present! You did not save_lapw a previous clculation.' 
            print 'You have '+str(sleeptime)+' seconds to kill this job ( ^C   or   kill '+str(os.getpid())+' )'
            print 'or the script will rm *.broyd* and continue (use -NI to avoid automatic rm)'
            time.sleep(sleeptime)
            cmd = 'rm -f *.broyd*'
            os.popen3(cmd)
            stdin, stdout, stderr = os.popen3(cmd)
            err,out = stderr.read().strip(), stdout.read().strip()
            print >> fday, case+'.broyd* files removed!'


    
    next='lapw0'
    print >> fday, '\n   start'+(' '*8)+ time.asctime()+' with '+next+' ('+str(1)+'/'+str(riter)+' to go)'
    #nohup echo in cycle $icycle "   ETEST: $etest[3]   CTEST: $ctest[3]"


    for icycle in range(iter):
        
        if ((iter-icycle)%riter == 0):
            cmd = 'rm -f *.broyd*'
            os.popen3(cmd)
            stdin, stdout, stderr = os.popen3(cmd)
            err,out = stderr.read().strip(), stdout.read().strip()
            print >> fday, case+'.broyd* files removed!'
            
        print >> fday, '\n   cycle %s \t%s %s/%s to go\n' % (icycle+1,time.asctime(),iter-icycle,(iter-icycle)%riter)
        lapw0(fday,case);       fday.flush()
        lapw1(fday,case);       fday.flush()
        lapwso(fday,case);      fday.flush()
        dEF = dmft2(fday,case); fday.flush()
        lcore(fday,case);       fday.flush()
        scf1(fday,case);        fday.flush()
        mixer(fday,case);       fday.flush()
        scf(fday,case);         fday.flush()
        drho, dene = Diff(fday, case, dEF)
        
    fday.close()
    
#os.popen3(cmd)
