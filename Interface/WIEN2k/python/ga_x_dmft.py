#!/usr/bin/env python
import sys, os, glob, re, numpy, shutil, subprocess
from scipy import *
from scipy import optimize, interpolate
import copy,optparse
# customer imports
import ga_utils,ga_indmffile

def PrepareDefinitionFile(idmf, mode, case, cixs, updn, dnup, so, para, scratch, cmplx, _band='', m_ext=''):

    fdef = open('dmft'+idmf+'.def', 'w')
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 3, "'"+case+".in2"+cmplx+"'",     "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".indmfl"+m_ext+"'",             "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf+updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 7, "'"+case+".in1c"+"'",               "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (13, "'"+case+".klist"+_band+"'",        "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (14, "'"+case+".kgen"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",           "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (19, "'"+case+".vsp"+dnup+"'",           "'unknown'","'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",             "'old'",    "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",              "'unknown'","'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (59, "'"+case+".energy"+so+dnup+para+"'", "'unknown'","'formatted'",0) #
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (60, "'"+case+".energy"+so+updn+para+"'", "'unknown'","'formatted'",0) # 

    fdef.close()


def PrepareDefinitionFile2(idmf, mode, case, cixs, updn, dnup, so, para, scratch, cmplx='', m_ext=''):

    fdef = open('dmft'+idmf+'.def', 'w')
    if so=='so':
        sodum = 'dum'
    else:
        sodum=dnup

    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 3, "'"+case+".in1c"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 4, "'"+case+".inso"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 5, "'"+case+".in2"+cmplx+"'",     "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 6, "'"+case+".outputdmf"+idmf+updn+"'","'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 7, "'"+case+'.indmfl'+m_ext+"'",        "'old'",     "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 8, "'"+case+".clmval"+updn+m_ext+"'",   "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % ( 9, "'"+scratch+"/"+case+".vector"+so+updn+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (10, "'"+scratch+"/"+case+".vector"+so+dnup+para+"'", "'unknown'","'unformatted'",9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (13, "'"+case+".recprlist"+"'",     "'unknown'", "'formatted'", 9000)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (14, "'"+case+".kgen"+"'",          "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (18, "'"+case+".vsp"+updn+"'",      "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (20, "'"+case+".struct"+"'",        "'old'",     "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (21, "'"+case+".scf2"+updn+m_ext+"'",     "'unknown'", "'formatted'", 0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (22, "'"+case+".rotlm"+"'",         "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (29, "'"+case+".energy"+sodum+"'",  "'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (30, "'"+case+".energy"+so+updn+para+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (31, "'"+case+".energy"+so+dnup+para+"'","'unknown'", "'formatted'",0)
    print >> fdef, "%3d, %-15s, %-10s, %-13s, %-4d" % (12, "'"+case+".norm"+so+"'",       "'unknown'", "'formatted'",0)

    fdef.close()

def runExternal(idmf, ROOT, MPI, dmft_exe='dmft', m_ext=''):

    deffile = 'dmft'+idmf+m_ext+'.def'

    if MPI: # In parallel mode needs to copy the executable to the current directory
        shutil.copy(ROOT+'/'+dmft_exe, os.getcwd() )
        exe = MPI + ' ./'+dmft_exe
    else:
        exe = dmfe.ROOT+'/'+dmft_exe

    cmd = exe + ' ' + deffile

    print '..... running:', cmd
    subprocess.call(cmd, bufsize = 1, shell=True)  # line-buffered

if __name__=='__main__':
    
    usage = """usage: %prog [ options ] mode

    Executes one of the dmft steps (dmft0|dmft1|dmft2|mu|dmftp|dmftu):
      dmft1  -- computes local green's function (case.gc?), hybridization (case.dlt?) and impurity levels (case.Eimp?)
      dmft2  -- computes the valence LDA+DMFT charge
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("--so", dest="so", action='store_true', default=False, help="switches on spin-orbit coupling")
    parser.add_option("--up", dest="up", action='store_true', default=False, help="magnetic LDA calculation with vector-up first")
    parser.add_option("--dn", dest="dn", action='store_true', default=False, help="magnetic LDA calculation with vector-dn first")
    parser.add_option("-p", dest="para", action='store_true', default=False, help="turns on parallel mode")
    parser.add_option("-c",  dest="c",  action='store_true', default=False, help="complex version")
    parser.add_option("-d", dest="def_only", action='store_true', default=False, help="prepared only def and input files, but do not run")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print verbose messages to stdout")
    parser.add_option("-l", "--lext", dest="m_ext",  default='',  help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-b", "--band", action="store_true", dest="_band", default=False, help="use case.klist_band instead of case.klist")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    if len(args)!=1:
        print 'Need exactly one argument. One of the dmft1|dmft2'
        sys.exit(1)
    
    dmfe= ga_utils.DmftEnvironment()  # DMFT paths
    w2k = ga_utils.W2kEnvironment()    # W2k filenames and paths
    
    if not options.so and os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
        print 'Found '+w2k.case+'.vectorso file, hence assuming so-coupling exists. Switching -so switch!'
        options.so = True

    if not options.c and os.path.isfile(w2k.case+".in1c") and os.path.getsize(w2k.case+".in1c")>0 :
        options.c = True

    # processing arguments
    updn=''
    dnup='dn'
    so=''
    cmplx=''
    para=''
    
    if options.para:
        para = '_x'
    if options.so:
        so = 'so'
        sodum = 'dum'
        cmplx = 'c'
    if options.up:
        spin=2
        updn = 'up'
        dnup = 'dn'
    if options.dn:
        spin=2
        updn = 'dn'
        dnup = 'up'
    if options.c:
        cmplx = 'c'
        
    # Processing 'case.indmfl' file
    inl = ga_indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read
    cixs = inl.siginds.keys()        # all columns for cix

    #print 'name=', args[0]
    print 'case=', w2k.case
    #print 'c=', cmplx
    #print 'scratch=', w2k.SCRATCH
    print 'root=', dmfe.ROOT
    print 'cixs=', cixs
    #print 'def_only=', options.def_only
    
    if args[0] == 'dmft1':
        idmf = '1'
        mode = 'g'               # mode for computing eigenvalues
        k_band = ''
        if (options._band):
            k_band = '_band'

        PrepareDefinitionFile(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, para,w2k.SCRATCH, cmplx, k_band, options.m_ext)
        if not options.def_only:
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, 'dmft', options.m_ext)
    elif args[0] == 'dmft2':
        idmf = '2'
        mode = 'c'
        PrepareDefinitionFile2(idmf+options.m_ext, mode, w2k.case, cixs, updn, dnup, so, para, w2k.SCRATCH, cmplx, options.m_ext)
        if not options.def_only:
            runExternal(idmf, dmfe.ROOT, dmfe.MPI, 'dmft2',options.m_ext)
