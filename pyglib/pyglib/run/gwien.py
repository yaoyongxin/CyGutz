from __future__ import print_function
from builtins import range, zip

import os, sys, glob, h5py, socket, shutil, time, re
import numpy as np
from collections import deque
from subprocess import Popen, PIPE
from pyglib.io.fio import file_exists


def get_file_info(fname, unit, idmf, case, scratch, so, para, cmplx, _band,
        updn, dnup):
    '''help function to setup informations in def file.
    '''
    if 'in2' == fname:
        return [unit, "'{}.in2{}'".format(case, cmplx), "'old'",
                "'formatted'", 0]
    elif 'inso' == fname:
        return [unit, "'{}.inso'".format(case), "'unknown'", "'formatted'", 0]
    elif 'indmfl' == fname:
        return [unit, "'{}.indmfl'".format(case), "'old'", "'formatted'", 0]
    elif 'outputdmf' == fname:
        return [unit, "'{}.outputdmf{}'".format(case, idmf), "'unknown'", \
                "'formatted'", 0]
    elif 'in1c' == fname:
        return [unit, "'{}.in1c'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vectorupdn' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                updn, para), "'unknown'","'unformatted'",9000]
    elif 'vectordnup' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                dnup, para), "'unknown'","'unformatted'",9000]
    elif 'klist' == fname:
        return [unit, "'{}.klist{}'".format(case, _band), "'old'", \
                "'formatted'", 0]
    elif 'kgen' == fname:
        return [unit, "'{}.kgen'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vspupdn' == fname:
        return [unit, "'{}.vsp{}'".format(case, updn), "'old'", \
                "'formatted'", 0]
    elif 'vspdnup' == fname:
        return [unit, "'{}.vsp{}'".format(case, dnup), "'unknown'", \
                "'formatted'", 0]
    elif 'struct' == fname:
        return [unit, "'{}.struct'".format(case), "'old'", "'formatted'", 0]
    elif 'rotlm' == fname:
        return [unit, "'{}.rotlm'".format(case), "'unknown'", "'formatted'", 0]
    elif 'energysodum' == fname:
        if so == 'so':
            sodum = 'dum'
        else:
            sodum = dnup
        return [unit, "'{}.energy{}'".format(case, sodum), \
               "'unknown'", "'formatted'", 0]
    elif 'energyupdn' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, updn, para), \
                "'unknown'", "'formatted'", 0]
    elif 'energydnup' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, dnup, para), \
                "'unknown'", "'formatted'", 0]
    elif 'clmval' == fname:
        return [unit, "'{}.clmval{}'".format(case, updn), "'unknown'", \
                "'formatted'", 0]
    elif 'recprlist' == fname:
        return [unit, "'{}.recprlist'".format(case), "'unknown'", \
                "'formatted'", 9000]
    elif 'scf2' == fname:
        return [unit, "'{}.scf2'".format(case), "'unknown'", "'formatted'", 0]
    elif 'norm' == fname:
        return [unit, "'{}.norm{}'".format(case, so), "'unknown'", \
                "'formatted'", 0]
    else:
        raise ValueError('No matching file name {}!'.format(fname))


def fcreate_def_dmft(case, scratch='.', so='', para='', idmf='1', cmplx='',
        _band='', updn='', dnup='dn'):
    '''create dmft1/2.def file.
    '''
    fdef = open('dmft{}.def'.format(idmf), 'w')
    if idmf == '1':
        fname_list = ['in2', 'inso', 'indmfl', 'outputdmf', \
                'in1c', 'vectorupdn', 'vectordnup', 'klist', \
                'kgen', 'vspupdn', 'vspdnup', 'struct', \
                'rotlm', 'energydnup', 'energyupdn']
        unit_list = [3, 4, 5, 6, 7, 9, 10, 13, 14, 18, 19, 20, 22, 59, 60]
    else:
        fname_list = ['in1c', 'inso', 'in2', 'outputdmf', 'indmfl',
                'clmval', 'vectorupdn', 'vectordnup', 'recprlist',
                'kgen', 'vspupdn', 'struct', 'scf2', 'rotlm', 'energysodum',
                'energyupdn', 'energydnup', 'norm']
        unit_list = [3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 18, 20, 21, 22, 29,
                30, 31, 12]

    for fname, unit in zip(fname_list, unit_list):
        fdef.write("{:3d}, {:<15s}, {:<10s}, {:<13s}, {:<4d}\n".format(\
                *get_file_info(fname, unit, idmf, case, scratch, so, \
                para, cmplx, _band, updn, dnup)))
    fdef.close()


def onestep(fday, case, exec_name, w_root, para):
    '''wien2k steps
    '''
    time_start = time.strftime("%H:%M:%S")
    cmd = ['{}/x'.format(w_root), para, '-f', case, exec_name]
    print(' '.join(x for x in cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    fday.write('>{:<10s} ({}) {}'.format(exec_name, time_start, out))
    fday.flush()
    for f in glob.glob('{}.error*'.format(exec_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def gonestep(fday, exec_name, mpi):
    '''dmft1, CyGutz and dmft2 steps.
    '''
    time_start = time.strftime("%H:%M:%S")
    with open(':log', 'a') as f:
        f.write('{}>   {}\n'.format(time.strftime("%a %b %d %H:%M:%S %Z %Y"), \
                exec_name))

    cmd = ['/usr/bin/time']
    if mpi != '':
        cmd.extend(mpi)
    cmd.append('./{}'.format(exec_name))
    if 'dmft' in exec_name:
        cmd.append('{}.def'.format(exec_name))

    print(' '.join(x for x in cmd))
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    with open('{}_info.out'.format(exec_name), 'w') as f:
        f.write(out)
    fday.write('>{:<10s} ({}) {}'.format(exec_name, time_start, \
            err.splitlines()[-2]))
    fday.flush()
    for f in glob.glob('{}.error*'.format(exec_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def get_file_content(fname):
    if os.path.exists(fname):
        data = '\n------- {} --------\n'.format(fname)
        with open(fname, 'r') as f:
            data += f.read()
        return data
    else:
        return ''


def scf(case):
    f_list = ['{}.scf{}'.format(case, i) for i in ['0', '1', 'so', '2']]
    f_list.append('Eorb.dat')
    f_list += ['{}.scf{}'.format(case, i) for i in ['1s', '2s', 'c']]
    data = ''.join(get_file_content(f) for f in f_list)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)

    for i in ['clmsum', 'vsp', 'vns', 'vrespsum']:
        name = '{}.{}'.format(case, i)
        if file_exists(name):
            shutil.copy2(name, '{}_old'.format(name))


def scfm(case):
    f_scf = '{}.scfm'.format(case)
    data = get_file_content(f_scf)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)


def diff(fday, case, mix_dc):
    e_que = deque([], 2)
    with open('{}.scf'.format(case), 'r')  as f:
        for line in f:
            if ':DIS' in line:
                d_rho = float(line.split()[-1])
            if ':ENE' in line:
                e_que.append(float(line.split()[-1]))
    if len(e_que) == 2:
        d_etot = e_que[1] - e_que[0]
    else:
        d_etot = 1.0

    with h5py.File("GPARAM.h5", 'a') as f:
        dc_mode = f["/dc_mode"][0]
        if dc_mode == 12:
            nelf_list_in = f["/dc_nelf_list"][()]
            with h5py.File("GDC_NELF_OUT.h5", 'r') as fp:
                nelf_list_out = fp["/dc_nelf_list_out"][()]
            nelf_diff_list = nelf_list_out - nelf_list_in
            dcv_err = np.max(np.abs(nelf_diff_list))
            nelf_list_mix = nelf_list_in + mix_dc*nelf_diff_list
            f["/dc_nelf_list"][()] = nelf_list_mix
        else:
            dcv_err = 0.

    fday.write(':ENERGY convergence: {}\n'.format(d_etot))
    fday.write(':CHARGE convergence: {}\n'.format(d_rho))
    fday.write(':VDC convergence: {}\n'.format(dcv_err))
    return d_rho, d_etot, dcv_err


def processes_convert(so):
    if not file_exists('.processes'):
        print('.processes file not present. It must be a serial run.')
        return
    lines = open('.processes').readlines()
    work = {}
    nkstart = 0
    for line in lines:
        data = line.split(':')
        if re.match('\s*\d+', data[0]):
            vecn = [None]*4
            i, nkp, nprc = map(int,data[::2])
            if not so:
                fdef = open('lapw1_{}.def'.format(i), 'r')
                for line in fdef:
                    data = line.split(',')
                    data0 = int(data[0])
                    if data0 == 10 or data0 == 11:
                        data0 = data0 % 10
                        m = re.search('.*[\'|\"](.*)_(\d+)', data[1])
                        assert m is not None, 'vector file to macth ' + \
                                ' lapw1.def not found!'
                        vecn[data0*2] = '{}_{}'.format(m.group(1), m.group(2))
                        vecn[data0*2+1] = '{}dn_{}'.format(m.group(1), \
                                 m.group(2))
                fdef.close()
            else:
                fdef = open('lapwso_{}.def'.format(i), 'r')
                for line in fdef:
                    data = line.split(',')
                    if int(data[0])==42: vecn[0]=data[1].split("'")[1]
                    if int(data[0])==41: vecn[1]=data[1].split("'")[1]
                    if int(data[0])==52: vecn[2]=data[1].split("'")[1]
                    if int(data[0])==51: vecn[3]=data[1].split("'")[1]
                fdef.close()

            if work.has_key(nprc):
                work[nprc].append((i, nkp, nkstart, vecn))
            else:
                work[nprc]=[(i, nkp, nkstart, vecn)]
            nkstart += nkp

    for prc in sorted(work.keys()):
        fo = open('_processes_{}'.format(prc-1), 'w')
        for (i, nkp, nkstart, vecn) in work[prc]:
            fo.write('{} {} {} "{}" "{}" "{}" "{}"\n'.format(\
                    i, nkp, nkstart, *vecn))


def create_gomp_file():
    '''
    Create GOMP.h5 file based on GMPI_X.h5 for openMP execution.
    '''
    with h5py.File('GMPI_0.h5', 'r') as f:
        num_procs = f["/nprocs"][0]
    nvec = 0
    kvec1 = []
    kvec2 = []
    for iproc in range(num_procs):
        with h5py.File('GMPI_' + str(iproc) + '.h5', 'r') as f:
            nvec += f["/nvec"][0]
            kvec = f["/KVEC"][()].T
            kvec1.append(kvec[0])
            if kvec.shape[1] == 2:
                kvec2.append(kvec[1])
    kvec = np.asarray(kvec1 + kvec2)

    with h5py.File('GOMP.h5', 'w') as f:
        f['/nvec'] = np.asarray([nvec])
        f['/KVEC'] = kvec.T


def run_gwien(nmaxiter=100, mix_dc=0.2, cc=1.e-3, ec=1.e-5, vc=1.e-2,
        startp='lapw0', endp='', band='', openmp=False, cygutz='CyGutz',
        pa_list=[], recycle_rl=False):
    '''Driver for Wien2k + Gutzwiller-Slave-boson job.
    '''
    if '-s' in sys.argv:
        startp = sys.argv[sys.argv.index('-s') + 1]
    if '-e' in sys.argv:
        endp = sys.argv[sys.argv.index('-e') + 1]
    if '-cc' in sys.argv:
        cc = float(sys.argv[sys.argv.index('-cc') + 1])
    if '-ec' in sys.argv:
        ec = float(sys.argv[sys.argv.index('-ec') + 1])
    if '-vc' in sys.argv:
        vc = float(sys.argv[sys.argv.index('-vc') + 1])
    if '-n' in sys.argv:
        nmaxiter = int(sys.argv[sys.argv.index('-n') + 1])
    if '-omp' in sys.argv:
        openmp = True
        print('Using Open-MP instead of MPI of CyGutz.')
    if '-amix' in sys.argv:
        mix_dc = float(sys.argv[sys.argv.index('-amix') + 1])
    if '-band' in sys.argv:
        band='-band'
    if '-rl' in sys.argv:
        recycle_rl = True
    if band == '-band':
        _band = '_band'
        nmaxiter = 1
    else:
        _band = ''

    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        help = '''
    The script is a wrapper to run Wien2k + Gutzwiller.
    It usually loops over the following steps:

        x lapw0   : computes LDA potential with current DFT+G-RISB charge
        x lapw1   : solves LDA eigenvalue equations
        [x lapwso]: second variational treatment of spin-orbit coupling
        x dmft1   : compute the local projector in the basis of DFT bands
        x cygutz  : solve the generic KS-Hubbard model using G-RISB
        x dmft2   : computes DFT+G-RISB valence charge
        x lcore   : computes DFT core charge
        x mixer   : mixes new charge density with the previous result

    The parameters with default values are as follows:

        name     default  inline-argument  help
        --------------------------------------------------------------------
        nmaxiter 100      -n 100           max charge mixing steps
        mix_dc   0.2      -amix            D.C. potential mxing param
        cc       1.e-3    -cc 1.e-3        charge density cutoff to exit
        ec       1.e-5    -ec 1.e-5        total energy cutoff to exit
        startp   'lapw0'  -s lapw0         start program
        endp     ''       -e ''            end program
        openmp   False    -omp             use openMP instead of openMPI
        rl       False    -rl              start from previous GA solutions
        '''
        print(help)
        sys.exit(0)

    para = ''
    _para = ''
    if file_exists('.machines') :
        para = ' -p'
        _para = '_x'

    toclean = glob.glob('*.scf') + glob.glob('*.error*') + \
            glob.glob('*.outputdmf?.*') + glob.glob('EMBED_HAMIL_RES*')
    for f in toclean:
        os.remove(f)

    struct_file = glob.glob('*.struct')
    if len(struct_file) != 1:
        raise ValueError('{} struct files present while only one must exist!'. \
                format(len(struct_file)))
    w_case = struct_file[0].split('.')[0]
    w_root = os.environ['WIENROOT']
    w_scratch = os.environ['SCRATCH']
    g_root = os.environ['WIEN_GUTZ_ROOT2']

    # infomation file
    fday = open(w_case + '.dayfile', 'w')
    fday.write('Calculating {} in {} \non {} with PID {}\n'.format(\
            w_case, os.getcwd(), socket.gethostname(), os.getpid()))

    # spin-orbit calculation?
    p_so = file_exists(w_case+'.inso')
    print('calculation with spin-orbit = {}'.format(p_so))
    if p_so:
        so='so'
        cmplx = 'c'
    else:
        so = ''
        cmplx = ''

    # In addition, check in1c file
    if file_exists(w_case+'.in1c'):
        cmplx = 'c'

    f_mpi = 'mpi_prefix.dat'
    if os.path.isfile(f_mpi):
        with open(f_mpi, 'r') as f:
            mpi = f.readline().split()
        print('{} exists -- running in parallel mode.'.format(f_mpi))
        print(' '.join(x for x in mpi))
    else:
        if para != '':
            raise ValueError('missing mpi_prefix.dat with .machines present!')
        mpi = ''
        print('{} not available -- running in serial mode.'.format(f_mpi))

    if openmp:
        _mpi = ''
    else:
        _mpi = mpi

    p_list = ['dmft1', 'dmft2', 'CyGutz', 'exe_spci', 'exe_spci_s2_mott']
    for pa in pa_list:
        if pa not in p_list:
            p_list.append(pa)
    for p in p_list:
        shutil.copy2(g_root+'/'+p, '.')

    # create dmft1/2.def files
    fcreate_def_dmft(w_case, scratch=w_scratch, so=so, para=_para,
            idmf='1', cmplx=cmplx, _band=_band)
    fcreate_def_dmft(w_case, scratch=w_scratch, so=so, para=_para,
            idmf='2', cmplx=cmplx, _band=_band)

    if nmaxiter > 0:
        if os.path.isfile('{}.clmsum_old'.format(w_case)):
            shutil.copy2('{}.clmsum_old'.format(w_case),
                    '{}.clmsum'.format(w_case))
        if not os.path.isfile('{}.clmsum'.format(w_case)):
            err_msg = 'no {}.clmsum(_old) file found--necessary for lapw0!'.\
                    format(w_case)
            print(err_msg)
            fday.print(err_msg+'\n')
            sys.exit(1)
        for f in glob.glob('*.broyd*'):
            os.remove(f)

        fday.write('   start at {} with {} \n    1/{} to go.\n'.format(
                time.asctime(), startp, nmaxiter))
        if startp in 'lapw0':
            onestep(fday, w_case, 'lapw0', w_root, para)
        if startp in 'lapw0 lapw1':
            onestep(fday, w_case, 'lapw1', w_root, para)
        if startp in 'lapw0 lapw1 lapwso' and p_so:
            onestep(fday, w_case, 'lapwso', w_root, para)

    if para != '':
        processes_convert(p_so)

    # Major charge density loop
    for icycle in range(nmaxiter):
        if icycle > 0 or (icycle == 0 and startp in \
                'lapw0 lapw1 lapwso dmft1'):
            gonestep(fday, 'dmft1', mpi)
            if openmp:
                create_gomp_file()
            elif os.path.isfile("GOMP.h5"):
                os.remove("GOMP.h5")
        if endp == 'dmft1':
            sys.exit(0)

        if icycle > 0 or (icycle == 0 and startp in \
                'lapw0 lapw1 lapwso dmft1 CyGutz'):
            gonestep(fday, 'CyGutz', _mpi)
        shutil.copy2('GUTZ.LOG', 'SAVE_GUTZ.LOG')
        if endp == 'CyGutz' or band == '-band':
            sys.exit(0)
        if recycle_rl:
            shutil.copy2('WH_RL_OUT.h5', 'WH_RL_INP.h5')

        gonestep(fday, 'dmft2', mpi)
        if endp == 'dmft2':
            sys.exit(0)

        onestep(fday, w_case, 'lcore', w_root, '')
        scf(w_case)
        onestep(fday, w_case, 'mixer', w_root, '')
        scfm(w_case)
        drho, dene, dvdc = diff(fday, w_case, mix_dc)

        with h5py.File('GLOG.h5', 'r') as f:
            gerr = f['/rl_maxerr'][0]

        print(('dc={:.1e}, cc={:.1e} -> {:.0e}, ec={:.1e} ' + \
                '-> {:.0e}, gc={:.1e} icycle={}').format(
                dvdc, drho, cc, dene, ec, gerr, icycle))
        if drho < cc and dene < ec and dvdc < vc:
            sys.exit(0)

        onestep(fday, w_case, 'lapw0', w_root, para)
        onestep(fday, w_case, 'lapw1', w_root, para)
        if p_so:
            onestep(fday, w_case, 'lapwso', w_root, para)



if __name__=='__main__':
    fcreate_def_dmft('FeSb2', scratch='.', so='', para='',
            idmf='1', cmplx='', _band='', updn='', dnup='dn')
