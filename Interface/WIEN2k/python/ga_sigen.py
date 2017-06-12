#!/usr/bin/env python


'''
sigen.py -- Sigma Index Generator
Also computes legends and crystal field transformations.

Takes as input case.indmf file
Writes case.indmfl file

'''


from numpy import asarray, zeros, diag, bmat, cumsum, where, identity
import optparse, os, traceback
import trafoso
import ga_indmffile
from ga_utils import eV2Ry, W2kEnvironment


# functions to handle nspins in generating sigind from raw lists given in qsplit_table
#
def no_op(sigind, nspins): return sigind

def dup(sigind, nspins):   return sigind*nspins

def dup_shift(sigind, mult):  # duplicate and shift
    ret = []
    for i in range(mult):
        ret += [x+i*max(sigind) for x in sigind]
    return ret

qsplit_table = {
    #    2spin?  action     transtype        L = 0  L = 1          L = 2                  L = 3
    #    ------  ------     ---------        -----  -----          -----                  -----
    0  : (False, dup,       "none",         ([1],   [1,1,1],       [1,1,1,1,1],           [1,1,1,1,1,1,1]              )),  # averaged
    1  : (True,  no_op,     "relativistic", ([1,1], [1,1,2,3,3,2], [1,2,2,1,3,4,5,5,4,3], [1,2,3,3,2,1,4,5,6,7,7,6,5,4])),  # |j,mj>
    -1 : (True,  no_op,     "relativistic", ([1,2], range(1,7),    range(1,11),           range(1,15)                  )),  # |j,mj> no time-reversal symmetry
    2  : (False, dup,       "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # no symmetry, real
    -2 : (False, dup_shift, "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # no symmetry, real spin-polarized
    3  : (False, dup_shift, "none",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # no symmetry, spherical
    4  : (True,  no_op,     "relativistic", ([1,1], [1,1,2,2,2,2], [1,1,1,1,2,2,2,2,2,2], [1,1,1,1,1,1,2,2,2,2,2,2,2,2])),  # |j,mj> jm's equivalent
    5  : (False, dup,       "real",         ([1],   [1,1,2],       [1,2,3,3,4],           [1,2,2,3,4,4,5]              )),  # axial
    6  : (False, dup,       "real",         ([1],   [1,1,2],       [1,2,3,3,2],           [1,2,3,4,2,3,1]              )),  # hexagonal
    7  : (False, dup,       "real",         ([1],   [1,1,1],       [1,1,2,2,2],           [1,2,2,2,3,3,3]              )),  # cubic
    8  : (True,  dup_shift, "real",         ([1],   [1,1,2],       [1,2,3,3,4],           [1,2,2,3,4,4,5]              )),  # axial spin-polarized
    9  : (True,  dup_shift, "real",         ([1],   [1,1,2],       [1,2,3,3,2],           [1,2,3,4,2,3,1]              )),  # hexagonal spin-polarized
    10 : (True,  dup_shift, "real",         ([1],   [1,1,1],       [1,1,2,2,2],           [1,2,2,2,3,3,3]              )),  # cubic spin-polarized
    11 : (True,  no_op,     "relativistic", ([1,1], [1,1,2,3,3,2], [1,2,2,1,3,4,5,5,4,3], [1,2,3,3,2,1,4,5,6,7,7,6,5,4])),  # |j,mj> + off-diagonal
    12 : (False, dup_shift, "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # real + off-diagonal
    }

def make_legend_table():
    def _j(jj):  return ' '.join(["(%d/2,%d/2)" % (2*jj+1, x) for x in range(2*jj+1, 0, -2)])
    def _jm(jj): return ' '.join([str(x)+'/2' for x in range(-2*jj-1, 2*jj+2, 2)])

    j_labels  = (_j(0),  _j(0)+' '+_j(1),  _j(1)+' '+_j(2), _j(2)+' '+_j(3))
    real      = ('r', 'x y z', 'z^2 x^2-y^2 yz xz xy', 'xyz x^3 y^3 z^3 x(y^2-z^2) y(x^2-z^2) z(x^2-y^2)')
    axial     = ('r', 'x+y z', 'z^2 x^2-y^2 yz+xz xy', 'xyz x^3+y^3 z^3 x(y^2-z^2)+y(x^2-z^2) z(x^2-y^2)')
    hexagonal = ('r', 'x+y z', 'z^2 x^2-y^2+xy yz+xz', 'xyz+zeta(T2) x(T1)+ksi(T2) y(T1)+eta(T2) z(T1)')
    cubic     = ('r', 'x+y+z', 'eg t2g', 'xyz T1 T2')

    table = {
        0  : ('s', 'p', 'd', 'f'),
        -1 : (_jm(0), _jm(0)+' '+_jm(1), _jm(1)+' '+_jm(2), _jm(2)+' '+_jm(3)),
        1  : j_labels,
        2  : real,
        -2 : real,
        3  : ('0', '-1 0 1', '-2 -1 0 1 2', '-3 -2 -1 0 1 2 3'),
        4  : ('1/2', '1/2 3/2', '3/2 5/2', '5/2 7/2'),
        5  : axial,
        6  : hexagonal,
        7  : cubic,
        8  : axial,
        9  : hexagonal,
        10 : cubic,
        11 : j_labels,
        12 : real,
        }
    return table

legend_table = make_legend_table()



def parse_cmdline_args():
    '''Instantiates command line argument parser, and parses the passed arguments.'''

    # First, instantiate the parser by defining usage of program,
    # and all possible options
    usage = """usage: %prog [options]
    Generates sigma index matrices and transformation matrices implementing desired qsplit.

    Input:  case.indmf

    Output: case.indmfl

    The possible values of qsplit are as follows:
    """ + '\n' + ga_indmffile.qsplit_doc

    parser = optparse.OptionParser(usage)

    parser.add_option("--so", "--spin-orbit", action="store_true", default=False, help="perform calculation with spin-orbit")
    parser.add_option("--sig", "--self-energy", default="sig.inp", help="filename for self-energy: -sig SIGNAME")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    # There should be no arguments passed (only options)
    if len(args) > 0:
        parser.error("Unexpected argument(s): " + " ".join(args))

    return options


def add_offdiag_sigind(sigind):
    """ Changes Sigind such that the off-diagonal matrix elements are present."""
    icounter = max(diag(sigind)) + 1
    for i in range(len(sigind)):
        for j in range(i+1,len(sigind)):
            sigind[i,j] = icounter
            sigind[j,i] = icounter
            icounter += 1

def add_offdiag_legend(legend):
    n = len(legend)
    icounter = n+1
    for i in range(n):
        for j in range(i+1,n):
            legend.append('('+str(i+1)+','+str(j+1)+')')
            icounter += 1


def make_legend(qsplit, L, nsites, nspins):
    func = qsplit_table[qsplit][1]
    legend_base = legend_table[qsplit][L].split()
    legend_one  = func(legend_base, nspins)
    legend = dup_shift(legend_one, nsites)
    if qsplit in [11, 12]:
        add_offdiag_legend(legend)
    return legend


def cmp_sigind(qsplit, L, nsites, nspins):
    """Computes indices to self-energy, which determines non-zero matrix elements and their symmetry."""
    func = qsplit_table[qsplit][1]
    sigind_base = qsplit_table[qsplit][3][L]
    sigind_one  = func(sigind_base, nspins)
    sigind = diag(dup_shift(sigind_one, nsites))
    if qsplit in [11, 12]:
        add_offdiag_sigind(sigind)
    return sigind


def cmp_cftrans(qsplit, L, nsites, nspins):
    """Computes transformation matrix from spherical harmonics to basis specified by qsplit."""
    transtype = qsplit_table[qsplit][2]
    if transtype == 'relativistic':
        T = trafoso.trafoso(L)
    elif transtype == 'none':
        T = identity(nspins*(2*L+1), dtype=complex)
    else:
        raise Exception, 'ERROR: Transformation `'+transtype+'` not yet implemented.'

    # make block matrix, one block for each site of (potentially) cluster problem
    Z = zeros((len(T), len(T)), dtype=complex)
    S = [[Z]*nsites for i in range(nsites)]
    for i in range(nsites):
        S[i][i] = T
    return asarray(bmat(S))


def check_nspins(qsplit, nspins):
    '''Make sure user-specified nspins is compatible with qsplit.'''
    require_twospins = qsplit_table[qsplit][0]
    if require_twospins and nspins != 2:
        raise Exception, 'ERROR: qsplit = %d requires spin-polarized/spin-orbit calculation.' % qsplit


def offset_sigind(sigind, offset):
    sigind += where(sigind == 0, 0, offset)


def cmp_sigind_cftrans(indmf, nspins):
    '''Loops over all correlated problems and constructs sigind, cftrans and legend.
    nspins - 1 or 2
    '''
    siginds = {}
    cftrans = {}
    legends = {}

    for icp,cp in indmf.cps.iteritems():
        iatoms = [iatom for iatom,L,qsplit in cp]
        Ls = set([L for iatom,L,qsplit in cp])
        qsplits = set([qsplit for iatom,L,qsplit in cp])

        # currently only support cps where all Ls and qsplits are the same
        if len(Ls) > 1:
            raise Exception, 'ERROR: each correlated problem must only have single L.'
        if len(qsplits) > 1:
            raise Exception, 'ERROR: each correlated problem must only have single qsplit.'

        L = Ls.pop()
        qsplit = qsplits.pop()
        nsites = len(iatoms)

        check_nspins(qsplit, nspins)
        siginds[icp] = cmp_sigind(qsplit, L, nsites, nspins)
        cftrans[icp] = cmp_cftrans(qsplit, L, nsites, nspins)
        legends[icp] = make_legend(qsplit, L, nsites, nspins)

    # offset the siginds so that the range spanned by each nonequivalent
    # correlated problem does not overlap with the range of any other ucp
    #
    maxind = {}   # mapping iucp -> max(sigind)
    for iucp,icps in indmf.ucps.iteritems():
        maxind[iucp] = max(siginds[icps[0]].flat)

    iucps = sorted(maxind.keys())
    maxinds = [maxind[iucp] for iucp in iucps]
    offsets = [0] + list(cumsum(maxinds)[:-1])
    offsets_dict = {}  # mapping iucp -> offset
    for iucp,offset in zip(iucps, offsets):
        offsets_dict[iucp] = offset

    for icp,cp in indmf.cps.iteritems():
        offset_sigind(siginds[icp], offsets_dict[indmf.iucps[icp]])

    return siginds, cftrans, legends


if __name__ == '__main__':
    # parse command line arguments
    options = parse_cmdline_args()
    w2kenv = W2kEnvironment()

    # automatically detect so/complex runs
    if not options.so and os.path.isfile(w2kenv.SCRATCH+"/"+w2kenv.case+".vectorso") and os.path.getsize(w2kenv.SCRATCH+"/"+w2kenv.case+".vectorso")>0 :
        print 'Found '+w2kenv.case+'.vectorso file, hence assuming so-coupling exists. Applying -so switch!'
        options.so = True

    # parse main input file case.indmf
    try:
        indmf = ga_indmffile.Indmf(w2kenv.case)
        indmf.read()
    except Exception, e:
        print "ERROR: reading file `%s.indmf` failed." % w2kenv.case
        print e
        exit(1)
    # instantiate case.indmfl object and initialize with values read from case.indmf
    inl = ga_indmffile.Indmfl(w2kenv.case)
    inl.copy_construct(indmf)

    # for each correlated problem, compute sigind, cftrans and create labels
    nspins = 2 if options.so else 1
    inl.siginds, inl.cftrans, inl.legends = cmp_sigind_cftrans(indmf, nspins)

    # write input file for x_dmft.py dmft0
    inl.write()
