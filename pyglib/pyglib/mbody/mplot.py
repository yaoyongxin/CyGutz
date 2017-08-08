from __future__ import print_function
import matplotlib.pyplot as plt
from builtins import zip
import numpy, h5py
from itertools import product
from pyglib.basic.splot import colors

def s_to_index(s):
    '''get integer index for s, which can either be integer or halves.
    '''
    if 0.49 < s % 1 < 0.51: # half s, 1/2, 3/2, ...
        return int(s*2 + 0.1) - 1
    elif 0 <= s % 1 < 0.01 or 0.99 < s % 1:
        return int(s + 0.1)
    else:
        raise ValueError(' s = {}'.format(s))


def get_group_weights(weights, n_labels, s_labels):
    '''group weights according to valence n and s quantum number.
    '''
    nmin = numpy.min(n_labels)
    nmax = numpy.max(n_labels)
    smax = numpy.max(s_labels)
    ns = s_to_index(smax)
    g_weights = numpy.zeros([nmax-nmin+1, ns+1])
    for wt, n, s in zip(weights, n_labels, s_labels):
        if wt < 1.e-7:
            continue
        try:
            i = s_to_index(s)
        except ValueError:
            print(' wt = {}: to be ignored!'.format(wt))
            continue
        g_weights[n-nmin, i] += wt
    return g_weights


def hist_ns(weights, n_labels, s_labels, amlabel='S', tol=1.e-3):
    '''generate histogram plot for the local reduced density matrix
    with label n and s(j).
    '''
    x_ticks = []
    label_ticks = []
    icolor = 0
    color_map = [{}, {}]
    nmin = numpy.min(n_labels)
    nmax = numpy.max(n_labels)
    g_weights = get_group_weights(weights, n_labels, s_labels)
    idx = -1
    fig, ax = plt.subplots(figsize=(5,4))
    for n in range(g_weights.shape[0]):
        idx_start = idx
        for s in range(g_weights.shape[1]):
            w = g_weights[n, s]
            jn = (n+nmin) % 2
            if w >= tol:
                idx += 1
                if s in color_map[jn]:
                    color = color_map[jn][s]
                    label = None
                else:
                    icolor += 1
                    color = color_map[jn][s] = colors[icolor]
                    if jn == 0:
                        label = '{}={}'.format(amlabel, s)
                    else:
                        label = '{}={}/2'.format(amlabel, 2*s+1)
                ax.hist([idx], 1, weights=[w], rwidth=0.5, align='mid',
                        color=color, label=label)
        if idx > idx_start:
            ax.axvspan(idx_start+0.5, idx+0.5, facecolor='g',
                    alpha=float(n)/g_weights.shape[0]*0.8, zorder=-1)
            x_ticks.append((idx_start+idx+1)/2.)
            label_ticks.append('N='+str(n+nmin))
    ax.set_xlim(-0.5, idx+0.5)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(label_ticks)
    ax.legend()
    ax.set_ylabel('Probability')
    fig.tight_layout()
    fig.savefig('hist_s.pdf')


def hist_ns2(weights, n_labels, s_labels, amlabel='S', tol=1.e-3):
    '''generate histogram plot for the multiple local reduced density matrix
    with label n and s(j).
    '''
    x_ticks = []
    label_ticks = []
    icolor = 0
    color_map = [{}, {}]
    nmin = numpy.min(n_labels)
    nmax = numpy.max(n_labels)
    g_weights = get_group_weights(weights, n_labels, s_labels)
    idx = -1
    fig, ax = plt.subplots(figsize=(5,4))
    for n in range(g_weights.shape[0]):
        idx_start = idx
        for s in range(g_weights.shape[1]):
            w = g_weights[n, s]
            jn = (n+nmin) % 2
            if w >= tol:
                idx += 1
                if s in color_map[jn]:
                    color = color_map[jn][s]
                    label = None
                else:
                    icolor += 1
                    color = color_map[jn][s] = colors[icolor]
                    if jn == 0:
                        label = '{}={}'.format(amlabel, s)
                    else:
                        label = '{}={}/2'.format(amlabel, 2*s+1)
                ax.hist([idx], 1, weights=[w], rwidth=0.5, align='mid',
                        color=color, label=label)
        if idx > idx_start:
            ax.axvspan(idx_start+0.5, idx+0.5, facecolor='g',
                    alpha=float(n)/g_weights.shape[0]*0.8, zorder=-1)
            x_ticks.append((idx_start+idx+1)/2.)
            label_ticks.append('N='+str(n+nmin))
    ax.set_xlim(-0.5, idx+0.5)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(label_ticks)
    ax.legend()
    ax.set_ylabel('Probability')
    fig.tight_layout()
    fig.savefig('hist_s.pdf')



if __name__=='__main__':
    with h5py.File('multiplets.h5', 'r') as f:
        weights, n_labels, s_labels = f['/impurity_1/weights'][()], \
                f['/impurity_1/n_labels'][()], f['/impurity_1/s_labels'][()]
    hist_ns(weights, n_labels, s_labels)
