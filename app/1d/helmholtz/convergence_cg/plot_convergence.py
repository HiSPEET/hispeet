#! /usr/bin/env python3
#===============================================================================
# 1D Helmholtz SEM convergence study: maximum error e_max vs. number of
# unknowns N = po * ne
#
# Reads series_p1.dat, series_p2.dat, series_ne1.dat (run_convergence.sh)
# and writes convergence.pdf
#===============================================================================

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl

mpl.rcParams["font.size"] = 13
mpl.rc('text', usetex=True)

#===============================================================================
# data

def load(name):
    dat = np.genfromtxt(name, skip_header=0, names=True)
    N = dat['po'] * dat['ne']
    return N, dat['e_max']

N_p1 , e_p1  = load('series_p1.dat')
N_p2 , e_p2  = load('series_p2.dat')
N_ne1, e_ne1 = load('series_ne1.dat')

#===============================================================================
# reference slopes for h-convergence
#
# The asymptotic rate is estimated from the last points above the round-off
# level and rounded to an integer k. The reference line ~ N^(-k) is drawn
# below the data over the asymptotic range.

def reference(N, e, n_fit=4, e_min=1e-12, shift=0.25):
    mask = e > e_min
    N, e = N[mask], e[mask]
    N, e = N[-n_fit:], e[-n_fit:]
    s = np.polyfit(np.log(N), np.log(e), 1)[0]
    k = int(round(-s))
    Nr = np.array([N[0], N[-1]])
    er = shift * e[-1] * (Nr / N[-1])**(-k)
    return Nr, er, k, s

def plot_reference(N, e, **kw):
    Nr, er, k, s = reference(N, e, **kw)
    plt.plot(Nr, er, linestyle='--', color='gray', linewidth=1)
    xm = np.sqrt(Nr[0] * Nr[1])
    ym = np.sqrt(er[0] * er[1])
    plt.text(xm, ym / 2.5, r'$N^{-%d}$' % k, color='gray',
             ha='center', va='top', fontsize=12)
    print('fitted slope: %6.3f  ->  reference N^-%d' % (s, k))

#===============================================================================
# plot e_max(N)

fig,ax = plt.subplots()
fig.set_figheight(4.8)
fig.set_figwidth(5.8)

plt.loglog(N_p1, e_p1, linestyle=':',
           marker='v', markersize=7, fillstyle='none',
           label=r'$p = 1$')

plt.loglog(N_p2, e_p2, linestyle=':',
           marker='o', markersize=7, fillstyle='none',
           label=r'$p = 2$')

plt.loglog(N_ne1, e_ne1, linestyle=':',
           marker='s', markersize=7, fillstyle='none',
           label=r'$n_{\textsc e} = 1$')

plot_reference(N_p1, e_p1)
plot_reference(N_p2, e_p2)

plt.xlabel(r'$N$', fontsize=14)
plt.ylabel(r'$e_{\max}$', fontsize=14)

plt.legend(alignment='right')

plt.suptitle(r'1D Helmholtz: SEM convergence')
plt.tight_layout(rect=[-0.03, -0.02, 1.028, 1.05])
plt.savefig('convergence.pdf')
