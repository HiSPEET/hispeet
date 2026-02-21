#! /usr/bin/env python3

import matplotlib as mpl
mpl.rc("xtick", labelsize=14)
mpl.rc("ytick", labelsize=14)

import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage[sfdefault]{notomath}"

import numpy as np
import math

#res = np.genfromtxt('test_hp_coarsening_mode1.dat', skip_header=0, names=True)
res = np.genfromtxt('test_hp_coarsening_mode2.dat', skip_header=0, names=True)

#  x     sampling points
#  u_f   fine (original) function
#  u_c   coarse (projected) function

fig, ax = plt.subplots()
fig.set_figheight(4.8)
fig.set_figwidth(5.8)

#plt.plot(res['x'], res['f']   , linestyle=':', label=r'$f$')
plt.plot(res['x'], res['u_f'], linestyle='--', label=r'$u_{\mathrm f}$')
plt.plot(res['x'], res['u_c'], linestyle=':' , label=r'$u_{\mathrm c}$')


plt.xlabel(r'$x$', size='16')
plt.ylabel(r'$u$', size='16')

plt.legend(loc='lower right', fontsize=14)
plt.tight_layout(rect=[-0.01, -0.02, 1.01, 1.01])
plt.savefig('test_hp_coarsening.pdf')
