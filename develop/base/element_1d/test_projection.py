#! /usr/bin/env python3

import matplotlib as mpl
mpl.rc("xtick", labelsize=14)
mpl.rc("ytick", labelsize=14)

import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage[sfdefault]{notomath}"

import numpy as np
import math

res = np.genfromtxt('test_projection.dat', skip_header=0, names=True)

#  x        sampling points
#  f        exact function
#  f_cs     interpolant of exact function given at collocation points
#  f_fs     filtered function
#  f_ps     projected function
#  e_i      interpolation error
#  e_p      projection error
#  e_f      filtering error
#  e_fp     filtering + projection error

fig, ax = plt.subplots()
fig.set_figheight(4.8)
fig.set_figwidth(5.8)

#plt.plot(res['x'], res['f']   , linestyle=':', label=r'$f$')
plt.plot(res['x'], res['f_cs'], linestyle='--', label=r'$If$')
plt.plot(res['x'], res['f_ps'], linestyle=':' , label=r'$PIf$')


plt.xlabel(r'$x$', size='16')
plt.ylabel(r'$f$', size='16')

plt.legend(loc='lower right', fontsize=14)
plt.tight_layout(rect=[-0.01, -0.02, 1.01, 1.01])
plt.savefig('test_projection.pdf')
