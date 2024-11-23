#! /usr/bin/env python3

import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 12

#===============================================================================
   
fig,ax = plt.subplots()
fig.set_figheight(4.8)
fig.set_figwidth(5.8)

res = np.genfromtxt('convergence_dx.dat', skip_header=0, names=True)

plt.plot(res['dx_min'], res['err_v'], linestyle='None', 
         markersize=7, fillstyle='none', marker='o', 
         label=r'$\varepsilon(\Delta x)$')

plt.xscale('log', base=10)
plt.xlabel(r'$\Delta x$', size='14')

plt.yscale('log', base=10)
plt.ylabel(r'$\varepsilon_v$', size='14')

plt.legend(loc=0, ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.0, 0.95])
plt.suptitle(r'Vortex TG -- convergence in space')
plt.savefig('convergence_dx.pdf')
