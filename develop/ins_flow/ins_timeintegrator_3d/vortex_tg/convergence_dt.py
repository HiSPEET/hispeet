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

# ------------------------------------------------------------------------------
# generic

if os.path.exists('convergence_dt.dat'):
    res = np.genfromtxt('convergence_dt.dat', skip_header=0, names=True)
    plt.plot(res['dt'], res['err_v'], linestyle='None', 
             markersize=7, fillstyle='none', marker='o', 
             label=r'$\varepsilon(\Delta t)$')

# ------------------------------------------------------------------------------
# comparison of results in files 'convergence_dt' + tag + '.dat'

tag    = [ 'euler', 'bdf2', 'ars443' ]
label  = [ 'Euler', 'BDF2', 'ARS(4,4,3)' ]
marker = [ 's', 'D', 'o']

for i in range(len(tag)):
    file = 'convergence_dt' + tag[i] + '.dat'
    if os.path.exists(file):
        res = np.genfromtxt(file[i], skip_header=0, names=True)
        plt.plot(res['dt'], res['err_v'], linestyle='None', 
                 markersize=7, fillstyle='none', marker=marker[i], 
                 label=label[i])

#=================================================================================

plt.xscale('log', base=10)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log', base=10)
plt.ylabel(r'$\varepsilon_v$', size='14')

plt.legend(loc=0, ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.0, 0.95])
plt.suptitle(r'Vortex TG -- convergence in time')
plt.savefig('convergence_dt.pdf')
