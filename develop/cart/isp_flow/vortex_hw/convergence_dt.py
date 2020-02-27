#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


res = np.genfromtxt('convergence_dt.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Vortex HW -- convergence in time', fontsize=14)

plt.plot(res['dt'], res['err_v_rms'], linestyle='None', marker='o')

# line indicating a slope of 2
n = res['dt'].size - 3
x = np.array([ res['dt'][n], res['dt'][0] ])
y = res['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**2 ])
plt.plot(x, y, linestyle=':', marker='None', color='C0', label=r'$\sim\Delta t^2$')

plt.xscale('log', basex=2)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_v$', size='14')

plt.legend(loc='best', ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('convergence_dt.pdf')
