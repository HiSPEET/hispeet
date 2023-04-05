#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


res_0 = np.genfromtxt('convergence_dt_0.dat', skip_header=0, names=True)
res_1 = np.genfromtxt('convergence_dt_1.dat', skip_header=0, names=True)
res_2 = np.genfromtxt('convergence_dt_2.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'convergence in time', fontsize=14)

plt.plot(res_0['dt'], res_0['err_2'], linestyle='None', marker='o')
plt.plot(res_1['dt'], res_1['err_2'], linestyle='None', marker='o')
plt.plot(res_2['dt'], res_2['err_2'], linestyle='None', marker='o')

# line indicating a slope of 2
n = res_0['dt'].size - 1
x = np.array([ res_0['dt'][n], res_0['dt'][0] ])
y = res_0['err_2'][n] * np.array([ 1.0, (x[1] / x[0])**3 ])
plt.plot(x, y, linestyle=':', marker='None', color='C0', label=r'$\sim\Delta t^3$')

plt.xscale('log', basex=2)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_2$', size='14')

plt.legend(loc='best', ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('convergence_dt.pdf')
