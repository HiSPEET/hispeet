#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


res = np.genfromtxt('convergence_dx.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Vortex HW -- $h$-convergence with $P=5$', fontsize=14)

plt.plot(res['dx'], res['err_v_rms'], linestyle='None', marker='o')

# line indicating a slope of 5
n = res['dt'].size - 1
x = np.array([ res['dx'][n], res['dx'][0] ])
y = res['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**5 ])
plt.plot(x, y, linestyle=':', marker='None', color='C0', label=r'$\sim\Delta x^5$')

plt.xscale('log', basex=2)
plt.xlabel(r'$\Delta x$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_v$', size='14')

plt.legend(loc='best', ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('convergence_dx.pdf')
