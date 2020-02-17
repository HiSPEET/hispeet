#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


sdc10 = np.genfromtxt('stability_dt-sdc-1-0.dat', skip_header=0, names=True)
sdc11 = np.genfromtxt('stability_dt-sdc-1-1.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Stokes DKM -- temporal stability', fontsize=14)

plt.plot(sdc10['dt'], sdc10['err_v_rms'], linestyle='None', marker='s', label='SDC(1,0)')
plt.plot(sdc11['dt'], sdc11['err_v_rms'], linestyle='None', marker='o', label='SDC(1,1)')

plt.xscale('log', basex=10)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_v$', size='14')

plt.legend(loc='best', ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('stability_dt.pdf')
