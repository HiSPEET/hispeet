#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


sdc_0 = np.genfromtxt('sdc_0-0.dat', skip_header=0, names=True)
sdc_1 = np.genfromtxt('sdc_1-3.dat', skip_header=0, names=True)
sdc_2 = np.genfromtxt('sdc_2-6.dat', skip_header=0, names=True)
sdc_3 = np.genfromtxt('sdc_3-9.dat', skip_header=0, names=True)


plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Vortex TG -- SDC(M,3M) convergence', fontsize=14)

plt.plot(sdc_0['dt'], sdc_0['err_v_rms'], linestyle='None', color='C0', marker='$0$', markersize=8)
plt.plot(sdc_1['dt'], sdc_1['err_v_rms'], linestyle='None', color='C1', marker='$1$', markersize=8)
plt.plot(sdc_2['dt'], sdc_2['err_v_rms'], linestyle='None', color='C2', marker='$2$', markersize=8)
plt.plot(sdc_3['dt'], sdc_3['err_v_rms'], linestyle='None', color='C3', marker='$3$', markersize=8)

plt.xscale('log', basex=2)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_{v}$', size='14')

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('sdc_m.pdf')
