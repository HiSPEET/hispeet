#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


sdc_3_0 = np.genfromtxt('sdc_3-0.dat', skip_header=0, names=True)
sdc_3_1 = np.genfromtxt('sdc_3-1.dat', skip_header=0, names=True)
sdc_3_2 = np.genfromtxt('sdc_3-2.dat', skip_header=0, names=True)
sdc_3_3 = np.genfromtxt('sdc_3-3.dat', skip_header=0, names=True)
sdc_3_4 = np.genfromtxt('sdc_3-4.dat', skip_header=0, names=True)
sdc_3_5 = np.genfromtxt('sdc_3-5.dat', skip_header=0, names=True)
sdc_3_6 = np.genfromtxt('sdc_3-6.dat', skip_header=0, names=True)
sdc_3_9 = np.genfromtxt('sdc_3-9.dat', skip_header=0, names=True)


plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Vortex TG -- SDC(3,k) convergence', fontsize=14)

plt.plot(sdc_3_0['dt'], sdc_3_0['err_v_rms'], linestyle='None', color='C0', marker='$0$', markersize=8)
plt.plot(sdc_3_1['dt'], sdc_3_1['err_v_rms'], linestyle='None', color='C1', marker='$1$', markersize=8)
plt.plot(sdc_3_2['dt'], sdc_3_2['err_v_rms'], linestyle='None', color='C2', marker='$2$', markersize=8)
plt.plot(sdc_3_3['dt'], sdc_3_3['err_v_rms'], linestyle='None', color='C3', marker='$3$', markersize=8)
plt.plot(sdc_3_4['dt'], sdc_3_4['err_v_rms'], linestyle='None', color='C4', marker='$4$', markersize=8)
plt.plot(sdc_3_5['dt'], sdc_3_5['err_v_rms'], linestyle='None', color='C5', marker='$5$', markersize=8)
plt.plot(sdc_3_6['dt'], sdc_3_6['err_v_rms'], linestyle='None', color='C6', marker='$6$', markersize=8)
plt.plot(sdc_3_9['dt'], sdc_3_9['err_v_rms'], linestyle='None', color='C7', marker='$9$', markersize=8)

plt.xscale('log', basex=2)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_{v}$', size='14')

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('sdc_k.pdf')
