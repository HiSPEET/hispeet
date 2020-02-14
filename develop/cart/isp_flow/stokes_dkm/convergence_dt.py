#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


sdc10 = np.genfromtxt('convergence_dt-sdc-1-0.dat', skip_header=0, names=True)
sdc11 = np.genfromtxt('convergence_dt-sdc-1-1.dat', skip_header=0, names=True)

plt.figure(figsize=(5.8, 4.8))

plt.plot(sdc10['dt'], sdc10['err_v_rms'], linestyle='None', marker='s', label='SDC(1,0)')
plt.plot(sdc11['dt'], sdc11['err_v_rms'], linestyle='None', marker='o', label='SDC(1,1)')

# line indicating order 1 for sdc10
n = sdc10['dt'].size
m = 1
c = 0.6 * sdc10['err_v_rms'][n-m]
tau = np.array([ sdc10['dt'][n-m], sdc10['dt'][3] ])
err = np.array([ c               , c * (tau[1]/tau[0]) ])
plt.plot(tau, err, linestyle='--', color='C0')
plt.text(8e-4, 9e-5, r'$\sim \Delta t$', color='C0', size='12')

# line indicating order 1 for CS
n = sdc11['dt'].size - 1
m = 2
c = 0.5 * sdc11['err_v_rms'][n-m]
tau = np.array([ sdc11['dt'][n-m], sdc11['dt'][2] ])
err = np.array([ c            , c * (tau[1]/tau[0])**2 ])
plt.plot(tau, err, linestyle='--', color='C1')
plt.text(2e-3, 5e-6, r'$\sim \Delta t^2$', color='C1', size='12')


plt.xscale('log', basex=10)
#plt.xlim(xmin=1e-6,xmax=1e-1)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_2$', size='14')
#plt.ylim(ymin=1e-5, ymax=3e-1)


plt.legend(loc='best')
plt.text(1e-3, 0.08, r'Velocity', color='black', size='16')


plt.tight_layout()


plt.savefig('convergence_dt.pdf')
