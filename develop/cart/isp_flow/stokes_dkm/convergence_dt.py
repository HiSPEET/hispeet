#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math


sdc0 = np.genfromtxt('convergence_dt-sdc-0.dat', skip_header=0, names=True)
sdc1 = np.genfromtxt('convergence_dt-sdc-1.dat', skip_header=0, names=True)
sdc2 = np.genfromtxt('convergence_dt-sdc-2.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Stokes DKM -- temporal convergence', fontsize=14)

plt.plot(sdc0['dt'], sdc0['err_v_rms'], linestyle='None', marker='v', label='SDC(1,0)')
plt.plot(sdc1['dt'], sdc1['err_v_rms'], linestyle='None', marker='s', label='SDC(1,3)')
plt.plot(sdc2['dt'], sdc2['err_v_rms'], linestyle='None', marker='o', label='SDC(2,6)')

# lines indicating order 

n = 6
m = 2
c = 0.75 * sdc0['err_v_rms'][n]
tau = np.array([ sdc0['dt'][n], sdc0['dt'][m] ])
err = np.array([ c            , c * (tau[1]/tau[0])**1 ])
plt.plot(tau, err, linestyle=':', color='C0', label=r'$\Delta t^{1}$')

n = 6
m = 1
c = 0.75 * sdc1['err_v_rms'][n]
tau = np.array([ sdc1['dt'][n], sdc1['dt'][m] ])
err = np.array([ c            , c * (tau[1]/tau[0])**2 ])
plt.plot(tau, err, linestyle='-.', color='C1', label=r'$\Delta t^{2}$')

n = 3
m = 1
c = 0.75 * sdc2['err_v_rms'][n]
tau = np.array([ sdc2['dt'][n], sdc2['dt'][m] ])
err = np.array([ c            , c * (tau[1]/tau[0])**3 ])
plt.plot(tau, err, linestyle='-.', color='C2', label=r'$\Delta t^{3}$')


plt.xscale('log', basex=10)
plt.xlabel(r'$\Delta t$', size='14')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_2$', size='14')

plt.legend(loc='best', ncol=2)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('convergence_dt.pdf')
