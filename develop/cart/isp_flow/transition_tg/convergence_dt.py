#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

#IMEX Euler
res= np.genfromtxt('convergence_imex_euler.dat', skip_header=0, names=True)

plt.figure(figsize=(5.22, 4.8))
plt.suptitle(r'Vortex TG--PPP--convergence in time', fontsize=12)

plt.plot(res['dt'], res['err_v_rms'], linestyle='None', marker='o',color='C1',label=r'IMEX Euler') #$\sim\Delta t$

# line indicating a slope of 1
n  = res['dt'].size - 1
x = np.array([ res['dt'][n], res['dt'][0] ])
y = res['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0]) ])
plt.plot(x, y, linestyle=':', marker='None', color='C1')



#RK (ns=3,method=1)
res2 = np.genfromtxt('convergence_rk_3.dat', skip_header=0, names=True)

plt.plot(res['dt'], res2['err_v_rms'], linestyle='None', marker='o',color='C3', label=r'RK(3,1)') #$\sim\Delta t^{2}$
# line indicating ......
n = res['dt'].size - 1
x = np.array([ res['dt'][n], res['dt'][0] ])
y2= res2['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**2 ])
plt.plot(x, y2, linestyle=':', marker='None', color='C3')

#RK (ns=4,method=1)
res41 = np.genfromtxt('convergence_rk_4_1.dat', skip_header=0, names=True)

plt.plot(res['dt'], res41['err_v_rms'], linestyle='None', marker='o',color='C4', label=r'RK(4,1)')
# line indicating ...
n = res41['dt'].size - 3
x = np.array([ res['dt'][n], res['dt'][0] ])
y41 = res41['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**3 ])
plt.plot(x, y41, linestyle=':', marker='None', color='C4')

#RK (ns=4,method=2)
res42 = np.genfromtxt('convergence_rk_4_2.dat', skip_header=0, names=True)

plt.plot(res['dt'], res42['err_v_rms'], linestyle='None', marker='o',color='C5', label=r'RK(4,2) ')
# line indicating ....
n = res42['dt'].size - 3
x = np.array([ res['dt'][n], res['dt'][0] ])
y42 = res42['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**3 ])
plt.plot(x, y42, linestyle=':', marker='None', color='C5')


#RK(ns=6,method=1)
res6 = np.genfromtxt('convergence_rk_6_1.dat', skip_header=0, names=True)

plt.plot(res['dt'], res6['err_v_rms'], linestyle='None', marker='o',color='C6', label=r'RK(6,1) ')
# line indicating a slope of 2
n = res['dt'].size - 5
x  = np.array([ res['dt'][n], res['dt'][0] ])
y6 = res6['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**3.5 ])
plt.plot(x, y6, linestyle=':', marker='None', color='C6')   

#SDC (2,4)
res1 = np.genfromtxt('convergence_sdc_2_4.dat', skip_header=0, names=True)

plt.plot(res['dt'], res1['err_v_rms'], linestyle='None', marker='o',color='C2', label=r'SDC(2,4) ')
# line indicating ....
n  = res['dt'].size - 4
x  = np.array([ res['dt'][n], res['dt'][0] ])
y1 = res1['err_v_rms'][n] * np.array([ 1.0, (x[1] / x[0])**3.3 ])
plt.plot(x, y1, linestyle=':', marker='None', color='C2')




#=================================================================================
plt.xscale('log', basex=10)
plt.xlabel(r'$\Delta t$', size='12')

plt.yscale('log')
plt.ylabel(r'$\varepsilon_v$', size='12')

plt.legend(loc=0, ncol=1)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.95])
plt.savefig('convergence_dt_1.pdf')
