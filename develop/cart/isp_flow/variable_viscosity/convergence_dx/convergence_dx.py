#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

################################################################################
# Data 

p06 = np.genfromtxt('p06.dat', names=True)
p11 = np.genfromtxt('p11.dat', names=True)
p16 = np.genfromtxt('p16.dat', names=True)

################################################################################
# Plot

plt.figure(figsize=(5.22, 4.8))

plt.suptitle(r' Nonlinear viscosity -- $h$-convergence', fontsize=14)

plt.xscale('log', basex=2)
plt.xlabel(r'$h$', size='14')
#plt.xlim(xmin=1./2**8, xmax=1/2**4)

plt.yscale('log')
plt.ylabel(r'$\varepsilon_{v}$', size='14')
plt.ylim(ymin=1e-11, ymax=1.20e-1)

# lines indicating order 

n = 12
m = 1
c = p06['err_v_rms'][n]
tau = np.array([ p06['dx'][n], 1.05 * p06['dx'][m] ])
err = np.array([ c           , c * (tau[1]/tau[0])**7 ])
l06, = plt.plot(tau, err, linestyle=':', color='C0')

n = 10
m = 1
c = 1.0 * p11['err_v_rms'][n]
tau = np.array([ 1.00 * p11['dx'][n], 1.00 * p11['dx'][m] ])
err = np.array([ c                  , c * (tau[1]/tau[0])**12 ])
l11, = plt.plot(tau, err, linestyle='--', color='C2')

n = 6
m = 1
c = 1.0 * p16['err_v_rms'][n]
tau = np.array([ 1.00 * p16['dx'][n], 1.00 * p16['dx'][m] ])
err = np.array([ c                  , c * (tau[1]/tau[0])**17 ])
l16, = plt.plot(tau, err, linestyle='-.', color='C3')

plt.plot(p06['dx'][0:], p06['err_v_rms'][0:], linestyle='None', marker='s', markersize=8, fillstyle=u'full', color='w')
plt.plot(p11['dx'][1:], p11['err_v_rms'][1:], linestyle='None', marker='D', markersize=8, fillstyle=u'full', color='w')
plt.plot(p16['dx'][1:], p16['err_v_rms'][1:], linestyle='None', marker='o', markersize=8, fillstyle=u'full', color='w')

m06, = plt.plot(p06['dx'][0:], p06['err_v_rms'][0:], linestyle='None', marker='s', markersize=8, fillstyle=u'none', color='C0')
m11, = plt.plot(p11['dx'][1:], p11['err_v_rms'][1:], linestyle='None', marker='D', markersize=8, fillstyle=u'none', color='C2')
m16, = plt.plot(p16['dx'][1:], p16['err_v_rms'][1:], linestyle='None', marker='o', markersize=8, fillstyle=u'none', color='C3')

plt.legend( ( m06        , m11         , m16         , l06      , l11       , l16       ), 
            ( '$P\!=\!6$', '$P\!=\!11$', '$P\!=\!16$', '$h^{7}$', '$h^{12}$', '$h^{17}$'), 
            loc='best', fontsize=12)

plt.tight_layout(rect=[-0.02, -0.02, 1.02, 0.96])
plt.savefig('convergence_dx.pdf')
 