#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica',
    'font.size': 12
})

print('\nStability diagram for standalone or SDC time-integration method\n')

method = input('method: ')
ns = int(input('num subintervals = '))
sp = int(input('predictor stages = '))
sc = int(input('corrector stages = '))
nc = int(input('corrector sweeps = '))
plot_file = input('plot file base name = ')

# work units
nw = max(ns,1) * (max(sp-1,1) + nc * max(sc-1,1))

print('\nnumber of works units = ', nw,'\n')

# unscaled ---------------------------------------------------------------------
                       
x = np.loadtxt('lambda_re.dat') 
y = np.loadtxt('lambda_im.dat') 
a = np.loadtxt('amplification.dat')

fig1, ax = plt.subplots()
fig1.suptitle(method)

diag = plt.contour( x, y, a, levels=[1], linestyles='-')

handles, labels = diag.legend_elements()

plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')

labels = [r'$|R| = 1$']
plt.legend(handles, labels, loc = 'upper left')

#plt.xlim(left=-10, right=2)
#plt.ylim(bottom=0, top=8)

plt.xlabel(r'Re($z$)')
plt.ylabel(r'Im($z$)')

#plt.grid(which='both', axis='y')

fig1.tight_layout(rect=[0, 0, 1, 1])
fig1.savefig(plot_file+'.pdf')

print('\nunscaled stability diagram saved to',plot_file+'.pdf','\n')

# scaled -----------------------------------------------------------------------

fig2, ax = plt.subplots()
fig2.suptitle(method)

diag = plt.contour( x/nw, y/nw, a, levels=[1], linestyles='-') 
# additional options
# linewidths= 0.75
# colors='red' 

handles, labels = diag.legend_elements()

plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')

labels = [r'$|R| = 1$']
plt.legend(handles, labels, loc = 'upper left')

#plt.xlim(left=-0.6, right=0.2)
#plt.ylim(bottom=0.0, top=1.6)

plt.xlabel(r'Re($z_{\mathrm{s}}$)')
plt.ylabel(r'Im($z_{\mathrm{s}}$)')

fig2.tight_layout(rect=[0, 0, 1, 1])
fig2.savefig(plot_file+'_scaled.pdf')

print('\nscaled stability diagram saved to',plot_file+'_scaled.pdf','\n')
