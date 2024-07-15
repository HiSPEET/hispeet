#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica',
    'font.size': 12
})

print('\nScaled accuracy diagram for standalone or SDC time-integration method\n')

method = input('method: ')
ns = int(input('num subintervals = '))
sp = int(input('predictor stages = '))
sc = int(input('corrector stages = '))
nc = int(input('corrector sweeps = '))
plot_file = input('plot file (*.pdf) = ')

# work units
nw = max(ns,1) * (max(sp-1,1) + nc * max(sc-1,1))

print('\nnumber of works units = ', nw,'\n')

x = np.loadtxt('lambda_re.dat') / nw
y = np.loadtxt('lambda_im.dat') / nw
e = np.loadtxt('error.dat')

print('ranges')
print('  zr_min =',np.min(x))
print('  zr_max =',np.max(x))
print('  zi_min =',np.min(y))
print('  zi_max =',np.max(y))
print('  e_min  =',np.min(e))
print('  e_max  =',np.max(e))

print('\nselect')
x_min = float(input('  zr_min = '))
x_max = float(input('  zr_max = '))
y_min = float(input('  zi_min = '))
y_max = float(input('  zi_max = '))
e_lvl = float(input('  e_lvl  = '))

fig1, ax = plt.subplots()

diag = plt.contour( x, y, e, levels=[e_lvl], colors='red' 
                  , linestyles='-', linewidths= 0.75 ) 

handles, labels = diag.legend_elements()

labels = [r'$\varepsilon = $'+str(e_lvl)]
plt.legend(handles, labels)

plt.xlim(left=x_min, right=x_max)
plt.ylim(bottom=y_min, top=y_max)

plt.xlabel(r'Re($z_{\mathrm{s}}$)')
plt.ylabel(r'Im($z_{\mathrm{s}}$)')

fig1.tight_layout(rect=[0, 0, 1, 1])
fig1.savefig(plot_file+'.pdf')

print('\nscaled accuracy diagram saved to',plot_file+'.pdf','\n')
