#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica',
    'font.size': 12
})

print('\nStability diagram for MLSDC time-integration method\n')

method = input('method: ')
nl = int(input('num levels = '))
plot_file = input('plot file  = ')

x = np.loadtxt('lambda_re.dat') 
y = np.loadtxt('lambda_im.dat') 

fig, ax = plt.subplots()
fig.suptitle(method)


for i in range(nl):
    l    = 1 + i
    cols = 'C' + str(l)
    file = 'amplification_level_' + str(l) + '.dat'
    print(file)
    a = np.loadtxt(file)
    line = ax.contour(x, y, a, levels=[1], colors=cols, linestyles='-') 
    handle, label = line.legend_elements()
    if i == 0:
        handles = [handle[0]]
        labels  = [r'$l =$' + ' ' + str(l)]
    else:
        handles.append(handle[0])
        labels.append(r'$l =$' + ' ' + str(l))


plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')
plt.legend(handles, labels, loc = 'upper left') # loc = 'upper left' 'upper right'

plt.xlabel(r'Re($z$)')
plt.ylabel(r'Im($z$)')

#plt.grid(which='both', axis='y')

plt.tight_layout(rect=[-0.02, -0.02, 1.03, 1.05])
plt.savefig(plot_file+'.pdf')
#plt.show()

print('\nstability diagram saved to',plot_file+'.pdf','\n')
