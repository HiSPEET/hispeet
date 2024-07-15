#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica',
    'font.size': 12
})

print('\nStability diagram for range of sweeps from K = K_min to K_max\n')
  
K_min      = int(input('K_min = '))
K_max      = int(input('K_max = '))
ampl_file  = input('data files (*-K.dat) = ')
plot_title = input('plot title = ')
plot_file  = input('plot file (*.pdf) = ')

x = np.loadtxt('lambda_re.dat') 
y = np.loadtxt('lambda_im.dat') 

fig, ax = plt.subplots()
fig.set_figheight(4.8)
fig.set_figwidth(5.8) 
fig.suptitle(plot_title)


for i in range(K_max+1 - K_min):
    k = i + K_min
    cols = 'C' + str(i+1)
    file = ampl_file + '-' + str(k) + '.dat'
    print(file)
    a = np.loadtxt(file)
    line = plt.contour(x, y, a, levels=[1], colors=cols, linestyles='-') 
    handle, label = line.legend_elements()
    if i == 0:
        handles = [handle[0]]
        labels  = [r'$K =$' + ' ' + str(k)]
    else:
        handles.append(handle[0])
        labels.append(r'$K =$' + ' ' + str(k))


plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')
plt.legend(handles, labels, ) # loc = 'upper left' 'upper right'

#plt.xlim(left=-5e-7, right=1e-7)
#plt.ylim(bottom=0, top=10)

plt.xlabel(r'Re($z$)')
plt.ylabel(r'Im($z$)')

#plt.grid(which='both', axis='y')

plt.tight_layout(rect=[-0.02, -0.02, 1.03, 1.05])
fig.savefig(plot_file+'.pdf')

print('\nstability diagram saved to',plot_file+'.pdf','\n')
