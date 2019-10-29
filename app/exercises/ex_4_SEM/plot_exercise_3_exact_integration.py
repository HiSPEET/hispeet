#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

# Create an A4-sized figures (sizes in inch)
plt.rc('figure', figsize=(8.27,11.69))

inpfile = 'output_exercise_3_exact_integration.dat'
outfile = 'plot_exercise_3_exact_integration.pdf'

# Extract the data .............................................................

A = np.loadtxt(inpfile);

po       = A[:,0];
int_appr = A[:,1];
int_ex   = A[:,2];
err_appr = A[:,3];
err_ex   = A[:,4];

# Plotting .....................................................................
fig = plt.figure()

# First plot: Integrals
plt.subplot(2, 1, 1)

plt.plot(po,int_appr, 'x', label = 'Approximative mass matrix')
plt.plot(po,int_ex,   'o', label = 'Exact mass matrix')

plt.legend()
plt.title('Integration of a polynomial of degree 2p')
plt.xlabel('p')
plt.ylabel('Integral')

# Second plot: Errors
plt.subplot(2, 1, 2)
plt.loglog(po,err_appr, 'x', label = 'Approximative mass matrix')
plt.loglog(po,err_ex,   'o', label = 'Exact mass matrix')

plt.legend()
plt.xlabel('p')
plt.ylabel('Error of integral')

# save figure
fig.savefig(outfile)
