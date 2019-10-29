#!/usr/bin/python2

import numpy             as np
import matplotlib.pyplot as plt

plt.rc('figure', figsize=(8.27,11.69))

# read input ...................................................................
plot_data  = np.loadtxt('output_e2_data')
nodal_data = np.loadtxt('output_e2_nodal')

# extract data .................................................................

po        = plot_data[0,0]
n_element = plot_data[0,1]
n_plot    = plot_data[0,2]

x_p       = plot_data[1::,0]
u_h_p     = plot_data[1::,1]
u_p       = plot_data[1::,2]

x_h = nodal_data[:,0]
u_h = nodal_data[:,1]
u   = nodal_data[:,2]

# plot .........................................................................

fig, (ax1,ax2) = plt.subplots(2, 1)

# First plot: Only nodal values
ax1.plot(x_h,u_h,  '-x', label = 'Approximation', color = "blue")
ax1.plot(x_h,u,    '-o', label = 'Exact'        , color = "green")

ax1.set_title ('Linearly interpolated nodal values')
ax1.set_xlabel('x [-]')
ax1.set_ylabel('Values [-]')
ax1.legend()

# Second plot: nodal values + interpolation
ax2.plot(x_p,u_h_p,  '-', label = 'Approximation', color = "blue")
ax2.plot(x_h,u_h,    'x',                          color = "blue")
ax2.plot(x_p,u_p,    '-', label = 'Exact',         color = "green")
ax2.plot(x_h,u,      'o',                          color = "green")

ax2.set_title ('Interpolated values')
ax2.set_xlabel('x [-]')
ax2.set_ylabel('Values [-]')
ax2.legend()


plt.tight_layout()
plt.savefig('plot_exercise_2.pdf')


