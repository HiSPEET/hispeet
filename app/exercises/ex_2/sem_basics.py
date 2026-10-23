#!/usr/bin/python3

import numpy             as np
import matplotlib.pyplot as plt

plt.rc('figure', figsize=(8.27,11.69))

# read input ...................................................................
data  = np.loadtxt("e1_data")
nodal = np.loadtxt("e1_nodal")

# extract data .................................................................

po     = data[0,0]
poly   = data[0,1]
n_plot = data[0,2] + 1

x       = data[1::,0]
u       = data[1::,1]
du      = data[1::,2]
ux      = data[1::,3]
dux     = data[1::,4]

# Select the polynomial family
if   (poly == 1):
    polynomial = 'GLL'
elif (poly == 2):
    polynomial = 'GL'
elif (poly == 3):
    polynomial = 'Lagrange'
else:
    polynomial = 'Unknown'

title_str = polynomial + '-interpolation using a polynomial degree of '        \
          + str(int(po)) + ' and ' + str(int(n_plot)) + ' plot nodes'

# plot .........................................................................

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

ax1.plot(x,u ,  'x-', label = 'Approximation')
ax1.plot(x,ux,  '-',  label = 'Exact')

ax1.set_title (title_str)
ax1.set_xlabel('x [-]')
ax1.set_ylabel('values [-]')
ax1.legend()

ax2.plot(x,du , 'x-', label = 'Approximation')
ax2.plot(x,dux, '-',  label = 'Exact')

ax2.set_title ('Differentiation')
ax2.set_xlabel('x [-]')
ax2.set_ylabel('values [-]')
ax2.legend()


for i in np.arange(1,po+2):
    ax3.plot(nodal[0,:],nodal[int(i),:])

ax3.set_title ('Basis functions')

plt.tight_layout()
plt.savefig('sem_basics.pdf')


