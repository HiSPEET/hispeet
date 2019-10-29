#!/usr/bin/python3

import numpy             as np
import matplotlib.pyplot as plt

plt.rc('figure', figsize=(8.27,11.69))

# read input ...................................................................
data = np.loadtxt('study_input')

# extract the data
all_po        = data[:,0]
all_n_element = data[:,1]
all_errors    = data[:,2]

polynomial_degrees = np.unique(all_po)

# plot .........................................................................

fig, (ax1) = plt.subplots(1, 1)

for po in polynomial_degrees:
    n_element = all_n_element[all_po == po]
    error     = all_errors   [all_po == po]
    plt.loglog(n_element,error,'x-',label = '$p = ' + str(po) + '$')

ax1.set_title ('Errors')
ax1.set_xlabel('$n_{e}$ [-]')
ax1.set_ylabel('$\\varepsilon$ [-]')
ax1.legend()


plt.tight_layout()
plt.savefig('plot_exercise_2_study.pdf')


