import numpy             as np
import matplotlib.pyplot as plt

data = np.loadtxt('jacobi_method.dat')

x = data[0::,0]
u = data[0::,1]

plt.plot(x,u)
plt.savefig("jacobi_method.pdf")
