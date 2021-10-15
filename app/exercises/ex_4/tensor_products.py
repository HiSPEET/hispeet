#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

inpfile = 'tensor_products.dat'
outfile = 'tensor_products.pdf'

A = np.loadtxt(inpfile);

po     = [int(A[0,0]  ), int(A[0,1]  )]
n_plot = [int(A[0,2]+1), int(A[0,3]+1)]

# Extract the data .............................................................
x      = np.reshape(A[1::,0],(n_plot[0],n_plot[1]))
y      = np.reshape(A[1::,1],(n_plot[0],n_plot[1]))
u_h    = np.reshape(A[1::,2],(n_plot[0],n_plot[1]))
du_h_1 = np.reshape(A[1::,3],(n_plot[0],n_plot[1]))
du_h_2 = np.reshape(A[1::,4],(n_plot[0],n_plot[1]))
err_h  = np.reshape(A[1::,5],(n_plot[0],n_plot[1]))


# generate figure ..............................................................
# Create an A4-sized figures with landscape mode (sizes in inch)
plt.rc('figure', figsize=(11.69,8.27))
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

#...............................................................................
# A short function for the plotting

def CreatePlot(ax,x,y,variable,variable_name):
    # elevation and rotation of the plots
    elevation = 30
    rotation  = 45
    ax.view_init(elev=elevation, azim=rotation)

    # Create a surface plot with the variable, use a 'hot' color map
    ax.plot_surface(x,y,variable,linewidth=1, cmap = plt.get_cmap("hot"))

    # Set the labels for the coordinates, use the title for the variable name
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.set_title (variable_name)

fig = plt.figure()

# Create a plot for the variable u in the axis 221
CreatePlot(fig.add_subplot(221, projection='3d'),x,y,u_h,  r'$u$')

# Create a plot for the error in the axis 221
CreatePlot(fig.add_subplot(222, projection='3d'),x,y,err_h,r'$\varepsilon$')

# Create a plot for the derivatives
CreatePlot(fig.add_subplot(223, projection='3d'),x,y,du_h_1,r'$\partial_1 u$')
CreatePlot(fig.add_subplot(224, projection='3d'),x,y,du_h_2,r'$\partial_2 u$')

fig.savefig(outfile, layout = 'tight')
