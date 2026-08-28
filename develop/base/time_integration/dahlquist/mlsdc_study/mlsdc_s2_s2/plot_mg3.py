#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

import matplotlib as mpl

mpl.rcParams["font.size"] = 12

# name of the method
method = "mlsdc_s2_s2_rr3_5_7"

# --------------------------------------------------------------------------------------
# stability

fig1, ax = plt.subplots()
fig1.set_figheight(4.8)
fig1.set_figwidth(5.8)

x = np.loadtxt("lambda_re_amp.dat") 
y = np.loadtxt("lambda_im_amp.dat") 

a0  = np.loadtxt("amplification_level_3_mg3_cyc0.dat")
a1  = np.loadtxt("amplification_level_3_mg3_cyc1.dat")
a2  = np.loadtxt("amplification_level_3_mg3_cyc2.dat")
a3  = np.loadtxt("amplification_level_3_mg3_cyc3.dat")
a4  = np.loadtxt("amplification_level_3_mg3_cyc4.dat")
a5  = np.loadtxt("amplification_level_3_mg3_cyc5.dat")
a6  = np.loadtxt("amplification_level_3_mg3_cyc6.dat")
a7  = np.loadtxt("amplification_level_3_mg3_cyc7.dat")

s0  = plt.contour(x, y, a0 , levels=[1], linestyles=['-' ], linewidths=2, colors=["C1"])
s1  = plt.contour(x, y, a1 , levels=[1], linestyles=['--'], linewidths=2, colors=["C2"])
s2  = plt.contour(x, y, a2 , levels=[1], linestyles=['-.'], linewidths=2, colors=["C3"])
s3  = plt.contour(x, y, a3 , levels=[1], linestyles=[':' ], linewidths=2, colors=["C4"])
s4  = plt.contour(x, y, a4 , levels=[1], linestyles=['-' ], linewidths=2, colors=["C5"])
s5  = plt.contour(x, y, a5 , levels=[1], linestyles=['--'], linewidths=2, colors=["C6"])
s6  = plt.contour(x, y, a6 , levels=[1], linestyles=['-.'], linewidths=2, colors=["C7"])
s7  = plt.contour(x, y, a7 , levels=[1], linestyles=[':' ], linewidths=2, colors=["C8"])

h0 , l0  = s0.legend_elements()
h1 , l1  = s1.legend_elements()
h2 , l2  = s2.legend_elements()
h3 , l3  = s3.legend_elements()
h4 , l4  = s4.legend_elements()
h5 , l5  = s5.legend_elements()
h6 , l6  = s6.legend_elements()
h7 , l7  = s7.legend_elements()

handles = [ h0[0],
            h1[0],
            h2[0],
            h3[0],
            h4[0],
            h5[0],
            h6[0],
            h7[0]]
labels  = [ "K = 0",
            "K = 1",
            "K = 2",
            "K = 3",
            "K = 4",
            "K = 5",
            "K = 6",
            "K = 7"]

plt.legend(handles, labels, loc = "upper right")

#plt.xlabel(r"$z_\mathrm{r}$", fontsize=14)
plt.xlabel(r"$z_r$", fontsize=14)
#plt.ylabel(r"$z_\mathrm{i}$", fontsize=14)
plt.ylabel(r"$z_i$", fontsize=14)

plt.xlim(left=-10, right=5)
plt.ylim(bottom=0, top=20)

plt.axvline(0, color='black', linewidth=0.5, linestyle='-')

fig1.suptitle(r"MLSDC-S2$_{3,5,7}^{K+1}$ stability domains level 3, FMG")
fig1.tight_layout(rect=[-0.02, -0.02, 1.0, 1.05])
fig1.savefig(method+"_amp_mg3.pdf")

# --------------------------------------------------------------------------------------
# accuracy


fig2, ax = plt.subplots()
fig2.set_figheight(4.8)
fig2.set_figwidth(5.8)

x = np.loadtxt("lambda_re_acc.dat") 
y = np.loadtxt("lambda_im_acc.dat") 

e0  = np.loadtxt("error_level_3_acc_mg3_cyc0.dat")
e1  = np.loadtxt("error_level_3_acc_mg3_cyc1.dat")
e2  = np.loadtxt("error_level_3_acc_mg3_cyc2.dat")
e3  = np.loadtxt("error_level_3_acc_mg3_cyc3.dat")
e4  = np.loadtxt("error_level_3_acc_mg3_cyc4.dat")
e5  = np.loadtxt("error_level_3_acc_mg3_cyc5.dat")
e6  = np.loadtxt("error_level_3_acc_mg3_cyc6.dat")
e7  = np.loadtxt("error_level_3_acc_mg3_cyc7.dat")

s0  = plt.contour(x, y, e0 , levels=[1e-6], linestyles=['-' ], linewidths=2, colors=["C1"])
s1  = plt.contour(x, y, e1 , levels=[1e-6], linestyles=['--'], linewidths=2, colors=["C2"])
s2  = plt.contour(x, y, e2 , levels=[1e-6], linestyles=['-.'], linewidths=2, colors=["C3"])
s3  = plt.contour(x, y, e3 , levels=[1e-6], linestyles=[':' ], linewidths=2, colors=["C4"])
s4  = plt.contour(x, y, e4 , levels=[1e-6], linestyles=['-' ], linewidths=2, colors=["C5"])
s5  = plt.contour(x, y, e5 , levels=[1e-6], linestyles=['--'], linewidths=2, colors=["C6"])
s6  = plt.contour(x, y, e6 , levels=[1e-6], linestyles=['-.'], linewidths=2, colors=["C7"])
s7  = plt.contour(x, y, e7 , levels=[1e-6], linestyles=[':' ], linewidths=2, colors=["C8"])

h0 , l0  = s0.legend_elements()
h1 , l1  = s1.legend_elements()
h2 , l2  = s2.legend_elements()
h3 , l3  = s3.legend_elements()
h4 , l4  = s4.legend_elements()
h5 , l5  = s5.legend_elements()
h6 , l6  = s6.legend_elements()
h7 , l7  = s7.legend_elements()

handles = [ h0[0],
            h1[0],
            h2[0],
            h3[0],
            h4[0],
            h5[0],
            h6[0],
            h7[0] ]
labels  = [ "K = 0",
            "K = 1",
            "K = 2",
            "K = 3",
            "K = 4",
            "K = 5",
            "K = 6",
            "K = 7" ]

plt.legend(handles, labels, loc = "upper right")

#plt.xlabel(r"$z_\mathrm{r}$", fontsize=14)
plt.xlabel(r"$z_r$", fontsize=14)
#plt.ylabel(r"$z_\mathrm{i}$", fontsize=14)
plt.ylabel(r"$z_i$", fontsize=14)

plt.xlim(left=-2, right=2)
plt.ylim(bottom=0, top=3)

plt.axvline(0, color='black', linewidth=0.5, linestyle='-')

fig2.suptitle(r"MLSDC-S2$_{3,5,7}^{K+1}$ error domain, level 3, FMG")
fig2.tight_layout(rect=[-0.02, -0.02, 1.0, 1.05])
fig2.savefig(method+"_acc_mg3.pdf")

print("ready")
