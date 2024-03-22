#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 12
})

# name of the method
method = "sdc_eu-eu_35"

# number of subintervals
ns = 3
# standalone or predictor method stages
sp = 2
# standalone or predictor method stages
sc = 2
# number of correction sweeps
nc = 5

# work units
nw = max(ns,1) * (max(sp-1,1) + nc * max(sc-1,1))

print("Accuracy evaluation of ", method,"\n")
print("  subintervals     =", ns)
print("  predictor stages =", sp)
print("  corrector stages =", sc)
print("  corrector sweeps =", nc)
print("  works units      =", nw, "\n")

# scaled 1e-05 -----------------------------------------------------------------
                       
x = np.loadtxt("lambda_re.dat") / nw
y = np.loadtxt("lambda_im.dat") / nw
e = np.loadtxt("error.dat")

fig1, ax = plt.subplots()

diag = plt.contour( x, y, e, levels=[1.e-5], colors='red' 
                  , linestyles='-', linewidths= 0.75 ) 

handles, labels = diag.legend_elements()

labels = [r"$\varepsilon = 10^{-5}$"]
plt.legend(handles, labels)

plt.xlim(left=-0.12, right=0.12)
plt.ylim(bottom=0.00, top=0.13)

plt.xlabel(r"Re($z_{\mathrm{s}}$)")
plt.ylabel(r"Im($z_{\mathrm{s}}$)")

fig1.tight_layout(rect=[0, 0.0, 1, 1.0])
fig1.savefig("accuracy-05_"+method+".pdf")


# scaled 1e-10 -----------------------------------------------------------------

fig2, ax = plt.subplots()

diag = plt.contour( x, y, e, levels=[1.e-10], colors='red' 
                  , linestyles='-', linewidths= 0.75 ) 

handles, labels = diag.legend_elements()

labels = [r"$\varepsilon = 10^{-10}$"]
plt.legend(handles, labels)

ax.set_xlim([-0.012, 0.012])
ax.set_ylim([ 0.000, 0.013])

plt.xlabel(r"Re($z_{\mathrm{s}}$)")
plt.ylabel(r"Im($z_{\mathrm{s}}$)")

fig2.tight_layout(rect=[0, 0.0, 1, 1.0])
fig2.savefig("accuracy-10_"+method+".pdf")


print("ready")
