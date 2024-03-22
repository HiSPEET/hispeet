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
ns = 4
# standalone or predictor method stages
sp = 2
# standalone or corrector method stages
sc = 2
# number of correction sweeps
nc = 6

# work units
nw = max(ns,1) * (max(sp-1,1) + nc * max(sc-1,1))

print("Stability evaluation of ", method,"\n")
print("  subintervals     =", ns)
print("  predictor stages =", sp)
print("  corrector stages =", sc)
print("  corrector sweeps =", nc)
print("  works units      =", nw, "\n")

# unscaled ---------------------------------------------------------------------
                       
x = np.loadtxt("lambda_re.dat") 
y = np.loadtxt("lambda_im.dat") 
a = np.loadtxt("amplification.dat")

fig1, ax = plt.subplots()

diag = plt.contour( x, y, a, levels=[1], colors='red' 
                  , linestyles='-', linewidths=0.5 ) 

#diag = plt.contourf( x, y, a, levels=[0,1], colors=['green','white'] 
#                  , linestyles='-', linewidths= 0.75, extend='both' ) 

handles, labels = diag.legend_elements()

plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')

labels = [r"$|R| = 1$"]
plt.legend(handles, labels, loc = "upper left")

plt.xlim(left=-10, right=2)
plt.ylim(bottom=0, top=8)

plt.xlabel(r"Re($z$)")
plt.ylabel(r"Im($z$)")

#plt.grid(which='both', axis='y')

fig1.tight_layout(rect=[0, 0.0, 1, 1.0])
fig1.savefig("stability_"+method+".pdf")

# scaled -----------------------------------------------------------------------

fig2, ax = plt.subplots()

diag = plt.contour( x/nw, y/nw, a, levels=[1], colors='red' 
                  , linestyles='-', linewidths= 0.75 ) 

handles, labels = diag.legend_elements()

plt.axvline(0, color='grey', linewidth=0.5, linestyle='--')

labels = [r"$|R| = 1$"]
plt.legend(handles, labels, loc = "upper left")

plt.xlim(left=-0.6, right=0.2)
plt.ylim(bottom=0.0, top=1.6)

plt.xlabel(r"Re($z_{\mathrm{s}}$)")
plt.ylabel(r"Im($z_{\mathrm{s}}$)")

fig2.tight_layout(rect=[0, 0.0, 1, 1.0])
fig2.savefig("stability-scaled_"+method+".pdf")

print("ready")
