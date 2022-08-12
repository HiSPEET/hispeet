#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

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
                  , linestyles='-', linewidths= 0.75 ) 

lines  = [ diag.collections[0] ]
labels = [ method ]

plt.legend(lines, labels)

ax.set_xlim([-10, 2])
ax.set_ylim([  0, 8])

plt.xlabel("Re($z$)")
plt.ylabel("Im($z$)")

fig1.savefig("stability_"+method+".pdf")

# scaled -----------------------------------------------------------------------

fig2, ax = plt.subplots()

diag = plt.contour( x/nw, y/nw, a, levels=[1], colors='red' 
                  , linestyles='-', linewidths= 0.75 ) 

lines  = [ diag.collections[0] ]
labels = [ method ]

plt.legend(lines, labels)

ax.set_xlim([-0.6, 0.2])
ax.set_ylim([ 0.0, 1.6])

plt.xlabel("Re($z_{\mathrm{s}}$)")
plt.ylabel("Im($z_{\mathrm{s}}$)")

fig2.savefig("stability-scaled_"+method+".pdf")

print("ready")
