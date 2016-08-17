#!/usr/bin/env python2
import numpy as np

gamma = [0.5, 1, 1.5, 2]
rc  = [12500]
rc_atan  = [5, 15, 25, 35]
gamma_atan = [1, 2, 3, 4]
Sigma_0 = [1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2]

gamma = np.asarray(gamma)
rc = np.asarray(rc)
gamma_atan = np.asarray(gamma_atan)
rc_atan = np.asarray(rc_atan)
Sigma_0 = np.asarray(Sigma_0)

chi2_all = np.genfromtxt("./grid_chi2.txt", dtype=float)

grid_S0 = chi2_all[:,1]
grid_rc = chi2_all[:,2]
grid_gamma = chi2_all[:,3]
grid_rc_atan = chi2_all[:,4]
grid_gamma_atan = chi2_all[:,5]

chi2_dust = chi2_all[:,6]
chi2_12co = chi2_all[:,7]
chi2_13co = chi2_all[:,8]
chi2_c18o = chi2_all[:,9]
chi2_sum = chi2_all[:,10]

print grid_S0.shape
print len(Sigma_0),len(gamma),len(gamma),len(rc_atan),len(gamma_atan)

grid_S0 = np.array(grid_S0).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
grid_rc = np.array(grid_rc).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
grid_gamma = np.array(grid_gamma).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
grid_rc_atan = np.array(grid_rc_atan).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
grid_gamma_atan = np.array(grid_gamma_atan).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))

chi2_dust = np.array(chi2_dust).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
chi2_12co = np.array(chi2_12co).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
chi2_13co = np.array(chi2_13co).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
chi2_c18o = np.array(chi2_c18o).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))
chi2_sum = np.array(chi2_sum).reshape(len(Sigma_0),len(rc),len(gamma),len(rc_atan),len(gamma_atan))


import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from matplotlib.colors import SymLogNorm


rc=0
rc_atan=NUM_I
gamma_atan=NUM_J


x=grid_gamma[:,rc,:,rc_atan,gamma_atan]
y=grid_S0[:,rc,:,rc_atan,gamma_atan]
x=np.array(x).reshape(len(Sigma_0)*len(gamma))
y=np.array(y).reshape(len(Sigma_0)*len(gamma))

color = chi2_12co[:,rc,:,rc_atan,gamma_atan]
color=np.array(color).reshape(len(Sigma_0)*len(gamma))
#plt.xlim(0,400)
plt.ylim(5.e-4,1e-1)
plt.xlabel("gamma")
plt.ylabel("Mdisk (Msun)")
plt.xticks(np.arange(0.0, 2.5, 0.5))
#plt.yticks(np.arange(-0.5, 1.5, 0.5))
plt.scatter( x, y, c=color, s=100)
plt.yscale('log')
print x
print y
print color
plt.colorbar(label='reduced-chi2 of 12co')
figname = "%d%s%d%s" % (rc_atan,"_",gamma_atan, ".png")
plt.savefig(figname)

plt.clf()

