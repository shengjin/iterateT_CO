#!/usr/bin/env python2
import sys
import os
import math
import time

import numpy as np

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window in plotting

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator

#sys.path.append('/turquoise/users/sjin/python/')


######### ######### ######### ######### ######### ######### 
###### the function that read one ..
def read_once(dir_name):
    import numpy as np
    au_jy_dust_23 = np.genfromtxt("./%s/Au_fluxJy_dust"%(dir_name), dtype=float)
    au_jy_12co = np.genfromtxt("./%s/Au_fluxJy_12CO"%(dir_name), dtype=float)
    au_jy_13co = np.genfromtxt("./%s/Au_fluxJy_13CO"%(dir_name), dtype=float)
    au_jy_c18o = np.genfromtxt("./%s/Au_fluxJy_C18O"%(dir_name), dtype=float)
    nrad = au_jy_12co.shape[0]
    intensity = np.zeros((nrad,4), dtype=float)
    intensity[0:23:1,0] = au_jy_dust_23[:,1]
    # NOTE: *3.0
    intensity[:,1] = au_jy_12co[:,1]*3.0
    intensity[:,2] = au_jy_13co[:,1]*3.0
    intensity[:,3] = au_jy_c18o[:,1]
    return intensity

# read the observed intensity and uncertainty
data = np.genfromtxt("./AU_fluxJy_noissJy_dust", dtype=float)
obs_jy_dust = np.zeros((data.shape[0],2), dtype=float)
obs_jy_dust[0:23:1,0] = data[:,1]
obs_jy_dust[0:23:1,1] = data[:,2]
data = np.genfromtxt("./AU_fluxJy_noissJy_12CO", dtype=float)
obs_jy_12co = np.zeros((data.shape[0],2), dtype=float)
obs_jy_12co[:,0] = data[:,1]
obs_jy_12co[:,1] = data[:,2]
data = np.genfromtxt("./AU_fluxJy_noissJy_13CO", dtype=float)
obs_jy_13co = np.zeros((data.shape[0],2), dtype=float)
obs_jy_13co[:,0] = data[:,1] 
obs_jy_13co[:,1] = data[:,2] 
data = np.genfromtxt("./AU_fluxJy_noissJy_C18O", dtype=float)
obs_jy_c18o = np.zeros((data.shape[0],2), dtype=float)
obs_jy_c18o[:,0] = data[:,1]
obs_jy_c18o[:,1] = data[:,2]
    
# read simulation grid 
run_grid = np.genfromtxt("./run_grid.txt", dtype=float)

# create array to store the chi2 of each model
# chi2 = Sum_over_r (I(r)^2 - Im(r)^2)/dI(r)^2
chi2_models = np.zeros((run_grid.shape[0],5), dtype=float)


for i in range(run_grid.shape[0]):
    #dir_name = '0fit'
    dir_name = '%04d'%(run_grid[i,0])
    model_intensity = read_once(dir_name)
    chi2_once = (obs_jy_dust[0:23:1,0]-model_intensity[0:23:1,0])**2.0/obs_jy_dust[:,1]**2.0
    for j in range(chi2_once.shape[0]):
        if j < 23:
            chi2_models[i,0] = chi2_models[i,0] + chi2_once[j]
            #print j, chi2_models[i,0]
    chi2_models[i,0] = chi2_models[i,0]/(23-5-1)

    chi2_once = (obs_jy_12co[:,0]-model_intensity[:,1])**2.0/obs_jy_12co[:,1]**2.0
    chi2_models[i,1] = chi2_once.sum(axis=0)/(60-5-1)

    chi2_once = (obs_jy_13co[:,0]-model_intensity[:,2])**2.0/obs_jy_13co[:,1]**2.0
    chi2_models[i,2] = chi2_once.sum(axis=0)/(60-5-1)

    chi2_once = (obs_jy_c18o[:,0]-model_intensity[:,3])**2.0/obs_jy_c18o[:,1]**2.0
    chi2_models[i,3] = chi2_once.sum(axis=0)/(60-5-1)

    chi2_models[i,4] = chi2_models[i,0]+chi2_models[i,1]+chi2_models[i,2]+chi2_models[i,3]
    print i, dir_name, chi2_models[i,:]


np.savetxt('chi2.txt', np.transpose([chi2_models[:,0], chi2_models[:,1], chi2_models[:,2], chi2_models[:,3], chi2_models[:,4]]))

